"""CID quantization tests.

Determinism contract: "quantized-equal". Every emitted float column is rounded
at the pipeline output boundary to one digit below its display precision (and
at least 4 dp). Rounding stays below the `.Nf` display format, so users see no
change, while the upcoming BioPython → numpy code-path swap's floating-point
drift is absorbed — a same-machine re-run produces byte-identical output and
the workflow's content-addressable id stays stable.

Per-family decimals (see `CID_QUANTIZE_DECIMALS_BY_PREFIX`): charge / chargeShift
/ pi at 3 dp; instability / mw / aliphatic at 4 dp; gravy / aromaticity at 5 dp;
eox / ered at 0 dp (integer-valued, exact). The aa_fraction frame's `value`
column rounds to 5 dp.

These tests guard the boundary behaviour. Internal property functions
(tested elsewhere) keep full precision so golden-value tests stay sharp.
"""

from __future__ import annotations

import math

import polars as pl
import pytest

from pipeline import (
    _decimals_for,
    _quantize_for_cid,
    run,
)


class TestQuantizeHelper:
    """Direct tests of the boundary helper."""

    # Quantization rounds every float family to its per-column decimals.
    # charge/chargeShift/pi at 3 dp; instability/mw/aliphatic at 4 dp;
    # gravy/aromaticity at 5 dp; eox/ered at 0 dp (integer-valued, exact).
    def test_all_float_families_rounded(self):
        df = pl.DataFrame(
            {
                "entity_key": ["x"],
                "charge_peptide": [0.1234567],
                "pi_peptide": [7.123456],
                "gravy_peptide": [0.4142857],
                "mw_peptide": [1234.56789],
                "instability_peptide": [38.7531234],
                "aliphatic_peptide": [61.296296],
                "aromaticity_peptide": [0.185185],
                "eox_peptide": [22460.0],
                "ered_peptide": [22460.0],
            }
        )
        out = _quantize_for_cid(df)
        # charge / pi — 3 dp:
        assert out["charge_peptide"][0] == round(0.1234567, 3)
        assert out["pi_peptide"][0] == round(7.123456, 3)
        # instability / mw / aliphatic — 4 dp:
        assert out["mw_peptide"][0] == round(1234.56789, 4)
        assert out["instability_peptide"][0] == round(38.7531234, 4)
        assert out["aliphatic_peptide"][0] == round(61.296296, 4)
        # gravy / aromaticity — 5 dp:
        assert out["gravy_peptide"][0] == round(0.4142857, 5)
        assert out["aromaticity_peptide"][0] == round(0.185185, 5)
        # eox / ered — 0 dp, integer-valued, exact:
        assert out["eox_peptide"][0] == 22460.0
        assert out["ered_peptide"][0] == 22460.0

    # Per-chain antibody columns share the same prefix matching — `charge_A_CDR3`,
    # `charge_A_VDJRegion`, `charge_Fv`, `pi_A_VDJRegion`, `pi_Fv` all round.
    def test_antibody_charge_and_pi_columns_match_prefix(self):
        df = pl.DataFrame(
            {
                "entity_key": ["c"],
                "charge_A_CDR3": [1.234567],
                "charge_B_CDR3": [-0.987654],
                "charge_A_VDJRegion": [3.141592],
                "charge_B_VDJRegion": [2.718281],
                "charge_Fv": [5.859873],
                "pi_A_VDJRegion": [7.018372],
                "pi_B_VDJRegion": [9.798889],
                "pi_Fv": [9.330627],
                # Other float families now round too — gravy 5 dp, mw 4 dp:
                "gravy_A_VDJRegion": [-0.11111156],
                "mw_A_VDJRegion": [6050.730289],
            }
        )
        out = _quantize_for_cid(df)
        for col in (
            "charge_A_CDR3",
            "charge_B_CDR3",
            "charge_A_VDJRegion",
            "charge_B_VDJRegion",
            "charge_Fv",
            "pi_A_VDJRegion",
            "pi_B_VDJRegion",
            "pi_Fv",
        ):
            v = out[col][0]
            assert v == pytest.approx(round(v, _decimals_for(col)), abs=0)
        # gravy / mw round to their per-column decimals (5 dp / 4 dp):
        assert out["gravy_A_VDJRegion"][0] == round(-0.11111156, 5)
        assert out["mw_A_VDJRegion"][0] == round(6050.730289, 4)

    # No-op when no column matches any known property prefix. Uses a fabricated
    # column name so the genuine no-match branch stays covered.
    def test_passthrough_when_no_matching_columns(self):
        df = pl.DataFrame({"entity_key": ["x"], "notaproperty_x": [0.123456789]})
        out = _quantize_for_cid(df)
        assert out["notaproperty_x"][0] == 0.123456789

    # Sanity — the per-prefix decimals table (via `_decimals_for`) maps each
    # property family to the decimals the docstring promises. This is the single
    # source of truth the pipeline rounds against; a column with no rule returns
    # None (left full-precision).
    def test_decimals_for_tracks_documented_values(self):
        documented = {
            "charge_peptide": 3,
            "chargeShift_peptide": 3,
            "pi_peptide": 3,
            "instability_peptide": 4,
            "mw_peptide": 4,
            "aliphatic_peptide": 4,
            "gravy_peptide": 5,
            "aromaticity_peptide": 5,
            "eox_peptide": 0,
            "ered_peptide": 0,
        }
        for col, dp in documented.items():
            assert _decimals_for(col) == dp, f"{col}: expected {dp} dp"
        # A non-property column has no quantization rule.
        assert _decimals_for("entity_key") is None

    # Signed-zero canonicalization: a `-0.0` input (FP-residual-sign drift on a
    # ~0 property) must emit as `+0.0`, so the TSV writer produces identical bytes
    # regardless of summation order. `-0.0 == 0.0` numerically, so the guard is
    # on the SIGN bit, not equality.
    def test_negative_zero_canonicalized_to_positive_zero(self):
        out = _quantize_for_cid(pl.DataFrame({"entity_key": ["x"], "gravy_x": [-0.0]}))
        v = out["gravy_x"][0]
        assert v == 0.0
        assert math.copysign(1.0, v) == 1.0, f"expected +0.0, got {v!r} (negative-zero bit set)"


class TestPipelineQuantizationApplied:
    """Quantization fires at the pipeline boundary, not just in the helper."""

    # Peptide pipeline output: charge / pi land at their 3-dp boundary.
    def test_peptide_run_rounds_charge_and_pi(self):
        reads = pl.DataFrame(
            {
                "entity_key": ["p1"],
                "sequence": ["ACDEFGHIKL"],  # 10 aa — every property defined
            }
        )
        out = run(reads, {"mode": "peptide"})
        row = out["properties"].row(0, named=True)

        # charge_peptide and pi_peptide are bisection / libm-derived; rounded.
        for c in ("charge_peptide", "pi_peptide"):
            v = row[c]
            assert v is not None
            assert v == pytest.approx(round(v, _decimals_for(c)), abs=0), f"{c}={v} not at 3-decimal precision"

    # Antibody full-coverage output: charge_*, pi_*, including Fv, all land at
    # their 3-dp boundary. (Other families round too — covered by the helper
    # tests; here we assert the 3-dp family specifically.)
    def test_antibody_run_rounds_all_charge_and_pi_columns(
        self, antibody_full_one_clone: pl.DataFrame, antibody_full_plan: dict
    ):
        out = run(antibody_full_one_clone, antibody_full_plan)
        row = out["properties"].row(0, named=True)

        # The 3-dp family (charge / chargeShift / pi) — every such column lands
        # exactly on its boundary.
        rounded_cols = [c for c in row if _decimals_for(c) == 3]
        # Sanity: at least one charge_ and one pi_ column present.
        assert any(c.startswith("charge_") for c in rounded_cols)
        assert any(c.startswith("pi_") for c in rounded_cols)
        for c in rounded_cols:
            v = row[c]
            if v is None:
                continue
            assert v == pytest.approx(round(v, _decimals_for(c)), abs=0), f"{c}={v} not at 3-decimal precision"


class TestQuantizationDoesNotPropagateInternally:
    """Internal property functions stay full-precision — golden values remain
    valid against `properties.py` direct calls. The pipeline boundary is the
    only place quantization happens.
    """

    def test_isoelectric_point_returns_unrounded_value(self):
        from pka_tables import IPC2_PROTEIN
        from properties import isoelectric_point

        # Same VH chain used in test_golden_values; pi pinned at 6.006653
        # under IPC 2.0 protein pKa. The internal function must keep digits
        # beyond the 3rd decimal — quantization is a boundary concern only.
        vh = "EVQLVESGFTFSSYAMSWVRQISGSGGSTYYAESVKGRFTICARDYWWGQGTLV"
        pi = isoelectric_point(vh, IPC2_PROTEIN, include_cys=False)
        assert pi == pytest.approx(6.006653, abs=1e-6)
        # The boundary rounds pi to its 3-dp family; the internal function must
        # keep more digits than that.
        assert pi != round(pi, _decimals_for("pi_peptide"))


class TestExtendedQuantizationBoundary:
    # Every float column rounds; integer-valued ε columns are unaffected by 0-dp rounding.
    def test_all_float_columns_rounded(self):
        from pipeline import _quantize_for_cid

        df = pl.DataFrame(
            {
                "entity_key": ["x"],
                "charge_peptide": [0.12345678],
                "gravy_peptide": [0.41428571],
                "mw_peptide": [1234.5678912],
                "instability_peptide": [38.75312345],
                "aliphatic_peptide": [61.29629629],
                "aromaticity_peptide": [0.18518518],
                "eox_peptide": [22460.0],
            }
        )
        out = _quantize_for_cid(df)
        assert out["charge_peptide"][0] == round(0.12345678, 3)
        assert out["gravy_peptide"][0] == round(0.41428571, 5)
        assert out["mw_peptide"][0] == round(1234.5678912, 4)
        assert out["instability_peptide"][0] == round(38.75312345, 4)
        assert out["aliphatic_peptide"][0] == round(61.29629629, 4)
        assert out["aromaticity_peptide"][0] == round(0.18518518, 5)
        assert out["eox_peptide"][0] == 22460.0

    # Rounding stays below display precision: the .Nf-formatted value is unchanged.
    @pytest.mark.parametrize(
        "col, value, display_dp",
        [
            ("gravy_peptide", 0.4142857, 3),
            ("mw_peptide", 1234.56789, 1),
            ("instability_peptide", 38.7531234, 2),
        ],
    )
    def test_display_precision_unchanged(self, col, value, display_dp):
        from pipeline import _quantize_for_cid

        out = _quantize_for_cid(pl.DataFrame({"entity_key": ["x"], col: [value]}))
        assert round(out[col][0], display_dp) == round(value, display_dp)


class TestAaFractionQuantization:
    """The aa_fraction frame's `value` column is rounded at the boundary too."""

    # Peptide run rounds aa_fraction value to 5 dp (one below the .3f display).
    # A length-7 sequence yields repeating-decimal fractions (1/7, 3/7) that are
    # NOT at 5 dp pre-quantization — so this exercises the rounding, not a no-op.
    def test_aa_fraction_value_rounded_to_5dp(self):
        reads = pl.DataFrame(
            {
                "entity_key": ["p1"],
                "sequence": ["AAACDEF"],  # 3/7, 1/7 — repeating decimals
            }
        )
        out = run(reads, {"mode": "peptide"})
        values = [v for v in out["aa_fraction"]["value"].to_list() if v is not None]
        assert values  # non-empty
        # Pre-quantization at least one value carries > 5 dp; all must end at 5 dp.
        for v in values:
            assert v == round(v, 5), f"aa_fraction value {v} not at 5-decimal precision"
