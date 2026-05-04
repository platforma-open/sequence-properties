"""CID quantization tests.

Charge and pI columns depend on `10**x` via libm — a real source of ULP
variance across library upgrades or Python→numpy code-path swaps. The
pipeline rounds these columns to bisection-tolerance precision (1e-3) at
the output boundary so a same-machine re-run produces byte-identical
output and the workflow's content-addressable id stays stable.

Every other property is closed-form integer / constant arithmetic and
bit-exact under IEEE-754; quantization there would be churn.

These tests guard the boundary behaviour. Internal property functions
(tested elsewhere) keep full precision so golden-value tests stay sharp.
"""

from __future__ import annotations

import polars as pl
import pytest

from pipeline import (
    CID_QUANTIZE_DECIMALS,
    CID_QUANTIZE_PREFIXES,
    _quantize_for_cid,
    run,
)


class TestQuantizeHelper:
    """Direct tests of the boundary helper."""

    # Quantization rounds only charge_* and pi_* columns; everything else
    # passes through bit-identical.
    def test_only_charge_and_pi_are_rounded(self):
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
        # Rounded:
        assert out["charge_peptide"][0] == pytest.approx(0.123, abs=0)
        assert out["pi_peptide"][0] == pytest.approx(7.123, abs=0)
        # Untouched (each value retains every digit it carried in):
        assert out["gravy_peptide"][0] == 0.4142857
        assert out["mw_peptide"][0] == 1234.56789
        assert out["instability_peptide"][0] == 38.7531234
        assert out["aliphatic_peptide"][0] == 61.296296
        assert out["aromaticity_peptide"][0] == 0.185185
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
                # Non-quantized — must pass through:
                "gravy_A_VDJRegion": [-0.111111],
                "mw_A_VDJRegion": [6050.7302],
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
            assert v == pytest.approx(round(v, CID_QUANTIZE_DECIMALS), abs=0)
        # Non-quantized columns retain full precision:
        assert out["gravy_A_VDJRegion"][0] == -0.111111
        assert out["mw_A_VDJRegion"][0] == 6050.7302

    # No-op when nothing matches the prefix list.
    def test_passthrough_when_no_matching_columns(self):
        df = pl.DataFrame({"entity_key": ["x"], "gravy_peptide": [0.123456789]})
        out = _quantize_for_cid(df)
        assert out["gravy_peptide"][0] == 0.123456789

    # Sanity — module-level constants match what the docstring promises.
    def test_constants_track_documented_values(self):
        assert CID_QUANTIZE_DECIMALS == 3
        assert CID_QUANTIZE_PREFIXES == ("charge_", "pi_")


class TestPipelineQuantizationApplied:
    """Quantization fires at the pipeline boundary, not just in the helper."""

    # Peptide pipeline output: charge / pi rounded, others not.
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
            assert v == pytest.approx(round(v, CID_QUANTIZE_DECIMALS), abs=0), f"{c}={v} not at 3-decimal precision"

    # Antibody full-coverage output: charge_*, pi_*, including Fv, all rounded.
    # Non-rounded columns (gravy / mw / instability / aliphatic / aromaticity / ε)
    # must keep full precision.
    def test_antibody_run_rounds_all_charge_and_pi_columns(
        self, antibody_full_one_clone: pl.DataFrame, antibody_full_plan: dict
    ):
        out = run(antibody_full_one_clone, antibody_full_plan)
        row = out["properties"].row(0, named=True)

        rounded_cols = [c for c in row if any(c.startswith(p) for p in CID_QUANTIZE_PREFIXES)]
        # Sanity: at least one charge_ and one pi_ column present.
        assert any(c.startswith("charge_") for c in rounded_cols)
        assert any(c.startswith("pi_") for c in rounded_cols)
        for c in rounded_cols:
            v = row[c]
            if v is None:
                continue
            assert v == pytest.approx(round(v, CID_QUANTIZE_DECIMALS), abs=0), f"{c}={v} not at 3-decimal precision"


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
        assert pi != round(pi, CID_QUANTIZE_DECIMALS)
