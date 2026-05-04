"""Behavioral tests for src/pipeline.py.

Verifies mode dispatch (peptide / antibody-tcr), expected output column sets,
NA propagation per-clone for missing regions, and Fv emission gating.

Run from blocks/sequence-properties/software/:
    uv sync
    uv run pytest tests/unit/test_pipeline.py
"""

from __future__ import annotations

import logging

import polars as pl
import pytest

from pipeline import (
    CDR3_PROPS,
    FULL_CHAIN_PROPS,
    FV_PROPS,
    PEPTIDE_PROPERTY_COLUMNS,
    run,
)

# ---------------------------------------------------------------------------
# Peptide mode
# ---------------------------------------------------------------------------


class TestPeptideMode:
    """Peptide pipeline: 9 scalar columns + AA-fraction long frame."""

    @pytest.fixture
    def reads(self) -> pl.DataFrame:
        return pl.DataFrame(
            {
                "entity_key": ["pep1", "pep2", "pep3"],
                "sequence": ["ACDEFGHIKL", "AAAAAAAAAA", ""],
            }
        )

    def test_emits_all_peptide_columns(self, reads: pl.DataFrame):
        out = run(reads, {"mode": "peptide"})
        assert "properties" in out
        properties = out["properties"]
        for col in PEPTIDE_PROPERTY_COLUMNS:
            assert col in properties.columns

    # Empty peptide → every property NA for that row.
    def test_empty_sequence_produces_all_na_row(self, reads: pl.DataFrame):
        out = run(reads, {"mode": "peptide"})
        row = out["properties"].filter(pl.col("entity_key") == "pep3").row(0, named=True)
        for col in PEPTIDE_PROPERTY_COLUMNS:
            assert row[col] is None

    # 10-residue all-A → instability index = 9.0 (closed form).
    def test_polyalanine_instability_value(self, reads: pl.DataFrame):
        out = run(reads, {"mode": "peptide"})
        row = out["properties"].filter(pl.col("entity_key") == "pep2").row(0, named=True)
        assert row["instability_peptide"] == pytest.approx(9.0)

    # AA fraction frame: 20 rows × N entities, sums to 1.0 per non-empty entity.
    def test_aa_fraction_long_format(self, reads: pl.DataFrame):
        out = run(reads, {"mode": "peptide"})
        aa = out["aa_fraction"]
        assert set(aa.columns) == {"entity_key", "aminoAcid", "value"}
        # 3 entities × 20 AAs = 60 rows.
        assert aa.height == 60
        # Sum of fractions for non-empty entity equals 1.0 (modulo float).
        pep1_sum = aa.filter(pl.col("entity_key") == "pep1")["value"].sum()
        assert pep1_sum == pytest.approx(1.0)
        # Empty entity has all-NA values.
        pep3_values = aa.filter(pl.col("entity_key") == "pep3")["value"].to_list()
        assert all(v is None for v in pep3_values)

    # R7 shape contract: an invalid-but-non-empty sequence (stop codon) must still
    # emit exactly 20 rows of None — same shape as an empty cell. The 2-axis
    # PColumn requires uniform width across entities.
    def test_aa_fraction_stop_codon_emits_20_na_rows(self):
        reads = pl.DataFrame(
            {
                "entity_key": ["valid", "stop"],
                "sequence": ["ACDEFGHIKL", "ACDE*FGH"],
            }
        )
        out = run(reads, {"mode": "peptide"})
        stop_rows = out["aa_fraction"].filter(pl.col("entity_key") == "stop")
        assert stop_rows.height == 20
        assert all(v is None for v in stop_rows["value"].to_list())
        # And every property NA for the stop-codon row.
        prop_row = out["properties"].filter(pl.col("entity_key") == "stop").row(0, named=True)
        for col in PEPTIDE_PROPERTY_COLUMNS:
            assert prop_row[col] is None

    # Spec L459: "Normalize all sequences to uppercase". Lowercase input must
    # produce identical output to the same TSV uppercased. Guards against the
    # silent failure mode where lowercase residues get filtered as non-standard.
    def test_lowercase_input_normalized(self):
        upper = pl.DataFrame({"entity_key": ["p"], "sequence": ["ACDEFGHIKL"]})
        lower = pl.DataFrame({"entity_key": ["p"], "sequence": ["acdefghikl"]})
        out_upper = run(upper, {"mode": "peptide"})
        out_lower = run(lower, {"mode": "peptide"})
        for col in PEPTIDE_PROPERTY_COLUMNS:
            up = out_upper["properties"].row(0, named=True)[col]
            lo = out_lower["properties"].row(0, named=True)[col]
            assert up == pytest.approx(lo) if up is not None else lo is None
        # AA fractions identical too.
        u_aa = out_upper["aa_fraction"].sort(["entity_key", "aminoAcid"])
        l_aa = out_lower["aa_fraction"].sort(["entity_key", "aminoAcid"])
        assert u_aa.equals(l_aa)


# ---------------------------------------------------------------------------
# Antibody/TCR mode
# ---------------------------------------------------------------------------


class TestAntibodyTcrCdr3Only:
    """CDR3-only input — only CDR3 columns emitted, no full-chain or Fv."""

    @pytest.fixture
    def reads(self) -> pl.DataFrame:
        return pl.DataFrame(
            {
                "entity_key": ["c1", "c2"],
                "A_CDR3": ["CARDYW", "CARGFW"],
                "B_CDR3": ["CQQYNS", "CQHFSS"],
            }
        )

    @pytest.fixture
    def plan(self) -> dict:
        return {
            "mode": "antibody_tcr_legacy_sc",
            "receptor": "IG",
            "chains": ["A", "B"],
            "fullChains": [],
            "hasFv": False,
        }

    def test_only_cdr3_columns_emitted(self, reads: pl.DataFrame, plan: dict):
        out = run(reads, plan)
        cols = set(out["properties"].columns)
        # entity_key + CDR3 props × {A, B}
        expected = {"entity_key"} | {f"{p}_{ch}_CDR3" for p in CDR3_PROPS for ch in "AB"}
        assert cols == expected

    def test_full_chain_columns_absent(self, reads: pl.DataFrame, plan: dict):
        out = run(reads, plan)
        for ch in "AB":
            for p in FULL_CHAIN_PROPS:
                assert f"{p}_{ch}_VDJRegion" not in out["properties"].columns

    def test_fv_columns_absent(self, reads: pl.DataFrame, plan: dict):
        out = run(reads, plan)
        for p in FV_PROPS:
            assert f"{p}_Fv" not in out["properties"].columns


class TestAntibodyTcrFullCoverage:
    """Both VH and VL fully covered + IG receptor → Fv columns emitted."""

    def test_emits_full_chain_and_fv_columns(self, antibody_full_two_clones: pl.DataFrame, antibody_full_plan: dict):
        out = run(antibody_full_two_clones, antibody_full_plan)
        cols = set(out["properties"].columns)
        for ch in "AB":
            for p in FULL_CHAIN_PROPS:
                assert f"{p}_{ch}_VDJRegion" in cols
        for p in FV_PROPS:
            assert f"{p}_Fv" in cols

    def test_full_chain_pi_in_range(self, antibody_full_two_clones: pl.DataFrame, antibody_full_plan: dict):
        out = run(antibody_full_two_clones, antibody_full_plan)
        row = out["properties"].filter(pl.col("entity_key") == "c1").row(0, named=True)
        # Variable-region pI should sit in [0, 14] when defined.
        assert 0.0 < row["pi_A_VDJRegion"] < 14.0
        assert 0.0 < row["pi_B_VDJRegion"] < 14.0


class TestAntibodyTcrPartialClone:
    """Per-clone partial coverage — full-chain columns NA for that clone only,
    not for the whole dataset."""

    def test_missing_region_yields_na_for_that_clone(self):
        # Clone c1 has all 7 regions for chain A. Clone c2 is missing FR1.
        regions = {
            "A_FR1": ["EVQLVES", ""],
            "A_CDR1": ["GFTFSSY", "GFTFSSY"],
            "A_FR2": ["AMSWVRQ", "AMSWVRQ"],
            "A_CDR2": ["ISGSGGS", "ISGSGGS"],
            "A_FR3": ["TYYAESVKGRFTI", "TYYAESVKGRFTI"],
            "A_CDR3": ["CARDYW", "CARGFW"],
            "A_FR4": ["WGQGTLV", "WGQGTLV"],
        }
        reads = pl.DataFrame({"entity_key": ["c1", "c2"], **regions})
        plan = {
            "mode": "antibody_tcr_legacy_bulk",
            "receptor": "IG",
            "chains": ["A"],
            "fullChains": ["A"],
            "hasFv": False,
        }
        out = run(reads, plan)
        rows = {r["entity_key"]: r for r in out["properties"].iter_rows(named=True)}
        # c1 — full chain reconstructable, all values present.
        assert rows["c1"]["mw_A_VDJRegion"] is not None
        # c2 — FR1 missing, all full-chain values NA.
        for p in FULL_CHAIN_PROPS:
            assert rows["c2"][f"{p}_A_VDJRegion"] is None
        # CDR3 still computed for c2 (CDR3 column present).
        assert rows["c2"]["charge_A_CDR3"] is not None

    # Spec edge cases table: stop codon `*` invalidates the whole sequence. When
    # `*` lands in any region (FR or CDR), the reconstructed full chain contains
    # `*` and `is_invalid_sequence` NAs every full-chain property for that clone.
    # CDR3 columns are independent — clean CDR3 must still produce CDR3 props.
    def test_stop_codon_in_fr_nas_full_chain_only(self):
        regions = {
            "A_FR1": ["EVQLVES", "EVQLVES"],
            "A_CDR1": ["GFTFSSY", "GFTFSSY"],
            "A_FR2": ["AMSWVRQ", "AMS*VRQ"],  # c2 has stop codon in FR2
            "A_CDR2": ["ISGSGGS", "ISGSGGS"],
            "A_FR3": ["TYYAESVKGRFTI", "TYYAESVKGRFTI"],
            "A_CDR3": ["CARDYW", "CARGFW"],
            "A_FR4": ["WGQGTLV", "WGQGTLV"],
        }
        reads = pl.DataFrame({"entity_key": ["c1", "c2"], **regions})
        plan = {
            "mode": "antibody_tcr_legacy_bulk",
            "receptor": "IG",
            "chains": ["A"],
            "fullChains": ["A"],
            "hasFv": False,
        }
        out = run(reads, plan)
        rows = {r["entity_key"]: r for r in out["properties"].iter_rows(named=True)}
        # c2 — `*` in FR2 invalidates every full-chain property.
        for p in FULL_CHAIN_PROPS:
            assert rows["c2"][f"{p}_A_VDJRegion"] is None
        # c2 CDR3 is clean ("CARGFW") — CDR3 props must still compute.
        assert rows["c2"]["charge_A_CDR3"] is not None
        assert rows["c2"]["gravy_A_CDR3"] is not None
        # c1 unaffected.
        assert rows["c1"]["mw_A_VDJRegion"] is not None


class TestSingleCellChainDropout:
    """Single-cell chain dropout — clone has all 7 regions for one chain but the
    other chain is entirely empty for that clone. Per-chain columns emit values
    for the present chain and NA for the absent one; Fv NA for that clone.

    Distinct from CDR3-only (R11a) and partial-region (R11b) — those are
    dataset-level shapes. This is per-clone dropout common in 10x sc data.
    """

    def test_chain_b_dropout_for_one_clone(self):
        regions = {
            # c1: both chains complete. c2: chain B every region empty.
            "A_FR1": ["EVQLVES", "EVQLVES"],
            "A_CDR1": ["GFTFSSY", "GFTFSSY"],
            "A_FR2": ["AMSWVRQ", "AMSWVRQ"],
            "A_CDR2": ["ISGSGGS", "ISGSGGS"],
            "A_FR3": ["TYYAESVKGRFTI", "TYYAESVKGRFTI"],
            "A_CDR3": ["CARDYW", "CARGFW"],
            "A_FR4": ["WGQGTLV", "WGQGTLV"],
            "B_FR1": ["DIQMTQS", ""],
            "B_CDR1": ["QSISSY", ""],
            "B_FR2": ["LNWYQQK", ""],
            "B_CDR2": ["AASSLQS", ""],
            "B_FR3": ["GVPSRFSGSG", ""],
            "B_CDR3": ["CQQYNS", ""],
            "B_FR4": ["FGQGTKV", ""],
        }
        reads = pl.DataFrame({"entity_key": ["c1", "c2"], **regions})
        plan = {
            "mode": "antibody_tcr_legacy_sc",
            "receptor": "IG",
            "chains": ["A", "B"],
            "fullChains": ["A", "B"],
            "hasFv": True,
        }
        out = run(reads, plan)
        rows = {r["entity_key"]: r for r in out["properties"].iter_rows(named=True)}

        # c1 paired and complete — every column populated.
        for p in CDR3_PROPS:
            assert rows["c1"][f"{p}_A_CDR3"] is not None
            assert rows["c1"][f"{p}_B_CDR3"] is not None
        for p in FULL_CHAIN_PROPS:
            assert rows["c1"][f"{p}_A_VDJRegion"] is not None
            assert rows["c1"][f"{p}_B_VDJRegion"] is not None
        for p in FV_PROPS:
            assert rows["c1"][f"{p}_Fv"] is not None

        # c2 chain A intact, chain B empty.
        for p in CDR3_PROPS:
            assert rows["c2"][f"{p}_A_CDR3"] is not None
            assert rows["c2"][f"{p}_B_CDR3"] is None
        for p in FULL_CHAIN_PROPS:
            assert rows["c2"][f"{p}_A_VDJRegion"] is not None
            assert rows["c2"][f"{p}_B_VDJRegion"] is None
        # Fv requires both chains — every Fv column NA for c2.
        for p in FV_PROPS:
            assert rows["c2"][f"{p}_Fv"] is None


class TestEmptyInput:
    """Zero-row inputs must produce well-formed outputs, not crash. The workflow
    panics earlier on missing columns, but the Python step must remain robust if
    the upstream filter selects zero entities.
    """

    # Peptide mode with zero rows — properties has the 9 scalar columns but no
    # rows; aa_fraction has the 3 schema columns but no rows; stats are empty.
    def test_peptide_mode_zero_rows(self):
        reads = pl.DataFrame(schema={"entity_key": pl.Utf8, "sequence": pl.Utf8})
        out = run(reads, {"mode": "peptide"})
        assert out["properties"].height == 0
        for col in PEPTIDE_PROPERTY_COLUMNS:
            assert col in out["properties"].columns
        assert out["aa_fraction"].height == 0
        assert set(out["aa_fraction"].columns) == {"entity_key", "aminoAcid", "value"}
        assert out["stats"] == {"medianCdr3Length": {}}

    # Antibody mode with zero rows — properties has entity_key + the columns
    # the plan asks for; medians dict is empty (no CDR3 lengths to compute).
    def test_antibody_mode_zero_rows(self):
        reads = pl.DataFrame(
            schema={
                "entity_key": pl.Utf8,
                "A_CDR3": pl.Utf8,
                "B_CDR3": pl.Utf8,
            }
        )
        plan = {
            "mode": "antibody_tcr_legacy_sc",
            "receptor": "IG",
            "chains": ["A", "B"],
            "fullChains": [],
            "hasFv": False,
        }
        out = run(reads, plan)
        assert out["properties"].height == 0
        cols = set(out["properties"].columns)
        assert "entity_key" in cols
        for p in CDR3_PROPS:
            assert f"{p}_A_CDR3" in cols
            assert f"{p}_B_CDR3" in cols
        assert out["stats"] == {"medianCdr3Length": {}}


class TestR11cStats:
    """`run()` returns a `stats` dict consumed by the workflow info layer for
    R11c VHH detection. Median CDR3 length per chain is the load-bearing field.
    """

    # Peptide mode has no chains — stats present, medianCdr3Length empty.
    def test_peptide_mode_emits_empty_medians(self):
        reads = pl.DataFrame({"entity_key": ["p1"], "sequence": ["ACDEFGHIKL"]})
        out = run(reads, {"mode": "peptide"})
        assert out["stats"] == {"medianCdr3Length": {}}

    # Antibody mode with only chain A CDR3 — only chain A appears in medians.
    # This is the VHH precondition shape the workflow checks.
    def test_chain_a_only_emits_only_chain_a_median(self):
        reads = pl.DataFrame(
            {
                "entity_key": ["c1", "c2", "c3"],
                "A_CDR3": ["CARDYW", "CARGFW", "CARWWWWWWWWWWWWWWWWW"],  # 6, 6, 21
            }
        )
        plan = {
            "mode": "antibody_tcr_legacy_sc",
            "receptor": "IG",
            "chains": ["A"],
            "fullChains": [],
            "hasFv": False,
        }
        out = run(reads, plan)
        medians = out["stats"]["medianCdr3Length"]
        assert "A" in medians
        assert "B" not in medians
        # Sorted lengths: [6, 6, 21] → odd count, middle is index 1 → 6.
        assert medians["A"] == pytest.approx(6.0)

    # Even-count median: average of the two middle values.
    def test_even_count_median_averages_middle_two(self):
        reads = pl.DataFrame(
            {
                "entity_key": ["c1", "c2", "c3", "c4"],
                "A_CDR3": ["CCCC", "CCCCCC", "CCCCCCCC", "CCCCCCCCCC"],  # 4, 6, 8, 10
            }
        )
        plan = {
            "mode": "antibody_tcr_legacy_sc",
            "receptor": "IG",
            "chains": ["A"],
            "fullChains": [],
            "hasFv": False,
        }
        out = run(reads, plan)
        # (6 + 8) / 2 = 7.0
        assert out["stats"]["medianCdr3Length"]["A"] == pytest.approx(7.0)

    # Effective-length convention: ambiguity codes excluded from the length
    # computation, matching all other property functions.
    def test_median_excludes_ambiguity_codes_from_length(self):
        reads = pl.DataFrame(
            {
                "entity_key": ["c1"],
                "A_CDR3": ["CARDYW" + "X" * 5],  # raw 11 chars, effective 6
            }
        )
        plan = {
            "mode": "antibody_tcr_legacy_sc",
            "receptor": "IG",
            "chains": ["A"],
            "fullChains": [],
            "hasFv": False,
        }
        out = run(reads, plan)
        assert out["stats"]["medianCdr3Length"]["A"] == pytest.approx(6.0)

    # Chain listed in plan but its CDR3 column missing from reads — chain absent
    # from medians. Defensive against plan/data mismatch.
    def test_chain_missing_cdr3_column_absent_from_medians(self):
        reads = pl.DataFrame({"entity_key": ["c1"], "A_CDR3": ["CARDYW"]})
        plan = {
            "mode": "antibody_tcr_legacy_sc",
            "receptor": "IG",
            "chains": ["A", "B"],  # B in plan but B_CDR3 not in reads
            "fullChains": [],
            "hasFv": False,
        }
        out = run(reads, plan)
        medians = out["stats"]["medianCdr3Length"]
        assert "A" in medians
        assert "B" not in medians

    # Chain CDR3 column present but every value empty — chain absent from
    # medians (no length data to compute over).
    def test_chain_with_all_empty_cdr3_absent_from_medians(self):
        reads = pl.DataFrame(
            {
                "entity_key": ["c1", "c2"],
                "A_CDR3": ["CARDYW", "CARGFW"],
                "B_CDR3": ["", ""],
            }
        )
        plan = {
            "mode": "antibody_tcr_legacy_sc",
            "receptor": "IG",
            "chains": ["A", "B"],
            "fullChains": [],
            "hasFv": False,
        }
        out = run(reads, plan)
        medians = out["stats"]["medianCdr3Length"]
        assert "A" in medians
        assert "B" not in medians


class TestTcrModeSkipsFv:
    """TCR receptor types must NOT emit Fv columns even with both chains full."""

    def test_tcr_no_fv(self):
        regions = {
            "A_FR1": ["EVQLVES"],
            "A_CDR1": ["GFTFSSY"],
            "A_FR2": ["AMSWVRQ"],
            "A_CDR2": ["ISGSGGS"],
            "A_FR3": ["TYYAESVKGRFTI"],
            "A_CDR3": ["CASSYW"],
            "A_FR4": ["WGQGTLV"],
            "B_FR1": ["DIQMTQS"],
            "B_CDR1": ["QSISSY"],
            "B_FR2": ["LNWYQQK"],
            "B_CDR2": ["AASSLQS"],
            "B_FR3": ["GVPSRFSGSG"],
            "B_CDR3": ["CASSF"],
            "B_FR4": ["FGQGTKV"],
        }
        reads = pl.DataFrame({"entity_key": ["c1"], **regions})
        plan = {
            "mode": "antibody_tcr_legacy_bulk",
            "receptor": "TCRAB",
            "chains": ["A", "B"],
            "fullChains": ["A", "B"],
            "hasFv": False,  # workflow plan sets this to false for TCR
        }
        out = run(reads, plan)
        for p in FV_PROPS:
            assert f"{p}_Fv" not in out["properties"].columns


# ---------------------------------------------------------------------------
# User-facing logging
# ---------------------------------------------------------------------------
#
# The block surfaces processing progress to users via stderr → workflow log
# stream → PlLogView. Every run-mode and major path emits an INFO-level line
# so users see what stage the block is in. Tests pin the canonical milestones;
# wording can evolve, but each path must keep emitting at least one line.


class TestPipelineLogging:
    """run() and its sub-paths emit INFO-level milestones on the pipeline logger."""

    def test_run_logs_peptide_mode_dispatch(self, caplog: pytest.LogCaptureFixture):
        reads = pl.DataFrame({"entity_key": ["p1"], "sequence": ["ACDEFGHIKL"]})
        with caplog.at_level(logging.INFO, logger="pipeline"):
            run(reads, {"mode": "peptide"})
        assert any("peptide" in r.message.lower() for r in caplog.records), caplog.text

    def test_peptide_path_logs_scalar_and_aa_milestones(self, caplog: pytest.LogCaptureFixture):
        reads = pl.DataFrame({"entity_key": ["p1"], "sequence": ["ACDEFGHIKL"]})
        with caplog.at_level(logging.INFO, logger="pipeline"):
            run(reads, {"mode": "peptide"})
        text = caplog.text.lower()
        # One milestone for scalar properties, one for AA fractions.
        assert "scalar" in text or "properties" in text, caplog.text
        assert "aa" in text or "amino" in text or "fraction" in text, caplog.text

    def test_run_logs_antibody_mode_dispatch(self, caplog: pytest.LogCaptureFixture):
        reads = pl.DataFrame(
            {
                "entity_key": ["c1"],
                "A_CDR3": ["CARDYW"],
                "B_CDR3": ["CQQYNS"],
            }
        )
        plan = {
            "mode": "antibody_tcr_legacy_sc",
            "receptor": "IG",
            "chains": ["A", "B"],
            "fullChains": [],
            "hasFv": False,
        }
        with caplog.at_level(logging.INFO, logger="pipeline"):
            run(reads, plan)
        assert any("antibody" in r.message.lower() or "tcr" in r.message.lower() for r in caplog.records), caplog.text

    def test_antibody_path_logs_cdr3_milestone(self, caplog: pytest.LogCaptureFixture):
        reads = pl.DataFrame(
            {
                "entity_key": ["c1"],
                "A_CDR3": ["CARDYW"],
                "B_CDR3": ["CQQYNS"],
            }
        )
        plan = {
            "mode": "antibody_tcr_legacy_sc",
            "receptor": "IG",
            "chains": ["A", "B"],
            "fullChains": [],
            "hasFv": False,
        }
        with caplog.at_level(logging.INFO, logger="pipeline"):
            run(reads, plan)
        assert any("cdr3" in r.message.lower() for r in caplog.records), caplog.text

    def test_antibody_path_logs_full_chain_and_fv_milestones(
        self,
        caplog: pytest.LogCaptureFixture,
        antibody_full_one_clone: pl.DataFrame,
        antibody_full_plan: dict,
    ):
        with caplog.at_level(logging.INFO, logger="pipeline"):
            run(antibody_full_one_clone, antibody_full_plan)
        text = caplog.text.lower()
        assert "full" in text or "vdj" in text or "vh" in text, caplog.text
        assert "fv" in text, caplog.text
