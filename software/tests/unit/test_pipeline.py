"""Behavioral tests for src/pipeline.py.

Verifies mode dispatch (peptide / antibody-tcr), expected output column sets,
NA propagation per-clone for missing regions, and Fv emission gating.

Run from blocks/sequence-properties/software/:
    uv sync
    uv run pytest tests/unit/test_pipeline.py
"""

from __future__ import annotations

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
                "peptide_seq": ["ACDEFGHIKL", "AAAAAAAAAA", ""],
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

    @pytest.fixture
    def reads(self) -> pl.DataFrame:
        # Synthetic 7-region clones for chains A and B. All 7 regions present
        # for both chains so reconstruction succeeds.
        regions = {
            "A_FR1": ["EVQLVES", "EVQLVES"],
            "A_CDR1": ["GFTFSSY", "GFTFSSY"],
            "A_FR2": ["AMSWVRQ", "AMSWVRQ"],
            "A_CDR2": ["ISGSGGS", "ISGSGGS"],
            "A_FR3": ["TYYAESVKGRFTI", "TYYAESVKGRFTI"],
            "A_CDR3": ["CARDYW", "CARGFW"],
            "A_FR4": ["WGQGTLV", "WGQGTLV"],
            "B_FR1": ["DIQMTQS", "DIQMTQS"],
            "B_CDR1": ["QSISSY", "QSISSY"],
            "B_FR2": ["LNWYQQK", "LNWYQQK"],
            "B_CDR2": ["AASSLQS", "AASSLQS"],
            "B_FR3": ["GVPSRFSGSG", "GVPSRFSGSG"],
            "B_CDR3": ["CQQYNS", "CQHFSS"],
            "B_FR4": ["FGQGTKV", "FGQGTKV"],
        }
        return pl.DataFrame({"entity_key": ["c1", "c2"], **regions})

    @pytest.fixture
    def plan(self) -> dict:
        return {
            "mode": "antibody_tcr_legacy_bulk",
            "receptor": "IG",
            "chains": ["A", "B"],
            "fullChains": ["A", "B"],
            "hasFv": True,
        }

    def test_emits_full_chain_and_fv_columns(self, reads: pl.DataFrame, plan: dict):
        out = run(reads, plan)
        cols = set(out["properties"].columns)
        for ch in "AB":
            for p in FULL_CHAIN_PROPS:
                assert f"{p}_{ch}_VDJRegion" in cols
        for p in FV_PROPS:
            assert f"{p}_Fv" in cols

    def test_full_chain_pi_in_range(self, reads: pl.DataFrame, plan: dict):
        out = run(reads, plan)
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
