"""Shared pytest fixtures for sequence-properties tests.

Hosts data shared across unit + integration suites. Keep individual test
files focused on assertions; put repeated synthetic-data construction here.
"""

from __future__ import annotations

import polars as pl
import pytest

# Seven IMGT regions per chain. The CDR3 column is the one R11c keys on,
# so it varies by test — but the framework regions are inert padding for
# anything other than full-chain reconstruction. Kept short (per-region
# 6–13 aa) so byte-stable tests don't churn on whitespace edits.
_BASE_REGIONS_TWO_CLONES: dict[str, list[str]] = {
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


@pytest.fixture
def antibody_full_two_clones() -> pl.DataFrame:
    """Two-clone antibody DataFrame with all 7 regions present on chains A and B.

    Use for tests that exercise full-chain reconstruction + Fv emission.
    """
    return pl.DataFrame({"entity_key": ["c1", "c2"], **_BASE_REGIONS_TWO_CLONES})


@pytest.fixture
def antibody_full_one_clone() -> pl.DataFrame:
    """Single-clone variant — every region defined for chain A and chain B."""
    return pl.DataFrame(
        {
            "entity_key": ["c1"],
            **{k: [v[0]] for k, v in _BASE_REGIONS_TWO_CLONES.items()},
        }
    )


@pytest.fixture
def antibody_full_plan() -> dict:
    """Plan dict matching `antibody_full_*` fixtures: IG, both chains full, Fv on."""
    return {
        "mode": "antibody_tcr_legacy_bulk",
        "receptor": "IG",
        "chains": ["A", "B"],
        "fullChains": ["A", "B"],
        "hasFv": True,
    }
