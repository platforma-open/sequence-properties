"""Golden-master case definitions — the single source of truth for both the
regeneration script and the byte-compare test. Each case is (input rows,
column order, plan dict). Inputs live in code; only the *outputs* are
committed under tests/data/golden/<name>/.

Edge cases are drawn from the spec's documented behaviour tables (NA
propagation, non-standard residues, stop codon, per-clone missing region,
TCR has no Fv).
"""

from __future__ import annotations

_AB_COLS = (
    ["entity_key"]
    + [f"A_{f}" for f in ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")]
    + [f"B_{f}" for f in ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")]
)

_AB_FULL_ROW = {
    "entity_key": "c1",
    "A_FR1": "EVQLVES",
    "A_CDR1": "GFTFSSY",
    "A_FR2": "AMSWVRQ",
    "A_CDR2": "ISGSGGS",
    "A_FR3": "TYYAESVKGRFTI",
    "A_CDR3": "CARDYW",
    "A_FR4": "WGQGTLV",
    "B_FR1": "DIQMTQS",
    "B_CDR1": "QSISSY",
    "B_FR2": "LNWYQQK",
    "B_CDR2": "AASSLQS",
    "B_FR3": "GVPSRFSGSG",
    "B_CDR3": "CQQYNS",
    "B_FR4": "FGQGTKV",
}
# Clone with one heavy region missing -> full-chain A is NA for this clone only.
_AB_MISSING_REGION = {**_AB_FULL_ROW, "entity_key": "c2", "A_FR3": "", "A_CDR3": "CARGFW", "B_CDR3": "CQHFSS"}
# Clone with empty heavy CDR3 -> CDR3-A NA for this clone only.
_AB_EMPTY_CDR3 = {**_AB_FULL_ROW, "entity_key": "c3", "A_CDR3": ""}

# CASES: name -> (rows, columns, plan)
CASES: dict[str, tuple[list[dict], list[str], dict]] = {
    "peptide": (
        [
            {"entity_key": "p_valid", "sequence": "ACDEFGHIKL"},
            {"entity_key": "p_basic", "sequence": "KKKKHHHHHH"},
            {"entity_key": "p_acidic", "sequence": "DDDDEEEEEE"},
            {"entity_key": "p_short", "sequence": "RPPGFSPF"},
            {"entity_key": "p_no_aromatic", "sequence": "AAAAAAAAAA"},
            {"entity_key": "p_paired_cys", "sequence": "CYIQNCPLG"},
            {"entity_key": "p_nonstd", "sequence": "ACDXEFGHIK"},
            {"entity_key": "p_all_nonstd", "sequence": "XXXXX"},
            {"entity_key": "p_stop", "sequence": "ACDE*GHIK"},
            {"entity_key": "p_empty", "sequence": ""},
        ],
        ["entity_key", "sequence"],
        {"mode": "peptide"},
    ),
    "antibody_full": (
        [_AB_FULL_ROW, _AB_MISSING_REGION, _AB_EMPTY_CDR3],
        _AB_COLS,
        {
            "mode": "antibody_tcr_legacy_bulk",
            "receptor": "IG",
            "chains": ["A", "B"],
            "fullChains": ["A", "B"],
            "hasFv": True,
        },
    ),
    "antibody_cdr3_only": (
        [
            {"entity_key": "c1", "A_CDR3": "CARDYW", "B_CDR3": "CQQYNS"},
            {"entity_key": "c2", "A_CDR3": "CARGFW", "B_CDR3": "CQHFSS"},
        ],
        ["entity_key", "A_CDR3", "B_CDR3"],
        {"mode": "antibody_tcr_legacy_sc", "receptor": "IG", "chains": ["A", "B"], "fullChains": [], "hasFv": False},
    ),
    "antibody_partial": (
        [_AB_FULL_ROW],
        _AB_COLS,
        {
            "mode": "antibody_tcr_legacy_bulk",
            "receptor": "IG",
            "chains": ["A", "B"],
            "fullChains": ["A"],
            "hasFv": False,
        },
    ),
    "tcr": (
        [_AB_FULL_ROW],
        _AB_COLS,
        {
            "mode": "antibody_tcr_legacy_bulk",
            "receptor": "TCRAB",
            "chains": ["A", "B"],
            "fullChains": ["A", "B"],
            "hasFv": False,
        },
    ),
}
