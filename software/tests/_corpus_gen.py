"""Deterministic generator for the characterization corpus + golden snapshot.

This is the keystone of the vectorization safety net. It builds a broad,
fixed corpus across every pipeline mode and edge case, then snapshots the
CURRENT quantized pipeline output to committed golden artifacts. The golden is
"whatever the current code produces" — NOT an independent correctness oracle.
The eventual vectorized engine must reproduce these exact golden bytes.

Re-runnable and idempotent: a fixed `random.Random(0)` seed means re-running
produces byte-identical inputs, plans, and golden. If a regen changes any
golden file, that is a real nondeterminism signal — investigate, do not commit.

Regenerate intentionally (writes inputs + plans + golden) with:

    uv run python -m tests._corpus_gen --write-golden

The snapshot test (`tests/integration/test_characterization_snapshot.py`) then
asserts `run()` reproduces the committed golden byte-for-byte.
"""

from __future__ import annotations

import json
import random
import sys
from pathlib import Path

# `pythonpath = ["src"]` in pyproject only applies under pytest. When run as a
# module (`python -m tests._corpus_gen`), wire up the same import root so
# `from pipeline import run` resolves exactly as the tests see it.
_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import polars as pl

from io_layer import read_input_tsv, read_plan, write_output_tsv
from pipeline import run

DATA = Path(__file__).resolve().parent / "data" / "characterization"
GOLDEN = DATA / "golden"

# Number of generated valid entities per case (edge-case rows are additional).
N_VALID = 50

# IMGT regions per chain, in reconstruction order.
REGIONS = ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")

# Realistic per-region length bands (short-ish, per corpus convention). The
# point is path coverage + a frozen snapshot, not biological realism.
_REGION_LEN = {
    "FR1": (7, 10),
    "CDR1": (6, 8),
    "FR2": (7, 9),
    "CDR2": (6, 8),
    "FR3": (10, 14),
    "CDR3": (5, 18),
    "FR4": (7, 9),
}

# Standard residues only — keep generated sequences valid so they exercise the
# compute path (edge cases below deliberately break this).
_STANDARD_AAS = "ACDEFGHIKLMNPQRSTVWY"
_NO_AROMATIC_AAS = "ACDEGHIKLMNPQRSTV"  # F, W, Y removed.

CASES: list[str] = [
    "peptide",
    "antibody_bulk_full",
    "antibody_cdr3_only",
    "antibody_partial",
    "sc_dropout",
    "tcr_ab",
    "tcr_gd",
]


def _rand_seq(rng: random.Random, lo: int, hi: int, alphabet: str = _STANDARD_AAS) -> str:
    n = rng.randint(lo, hi)
    return "".join(rng.choice(alphabet) for _ in range(n))


# ---------------------------------------------------------------------------
# Peptide case
# ---------------------------------------------------------------------------


def _peptide_rows(rng: random.Random) -> list[dict[str, str]]:
    """~50 generated valid peptides PLUS the required edge-case rows.

    Edge cases (per task): empty cell, stop-codon, non-standard-only, single
    residue, homopolymer, below-instability-floor (<10 aa), no-aromatic.
    """
    rows: list[dict[str, str]] = []
    for i in range(N_VALID):
        rows.append({"entity_key": f"pep_{i:03d}", "sequence": _rand_seq(rng, 8, 30)})

    edge = [
        ("pep_empty", ""),  # invalid -> NA row
        ("pep_stop_codon", "ACDE*FGHI"),  # stop codon -> NA row
        ("pep_nonstandard_only", "BZXJ"),  # cleans to empty -> NA row
        ("pep_single_residue", "A"),
        ("pep_homopolymer", "AAAAAAAAAA"),
        ("pep_below_floor", "ACDEFG"),  # 6 aa < instability floor (10)
        ("pep_no_aromatic", _NO_AROMATIC_AAS),  # no F/W/Y
    ]
    for key, seq in edge:
        rows.append({"entity_key": key, "sequence": seq})
    return rows


def _write_peptide_input(rng: random.Random) -> None:
    rows = _peptide_rows(rng)
    df = pl.DataFrame(
        {
            "entity_key": [r["entity_key"] for r in rows],
            "sequence": [r["sequence"] for r in rows],
        },
        schema={"entity_key": pl.Utf8, "sequence": pl.Utf8},
    )
    write_output_tsv(df, DATA / "peptide_input.tsv")
    _write_plan("peptide", {"mode": "peptide"})


# ---------------------------------------------------------------------------
# Antibody / TCR cases
# ---------------------------------------------------------------------------


def _full_clone(rng: random.Random, key: str, chains: tuple[str, ...]) -> dict[str, str]:
    """A clone with every region present on the given chains."""
    rec: dict[str, str] = {"entity_key": key}
    for ch in chains:
        for feat in REGIONS:
            lo, hi = _REGION_LEN[feat]
            rec[f"{ch}_{feat}"] = _rand_seq(rng, lo, hi)
    return rec


def _empty_regions(chains: tuple[str, ...]) -> dict[str, str]:
    return {f"{ch}_{feat}": "" for ch in chains for feat in REGIONS}


def _antibody_columns(chains: tuple[str, ...]) -> list[str]:
    return ["entity_key"] + [f"{ch}_{feat}" for ch in chains for feat in REGIONS]


def _rows_to_df(rows: list[dict[str, str]], chains: tuple[str, ...]) -> pl.DataFrame:
    cols = _antibody_columns(chains)
    data = {c: [r.get(c, "") for r in rows] for c in cols}
    return pl.DataFrame(data, schema={c: pl.Utf8 for c in cols})


def _write_plan(case: str, plan: dict) -> None:
    (DATA / f"{case}_plan.json").write_text(json.dumps(plan, indent=2) + "\n")


def _write_antibody_bulk_full(rng: random.Random) -> None:
    """All 7 regions on both chains; IG; full chains + Fv."""
    chains = ("A", "B")
    rows = [_full_clone(rng, f"clone_{i:03d}", chains) for i in range(N_VALID)]
    _rows_to_df(rows, chains).pipe(write_output_tsv, DATA / "antibody_bulk_full_input.tsv")
    _write_plan(
        "antibody_bulk_full",
        {
            "mode": "antibody_tcr_legacy_bulk",
            "receptor": "IG",
            "chains": ["A", "B"],
            "fullChains": ["A", "B"],
            "hasFv": True,
        },
    )


def _write_antibody_cdr3_only(rng: random.Random) -> None:
    """Only CDR3 present per chain; other regions empty. No full chains, no Fv."""
    chains = ("A", "B")
    rows: list[dict[str, str]] = []
    for i in range(N_VALID):
        rec = {"entity_key": f"clone_{i:03d}", **_empty_regions(chains)}
        for ch in chains:
            lo, hi = _REGION_LEN["CDR3"]
            rec[f"{ch}_CDR3"] = _rand_seq(rng, lo, hi)
        rows.append(rec)
    _rows_to_df(rows, chains).pipe(write_output_tsv, DATA / "antibody_cdr3_only_input.tsv")
    _write_plan(
        "antibody_cdr3_only",
        {
            "mode": "antibody_tcr_legacy_bulk",
            "receptor": "IG",
            "chains": ["A", "B"],
            "fullChains": [],
            "hasFv": False,
        },
    )


def _write_antibody_partial(rng: random.Random) -> None:
    """Full coverage on most clones, but a deterministic subset is missing one
    region on a chain (empty cell) -> full-chain is NA for those clones, and Fv
    is NA when the missing region is on A or B. fullChains + Fv still requested.
    """
    chains = ("A", "B")
    rows: list[dict[str, str]] = []
    for i in range(N_VALID):
        rec = _full_clone(rng, f"clone_{i:03d}", chains)
        # Every 3rd clone drops one region; rotate which chain/region is dropped
        # so both A-full-NA and B-full-NA paths are exercised.
        if i % 3 == 0:
            ch = "A" if (i // 3) % 2 == 0 else "B"
            feat = REGIONS[(i // 3) % len(REGIONS)]
            rec[f"{ch}_{feat}"] = ""
        rows.append(rec)
    _rows_to_df(rows, chains).pipe(write_output_tsv, DATA / "antibody_partial_input.tsv")
    _write_plan(
        "antibody_partial",
        {
            "mode": "antibody_tcr_legacy_bulk",
            "receptor": "IG",
            "chains": ["A", "B"],
            "fullChains": ["A", "B"],
            "hasFv": True,
        },
    )


def _write_sc_dropout(rng: random.Random) -> None:
    """Single-cell mode; a deterministic subset of clones is missing chain B
    entirely (all B regions empty) -> chain-B CDR3, full-chain, and Fv all NA.
    """
    chains = ("A", "B")
    rows: list[dict[str, str]] = []
    for i in range(N_VALID):
        rec = _full_clone(rng, f"clone_{i:03d}", chains)
        if i % 4 == 0:  # ~1/4 of clones drop chain B entirely.
            for feat in REGIONS:
                rec[f"B_{feat}"] = ""
        rows.append(rec)
    _rows_to_df(rows, chains).pipe(write_output_tsv, DATA / "sc_dropout_input.tsv")
    _write_plan(
        "sc_dropout",
        {
            "mode": "antibody_tcr_legacy_sc",
            "receptor": "IG",
            "chains": ["A", "B"],
            "fullChains": ["A", "B"],
            "hasFv": True,
        },
    )


def _write_tcr_ab(rng: random.Random) -> None:
    """TCRAB receptor; no Fv (TCR has no Fv). Full chains on both chains."""
    chains = ("A", "B")
    rows = [_full_clone(rng, f"clone_{i:03d}", chains) for i in range(N_VALID)]
    _rows_to_df(rows, chains).pipe(write_output_tsv, DATA / "tcr_ab_input.tsv")
    _write_plan(
        "tcr_ab",
        {
            "mode": "antibody_tcr_legacy_bulk",
            "receptor": "TCRAB",
            "chains": ["A", "B"],
            "fullChains": ["A", "B"],
            "hasFv": False,
        },
    )


def _write_tcr_gd(rng: random.Random) -> None:
    """TCRGD receptor; no Fv. Full chains on both chains."""
    chains = ("A", "B")
    rows = [_full_clone(rng, f"clone_{i:03d}", chains) for i in range(N_VALID)]
    _rows_to_df(rows, chains).pipe(write_output_tsv, DATA / "tcr_gd_input.tsv")
    _write_plan(
        "tcr_gd",
        {
            "mode": "antibody_tcr_legacy_bulk",
            "receptor": "TCRGD",
            "chains": ["A", "B"],
            "fullChains": ["A", "B"],
            "hasFv": False,
        },
    )


# ---------------------------------------------------------------------------
# Corpus + golden orchestration
# ---------------------------------------------------------------------------

# Each writer consumes the SAME rng in CASES order, so the byte output is fixed
# by the seed and the call order alone.
_WRITERS = {
    "peptide": _write_peptide_input,
    "antibody_bulk_full": _write_antibody_bulk_full,
    "antibody_cdr3_only": _write_antibody_cdr3_only,
    "antibody_partial": _write_antibody_partial,
    "sc_dropout": _write_sc_dropout,
    "tcr_ab": _write_tcr_ab,
    "tcr_gd": _write_tcr_gd,
}


def write_inputs_and_plans() -> None:
    """(Re)write every case's input TSV + plan JSON deterministically."""
    DATA.mkdir(parents=True, exist_ok=True)
    rng = random.Random(0)
    for case in CASES:
        _WRITERS[case](rng)


def regenerate_golden() -> None:
    """Run the CURRENT pipeline per case and freeze its quantized output to
    golden artifacts, using the exact same calls `main.py` makes on the wire.
    """
    GOLDEN.mkdir(parents=True, exist_ok=True)
    for case in CASES:
        reads = read_input_tsv(DATA / f"{case}_input.tsv")
        plan = read_plan(DATA / f"{case}_plan.json")
        out = run(reads, plan)
        write_output_tsv(out["properties"], GOLDEN / f"{case}.properties.tsv", sort_keys=["entity_key"])
        if plan["mode"] == "peptide":
            write_output_tsv(
                out["aa_fraction"],
                GOLDEN / f"{case}.aa_fraction.tsv",
                sort_keys=["entity_key", "aminoAcid"],
            )
        (GOLDEN / f"{case}.stats.json").write_text(json.dumps(out["stats"], sort_keys=True, separators=(",", ":")))


def main() -> None:
    if "--write-golden" not in sys.argv[1:]:
        raise SystemExit("usage: python -m tests._corpus_gen --write-golden")
    write_inputs_and_plans()
    regenerate_golden()
    print(f"Wrote inputs + plans + golden for {len(CASES)} cases under {DATA}")


if __name__ == "__main__":
    main()
