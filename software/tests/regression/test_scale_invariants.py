"""Structural invariants at scale — the per-row-map contract the parallel and
streaming refactors must preserve: one output row per input entity, same key
set, AA-fraction 20 rows/peptide, byte-stable across runs.

Synthetic data is built from a SEEDED RNG (random.Random(0)) so the suite is
deterministic. Marked slow: 5000 rows is enough to exercise >1 worker chunk
once Phase B lands, but too slow for every iteration.
"""

from __future__ import annotations

import hashlib
import random
from pathlib import Path

import polars as pl
import pytest

from io_layer import write_output_tsv
from pipeline import run

_AAS = "ACDEFGHIKLMNPQRSTVWY"
_REGIONS = ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")


def _rand_seq(rng: random.Random, lo: int, hi: int) -> str:
    return "".join(rng.choice(_AAS) for _ in range(rng.randint(lo, hi)))


def _peptide_frame(n: int) -> pl.DataFrame:
    rng = random.Random(0)
    return pl.DataFrame(
        {"entity_key": [f"p{i}" for i in range(n)], "sequence": [_rand_seq(rng, 5, 25) for _ in range(n)]},
        schema={"entity_key": pl.Utf8, "sequence": pl.Utf8},
    )


def _antibody_frame(n: int) -> pl.DataFrame:
    rng = random.Random(1)
    cols: dict[str, list[str]] = {"entity_key": [f"c{i}" for i in range(n)]}
    for ch in ("A", "B"):
        for feat in _REGIONS:
            cols[f"{ch}_{feat}"] = [_rand_seq(rng, 6, 14) for _ in range(n)]
    schema = {k: pl.Utf8 for k in cols}
    return pl.DataFrame(cols, schema=schema)


_AB_PLAN = {
    "mode": "antibody_tcr_legacy_bulk",
    "receptor": "IG",
    "chains": ["A", "B"],
    "fullChains": ["A", "B"],
    "hasFv": True,
}

N = 5000


@pytest.mark.slow
def test_peptide_row_and_key_preservation():
    reads = _peptide_frame(N)
    out = run(reads, {"mode": "peptide"})
    props = out["properties"]
    assert props.height == N
    assert set(props["entity_key"].to_list()) == set(reads["entity_key"].to_list())
    assert out["aa_fraction"].height == 20 * N


@pytest.mark.slow
def test_antibody_row_and_key_preservation():
    reads = _antibody_frame(N)
    out = run(reads, _AB_PLAN)
    props = out["properties"]
    assert props.height == N
    assert set(props["entity_key"].to_list()) == set(reads["entity_key"].to_list())
    assert "charge_Fv" in props.columns


@pytest.mark.slow
def test_peptide_byte_stable_two_runs(tmp_path: Path):
    reads = _peptide_frame(N)
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    write_output_tsv(run(reads, {"mode": "peptide"})["properties"], a, sort_keys=["entity_key"])
    write_output_tsv(run(reads, {"mode": "peptide"})["properties"], b, sort_keys=["entity_key"])
    assert hashlib.sha256(a.read_bytes()).hexdigest() == hashlib.sha256(b.read_bytes()).hexdigest()
