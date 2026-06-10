"""Worker-count invariance — the core safety property for parallelism.

The output of run() must not depend on how many workers compute it. We
compare workers=1 against workers=4 (and a larger frame that actually spills
into the pool) and assert byte-identical serialized output + identical stats.
"""

from __future__ import annotations

import random
from pathlib import Path

import polars as pl
import pytest

from io_layer import write_output_tsv
from pipeline import run

_AAS = "ACDEFGHIKLMNPQRSTVWY"
_REGIONS = ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")
_AB_PLAN = {
    "mode": "antibody_tcr_legacy_bulk",
    "receptor": "IG",
    "chains": ["A", "B"],
    "fullChains": ["A", "B"],
    "hasFv": True,
}


def _peptide_frame(n: int, seed: int) -> pl.DataFrame:
    rng = random.Random(seed)
    seqs = ["".join(rng.choice(_AAS) for _ in range(rng.randint(5, 25))) for _ in range(n)]
    return pl.DataFrame(
        {"entity_key": [f"p{i}" for i in range(n)], "sequence": seqs},
        schema={"entity_key": pl.Utf8, "sequence": pl.Utf8},
    )


def _antibody_frame(n: int, seed: int) -> pl.DataFrame:
    rng = random.Random(seed)
    cols = {"entity_key": [f"c{i}" for i in range(n)]}
    for ch in ("A", "B"):
        for feat in _REGIONS:
            cols[f"{ch}_{feat}"] = ["".join(rng.choice(_AAS) for _ in range(rng.randint(6, 14))) for _ in range(n)]
    return pl.DataFrame(cols, schema={k: pl.Utf8 for k in cols})


def _serialize(out: dict, tmp: Path, tag: str) -> tuple[bytes, bytes, dict]:
    p = tmp / f"props_{tag}.tsv"
    a = tmp / f"aa_{tag}.tsv"
    write_output_tsv(out["properties"], p, sort_keys=["entity_key"])
    write_output_tsv(out["aa_fraction"], a, sort_keys=["entity_key", "aminoAcid"])
    return p.read_bytes(), a.read_bytes(), out["stats"]


@pytest.mark.parametrize("n", [3, 3000], ids=["sequential-path", "pool-path"])
def test_peptide_invariant_to_workers(tmp_path: Path, n: int):
    reads = _peptide_frame(n, seed=0)
    one = _serialize(run(reads, {"mode": "peptide"}, workers=1), tmp_path, "1")
    four = _serialize(run(reads, {"mode": "peptide"}, workers=4), tmp_path, "4")
    assert one == four


@pytest.mark.parametrize("n", [3, 3000], ids=["sequential-path", "pool-path"])
def test_antibody_invariant_to_workers(tmp_path: Path, n: int):
    reads = _antibody_frame(n, seed=1)
    one = _serialize(run(reads, _AB_PLAN, workers=1), tmp_path, "1")
    four = _serialize(run(reads, _AB_PLAN, workers=4), tmp_path, "4")
    assert one == four
