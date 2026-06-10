"""Whole-frame characterization snapshot — the vectorization safety net.

Freezes the CURRENT quantized pipeline output over a broad fixed corpus. The
vectorized engine (later) must reproduce these byte-for-byte after the same
quantization. Regenerate intentionally with:
    uv run python -m tests._corpus_gen --write-golden
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from io_layer import read_input_tsv, read_plan, write_output_tsv
from pipeline import run

DATA = Path(__file__).resolve().parents[1] / "data" / "characterization"
GOLDEN = DATA / "golden"
CASES = sorted(p.name.removesuffix("_input.tsv") for p in DATA.glob("*_input.tsv"))


@pytest.mark.parametrize("case", CASES)
def test_quantized_output_matches_golden(case, tmp_path):
    reads = read_input_tsv(DATA / f"{case}_input.tsv")
    plan = read_plan(DATA / f"{case}_plan.json")
    out = run(reads, plan)

    props = tmp_path / "p.tsv"
    write_output_tsv(out["properties"], props, sort_keys=["entity_key"])
    assert props.read_bytes() == (GOLDEN / f"{case}.properties.tsv").read_bytes()

    if (GOLDEN / f"{case}.aa_fraction.tsv").exists():
        aa = tmp_path / "aa.tsv"
        write_output_tsv(out["aa_fraction"], aa, sort_keys=["entity_key", "aminoAcid"])
        assert aa.read_bytes() == (GOLDEN / f"{case}.aa_fraction.tsv").read_bytes()

    got_stats = json.dumps(out["stats"], sort_keys=True, separators=(",", ":"))
    assert got_stats == (GOLDEN / f"{case}.stats.json").read_text()
