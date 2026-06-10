"""Whole-frame characterization snapshot — the vectorization safety net.

Freezes the CURRENT quantized pipeline output over a broad fixed corpus. The
vectorized engine (later) must reproduce these byte-for-byte after the same
quantization. Regenerate intentionally with:
    uv run python -m tests._corpus_gen --write-golden
"""

from __future__ import annotations

import difflib
import json
from pathlib import Path

import pytest

from io_layer import read_input_tsv, read_plan, write_output_tsv
from pipeline import run

DATA = Path(__file__).resolve().parents[1] / "data" / "characterization"
GOLDEN = DATA / "golden"
CASES = sorted(p.name.removesuffix("_input.tsv") for p in DATA.glob("*_input.tsv"))


def _assert_file_bytes_equal(got: Path, golden: Path, label: str) -> None:
    """Byte-exact comparison; on mismatch, raise a case/artifact-attributed
    AssertionError carrying a unified line diff instead of an opaque byte blob."""
    got_b, gold_b = got.read_bytes(), golden.read_bytes()
    if got_b == gold_b:
        return
    diff = "\n".join(
        difflib.unified_diff(
            gold_b.decode().splitlines(),
            got_b.decode().splitlines(),
            fromfile=f"golden/{label}",
            tofile=f"got/{label}",
            lineterm="",
        )
    )
    raise AssertionError(f"{label}: byte mismatch vs golden\n{diff}")


@pytest.mark.parametrize("case", CASES)
def test_quantized_output_matches_golden(case, tmp_path):
    reads = read_input_tsv(DATA / f"{case}_input.tsv")
    plan = read_plan(DATA / f"{case}_plan.json")
    out = run(reads, plan)

    props = tmp_path / "p.tsv"
    write_output_tsv(out["properties"], props, sort_keys=["entity_key"])
    _assert_file_bytes_equal(props, GOLDEN / f"{case}.properties.tsv", f"{case}.properties.tsv")

    if plan["mode"] == "peptide":
        aa = tmp_path / "aa.tsv"
        write_output_tsv(out["aa_fraction"], aa, sort_keys=["entity_key", "aminoAcid"])
        aa_golden = GOLDEN / f"{case}.aa_fraction.tsv"
        assert aa_golden.exists(), f"{case}.aa_fraction.tsv: golden missing at {aa_golden}"
        _assert_file_bytes_equal(aa, aa_golden, f"{case}.aa_fraction.tsv")

    got_stats = json.dumps(out["stats"], sort_keys=True, separators=(",", ":"))
    golden_text = (GOLDEN / f"{case}.stats.json").read_text()
    assert got_stats == golden_text, f"{case}.stats.json mismatch"
