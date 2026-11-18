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

from io_layer import TsvChunkWriter, read_input_tsv, read_plan
from pipeline import DEFAULT_CHUNK_ROWS, run_streaming

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


# Tiny chunk sizes split every case into many chunks (1 -> single-row chunks,
# including all-NA ones); output must equal the single-chunk golden byte-for-byte.
@pytest.mark.parametrize("chunk_rows", [DEFAULT_CHUNK_ROWS, 7, 1])
@pytest.mark.parametrize("case", CASES)
def test_quantized_output_matches_golden(case, chunk_rows, tmp_path):
    # The golden is in key order; the peptide fixture is not.
    reads = read_input_tsv(DATA / f"{case}_input.tsv").sort("entity_key", maintain_order=True)
    plan = read_plan(DATA / f"{case}_plan.json")
    # Same calls as main.py: chunks are written as they are computed. The
    # golden is written from whole frames with a global sort (_corpus_gen).
    props = tmp_path / "p.tsv"
    aa = tmp_path / "aa.tsv"
    with TsvChunkWriter(props) as props_writer, TsvChunkWriter(aa) as aa_writer:

        def emit(props_chunk, aa_chunk):
            props_writer.write(props_chunk)
            aa_writer.write(aa_chunk)

        stats = run_streaming(reads, plan, emit, chunk_rows=chunk_rows)
    _assert_file_bytes_equal(props, GOLDEN / f"{case}.properties.tsv", f"{case}.properties.tsv")

    if plan["mode"] == "peptide":
        aa_golden = GOLDEN / f"{case}.aa_fraction.tsv"
        assert aa_golden.exists(), f"{case}.aa_fraction.tsv: golden missing at {aa_golden}"
        _assert_file_bytes_equal(aa, aa_golden, f"{case}.aa_fraction.tsv")
    else:
        assert aa.read_text() == "entity_key\taminoAcid\tvalue\n"

    got_stats = json.dumps(stats, sort_keys=True, separators=(",", ":"))
    golden_text = (GOLDEN / f"{case}.stats.json").read_text()
    assert got_stats == golden_text, f"{case}.stats.json mismatch"
