"""Behavioral tests for src/io_layer.py.

Round-trip TSV in → DataFrame → TSV out preserves entity keys and per-cell
NA semantics (empty string in input == empty string in output).
"""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl
import pytest

from io_layer import read_input_batches, read_input_tsv, read_plan, write_output_tsv


# Empty cell in input must survive parse + write as empty cell — not "null"
# or "NaN" — the Tengo workflow's xsv.importFile treats "" as missing.
# Polars reads empty as None and writes None as "" (with null_value="").
def test_round_trip_preserves_empty_cells(tmp_path: Path):
    raw = "entity_key\tsequence\nA\tACDE\nB\t\n"
    src = tmp_path / "in.tsv"
    src.write_text(raw)
    df = read_input_tsv(src)
    assert df["sequence"].to_list() == ["ACDE", None]
    out = tmp_path / "out.tsv"
    write_output_tsv(df, out)
    written = out.read_text()
    # Output row for entity B has its sequence cell empty.
    assert "B\t\n" in written or written.endswith("B\t\n")


def test_read_plan_returns_dict(tmp_path: Path):
    plan = {"mode": "peptide", "chains": ["A"]}
    p = tmp_path / "plan.json"
    p.write_text(json.dumps(plan))
    assert read_plan(p) == plan


def test_write_output_with_none_values(tmp_path: Path):
    df = pl.DataFrame({"k": ["a", "b"], "v": [1.0, None]}, schema={"k": pl.Utf8, "v": pl.Float64})
    out = tmp_path / "x.tsv"
    write_output_tsv(df, out)
    text = out.read_text()
    # Polars writes None as empty cell when null_value="".
    assert "b\t\n" in text


# ---------------------------------------------------------------------------
# read_input_batches — the file as consecutive frames of whole lines
# ---------------------------------------------------------------------------

_BATCH_TSV = "entity_key\tsequence\np1\tACDEF\np2\t\np3\tACDEFGHIKLMNPQRSTVWY\np4\tAC\np5\tACD\n"


@pytest.mark.parametrize("block_bytes", [1, 5, 9, 16, 30, 1000])
@pytest.mark.parametrize("trailing_newline", [True, False])
def test_batches_concatenate_to_the_whole_table(tmp_path: Path, block_bytes: int, trailing_newline: bool):
    src = tmp_path / "in.tsv"
    src.write_text(_BATCH_TSV if trailing_newline else _BATCH_TSV.rstrip("\n"))
    batches = list(read_input_batches(src, block_bytes))
    assert all(b.height >= 1 for b in batches)
    assert pl.concat(batches).equals(read_input_tsv(src), null_equal=True)


def test_batches_hold_whole_lines_within_the_block(tmp_path: Path):
    src = tmp_path / "in.tsv"
    src.write_text(_BATCH_TSV)
    keys = [b["entity_key"].to_list() for b in read_input_batches(src, 16)]
    # A batch ends at the last line end inside its 16-byte block; p3 (24 bytes)
    # spans two blocks and comes out when its line end is read.
    assert keys == [["p1", "p2"], ["p3", "p4"], ["p5"]]


@pytest.mark.parametrize("text", ["entity_key\tsequence\n", "entity_key\tsequence"])
def test_header_only_file_yields_one_empty_frame(tmp_path: Path, text: str):
    src = tmp_path / "in.tsv"
    src.write_text(text)
    (frame,) = read_input_batches(src)
    assert frame.height == 0
    assert frame.schema == {"entity_key": pl.Utf8, "sequence": pl.Utf8}
