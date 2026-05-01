"""Behavioral tests for src/io_layer.py.

Round-trip TSV in → DataFrame → TSV out preserves entity keys and per-cell
NA semantics (empty string in input == empty string in output).
"""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl

from io_layer import read_input_tsv, read_plan, write_output_tsv


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
