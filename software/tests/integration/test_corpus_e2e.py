"""End-to-end corpus tests against committed input TSVs + manifest.

Run from blocks/sequence-properties/software/:
    uv sync
    uv run pytest tests/integration/test_corpus_e2e.py -v

Each clonotype/peptide is one parametrised case. When a case fails, the
pytest ID identifies exactly which sequence regressed — instead of a single
giant failure listing every bad row.

The manifest uses three assertion shapes per cell:

- `"na"`         → cell must be None
- `"not_na"`     → cell must be defined (anything other than None)
- bare number    → exact match (use only for closed-form values)
- `{min, max}`   → numeric range — guards drift without locking the value

For an entire case, `expected_all_na: true` shortcuts to "every property
column for this entity is None".
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import polars as pl
import pytest

from io_layer import read_input_tsv, read_plan
from pipeline import run

CORPUS_DIR = Path(__file__).resolve().parent.parent / "data" / "corpus"
MANIFEST = json.loads((CORPUS_DIR / "manifest.json").read_text())


# ---------------------------------------------------------------------------
# Pipeline output fixtures (computed once per session — pure functions).
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def peptide_outputs() -> dict[str, Any]:
    reads = read_input_tsv(CORPUS_DIR / "peptide_input.tsv")
    plan = read_plan(CORPUS_DIR / "peptide_plan.json")
    return run(reads, plan)


@pytest.fixture(scope="module")
def antibody_outputs() -> dict[str, Any]:
    reads = read_input_tsv(CORPUS_DIR / "antibody_input.tsv")
    plan = read_plan(CORPUS_DIR / "antibody_plan.json")
    return run(reads, plan)


@pytest.fixture(scope="module")
def tcr_outputs() -> dict[str, Any]:
    reads = read_input_tsv(CORPUS_DIR / "antibody_input.tsv")
    plan = read_plan(CORPUS_DIR / "tcr_plan.json")
    return run(reads, plan)


# ---------------------------------------------------------------------------
# Assertion helpers
# ---------------------------------------------------------------------------


def _row_for(properties_df: pl.DataFrame, entity_key: str) -> dict[str, Any]:
    sub = properties_df.filter(pl.col("entity_key") == entity_key)
    assert sub.height == 1, f"expected 1 row for {entity_key}, found {sub.height}"
    return sub.row(0, named=True)


def _check_expectation(row: dict[str, Any], column: str, expected: Any) -> None:
    """One cell vs one expectation. Pytest assertion errors carry context."""
    assert column in row, f"column '{column}' missing from output schema"
    actual = row[column]
    if expected == "na":
        assert actual is None, f"expected NA for {column}, got {actual!r}"
    elif expected == "not_na":
        assert actual is not None, f"expected defined value for {column}, got None"
    elif isinstance(expected, (int, float)):
        assert actual is not None, f"expected ~{expected} for {column}, got None"
        assert actual == pytest.approx(expected, abs=1e-3), f"expected ~{expected} for {column}, got {actual}"
    elif isinstance(expected, dict):
        lo, hi = expected.get("min"), expected.get("max")
        assert actual is not None, f"expected value in [{lo}, {hi}] for {column}, got None"
        if lo is not None:
            assert actual >= lo, f"{column}={actual} below min {lo}"
        if hi is not None:
            assert actual <= hi, f"{column}={actual} above max {hi}"
    else:
        raise AssertionError(f"unrecognised expectation form for {column}: {expected!r}")


# ---------------------------------------------------------------------------
# Peptide section
# ---------------------------------------------------------------------------


PEPTIDE_CASES = sorted(MANIFEST["peptide"]["cases"].items())


@pytest.mark.parametrize("entity_key, entry", PEPTIDE_CASES, ids=[k for k, _ in PEPTIDE_CASES])
def test_peptide_corpus_case(peptide_outputs: dict, entity_key: str, entry: dict):
    """One pytest case per peptide entity — surfaces regressions cleanly."""
    properties = peptide_outputs["properties"]
    row = _row_for(properties, entity_key)

    if entry.get("expected_all_na"):
        for col in properties.columns:
            if col == "entity_key":
                continue
            assert row[col] is None, f"expected_all_na: {col} = {row[col]!r}"
        return

    for column, expected in entry.get("expected", {}).items():
        _check_expectation(row, column, expected)


# Peptide AA-fraction frame: every entity has 20 rows; fractions either sum
# to 1.0 (defined sequence) or are all None (invalid).
def test_peptide_aa_fraction_shape(peptide_outputs: dict):
    aa = peptide_outputs["aa_fraction"]
    expected_entities = {k for k, _ in PEPTIDE_CASES}
    actual_entities = set(aa["entity_key"].unique().to_list())
    assert actual_entities == expected_entities
    # 20 rows per entity.
    for k in expected_entities:
        sub = aa.filter(aa["entity_key"] == k)
        assert sub.height == 20


# ---------------------------------------------------------------------------
# Antibody section
# ---------------------------------------------------------------------------


ANTIBODY_CASES = sorted(MANIFEST["antibody"]["cases"].items())


@pytest.mark.parametrize("entity_key, entry", ANTIBODY_CASES, ids=[k for k, _ in ANTIBODY_CASES])
def test_antibody_corpus_case(antibody_outputs: dict, entity_key: str, entry: dict):
    properties = antibody_outputs["properties"]
    row = _row_for(properties, entity_key)
    for column, expected in entry.get("expected", {}).items():
        _check_expectation(row, column, expected)


# ---------------------------------------------------------------------------
# TCR section — Fv columns must be absent from schema (R12).
# ---------------------------------------------------------------------------


TCR_CASES = sorted(MANIFEST["tcr"]["cases"].items())


@pytest.mark.parametrize("entity_key, entry", TCR_CASES, ids=[k for k, _ in TCR_CASES])
def test_tcr_corpus_case(tcr_outputs: dict, entity_key: str, entry: dict):
    properties = tcr_outputs["properties"]
    for absent_col in entry.get("expected_columns_absent", []):
        assert absent_col not in properties.columns, f"R12 violation: TCR plan must not emit {absent_col}"
    row = _row_for(properties, entity_key)
    for column, expected in entry.get("expected", {}).items():
        _check_expectation(row, column, expected)
