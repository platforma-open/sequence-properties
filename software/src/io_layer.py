"""TSV in / TSV out. Defined by the Tengo workflow's contract.

Input TSV (entity per row):
  Peptide mode columns: `entity_key`, `peptide_seq`
  Antibody/TCR mode columns: `entity_key`, `<chain>_<feature>` for each
    (chain ∈ {A, B}) × (feature ∈ {FR1, CDR1, FR2, CDR2, FR3, CDR3, FR4})
    that the upstream block actually emitted. Missing region for a clone is
    represented by an empty cell — never the literal "NA".

Output properties TSV (entity per row): one column per requested property,
named per `process.tpl.tengo`'s expectations. Empty cell == NA.

Output AA-fraction TSV (long format): three columns — `entity_key`,
`aminoAcid`, `value`. One row per (entity, aminoAcid) pair where value is
the mole fraction. Empty body in non-peptide mode.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import polars as pl


def read_input_tsv(path: str | Path) -> pl.DataFrame:
    """Read the entity TSV. Empty cells come back as None — the pipeline
    treats both empty-string and None as "missing region" via `if not v`.
    """
    return pl.read_csv(path, separator="\t", infer_schema_length=0, has_header=True)


def read_plan(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text())


def write_output_tsv(
    df: pl.DataFrame,
    path: str | Path,
    sort_keys: list[str] | None = None,
) -> None:
    """Write a TSV. Polars writes None as empty by default with `null_value=""`.

    `sort_keys`, when provided, fixes row order — required for byte-stable
    output so the resulting resource has a stable CID across runs (dedup).
    """
    if sort_keys:
        df = df.sort(sort_keys)
    df.write_csv(path, separator="\t", null_value="")
