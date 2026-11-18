"""TSV in / TSV out. Defined by the Tengo workflow's contract.

Input TSV (entity per row):
  Peptide mode columns: `entity_key`, `sequence`
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
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import polars as pl


def read_input_tsv(path: str | Path) -> pl.DataFrame:
    """Read the entity TSV. Empty cells come back as None — the pipeline
    treats both empty-string and None as "missing region" via `if not v`.
    """
    return pl.read_csv(path, separator="\t", infer_schema_length=0, has_header=True)


# File bytes parsed per input batch. Bounds the sequence bytes of a compute
# chunk, so its per-residue buffers stay the same size for long sequences.
READ_BLOCK_BYTES = 25_000_000


def read_input_batches(path: str | Path, block_bytes: int = READ_BLOCK_BYTES) -> Iterator[pl.DataFrame]:
    """Yield the entity TSV as consecutive frames in file order, each parsed
    like `read_input_tsv` from about `block_bytes` of whole lines (one line at
    least), so the file is never held whole. A record is one line. A file with
    no data rows yields one empty frame with the header's columns.
    """
    with open(path, "rb") as file:
        header = file.readline()

        def parse(lines: bytes) -> pl.DataFrame:
            return pl.read_csv(header + lines, separator="\t", infer_schema_length=0, has_header=True)

        yielded = False
        rest = b""
        while block := file.read(block_bytes):
            data = rest + block
            end = data.rfind(b"\n") + 1  # 0 when the block has no line end yet
            rest = data[end:]
            if end:
                yield parse(data[:end])
                yielded = True
        if rest or not yielded:
            yield parse(rest)


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


class TsvChunkWriter:
    """Write a TSV from consecutive frames of one schema, in the given order:
    header from the first frame, then each frame's rows. Byte-identical to
    `write_output_tsv` on the concatenated frame, without materializing it.
    """

    def __init__(self, path: str | Path) -> None:
        self._file = open(path, "wb")
        self._header = True

    def write(self, frame: pl.DataFrame) -> None:
        frame.write_csv(self._file, separator="\t", null_value="", include_header=self._header)
        self._header = False

    def __enter__(self) -> TsvChunkWriter:
        return self

    def __exit__(self, *exc: object) -> None:
        self._file.close()
