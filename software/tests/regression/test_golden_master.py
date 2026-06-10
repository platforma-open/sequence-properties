"""Full-output byte snapshot — the refactor safety net.

Each case runs main() into a tmp dir and asserts every output file is
byte-identical (sha256) to the committed golden produced by regen_golden.py
from pre-refactor code. If any later refactor changes a byte, this fails —
which is the entire point.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest
from golden_cases import CASES

from main import main

GOLDEN_DIR = Path(__file__).resolve().parents[1] / "data" / "golden"
_OUTPUTS = ("properties.tsv", "aa_fraction.tsv", "stats.json")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join(row.get(c, "") for c in columns))
    path.write_text("\n".join(lines) + "\n")


@pytest.mark.parametrize("name", sorted(CASES), ids=sorted(CASES))
def test_output_matches_golden(tmp_path: Path, name: str):
    rows, columns, plan = CASES[name]
    in_tsv = tmp_path / "input.tsv"
    plan_json = tmp_path / "plan.json"
    _write_tsv(in_tsv, rows, columns)
    plan_json.write_text(json.dumps(plan))

    out = {n: tmp_path / n for n in _OUTPUTS}
    rc = main(
        [
            "--input",
            str(in_tsv),
            "--plan",
            str(plan_json),
            "--output",
            str(out["properties.tsv"]),
            "--aa-fraction",
            str(out["aa_fraction.tsv"]),
            "--stats",
            str(out["stats.json"]),
        ]
    )
    assert rc == 0

    golden = GOLDEN_DIR / name
    assert golden.is_dir(), f"missing golden dir for {name}; run regen_golden.py"
    for n in _OUTPUTS:
        assert _sha256(out[n]) == _sha256(golden / n), f"{name}/{n} diverged from golden — refactor changed behaviour"


# Non-vacuous guard: the antibody_full golden actually carries the Fv columns,
# and the tcr golden does NOT (R12). Proves the snapshot is meaningful, not empty.
def test_golden_column_presence():
    import polars as pl

    ab = pl.read_csv(GOLDEN_DIR / "antibody_full" / "properties.tsv", separator="\t")
    assert {"charge_Fv", "pi_Fv", "charge_A_VDJRegion"} <= set(ab.columns)
    tcr = pl.read_csv(GOLDEN_DIR / "tcr" / "properties.tsv", separator="\t")
    assert not any(c.endswith("_Fv") for c in tcr.columns)
