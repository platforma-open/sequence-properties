"""Regenerate committed golden outputs from the CURRENT code.

Run ONLY deliberately:
    uv run python tests/regression/regen_golden.py

DO NOT run this during the parallelism/streaming refactor — if the golden
needs regenerating to make test_golden_master pass, the refactor changed
observable behaviour and that is a bug, not a snapshot update. Regenerate
only on an intentional behaviour change (BioPython bump, pKa update, spec
change) and call it out in the PR.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from golden_cases import CASES  # noqa: E402

from main import main  # noqa: E402

GOLDEN_DIR = Path(__file__).resolve().parents[1] / "data" / "golden"


def _write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join(row.get(c, "") for c in columns))
    path.write_text("\n".join(lines) + "\n")


def regen() -> None:
    for name, (rows, columns, plan) in CASES.items():
        case_dir = GOLDEN_DIR / name
        case_dir.mkdir(parents=True, exist_ok=True)
        in_tsv = case_dir / "input.tsv"
        plan_json = case_dir / "plan.json"
        _write_tsv(in_tsv, rows, columns)
        plan_json.write_text(json.dumps(plan))
        rc = main(
            [
                "--input",
                str(in_tsv),
                "--plan",
                str(plan_json),
                "--output",
                str(case_dir / "properties.tsv"),
                "--aa-fraction",
                str(case_dir / "aa_fraction.tsv"),
                "--stats",
                str(case_dir / "stats.json"),
            ]
        )
        assert rc == 0, f"regen failed for case {name}"
        print(f"regenerated golden/{name}")


if __name__ == "__main__":
    regen()
