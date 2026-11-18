"""Peak memory of the CLI does not grow with the row count (marked slow).

Runs the CLI in a child process on two full-paired antibody inputs
(all 7 regions on both chains, Fv requested), the second with 4x the rows of
the first, and compares the peak resident set size of the two children. Both
inputs span more than one read batch and more than one compute chunk.

Run explicitly:
    uv run pytest tests/integration/test_peak_memory.py -m slow -v
"""

from __future__ import annotations

import json
import random
import subprocess
import sys
from pathlib import Path

import pytest

_SRC = Path(__file__).resolve().parents[2] / "src"

_REGIONS = ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")
_REGION_LEN = {"FR1": 25, "CDR1": 8, "FR2": 17, "CDR2": 8, "FR3": 38, "CDR3": 14, "FR4": 11}
_STANDARD_AAS = "ACDEFGHIKLMNPQRSTVWY"
_COLUMNS = [f"{chain}_{region}" for chain in "AB" for region in _REGIONS]

# The peak rises over the first chunks and is level from about 300k rows
# (measured 581 MiB at 100k, 803 MiB at 300k, 827 MiB at 2.4M on macOS arm64),
# so the smaller input is above that point.
_SMALL_ROWS = 400_000
_LARGE_ROWS = 4 * _SMALL_ROWS
# Allocator noise between two runs of the same shape. A peak that follows the
# row count would be near 4x.
_MAX_GROWTH = 1.25
_MAX_PEAK_BYTES = 1.5 * 1024**3

# Runs the CLI, then prints the peak resident set size of this process in bytes.
# ru_maxrss is in bytes on macOS and in KiB on Linux.
_CHILD = """
import resource, sys
sys.path.insert(0, sys.argv[1])
from main import main
assert main(sys.argv[2:]) == 0
peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
print(peak if sys.platform == "darwin" else peak * 1024)
"""


def _write_input(path: Path, rows: int) -> None:
    rng = random.Random(20260617)
    distinct = [
        "\t".join("".join(rng.choices(_STANDARD_AAS, k=_REGION_LEN[region])) for _ in "AB" for region in _REGIONS)
        for _ in range(1000)
    ]
    with path.open("w") as file:
        file.write("entity_key\t" + "\t".join(_COLUMNS) + "\n")
        for key in range(rows):
            file.write(f"{key}\t{distinct[key % len(distinct)]}\n")


def _peak_bytes(tmp_path: Path, name: str, rows: int) -> int:
    work = tmp_path / name
    work.mkdir()
    _write_input(work / "input.tsv", rows)
    plan = {
        "mode": "antibody_tcr",
        "receptor": "IG",
        "chains": ["A", "B"],
        "fullChains": ["A", "B"],
        "hasFv": True,
    }
    (work / "plan.json").write_text(json.dumps(plan))
    args = [
        "--input",
        str(work / "input.tsv"),
        "--plan",
        str(work / "plan.json"),
        "--output",
        str(work / "properties.tsv"),
        "--aa-fraction",
        str(work / "aa_fraction.tsv"),
        "--stats",
        str(work / "stats.json"),
    ]
    child = subprocess.run([sys.executable, "-c", _CHILD, str(_SRC), *args], capture_output=True, text=True, check=True)
    return int(child.stdout.strip().splitlines()[-1])


@pytest.mark.slow
def test_peak_memory_does_not_grow_with_row_count(tmp_path: Path):
    small = _peak_bytes(tmp_path, "small", _SMALL_ROWS)
    large = _peak_bytes(tmp_path, "large", _LARGE_ROWS)
    print(f"\n[peak] {_SMALL_ROWS} rows: {small / 1024**2:.0f} MiB, {_LARGE_ROWS} rows: {large / 1024**2:.0f} MiB")
    assert large < _MAX_GROWTH * small, f"peak grew {large / small:.2f}x for 4x the rows"
    assert large < _MAX_PEAK_BYTES
