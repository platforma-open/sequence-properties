"""Cross-process byte-stability regression test.

`test_cli.py` runs byte-comparison checks via in-process `main()` calls. Those
share the Python process's randomized ahash seed across runs, so a future
change that re-introduces hash-based ordering (Polars `group_by` without
`maintain_order=True`, `set()` iteration into output bytes, etc.) could pass
the in-process checks while still breaking production — two block instances
run in separate processes with independent hash seeds.

This file spawns the CLI via `subprocess.run` for each run. Fresh process,
fresh hash seed. Closes the subprocess-isolation gap titeseq-analysis PR #13's
reviewer flagged on its determinism test.

Every output the tool writes — `properties.tsv`, `aa_fraction.tsv`,
`stats.json` — gets its own byte-compare. A determinism fix on one file does
not protect the others; siblings inherit the same upstream non-determinism
sources.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pytest

# main.py is loaded via pythonpath = ["src"] (see software/pyproject.toml).
# That import-path setup applies to in-process imports, not subprocess
# invocations — for those we point Python at the script directly.
_MAIN_PY = Path(__file__).resolve().parents[2] / "src" / "main.py"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join(row.get(c, "") for c in columns))
    path.write_text("\n".join(lines) + "\n")


def _run_cli_subprocess(
    *,
    input_tsv: Path,
    plan_json: Path,
    out_tsv: Path,
    aa_tsv: Path,
    stats_json: Path,
) -> None:
    """Invoke main.py in a fresh subprocess. Independent hash seed per call."""
    result = subprocess.run(
        [
            sys.executable,
            str(_MAIN_PY),
            "--input",
            str(input_tsv),
            "--plan",
            str(plan_json),
            "--output",
            str(out_tsv),
            "--aa-fraction",
            str(aa_tsv),
            "--stats",
            str(stats_json),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, (
        f"main.py failed (rc={result.returncode}); stderr=\n{result.stderr}"
    )


@pytest.fixture()
def _peptide_inputs(tmp_path: Path) -> tuple[Path, Path]:
    """Shared peptide-mode fixture for cross-process byte-stability."""
    in_tsv = tmp_path / "input.tsv"
    plan_json = tmp_path / "plan.json"
    _write_tsv(
        in_tsv,
        [
            {"entity_key": "p1", "sequence": "ACDEFGHIKL"},
            {"entity_key": "p2", "sequence": "MNPQRSTVWY"},
            {"entity_key": "p3", "sequence": "GFTFSSYAMS"},
            {"entity_key": "p4", "sequence": "KKKKHHHHHH"},
            {"entity_key": "p5", "sequence": "DDDDEEEEEE"},
        ],
        ["entity_key", "sequence"],
    )
    plan_json.write_text(json.dumps({"mode": "peptide"}))
    return in_tsv, plan_json


_OUTPUT_FILE_NAMES = ("properties.tsv", "aa_fraction.tsv", "stats.json")


def _run_paths(tmp_path: Path, suffix: str) -> dict[str, Path]:
    return {
        "out_tsv": tmp_path / f"properties{suffix}.tsv",
        "aa_tsv": tmp_path / f"aa_fraction{suffix}.tsv",
        "stats_json": tmp_path / f"stats{suffix}.json",
    }


def _assert_all_three_byte_identical(a: dict[str, Path], b: dict[str, Path]) -> None:
    """Sibling-output rule: every file the binary writes gets its own check."""
    hashes_a = {
        "properties.tsv": _sha256(a["out_tsv"]),
        "aa_fraction.tsv": _sha256(a["aa_tsv"]),
        "stats.json": _sha256(a["stats_json"]),
    }
    hashes_b = {
        "properties.tsv": _sha256(b["out_tsv"]),
        "aa_fraction.tsv": _sha256(b["aa_tsv"]),
        "stats.json": _sha256(b["stats_json"]),
    }
    for name in _OUTPUT_FILE_NAMES:
        assert hashes_a[name] == hashes_b[name], (
            f"{name} diverged across subprocess runs: "
            f"sha256(A)={hashes_a[name]} sha256(B)={hashes_b[name]}"
        )


# Peptide mode — exercises the scalar-properties + AA-fraction + stats paths.
def test_peptide_outputs_byte_stable_across_subprocess_runs(
    tmp_path: Path, _peptide_inputs: tuple[Path, Path]
) -> None:
    in_tsv, plan_json = _peptide_inputs

    a = _run_paths(tmp_path, "_a")
    b = _run_paths(tmp_path, "_b")

    _run_cli_subprocess(input_tsv=in_tsv, plan_json=plan_json, **a)
    _run_cli_subprocess(input_tsv=in_tsv, plan_json=plan_json, **b)

    _assert_all_three_byte_identical(a, b)


# Antibody mode — exercises per-chain CDR3 + full-chain + Fv computation paths
# that peptide mode does not touch. Same sibling-rule coverage.
_ANTIBODY_COLUMNS = (
    ["entity_key"]
    + [f"A_{f}" for f in ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")]
    + [f"B_{f}" for f in ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")]
)
_ANTIBODY_ROWS = [
    {
        "entity_key": "c1",
        "A_FR1": "EVQLVES", "A_CDR1": "GFTFSSY", "A_FR2": "AMSWVRQ",
        "A_CDR2": "ISGSGGS", "A_FR3": "TYYAESVKGRFTI", "A_CDR3": "CARDYW",
        "A_FR4": "WGQGTLV",
        "B_FR1": "DIQMTQS", "B_CDR1": "QSISSY", "B_FR2": "LNWYQQK",
        "B_CDR2": "AASSLQS", "B_FR3": "GVPSRFSGSG", "B_CDR3": "CQQYNS",
        "B_FR4": "FGQGTKV",
    },
    {
        "entity_key": "c2",
        "A_FR1": "EVQLVES", "A_CDR1": "GFTFSSY", "A_FR2": "AMSWVRQ",
        "A_CDR2": "ISGSGGS", "A_FR3": "TYYAESVKGRFTI", "A_CDR3": "CARGFW",
        "A_FR4": "WGQGTLV",
        "B_FR1": "DIQMTQS", "B_CDR1": "QSISSY", "B_FR2": "LNWYQQK",
        "B_CDR2": "AASSLQS", "B_FR3": "GVPSRFSGSG", "B_CDR3": "CQHFSS",
        "B_FR4": "FGQGTKV",
    },
]
_ANTIBODY_PLAN = {
    "mode": "antibody_tcr_legacy_bulk",
    "receptor": "IG",
    "chains": ["A", "B"],
    "fullChains": ["A", "B"],
    "hasFv": True,
}


def test_antibody_outputs_byte_stable_across_subprocess_runs(tmp_path: Path) -> None:
    in_tsv = tmp_path / "input.tsv"
    plan_json = tmp_path / "plan.json"
    _write_tsv(in_tsv, _ANTIBODY_ROWS, _ANTIBODY_COLUMNS)
    plan_json.write_text(json.dumps(_ANTIBODY_PLAN))

    a = _run_paths(tmp_path, "_a")
    b = _run_paths(tmp_path, "_b")

    _run_cli_subprocess(input_tsv=in_tsv, plan_json=plan_json, **a)
    _run_cli_subprocess(input_tsv=in_tsv, plan_json=plan_json, **b)

    _assert_all_three_byte_identical(a, b)
