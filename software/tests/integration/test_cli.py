"""End-to-end CLI test — invokes main.main() with file args, verifies outputs.

Mirrors the contract that the Tengo workflow uses to call the script.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import polars as pl

from main import main


def _write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join(row.get(c, "") for c in columns))
    path.write_text("\n".join(lines) + "\n")


# Smoke test: peptide mode round-trip via CLI.
def test_cli_peptide_mode(tmp_path: Path):
    in_tsv = tmp_path / "input.tsv"
    plan_json = tmp_path / "plan.json"
    out_tsv = tmp_path / "out.tsv"
    aa_tsv = tmp_path / "aa.tsv"

    _write_tsv(
        in_tsv,
        [{"entity_key": "p1", "peptide_seq": "ACDEFGHIKL"}],
        ["entity_key", "peptide_seq"],
    )
    plan_json.write_text(json.dumps({"mode": "peptide"}))

    rc = main(
        [
            "--input",
            str(in_tsv),
            "--plan",
            str(plan_json),
            "--output",
            str(out_tsv),
            "--aa-fraction",
            str(aa_tsv),
        ]
    )
    assert rc == 0

    out = pl.read_csv(out_tsv, separator="\t")
    assert out.height == 1
    assert "charge_peptide" in out.columns
    assert "aromaticity_peptide" in out.columns

    aa = pl.read_csv(aa_tsv, separator="\t")
    assert aa.height == 20  # 1 entity × 20 standard AAs


# Antibody mode end-to-end — confirms full-chain + Fv emission via CLI.
def test_cli_antibody_full_coverage(tmp_path: Path):
    in_tsv = tmp_path / "input.tsv"
    plan_json = tmp_path / "plan.json"
    out_tsv = tmp_path / "out.tsv"
    aa_tsv = tmp_path / "aa.tsv"

    columns = (
        ["entity_key"]
        + [f"A_{f}" for f in ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")]
        + [f"B_{f}" for f in ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")]
    )
    row = {
        "entity_key": "c1",
        "A_FR1": "EVQLVES",
        "A_CDR1": "GFTFSSY",
        "A_FR2": "AMSWVRQ",
        "A_CDR2": "ISGSGGS",
        "A_FR3": "TYYAESVKGRFTI",
        "A_CDR3": "CARDYW",
        "A_FR4": "WGQGTLV",
        "B_FR1": "DIQMTQS",
        "B_CDR1": "QSISSY",
        "B_FR2": "LNWYQQK",
        "B_CDR2": "AASSLQS",
        "B_FR3": "GVPSRFSGSG",
        "B_CDR3": "CQQYNS",
        "B_FR4": "FGQGTKV",
    }
    _write_tsv(in_tsv, [row], columns)
    plan_json.write_text(
        json.dumps(
            {
                "mode": "antibody_tcr_legacy_bulk",
                "receptor": "IG",
                "chains": ["A", "B"],
                "fullChains": ["A", "B"],
                "hasFv": True,
            }
        )
    )

    rc = main(
        [
            "--input",
            str(in_tsv),
            "--plan",
            str(plan_json),
            "--output",
            str(out_tsv),
            "--aa-fraction",
            str(aa_tsv),
        ]
    )
    assert rc == 0

    out = pl.read_csv(out_tsv, separator="\t")
    for col in ("charge_A_CDR3", "charge_A_VDJRegion", "charge_B_VDJRegion", "charge_Fv", "pi_Fv"):
        assert col in out.columns


# ---------------------------------------------------------------------------
# Byte-stability — properties.tsv / aa_fraction.tsv must hash identical across
# runs of the same input so the resource lands on the dedup path.
# ---------------------------------------------------------------------------


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _run_peptide(tmp_path: Path, suffix: str, rows: list[dict[str, str]]) -> tuple[Path, Path]:
    """Run the CLI in peptide mode, returning (properties_path, aa_path)."""
    in_tsv = tmp_path / f"input{suffix}.tsv"
    plan_json = tmp_path / f"plan{suffix}.json"
    out_tsv = tmp_path / f"out{suffix}.tsv"
    aa_tsv = tmp_path / f"aa{suffix}.tsv"
    _write_tsv(in_tsv, rows, ["entity_key", "peptide_seq"])
    plan_json.write_text(json.dumps({"mode": "peptide"}))
    rc = main(
        [
            "--input", str(in_tsv),
            "--plan", str(plan_json),
            "--output", str(out_tsv),
            "--aa-fraction", str(aa_tsv),
        ]
    )
    assert rc == 0
    return out_tsv, aa_tsv


# Same input, two runs → same bytes. Guards against ULP drift in transcendentals
# (caught by the quantization step in pipeline._quantize_for_cid) and against
# row-order drift on the way out (caught by the sort_keys path in
# io_layer.write_output_tsv).
def test_byte_stable_across_two_runs(tmp_path: Path):
    rows = [
        {"entity_key": "p1", "peptide_seq": "ACDEFGHIKL"},
        {"entity_key": "p2", "peptide_seq": "MNPQRSTVWY"},
        {"entity_key": "p3", "peptide_seq": "GFTFSSYAMS"},
    ]
    out_a, aa_a = _run_peptide(tmp_path, "_a", rows)
    out_b, aa_b = _run_peptide(tmp_path, "_b", rows)
    assert _sha256(out_a) == _sha256(out_b)
    assert _sha256(aa_a) == _sha256(aa_b)


# Same content, different input row order → same output bytes. Proves the
# write-side sort actually normalises ordering rather than passing through.
def test_byte_stable_under_row_permutation(tmp_path: Path):
    rows = [
        {"entity_key": "p1", "peptide_seq": "ACDEFGHIKL"},
        {"entity_key": "p2", "peptide_seq": "MNPQRSTVWY"},
        {"entity_key": "p3", "peptide_seq": "GFTFSSYAMS"},
    ]
    out_sorted, aa_sorted = _run_peptide(tmp_path, "_sorted", rows)
    out_shuffled, aa_shuffled = _run_peptide(tmp_path, "_shuffled", list(reversed(rows)))
    assert _sha256(out_sorted) == _sha256(out_shuffled)
    assert _sha256(aa_sorted) == _sha256(aa_shuffled)


# Negative test: different content → different bytes. Proves we're actually
# distinguishing inputs and not just always cache-hitting on a constant.
def test_byte_changes_when_input_changes(tmp_path: Path):
    rows_a = [{"entity_key": "p1", "peptide_seq": "ACDEFGHIKL"}]
    rows_b = [{"entity_key": "p1", "peptide_seq": "ACDEFGHIKM"}]  # last residue differs
    out_a, _ = _run_peptide(tmp_path, "_a", rows_a)
    out_b, _ = _run_peptide(tmp_path, "_b", rows_b)
    assert _sha256(out_a) != _sha256(out_b)


# ---------------------------------------------------------------------------
# CLI logging — pipeline milestones reach stderr so the workflow log stream
# captures them and the UI surfaces them via PlLogView.
# ---------------------------------------------------------------------------


def test_cli_writes_progress_to_stderr(tmp_path: Path, capsys):
    in_tsv = tmp_path / "input.tsv"
    plan_json = tmp_path / "plan.json"
    out_tsv = tmp_path / "out.tsv"
    aa_tsv = tmp_path / "aa.tsv"
    _write_tsv(
        in_tsv,
        [{"entity_key": "p1", "peptide_seq": "ACDEFGHIKL"}],
        ["entity_key", "peptide_seq"],
    )
    plan_json.write_text(json.dumps({"mode": "peptide"}))

    rc = main(
        [
            "--input", str(in_tsv),
            "--plan", str(plan_json),
            "--output", str(out_tsv),
            "--aa-fraction", str(aa_tsv),
        ]
    )
    assert rc == 0
    captured = capsys.readouterr()
    assert "peptide" in captured.err.lower(), captured.err
    assert "scalar" in captured.err.lower() or "properties" in captured.err.lower(), captured.err
