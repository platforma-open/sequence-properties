"""End-to-end CLI test — invokes main.main() with file args, verifies outputs.

Mirrors the contract that the Tengo workflow uses to call the script.
"""

from __future__ import annotations

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
