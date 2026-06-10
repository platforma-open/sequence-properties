"""Pipeline-level parity: the vectorized antibody/TCR pipeline matches the
scalar BioPython oracle (properties.py) cell-for-cell, within the CID
quantization.

The per-function unit parity tests (test_vectorized_*.py) pin each vectorized
math function to the oracle on random sequences, and the characterization
snapshot freezes the pipeline's own output. Neither pins the pipeline
COMPOSITION to BioPython: which pKa set + include_cys each column uses, the
A=VH / B=VL Fv pairing, full-chain reconstruction, and column naming. This test
runs the committed characterization corpus through both the vectorized pipeline
and an independent oracle restatement and asserts they agree on every emitted
cell — so a wiring regression (wrong pKa set, include_cys flip, chain swap)
fails here even though every column is present and the self-generated golden
would still pass.

Tolerance is one display quantum (10**-decimals) per column: the pipeline value
is the oracle value rounded at the boundary, so they agree to display precision,
while a wiring error shifts values by far more than a quantum.

Run from blocks/sequence-properties/software/:
    uv sync
    uv run pytest tests/integration/test_pipeline_oracle_parity.py
"""

from __future__ import annotations

from pathlib import Path

import pytest

import properties as scalar
from io_layer import read_input_tsv, read_plan
from pipeline import _decimals_for, run
from pka_tables import IPC2_PEPTIDE, IPC2_PROTEIN

DATA = Path(__file__).resolve().parents[1] / "data" / "characterization"

# Cases with full-chain and/or Fv columns. peptide is fully covered by the
# per-function unit parity tests; these exercise the reconstruction + Fv wiring.
ANTIBODY_CASES = sorted(
    p.name.removesuffix("_input.tsv") for p in DATA.glob("*_input.tsv") if p.name != "peptide_input.tsv"
)

REQUIRED_FEATURES = ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")


def _reconstruct(record: dict, chain: str) -> str | None:
    """FR1+CDR1+...+FR4 for one clone; None if any region is empty or missing.

    An independent restatement of the pipeline's reconstruction contract (the
    old scalar `_reconstruct_chain` this PR removed), so the comparison does not
    lean on the pipeline's own polars implementation.
    """
    parts = []
    for feat in REQUIRED_FEATURES:
        v = record.get(f"{chain}_{feat}")
        if not v:  # None or "" -> missing region -> whole chain is NA
            return None
        parts.append(v)
    return "".join(parts)


def _oracle_expected(record: dict, plan: dict) -> dict[str, float | None]:
    """Every property column the plan emits, computed via the scalar oracle with
    the SAME pKa set / include_cys rule the pipeline uses per context. The oracle
    functions return None for invalid/missing input, so reconstructed-None chains
    and empty CDR3 cells flow through to None naturally.
    """
    exp: dict[str, float | None] = {}

    # CDR3 — IPC 2.0 peptide pKa, Cys included (free-thiol rule).
    for ch in plan.get("chains", []):
        cdr3 = record.get(f"{ch}_CDR3")
        exp[f"charge_{ch}_CDR3"] = scalar.charge_at_ph(cdr3, 7.0, IPC2_PEPTIDE, include_cys=True)
        exp[f"chargeShift_{ch}_CDR3"] = scalar.charge_shift(cdr3, IPC2_PEPTIDE, include_cys=True)
        exp[f"gravy_{ch}_CDR3"] = scalar.gravy(cdr3)

    # Full chain — IPC 2.0 protein pKa, Cys excluded (disulfide-bonded rule).
    for ch in plan.get("fullChains", []):
        seq = _reconstruct(record, ch)
        ox, red = scalar.extinction_coefficients(seq)
        exp[f"charge_{ch}_VDJRegion"] = scalar.charge_at_ph(seq, 7.0, IPC2_PROTEIN, include_cys=False)
        exp[f"pi_{ch}_VDJRegion"] = scalar.isoelectric_point(seq, IPC2_PROTEIN, include_cys=False)
        exp[f"gravy_{ch}_VDJRegion"] = scalar.gravy(seq)
        exp[f"mw_{ch}_VDJRegion"] = scalar.molecular_weight(seq)
        exp[f"eox_{ch}_VDJRegion"] = ox
        exp[f"ered_{ch}_VDJRegion"] = red
        exp[f"instability_{ch}_VDJRegion"] = scalar.instability_index(seq)
        exp[f"aliphatic_{ch}_VDJRegion"] = scalar.aliphatic_index(seq)
        exp[f"aromaticity_{ch}_VDJRegion"] = scalar.aromaticity(seq)

    # Fv — additive over VH=A, VL=B with the full-chain rule (Cys excluded).
    if plan.get("hasFv"):
        vh = _reconstruct(record, "A")
        vl = _reconstruct(record, "B")
        ox, red = scalar.fv_extinction_coefficients(vh, vl)
        exp["charge_Fv"] = scalar.fv_charge(vh, vl, 7.0, IPC2_PROTEIN)
        exp["chargeShift_Fv"] = scalar.fv_charge_shift(vh, vl, IPC2_PROTEIN)
        exp["pi_Fv"] = scalar.fv_isoelectric_point(vh, vl, IPC2_PROTEIN)
        exp["eox_Fv"] = ox
        exp["ered_Fv"] = red
        exp["mw_Fv"] = scalar.fv_molecular_weight(vh, vl)

    return exp


@pytest.mark.parametrize("case", ANTIBODY_CASES)
def test_pipeline_matches_oracle(case):
    reads = read_input_tsv(DATA / f"{case}_input.tsv")
    plan = read_plan(DATA / f"{case}_plan.json")
    out = run(reads, plan)
    actual_rows = {r["entity_key"]: r for r in out["properties"].iter_rows(named=True)}

    compared = 0
    for record in reads.iter_rows(named=True):
        actual = actual_rows[record["entity_key"]]
        for col, expected in _oracle_expected(record, plan).items():
            got = actual[col]
            if expected is None:
                assert got is None, f"{case}/{record['entity_key']}/{col}: oracle NA but pipeline {got}"
                continue
            assert got is not None, f"{case}/{record['entity_key']}/{col}: oracle {expected} but pipeline NA"
            # One display quantum below the column's rounding decimals: covers the
            # boundary rounding (<=0.5 quantum) + per-function FP slack; a wiring
            # error shifts values by far more.
            tol = 10.0 ** (-_decimals_for(col))
            assert got == pytest.approx(expected, abs=tol), (
                f"{case}/{record['entity_key']}/{col}: pipeline {got} vs oracle {expected} (tol {tol})"
            )
            compared += 1
    # Guard against the comparison silently doing nothing (column-name drift, a
    # case whose rows are all-NA) reading as a pass.
    assert compared > 0, f"{case}: no non-NA cells compared"
