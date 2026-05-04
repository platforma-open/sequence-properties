"""Characterization tests — pin the BioPython-backed numeric output for
representative sequences so refactoring drift is caught.

These are NOT correctness checks against an independent reference (those live
in test_properties.py with closed-form expectations like instability=9.0 for
poly-A-10). They snapshot the *current* BioPython-backed output to pin
behaviour during refactoring. The IPC 2.0 pKa values are spec-blessed; the
mass / hydropathy / extinction / instability tables come from BioPython
ProtParam (also spec-blessed per M1 strategy direction).

The reconstructed VH / VL pair below matches the `ab_full_paired` row in
the e2e corpus, so every refactor exercises the same sequence path the
corpus does — but at numeric precision the corpus's range bounds don't
catch.

If a value changes intentionally (BioPython version update, pKa update,
spec change), update the golden number here as part of the change.
"""

from __future__ import annotations

import pytest

from pka_tables import IPC2_PEPTIDE, IPC2_PROTEIN
from properties import (
    aliphatic_index,
    aromaticity,
    charge_at_ph,
    extinction_coefficients,
    fv_charge,
    fv_extinction_coefficients,
    fv_isoelectric_point,
    fv_molecular_weight,
    gravy,
    instability_index,
    isoelectric_point,
    molecular_weight,
)

# Reconstructed antibody chains — concatenation order from the spec:
# FR1 + CDR1 + FR2 + CDR2 + FR3 + CDR3 + FR4. Same chains used in the
# e2e corpus's ab_full_paired entry.
VH = "EVQLVESGFTFSSYAMSWVRQISGSGGSTYYAESVKGRFTICARDYWWGQGTLV"
VL = "DIQMTQSQSISSYLNWYQQKAASSLQSGVPSRFSGSGCQQYNSFGQGTKV"


# Tolerance: 1e-6 catches drift smaller than the spec's reporting precision
# (.2f / .3f formats) without being so tight that numpy/scipy version
# changes break the suite. Update the golden value if the change is
# intentional.
ABS_TOL = 1e-6


class TestVHGoldenValues:
    """VH (chain A, full reconstruction). Protein pKa set, Cys excluded."""

    @pytest.mark.parametrize(
        "fn_name, expected",
        [
            ("charge", -0.836811),
            ("pi", 6.006653),
            ("gravy", -0.111111),
            ("mw", 6050.6567),
            ("eox", 22460.0),
            ("ered", 22460.0),
            ("instability", 38.753704),
            ("aliphatic", 61.296296),
            ("aromaticity", 0.185185),
        ],
    )
    def test_vh_property(self, fn_name: str, expected: float):
        if fn_name == "charge":
            actual = charge_at_ph(VH, 7.0, IPC2_PROTEIN, include_cys=False)
        elif fn_name == "pi":
            actual = isoelectric_point(VH, IPC2_PROTEIN, include_cys=False)
        elif fn_name == "gravy":
            actual = gravy(VH)
        elif fn_name == "mw":
            actual = molecular_weight(VH)
        elif fn_name == "eox":
            actual = extinction_coefficients(VH)[0]
        elif fn_name == "ered":
            actual = extinction_coefficients(VH)[1]
        elif fn_name == "instability":
            actual = instability_index(VH)
        elif fn_name == "aliphatic":
            actual = aliphatic_index(VH)
        elif fn_name == "aromaticity":
            actual = aromaticity(VH)
        else:
            raise AssertionError(f"unknown fn: {fn_name}")
        assert actual == pytest.approx(expected, abs=ABS_TOL)


class TestVLGoldenValues:
    """VL (chain B). Same pKa context as VH."""

    def test_vl_charge(self):
        assert charge_at_ph(VL, 7.0, IPC2_PROTEIN, include_cys=False) == pytest.approx(1.149394, abs=ABS_TOL)

    def test_vl_pi(self):
        assert isoelectric_point(VL, IPC2_PROTEIN, include_cys=False) == pytest.approx(9.16571, abs=ABS_TOL)


class TestFvGoldenValues:
    """Fv (paired VH + VL). pI bisects the per-chain charge sum."""

    def test_fv_charge(self):
        # Additive: must equal charge(VH) + charge(VL).
        assert fv_charge(VH, VL, 7.0, IPC2_PROTEIN) == pytest.approx(0.312583, abs=ABS_TOL)

    def test_fv_pi(self):
        # Per-chain charge sum bisection — distinct from concatenated-chain pI.
        assert fv_isoelectric_point(VH, VL, IPC2_PROTEIN) == pytest.approx(7.633606, abs=ABS_TOL)

    def test_fv_extinction_oxidised(self):
        ox, _ = fv_extinction_coefficients(VH, VL)
        assert ox == pytest.approx(32430.0, abs=ABS_TOL)

    def test_fv_extinction_reduced(self):
        _, red = fv_extinction_coefficients(VH, VL)
        assert red == pytest.approx(32430.0, abs=ABS_TOL)

    def test_fv_mw(self):
        assert fv_molecular_weight(VH, VL) == pytest.approx(11517.5517, abs=ABS_TOL)


# Three peptides chosen for varied composition — covers the spec M1 acceptance
# criterion explicitly (charge / GRAVY / MW / ε pinned on ≥3 peptide sequences).
# Peptide pKa set, Cys included as ionizable per spec L535.
#
# bradykinin   — RPPGFSPFR, 9 aa, basic, no aromatic, no Cys, <10 aa (II = NA)
# scan12       — ACDEFGHIKLMN, 12 aa mixed-property scan
# oxytocin     — CYIQNCPLG, 9 aa, paired Cys + 1 Tyr (ε exercises floor(C/2)·125)
PEPTIDE_GOLDENS = [
    {
        "name": "bradykinin",
        "seq": "RPPGFSPFR",
        "charge": 1.898520,
        "pi": 11.493347,
        "gravy": -1.044444,
        "mw": 1060.2085,
        "eox": 0.0,
        "ered": 0.0,
        "instability": None,  # length < 10 → NA
        "aliphatic": 0.0,
        "aromaticity": 0.222222,
    },
    {
        "name": "scan12",
        "seq": "ACDEFGHIKLMN",
        "charge": -0.949365,
        "pi": 5.530701,
        "gravy": -0.058333,
        "mw": 1377.5880,
        "eox": 0.0,
        "ered": 0.0,
        "instability": 84.300000,
        "aliphatic": 73.333333,
        "aromaticity": 0.083333,
    },
    {
        "name": "oxytocin",
        "seq": "CYIQNCPLG",
        "charge": -0.115650,
        "pi": 5.435852,
        "gravy": 0.333333,
        "mw": 1010.1878,
        "eox": 1615.0,  # 1490 (Y) + 125 (one disulfide from 2 Cys)
        "ered": 1490.0,
        "instability": None,
        "aliphatic": 86.666667,
        "aromaticity": 0.111111,
    },
]


class TestPeptideGoldenValues:
    """Pinned numeric output for representative peptides — discharges the M1
    acceptance criterion of cross-validating ≥3 peptide sequences against the
    BioPython-backed pipeline. Charge / pI use IPC 2.0 peptide pKa with Cys
    included; ε / MW / GRAVY / II / AI / aromaticity use BioPython tables.
    """

    @pytest.mark.parametrize(
        "case",
        PEPTIDE_GOLDENS,
        ids=[c["name"] for c in PEPTIDE_GOLDENS],
    )
    def test_peptide_property(self, case: dict):
        seq = case["seq"]
        assert charge_at_ph(seq, 7.0, IPC2_PEPTIDE, include_cys=True) == pytest.approx(
            case["charge"], abs=ABS_TOL
        )
        assert isoelectric_point(seq, IPC2_PEPTIDE, include_cys=True) == pytest.approx(
            case["pi"], abs=ABS_TOL
        )
        assert gravy(seq) == pytest.approx(case["gravy"], abs=ABS_TOL)
        assert molecular_weight(seq) == pytest.approx(case["mw"], abs=ABS_TOL)
        eox, ered = extinction_coefficients(seq)
        assert eox == pytest.approx(case["eox"], abs=ABS_TOL)
        assert ered == pytest.approx(case["ered"], abs=ABS_TOL)
        ii = instability_index(seq)
        if case["instability"] is None:
            assert ii is None
        else:
            assert ii == pytest.approx(case["instability"], abs=ABS_TOL)
        assert aliphatic_index(seq) == pytest.approx(case["aliphatic"], abs=ABS_TOL)
        assert aromaticity(seq) == pytest.approx(case["aromaticity"], abs=ABS_TOL)
