"""Characterization tests — pin the exact numeric output for representative
sequences so refactoring drift is caught.

These are NOT correctness checks against an external reference (those live in
test_properties.py with closed-form expectations like instability=9.0 for
poly-A-10). They snapshot the *current* output to pin behaviour during
refactoring. If a value changes intentionally (formula fix, pKa update),
update the golden number here as part of the change.

The reconstructed VH / VL pair below matches the `ab_full_paired` row in
the e2e corpus, so every refactor exercises the same sequence path the
corpus does — but at numeric precision the corpus's range bounds don't
catch.
"""

from __future__ import annotations

import pytest

from pka_tables import IPC2_PROTEIN
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
            ("charge", 0.000578),
            ("pi", 7.018372),
            ("gravy", -0.111111),
            ("mw", 6050.730200),
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
        assert charge_at_ph(VL, 7.0, IPC2_PROTEIN, include_cys=False) == pytest.approx(1.994916, abs=ABS_TOL)

    def test_vl_pi(self):
        assert isoelectric_point(VL, IPC2_PROTEIN, include_cys=False) == pytest.approx(9.798889, abs=ABS_TOL)


class TestFvGoldenValues:
    """Fv (paired VH + VL). pI bisects the per-chain charge sum."""

    def test_fv_charge(self):
        # Additive: must equal charge(VH) + charge(VL).
        assert fv_charge(VH, VL, 7.0, IPC2_PROTEIN) == pytest.approx(1.995494, abs=ABS_TOL)

    def test_fv_pi(self):
        # Per-chain charge sum bisection — distinct from concatenated-chain pI.
        assert fv_isoelectric_point(VH, VL, IPC2_PROTEIN) == pytest.approx(9.330627, abs=ABS_TOL)

    def test_fv_extinction_oxidised(self):
        ox, _ = fv_extinction_coefficients(VH, VL)
        assert ox == pytest.approx(32430.0, abs=ABS_TOL)

    def test_fv_extinction_reduced(self):
        _, red = fv_extinction_coefficients(VH, VL)
        assert red == pytest.approx(32430.0, abs=ABS_TOL)

    def test_fv_mw(self):
        assert fv_molecular_weight(VH, VL) == pytest.approx(11517.686500, abs=ABS_TOL)
