"""Behavioral tests for src/properties.py.

Each test verifies an observable input → output relationship. Reference
values come from formula-based hand computations or well-established
literature; sources are noted in test docstrings/comments.

Run:
    cd blocks/sequence-properties/software
    uv sync
    uv run pytest tests/unit/test_properties.py
"""

from __future__ import annotations

import math

import pytest

from pka_tables import IPC2_PEPTIDE, IPC2_PROTEIN
from properties import (
    aa_counts,
    aa_fractions,
    aliphatic_index,
    aromaticity,
    charge_at_ph,
    charge_shift,
    clean_sequence,
    effective_length,
    extinction_coefficients,
    fv_charge,
    fv_charge_shift,
    fv_extinction_coefficients,
    fv_isoelectric_point,
    fv_molecular_weight,
    gravy,
    instability_index,
    is_invalid_sequence,
    isoelectric_point,
    molecular_weight,
)

# ---------------------------------------------------------------------------
# Sequence cleanup
# ---------------------------------------------------------------------------


class TestSequenceCleanup:
    """Edge cases at the input boundary: stop codons, gaps, ambiguity codes."""

    @pytest.mark.parametrize(
        "seq, expected",
        [
            (None, True),
            ("", True),
            ("ACDE*FG", True),  # stop codon → invalid as a whole, per spec
            ("*", True),
            ("ACDE", False),
            ("acde", False),  # case-insensitive — lower-case still valid
            ("ACDXBZ", False),  # non-stop ambiguities don't invalidate the whole seq
        ],
    )
    def test_is_invalid_sequence(self, seq: str | None, expected: bool):
        assert is_invalid_sequence(seq) is expected

    # X / B / Z / U / J / - all dropped; case folded to upper.
    @pytest.mark.parametrize(
        "seq, expected",
        [
            ("ACDE", "ACDE"),
            ("acde", "ACDE"),
            ("AC-DE", "ACDE"),
            ("ACXBZUJ", "AC"),
            # X / B / Z / U / J / - dropped; e (lower-case) → E (standard).
            ("aXbCdZeFg", "ACDEFG"),
        ],
    )
    def test_clean_sequence(self, seq: str, expected: str):
        assert clean_sequence(seq) == expected

    # Effective length excludes ambiguity codes — denominator for GRAVY etc.
    def test_effective_length_excludes_nonstandard(self):
        assert effective_length("ACX-BZdgY") == 5  # A C D G Y survive

    # AA-count over ambiguity-laden input — keeps only the 20 standard codes.
    def test_aa_counts_filters_nonstandard(self):
        counts = aa_counts("AAXC")
        assert counts["A"] == 2
        assert counts["C"] == 1
        # All other 18 residues == 0 (sum of values = effective length).
        assert sum(counts.values()) == 3


# ---------------------------------------------------------------------------
# Extinction coefficients — closed form, easy to verify
# ---------------------------------------------------------------------------


class TestExtinctionCoefficients:
    """ε = Y·1490 + W·5500 (+ floor(C/2)·125 when oxidised). Pace et al."""

    # Single Y, single W: ε_red = 1490 + 5500 = 6990; no Cys ⇒ ε_ox = ε_red.
    def test_yw_only(self):
        ox, red = extinction_coefficients("YW")
        assert red == pytest.approx(6990.0)
        assert ox == pytest.approx(6990.0)

    # Even Cys count: floor(2/2)·125 = 125 → ε_ox − ε_red = 125.
    def test_yw_with_paired_cys(self):
        ox, red = extinction_coefficients("YWCC")
        assert red == pytest.approx(6990.0)
        assert ox == pytest.approx(7115.0)

    # Odd Cys count: floor(1/2) = 0 — single Cys contributes no disulfide.
    def test_single_cys_does_not_contribute(self):
        ox, red = extinction_coefficients("YWC")
        assert ox == red

    # No Y, no W: ε = 0, NOT NA — A280 is not informative but the value is defined.
    def test_no_aromatic_yields_zero_not_na(self):
        ox, red = extinction_coefficients("AAA")
        assert ox == 0.0
        assert red == 0.0

    def test_invalid_sequence_returns_nones(self):
        assert extinction_coefficients("") == (None, None)
        assert extinction_coefficients("*") == (None, None)


# ---------------------------------------------------------------------------
# Molecular weight
# ---------------------------------------------------------------------------


class TestMolecularWeight:
    """MW = Σ residue masses + one H₂O (BioPython ProtParam mass table)."""

    # GGG: 3 × Gly + one H₂O. BioPython ProtParam value 189.1692 — pinned per
    # spec direction that BioPython is the reference. Tolerance 0.01 matches
    # the .1f display precision.
    def test_three_glycines(self):
        mw = molecular_weight("GGG")
        assert mw == pytest.approx(189.1692, abs=0.01)

    def test_invalid_sequence_returns_none(self):
        assert molecular_weight("") is None
        assert molecular_weight("***") is None

    # Non-standard residues skipped — denominator changes for GRAVY but MW
    # just sums what's present.
    def test_nonstandard_residues_skipped(self):
        mw_clean = molecular_weight("AGA")
        mw_with_x = molecular_weight("AGXA")  # X dropped — same residues remain
        assert mw_clean == pytest.approx(mw_with_x)


# ---------------------------------------------------------------------------
# GRAVY (Kyte-Doolittle hydropathy mean)
# ---------------------------------------------------------------------------


class TestGravy:
    """GRAVY = mean Kyte-Doolittle score over standard residues."""

    # All-A: every residue is +1.8 ⇒ GRAVY = +1.8.
    def test_homopolymer_alanine(self):
        assert gravy("AAA") == pytest.approx(1.8)

    # Mixed with non-standard: X excluded from both sum and denominator.
    def test_nonstandard_excluded_from_denominator(self):
        assert gravy("AXA") == pytest.approx(1.8)

    # Empty / invalid → NA.
    def test_invalid_yields_none(self):
        assert gravy("") is None
        assert gravy("XX") is None  # only non-standard ⇒ effective length 0
        assert gravy("ACG*") is None


# ---------------------------------------------------------------------------
# Charge & pI
# ---------------------------------------------------------------------------


class TestChargeAtPh:
    """Charge formula: termini contribute, side chains contribute per pKa.

    All tests use IPC 2.0 sets — values were transcribed from the 2021 paper.
    """

    # At a pH well below all acidic pKa's, the molecule is dominated by
    # the protonated N-terminus and protonated bases — strong positive.
    def test_low_ph_polylysine_strongly_positive(self):
        c = charge_at_ph("KKKKK", 1.0, IPC2_PEPTIDE)
        assert c > 4.0

    # At a pH well above all basic pKa's, the molecule is dominated by
    # deprotonated acids and the deprotonated C-terminus — strong negative.
    def test_high_ph_polyaspartate_strongly_negative(self):
        c = charge_at_ph("DDDDD", 13.0, IPC2_PEPTIDE)
        assert c < -4.0

    def test_invalid_sequence_returns_none(self):
        assert charge_at_ph("", 7.0, IPC2_PEPTIDE) is None
        assert charge_at_ph("ACG*", 7.0, IPC2_PEPTIDE) is None

    # Cys-include vs Cys-exclude must differ on a sequence containing Cys.
    # At pH 7 with peptide-set pKa_C=7.555, free Cys contributes ~−0.22 charge units.
    def test_cys_include_versus_exclude(self):
        c_with = charge_at_ph("C", 7.0, IPC2_PEPTIDE, include_cys=True)
        c_without = charge_at_ph("C", 7.0, IPC2_PEPTIDE, include_cys=False)
        assert c_with is not None and c_without is not None
        assert c_with < c_without
        # The difference is a single Cys side-chain charge contribution.
        diff = c_without - c_with
        expected = 1.0 / (1.0 + 10.0 ** (IPC2_PEPTIDE.side_chain["C"] - 7.0))
        assert diff == pytest.approx(expected, abs=1e-9)


class TestChargeShift:
    """ΔCharge (pH 7.4 → 6.0) — captures pH-switching capacity, dominated by His.

    All tests use IPC 2.0 sets. Histidine pKa ~6.0 sits inside the window;
    Asp / Glu / Lys / Arg are far outside and contribute ~0.
    """

    # Identity: charge_shift = charge(ph_from) − charge(ph_to). Property of the
    # implementation that should hold for any sequence, regardless of His count.
    @pytest.mark.parametrize("seq", ["ACDEFGHIK", "RRRR", "DDDD", "PEPTIDE", "MAGICK"])
    def test_matches_two_point_subtraction(self, seq: str):
        ds = charge_shift(seq, IPC2_PEPTIDE, include_cys=True)
        c_from = charge_at_ph(seq, 7.4, IPC2_PEPTIDE, include_cys=True)
        c_to = charge_at_ph(seq, 6.0, IPC2_PEPTIDE, include_cys=True)
        assert ds is not None and c_from is not None and c_to is not None
        assert ds == pytest.approx(c_from - c_to, abs=1e-9)

    # No His and no other titrators in [6.0, 7.4] window: |ΔCharge| stays small,
    # driven only by the IPC 2.0 N-terminus pKa (7.947) edge contribution.
    def test_no_histidine_small_magnitude(self):
        ds = charge_shift("AAAAAAAAAA", IPC2_PEPTIDE, include_cys=True)
        assert ds is not None
        assert abs(ds) < 0.5

    # Adding histidines monotonically grows |ΔCharge| — His is the metric's
    # primary driver per the spec.
    def test_histidine_count_dominates_magnitude(self):
        ds_0 = charge_shift("AAAAAAAAAA", IPC2_PEPTIDE, include_cys=True)
        ds_1 = charge_shift("AAAAAAAAAH", IPC2_PEPTIDE, include_cys=True)
        ds_3 = charge_shift("AAAAAAAHHH", IPC2_PEPTIDE, include_cys=True)
        ds_5 = charge_shift("AAAAAHHHHH", IPC2_PEPTIDE, include_cys=True)
        assert ds_0 is not None and ds_1 is not None and ds_3 is not None and ds_5 is not None
        # All negative (acidification gains positive charge → from − to is negative).
        assert ds_5 < ds_3 < ds_1 < ds_0
        # Each added His shifts ΔCharge measurably more negative (~−0.5 to
        # −0.7 per His at IPC 2.0 peptide pKa_H ~6.04, BioPython HH model).
        assert (ds_1 - ds_0) < -0.4
        assert (ds_5 - ds_0) < -2.0

    # NA propagation: invalid sequence at either pH point yields None.
    def test_invalid_sequence_returns_none(self):
        assert charge_shift("", IPC2_PEPTIDE) is None
        assert charge_shift("ACG*", IPC2_PEPTIDE) is None

    # Anti-symmetry: swapping the pH endpoints flips the sign.
    def test_swapped_endpoints_flip_sign(self):
        seq = "ACDEFGHIK"
        forward = charge_shift(seq, IPC2_PEPTIDE, include_cys=True, ph_from=7.4, ph_to=6.0)
        reverse = charge_shift(seq, IPC2_PEPTIDE, include_cys=True, ph_from=6.0, ph_to=7.4)
        assert forward is not None and reverse is not None
        assert forward == pytest.approx(-reverse, abs=1e-9)

    # Cys-include vs Cys-exclude affects ΔCharge for Cys-bearing sequences.
    # IPC 2.0 peptide pKa_C ~7.555 sits inside the 6.0-7.4 window, so Cys
    # contributes a non-trivial component when included.
    def test_cys_include_versus_exclude_differs_on_cys_sequence(self):
        with_cys = charge_shift("ACDEFG", IPC2_PEPTIDE, include_cys=True)
        without_cys = charge_shift("ACDEFG", IPC2_PEPTIDE, include_cys=False)
        assert with_cys is not None and without_cys is not None
        assert with_cys != pytest.approx(without_cys, abs=1e-3)


class TestFvChargeShift:
    """Fv ΔCharge — additivity and NA propagation across paired chains."""

    # Fv ΔCharge equals the sum of per-chain ΔCharge values (Cys-excluded,
    # protein pKa). Equivalent to (charge(VH, 7.4) + charge(VL, 7.4)) −
    # (charge(VH, 6.0) + charge(VL, 6.0)).
    @pytest.mark.parametrize(
        "vh, vl",
        [
            ("EVQLVQSGGGLVQPGGSLRLSCAAS", "DIQMTQSPSSLSASVGDRVTITC"),
            ("EVQLVQSGGGHHHLVQPGGSLRLSCAAS", "DIQMTHHHQSPSSLSASVGDRVTITC"),
        ],
        ids=["natural", "his_rich"],
    )
    def test_additivity_across_chains(self, vh: str, vl: str):
        s_fv = fv_charge_shift(vh, vl, IPC2_PROTEIN)
        s_vh = charge_shift(vh, IPC2_PROTEIN, include_cys=False)
        s_vl = charge_shift(vl, IPC2_PROTEIN, include_cys=False)
        assert s_fv is not None and s_vh is not None and s_vl is not None
        assert s_fv == pytest.approx(s_vh + s_vl, abs=1e-9)

    # Either chain invalid → Fv ΔCharge None.
    def test_invalid_chain_returns_none(self):
        assert fv_charge_shift("", "DIQMTQ", IPC2_PROTEIN) is None
        assert fv_charge_shift("EVQLVQ", "", IPC2_PROTEIN) is None
        assert fv_charge_shift("EVQLV*Q", "DIQMTQ", IPC2_PROTEIN) is None


class TestIsoelectricPoint:
    """pI from bisection over [0, 14]."""

    # Glycine peptide — only N- and C-termini ionise. Theoretical pI sits
    # between the two terminal pKa's (≈ 9.094 and ≈ 2.869 → ~5.98 for the peptide set).
    def test_pi_polyglycine_between_termini(self):
        pi = isoelectric_point("GGGG", IPC2_PEPTIDE)
        assert pi is not None
        assert 4.5 < pi < 7.5

    # No zero crossing → NA. Spec §Defaults: "Possible in extreme all-basic
    # or all-acidic synthetic sequences" — the same-sign-endpoint guard must
    # emit NA, not loop forever or clamp to a boundary.
    #
    # Construct a synthetic pKa set that's pure-base at every group; under
    # this set the molecule charge is positive at every pH in [0, 14] and
    # the bisection has no zero to find. The polybasic R*50 example used to
    # trigger this under IPC 1.0 pKa_R = 12.503 but doesn't under IPC 2.0
    # (pKa_R = 10.223 with C-term acid pKa = 6.065 produces a real zero
    # crossing around pH 12 even for pure-Arg sequences) — switching the
    # test to a direct demonstration of the no-crossing branch.
    def test_pi_no_zero_crossing_returns_none(self):
        from pka_tables import PKaSet

        all_base = PKaSet(
            name="synthetic_all_base",
            side_chain={"K": 100.0, "R": 100.0, "H": 100.0, "C": 100.0, "D": 100.0, "E": 100.0, "Y": 100.0},
            n_terminus=100.0,
            c_terminus=100.0,
        )
        pi = isoelectric_point("KRHKR", all_base, include_cys=True)
        assert pi is None

    # Polyacidic always crosses: even at pH 0, the protonated N-terminus
    # contributes ~+1, dominating the near-neutral acid side chains. Verifying
    # this guards against an over-eager NA shortcut on legitimate sequences.
    def test_pi_polyacidic_has_crossing(self):
        pi = isoelectric_point("DDDDDDDDDD", IPC2_PEPTIDE)
        assert pi is not None
        assert 0.0 < pi < 14.0

    # The result of bisection should give a charge near zero at the returned pI.
    def test_pi_returns_charge_near_zero(self):
        seq = "KEDARC"
        pi = isoelectric_point(seq, IPC2_PEPTIDE)
        assert pi is not None
        c_at_pi = charge_at_ph(seq, pi, IPC2_PEPTIDE)
        assert abs(c_at_pi) < 1e-2

    # Different pKa sets give different pI for the same sequence — confirms
    # the protein-vs-peptide distinction is wired correctly.
    def test_peptide_vs_protein_set_differ(self):
        seq = "ACDEHKR"
        peptide_pi = isoelectric_point(seq, IPC2_PEPTIDE, include_cys=True)
        protein_pi = isoelectric_point(seq, IPC2_PROTEIN, include_cys=False)
        assert peptide_pi is not None
        assert protein_pi is not None
        assert peptide_pi != protein_pi


# ---------------------------------------------------------------------------
# Instability index
# ---------------------------------------------------------------------------


class TestInstabilityIndex:
    """Guruprasad index — needs effective length ≥ 10."""

    # 9 residues — below the 10-residue floor, NA per spec R9.
    def test_short_sequence_returns_none(self):
        assert instability_index("AAAAAAAAA") is None

    # Effective length floor uses cleaned length, not raw — XAAAA + XXXXX = 5 effective.
    def test_floor_uses_effective_length(self):
        assert instability_index("AAAAAAAAA" + "X" * 20) is None

    # 10 alanines — every dipeptide is AA → DIWV = 1.0 → II = 10/10 × 9 = 9.0.
    def test_polyalanine_known_value(self):
        ii = instability_index("AAAAAAAAAA")
        assert ii == pytest.approx(9.0)

    def test_invalid_sequence_returns_none(self):
        assert instability_index("") is None
        assert instability_index("AAAAAAAAA*") is None  # stop codon ⇒ NA


# ---------------------------------------------------------------------------
# Aliphatic index
# ---------------------------------------------------------------------------


class TestAliphaticIndex:
    """AI = 100 × (X_A + 2.9·X_V + 3.9·(X_I + X_L))."""

    # Equal-mole quartet — closed-form expectation.
    def test_alvi_equal_fractions(self):
        ai = aliphatic_index("ALVI")
        # X_A = X_V = X_I = X_L = 0.25
        # AI = 100 × (0.25 + 2.9·0.25 + 3.9·0.5) = 100 × 2.925 = 292.5
        assert ai == pytest.approx(292.5)

    # Sequence with only non-aliphatic residues yields zero.
    def test_no_aliphatic_yields_zero(self):
        assert aliphatic_index("FFFF") == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# Aromaticity
# ---------------------------------------------------------------------------


class TestAromaticity:
    """Aromaticity = (F + W + Y) / effective_length."""

    @pytest.mark.parametrize(
        "seq, expected",
        [
            ("FFF", 1.0),
            ("FYW", 1.0),
            ("AAA", 0.0),
            ("FA", 0.5),
            # Non-standard excluded from denominator — F over (F+A) = 0.5, X dropped.
            ("FAX", 0.5),
        ],
    )
    def test_aromaticity_values(self, seq: str, expected: float):
        result = aromaticity(seq)
        assert result == pytest.approx(expected)


# ---------------------------------------------------------------------------
# AA fractions
# ---------------------------------------------------------------------------


class TestAaFractions:
    """20 standard AAs always present; fractions sum to 1.0."""

    def test_homopolymer_a_fraction_one(self):
        f = aa_fractions("AAAAA")
        assert f is not None
        assert f["A"] == pytest.approx(1.0)
        # Every other AA is zero.
        assert sum(v for k, v in f.items() if k != "A") == 0.0

    def test_fractions_sum_to_one(self):
        f = aa_fractions("ACDEFGHIKLMNPQRSTVWY")
        assert f is not None
        assert math.isclose(sum(f.values()), 1.0, abs_tol=1e-9)

    # Non-standard residues dropped from both numerator and denominator.
    def test_nonstandard_dropped(self):
        f_clean = aa_fractions("AC")
        f_with_x = aa_fractions("AXC")
        assert f_clean == f_with_x

    def test_invalid_returns_none(self):
        assert aa_fractions("") is None
        assert aa_fractions("XX") is None
        assert aa_fractions("AC*") is None


# ---------------------------------------------------------------------------
# Fv (paired-chain) properties
# ---------------------------------------------------------------------------


class TestFvProperties:
    """Fv columns must use the per-chain sum, NOT a concatenated string."""

    # MW additive: simple sum of two chain weights.
    def test_fv_mw_additive(self):
        vh = "AAAA"
        vl = "GGGG"
        mw_h = molecular_weight(vh)
        mw_l = molecular_weight(vl)
        mw_fv = fv_molecular_weight(vh, vl)
        assert mw_fv == pytest.approx(mw_h + mw_l)

    # Charge additive at any pH.
    def test_fv_charge_additive(self):
        vh = "ACDEFGHIKLMNPQRSTVWY"
        vl = "ACDEFGHIKLMNPQRSTVWY"
        c_fv = fv_charge(vh, vl, 7.0, IPC2_PROTEIN)
        c_h = charge_at_ph(vh, 7.0, IPC2_PROTEIN, include_cys=False)
        c_l = charge_at_ph(vl, 7.0, IPC2_PROTEIN, include_cys=False)
        assert c_fv == pytest.approx(c_h + c_l)

    # ε_ox additive — including disulfide bonds within each chain.
    def test_fv_extinction_additive(self):
        vh = "WYCC"  # 1W + 1Y + 1 disulfide bond
        vl = "WYCC"
        ox_fv, red_fv = fv_extinction_coefficients(vh, vl)
        ox_h, red_h = extinction_coefficients(vh)
        ox_l, red_l = extinction_coefficients(vl)
        assert ox_fv == pytest.approx(ox_h + ox_l)
        assert red_fv == pytest.approx(red_h + red_l)

    # Fv pI: by spec, bisect on the per-chain sum, not on a concatenated string.
    # The two approaches differ when the chains have different per-chain
    # charge profiles. We verify the correct (per-chain-sum) result by
    # checking that the function evaluates to ~0 at the returned pI.
    def test_fv_pi_zero_crossing(self):
        vh = "ACDEHKR"
        vl = "ACDEHKR"
        pi = fv_isoelectric_point(vh, vl, IPC2_PROTEIN)
        assert pi is not None
        c_at_pi = fv_charge(vh, vl, pi, IPC2_PROTEIN)
        assert abs(c_at_pi) < 1e-2

    # Either chain invalid → all Fv outputs are NA.
    def test_invalid_chain_propagates_na(self):
        assert fv_charge("ACDE", "", 7.0, IPC2_PROTEIN) is None
        assert fv_isoelectric_point("ACDE", "*", IPC2_PROTEIN) is None
        assert fv_molecular_weight("", "ACDE") is None
        assert fv_extinction_coefficients("ACDE*", "GGGG") == (None, None)
