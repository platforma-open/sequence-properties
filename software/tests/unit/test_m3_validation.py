"""M3 external-validation cross-checks.

Discharges the M3 (Validation) acceptance criteria from the spec:

* pI for >=5 VH sequences within 0.1 pH of the IPC 2.0 reference.
* VL pI on >=2 sequences.
* Fv charge / Fv pI manually verified on >=2 paired sequences.
* CDR-L3 charge by manual Henderson-Hasselbalch on >=3 sequences.
* CDR-H3 charge by manual computation on >=10 sequences.
* Aliphatic for >=3 VH sequences with closed-form check.

Two cross-check shapes:

1. **Pinned IPC 2.0 webserver values** (`TestVHPIWebserverCrosscheck`,
   `TestVLPIWebserverCrosscheck`). Submitted once to
   https://ipc2.mimuw.edu.pl/ and pinned here. Spec L518 says the
   webserver may be unavailable in CI; pinning avoids the dependency.

   The webserver includes Cys as ionizable; production code excludes Cys
   for full-chain VH/VL per spec L535-540 (disulfide-bonded). The
   cross-check calls `isoelectric_point(..., include_cys=True)` to
   isolate pKa transcription and HH formula from the Cys policy.

2. **Independent closed-form Henderson-Hasselbalch** (`_ref_charge_hh`).
   Pure-Python textbook HH with our IPC 2.0 pKa values, no BioPython.
   Production code IS BioPython's `IsoelectricPoint` under the hood, so
   this reference catches a BioPython formula bug.
"""

from __future__ import annotations

from collections import Counter

import pytest

from aa_tables import STANDARD_AA_SET, STANDARD_AAS
from pka_tables import IPC2_PEPTIDE, IPC2_PROTEIN, PKaSet
from properties import (
    aliphatic_index,
    charge_at_ph,
    fv_charge,
    fv_isoelectric_point,
    isoelectric_point,
)


# ---------------------------------------------------------------------------
# Reference sequences
# ---------------------------------------------------------------------------

# VH and VL pairs. VH1/VL1 is the synthetic pair from the e2e corpus; the
# others are well-known therapeutic-antibody V-region sequences (free in
# DrugBank). VL2 is paired with VH2 (trastuzumab).
VH_SYNTHETIC = "EVQLVESGFTFSSYAMSWVRQISGSGGSTYYAESVKGRFTICARDYWWGQGTLV"
VH_TRASTUZUMAB = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGR"
    "FTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
VH_CETUXIMAB = (
    "QVQLKQSGPGLVQPSQSLSITCTVSGFSLTNYGVHWVRQSPGKGLEWLGVIWSGGNTDYNTPFTSRL"
    "SINKDNSKSQVFFKMNSLQSNDTAIYYCARALTYYDYEFAYWGQGTLVTVSA"
)
VH_ADALIMUMAB = (
    "EVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSAITWNSGHIDYADSVEGR"
    "FTISRDNAKNSLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSS"
)
VH_PEMBROLIZUMAB = (
    "QVQLVQSGVEVKKPGASVKVSCKASGYTFTNYYMYWVRQAPGQGLEWMGGINPSNGGTNFNEKFKNR"
    "VTLTTDSSTTTAYMELKSLQFDDTAVYYCARRDYRFDMGFDYWGQGTTVTVSS"
)

VL_SYNTHETIC = "DIQMTQSQSISSYLNWYQQKAASSLQSGVPSRFSGSGCQQYNSFGQGTKV"
VL_TRASTUZUMAB = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRS"
    "GTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)


ABS_TOL = 1e-9  # impl vs reference HH agrees to floating-point noise


# ---------------------------------------------------------------------------
# Independent reference Henderson-Hasselbalch (pure-Python, no BioPython)
# ---------------------------------------------------------------------------


def _ref_charge_hh(seq: str, ph: float, pka_set: PKaSet, include_cys: bool) -> float:
    """Net charge by textbook Henderson-Hasselbalch with the given IPC 2.0
    pKa set. Pure-Python; no BioPython dependency.

    Basic group at pH:  f_i = +1 / (1 + 10**(pH - pKa_i))
    Acidic group at pH: f_j = -1 / (1 + 10**(pKa_j - pH))
    Net charge = sum(basic fractions) - sum(acidic fractions).

    Cys is in `pka_set.side_chain` for both IPC2_PEPTIDE and IPC2_PROTEIN
    but is included as an acid only when `include_cys=True` (peptide / CDR3
    rules). Full-chain VH/VL excludes Cys (disulfide-bonded).
    """
    cleaned = "".join(c for c in seq.upper() if c in STANDARD_AA_SET)
    counts = Counter(cleaned)
    pos = 1.0 / (1.0 + 10.0 ** (ph - pka_set.n_terminus))
    for aa in ("K", "R", "H"):
        pos += counts[aa] / (1.0 + 10.0 ** (ph - pka_set.side_chain[aa]))
    neg = 1.0 / (1.0 + 10.0 ** (pka_set.c_terminus - ph))
    for aa in ("D", "E", "Y"):
        neg += counts[aa] / (1.0 + 10.0 ** (pka_set.side_chain[aa] - ph))
    if include_cys:
        neg += counts["C"] / (1.0 + 10.0 ** (pka_set.side_chain["C"] - ph))
    return pos - neg


# ---------------------------------------------------------------------------
# 1. VH pI vs IPC 2.0 webserver (>=5 sequences)
# ---------------------------------------------------------------------------

# pI values reported by the IPC 2.0 web server (IPC2_protein consensus model)
# for each VH sequence below. The webserver includes Cys as ionizable; we
# therefore call `isoelectric_point` with `include_cys=True` to cross-check
# the pKa transcription and HH formula independently of the spec's Cys
# policy. The webserver's reported precision is 2 decimal places.
VH_WEBSERVER_PI = [
    ("VH_synthetic",     VH_SYNTHETIC,     6.00),
    ("VH_trastuzumab",   VH_TRASTUZUMAB,   6.95),
    ("VH_cetuximab",     VH_CETUXIMAB,     7.77),
    ("VH_adalimumab",    VH_ADALIMUMAB,    5.22),
    ("VH_pembrolizumab", VH_PEMBROLIZUMAB, 7.73),
]


class TestVHPIWebserverCrosscheck:
    """Pins IPC 2.0 webserver pI for >=5 VH sequences. Discharges M3:
    'pI for >=5 VH within 0.1 pH of IPC 2.0 reference'. Tolerance 0.05
    pH — tighter than the spec's 0.1 — catches drift from the webserver's
    2-dp reported precision."""

    @pytest.mark.parametrize("name,seq,webserver_pi", VH_WEBSERVER_PI)
    def test_vh_pi_matches_webserver(self, name: str, seq: str, webserver_pi: float):
        ours = isoelectric_point(seq, IPC2_PROTEIN, include_cys=True)
        assert ours == pytest.approx(webserver_pi, abs=0.05), (
            f"{name}: ours={ours:.4f} vs webserver IPC2_protein={webserver_pi}"
        )


# ---------------------------------------------------------------------------
# 2. VL pI vs IPC 2.0 webserver (>=2 sequences)
# ---------------------------------------------------------------------------

VL_WEBSERVER_PI = [
    ("VL_synthetic",   VL_SYNTHETIC,   8.38),
    ("VL_trastuzumab", VL_TRASTUZUMAB, 7.77),
]


class TestVLPIWebserverCrosscheck:
    """VL pI on >=2 sequences vs IPC 2.0 webserver."""

    @pytest.mark.parametrize("name,seq,webserver_pi", VL_WEBSERVER_PI)
    def test_vl_pi_matches_webserver(self, name: str, seq: str, webserver_pi: float):
        ours = isoelectric_point(seq, IPC2_PROTEIN, include_cys=True)
        assert ours == pytest.approx(webserver_pi, abs=0.05), (
            f"{name}: ours={ours:.4f} vs webserver IPC2_protein={webserver_pi}"
        )


# ---------------------------------------------------------------------------
# 3. VH/VL pI snapshot under spec convention (Cys excluded for full chains)
# ---------------------------------------------------------------------------


class TestVHVLPISpecConvention:
    """Pins production VH/VL pI (spec convention: `include_cys=False`,
    spec L535-540). Webserver values above use `include_cys=True`; this
    class catches Cys-handling regressions."""

    @pytest.mark.parametrize(
        "name,seq,expected",
        [
            ("VH_synthetic",     VH_SYNTHETIC,     6.006653),
            ("VH_trastuzumab",   VH_TRASTUZUMAB,   7.455017),
            ("VH_cetuximab",     VH_CETUXIMAB,     8.602600),
            ("VH_adalimumab",    VH_ADALIMUMAB,    5.218811),
            ("VH_pembrolizumab", VH_PEMBROLIZUMAB, 8.427429),
            ("VL_synthetic",     VL_SYNTHETIC,     9.165710),
            ("VL_trastuzumab",   VL_TRASTUZUMAB,   8.592346),
        ],
    )
    def test_pi_no_cys(self, name: str, seq: str, expected: float):
        actual = isoelectric_point(seq, IPC2_PROTEIN, include_cys=False)
        assert actual == pytest.approx(expected, abs=1e-5), (
            f"{name}: pi(no_Cys)={actual} expected~{expected}"
        )


# ---------------------------------------------------------------------------
# 4. CDR-H3 charge by manual Henderson-Hasselbalch (>=10 sequences)
# ---------------------------------------------------------------------------

# Mix of corpus (CARDYW), therapeutic CDR-H3 fragments, and synthetic
# sequences spanning length range and charge polarity.
CDRH3_SEQUENCES = [
    "CARDYW",                    # corpus synthetic
    "CSRWGGDGFY",                # trastuzumab CDR-H3 fragment
    "CARALTYYDYEFAY",            # cetuximab CDR-H3
    "CAKVSYLSTASSLDY",           # adalimumab CDR-H3
    "CARRDYRFDMGFDY",            # pembrolizumab CDR-H3
    "CARGSDPHGFDY",              # synthetic, mid-length
    "CSAQGGYPVT",                # synthetic, neutral
    "CARWHRLDFY",                # synthetic, His-bearing
    "CQHRGGTPLD",                # synthetic, basic
    "CARDDGYY",                  # synthetic, short / acidic
    "CARGSGGFDPMGTDIWGQGTLVT",   # synthetic, long
]


class TestCDRH3ChargeManualHH:
    """CDR-H3 net charge at pH 7 vs a textbook HH reference. >=10
    sequences, peptide pKa, Cys included per spec (CDR-H3 Cys = free
    thiol; intra-CDR3 disulfides not detected)."""

    @pytest.mark.parametrize("seq", CDRH3_SEQUENCES)
    def test_cdrh3_charge_at_ph7(self, seq: str):
        impl = charge_at_ph(seq, 7.0, IPC2_PEPTIDE, include_cys=True)
        ref = _ref_charge_hh(seq, 7.0, IPC2_PEPTIDE, include_cys=True)
        assert impl == pytest.approx(ref, abs=ABS_TOL)


# ---------------------------------------------------------------------------
# 5. CDR-L3 charge by manual Henderson-Hasselbalch (>=3 sequences)
# ---------------------------------------------------------------------------

CDRL3_SEQUENCES = [
    "CQQYNS",       # corpus synthetic
    "CQQHYTTPPT",   # trastuzumab CDR-L3
    "CSAQGGYPVT",   # synthetic, neutral
    "CQQYNSLPYT",   # synthetic, neutral / Tyr-rich
    "CHQYHRSPT",    # synthetic, two His + Arg
]


class TestCDRL3ChargeManualHH:
    """CDR-L3 net charge at pH 7 vs the same HH reference."""

    @pytest.mark.parametrize("seq", CDRL3_SEQUENCES)
    def test_cdrl3_charge_at_ph7(self, seq: str):
        impl = charge_at_ph(seq, 7.0, IPC2_PEPTIDE, include_cys=True)
        ref = _ref_charge_hh(seq, 7.0, IPC2_PEPTIDE, include_cys=True)
        assert impl == pytest.approx(ref, abs=ABS_TOL)


# ---------------------------------------------------------------------------
# 6. Fv charge / pI manually verified on >=2 paired sequences
# ---------------------------------------------------------------------------


def _manual_fv_charge(vh: str, vl: str, ph: float, pka_set: PKaSet) -> float:
    """Reference Fv net charge = HH(VH) + HH(VL), Cys excluded per
    full-chain rule. Independent of `fv_charge` and BioPython."""
    return _ref_charge_hh(vh, ph, pka_set, include_cys=False) + _ref_charge_hh(
        vl, ph, pka_set, include_cys=False
    )


def _manual_fv_pi(vh: str, vl: str, pka_set: PKaSet) -> float:
    """Reference Fv pI = bisect [0, 14] for pH where HH(VH)+HH(VL) = 0.
    Tolerance 1e-4 — tighter than impl's 1e-3, so an impl result within
    its own tolerance still satisfies the 1e-3 assertion."""
    f = lambda ph: _manual_fv_charge(vh, vl, ph, pka_set)
    lo, hi = 0.0, 14.0
    f_lo = f(lo)
    while hi - lo > 1e-4:
        mid = 0.5 * (lo + hi)
        f_mid = f(mid)
        if (f_mid > 0) == (f_lo > 0):
            lo, f_lo = mid, f_mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


class TestFvManualVerification:
    """Fv charge and pI on >=2 paired chains vs the HH reference
    (`_manual_fv_charge`, `_manual_fv_pi`)."""

    @pytest.mark.parametrize(
        "name,vh,vl",
        [
            ("synthetic",   VH_SYNTHETIC,   VL_SYNTHETIC),
            ("trastuzumab", VH_TRASTUZUMAB, VL_TRASTUZUMAB),
        ],
    )
    def test_fv_charge_at_ph7(self, name: str, vh: str, vl: str):
        impl = fv_charge(vh, vl, 7.0, IPC2_PROTEIN)
        ref = _manual_fv_charge(vh, vl, 7.0, IPC2_PROTEIN)
        assert impl == pytest.approx(ref, abs=ABS_TOL)

    @pytest.mark.parametrize(
        "name,vh,vl",
        [
            ("synthetic",   VH_SYNTHETIC,   VL_SYNTHETIC),
            ("trastuzumab", VH_TRASTUZUMAB, VL_TRASTUZUMAB),
        ],
    )
    def test_fv_pi(self, name: str, vh: str, vl: str):
        impl = fv_isoelectric_point(vh, vl, IPC2_PROTEIN)
        ref = _manual_fv_pi(vh, vl, IPC2_PROTEIN)
        # Both bisections converge to within their respective tolerances
        # (impl 1e-3, ref 1e-4); accept agreement within 1e-3.
        assert impl == pytest.approx(ref, abs=1e-3)


# ---------------------------------------------------------------------------
# 7. Aliphatic index closed-form check (>=3 VH)
# ---------------------------------------------------------------------------


def _closed_form_aliphatic(seq: str) -> float:
    """Ikai 1980: AI = 100 * (X_A + 2.9*X_V + 3.9*(X_I + X_L)).
    Pure-Python; no impl dependency."""
    cleaned = "".join(c for c in seq.upper() if c in STANDARD_AA_SET)
    counts = Counter(cleaned)
    n = sum(counts[a] for a in STANDARD_AAS)
    return 100.0 * (
        counts["A"] / n + 2.9 * counts["V"] / n + 3.9 * (counts["I"] + counts["L"]) / n
    )


class TestAliphaticIndexClosedForm:
    """Aliphatic index for >=3 VH vs the Ikai 1980 closed-form."""

    @pytest.mark.parametrize(
        "name,seq",
        [
            ("VH_synthetic",   VH_SYNTHETIC),
            ("VH_trastuzumab", VH_TRASTUZUMAB),
            ("VH_cetuximab",   VH_CETUXIMAB),
            ("VH_adalimumab",  VH_ADALIMUMAB),
        ],
    )
    def test_aliphatic_matches_closed_form(self, name: str, seq: str):
        impl = aliphatic_index(seq)
        ref = _closed_form_aliphatic(seq)
        assert impl == pytest.approx(ref, abs=1e-9)
