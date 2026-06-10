"""Vectorized compute substrate. A column of sequences -> per-residue count
matrix + length + validity, matching properties.py's scalar cleaning exactly.
Single-threaded by construction (pure numpy elementwise + per-seq counting)."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from Bio.Data.IUPACData import protein_weights
from Bio.SeqUtils.ProtParam import ProtParamData

from aa_tables import STANDARD_AAS
from pka_tables import PKaSet

_AA_INDEX = {aa: i for i, aa in enumerate(STANDARD_AAS)}

# --------------------------------------------------------------------------- #
# Constant vectors in STANDARD_AAS order, sourced verbatim from BioPython so
# the vectorized math uses the SAME numbers the scalar oracle (properties.py)
# resolves through ProteinAnalysis. Hand-retyping any of these would silently
# break per-property parity.
# --------------------------------------------------------------------------- #

# Kyte-Doolittle hydropathy (Bio.SeqUtils.ProtParam.ProtParamData.kd).
_KD = np.array([ProtParamData.kd[aa] for aa in STANDARD_AAS], dtype=np.float64)

# Average per-residue masses (Bio.Data.IUPACData.protein_weights).
_MASS = np.array([protein_weights[aa] for aa in STANDARD_AAS], dtype=np.float64)

# Average water mass subtracted per peptide bond. Bio/SeqUtils/__init__.py
# molecular_weight() uses water = 18.0153 for the non-monoisotopic protein path,
# and ProteinAnalysis.molecular_weight() calls it with monoisotopic=False.
_WATER = 18.0153

# Aromatic indicator (F, W, Y) — ProteinAnalysis.aromaticity() = sum of the
# F/W/Y mole fractions.
_AROMATIC = np.array([1.0 if aa in "FWY" else 0.0 for aa in STANDARD_AAS], dtype=np.float64)

# Reduced extinction contributions: Trp 5500, Tyr 1490 (per
# ProteinAnalysis.molar_extinction_coefficient(): W*5500 + Y*1490).
_EXT_RED = np.array(
    [5500.0 if aa == "W" else 1490.0 if aa == "Y" else 0.0 for aa in STANDARD_AAS],
    dtype=np.float64,
)

_C_INDEX = STANDARD_AAS.index("C")
_A_INDEX = _AA_INDEX["A"]
_V_INDEX = _AA_INDEX["V"]
_I_INDEX = _AA_INDEX["I"]
_L_INDEX = _AA_INDEX["L"]


def _mask(arr: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Return a float64 copy of `arr` with NaN wherever `valid` is False."""
    out = np.asarray(arr, dtype=np.float64).copy()
    out[~valid] = np.nan
    return out


def _safe_div(num: np.ndarray, length: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """`num / length` as float64, NaN where invalid or length == 0."""
    out = np.full(len(length), np.nan, dtype=np.float64)
    ok = valid & (length > 0)
    out[ok] = np.asarray(num, dtype=np.float64)[ok] / length[ok]
    return out


@dataclass(frozen=True)
class Substrate:
    counts: np.ndarray  # (N, 20) int64, STANDARD_AAS order
    length: np.ndarray  # (N,) int64 effective length
    valid: np.ndarray  # (N,) bool


def build_counts(seqs: list[str | None]) -> Substrate:
    n = len(seqs)
    counts = np.zeros((n, 20), dtype=np.int64)
    valid = np.zeros(n, dtype=bool)
    for i, s in enumerate(seqs):
        if s is None or s == "" or "*" in s:
            continue
        row = counts[i]
        any_std = False
        for c in s.upper():
            j = _AA_INDEX.get(c)
            if j is not None:
                row[j] += 1
                any_std = True
        valid[i] = any_std
    length = counts.sum(axis=1)
    return Substrate(counts=counts, length=length, valid=valid)


# --------------------------------------------------------------------------- #
# Linear properties — pure array ops over the count matrix. Each returns a
# float64 array (or pair of arrays) with NaN for invalid rows, matching the
# scalar properties.py oracle within 1e-6.
# --------------------------------------------------------------------------- #


def gravy(sub: Substrate) -> np.ndarray:
    """Kyte-Doolittle GRAVY: (sum of per-residue hydropathy) / length.

    Mirrors ProteinAnalysis.gravy(), which returns total_gravy / length.
    """
    return _safe_div(sub.counts @ _KD, sub.length, sub.valid)


def molecular_weight(sub: Substrate) -> np.ndarray:
    """Average mass (Da): sum(residue masses) - (length - 1) * water.

    Mirrors Bio.SeqUtils.molecular_weight(seq, "protein"); water = 18.0153.
    """
    mw = sub.counts @ _MASS - (sub.length - 1) * _WATER
    return _mask(mw, sub.valid & (sub.length > 0))


def aromaticity(sub: Substrate) -> np.ndarray:
    """Aromatic mole fraction (F + W + Y) / length."""
    return _safe_div(sub.counts @ _AROMATIC, sub.length, sub.valid)


def aliphatic_index(sub: Substrate) -> np.ndarray:
    """Ikai aliphatic index: 100 * (X_A + 2.9*X_V + 3.9*(X_I + X_L)).

    X_aa is the mole fraction; the shared `/length` is folded into _safe_div.
    """
    c = sub.counts
    num = 100.0 * (c[:, _A_INDEX] + 2.9 * c[:, _V_INDEX] + 3.9 * (c[:, _I_INDEX] + c[:, _L_INDEX]))
    return _safe_div(num, sub.length, sub.valid)


def extinction(sub: Substrate) -> tuple[np.ndarray, np.ndarray]:
    """Pace extinction coefficients at 280 nm, returned as (oxidized, reduced).

    reduced  = nW*5500 + nY*1490
    oxidized = reduced + (nC // 2) * 125   (cystine Cys-Cys bonds)

    Matches the scalar extinction_coefficients() column order (oxidized first);
    BioPython's molar_extinction_coefficient() returns (reduced, oxidized).
    """
    red = sub.counts @ _EXT_RED
    ox = red + (sub.counts[:, _C_INDEX] // 2) * 125.0
    return _mask(ox, sub.valid), _mask(red, sub.valid)


def aa_fractions(sub: Substrate) -> np.ndarray:
    """Per-AA mole fractions, (N, 20) float64 in STANDARD_AAS order.

    Mirrors ProteinAnalysis.amino_acids_percent / 100 = count / length.
    Invalid rows (and length == 0) are all-NaN.
    """
    frac = np.divide(
        sub.counts,
        sub.length[:, None],
        out=np.full(sub.counts.shape, np.nan, dtype=np.float64),
        where=(sub.length[:, None] > 0),
    )
    frac[~sub.valid] = np.nan
    return frac.astype(np.float64)


# --------------------------------------------------------------------------- #
# Charge — vectorized Henderson-Hasselbalch, byte-for-byte with BioPython's
# Bio.SeqUtils.IsoelectricPoint.charge_at_pH:
#
#   positive_charge = sum over pos_pKs of  content[aa] / (10**(pH - pK) + 1.0)
#   negative_charge = sum over neg_pKs of  content[aa] / (10**(pK - pH) + 1.0)
#   charge          = positive_charge - negative_charge
#
# content[Nterm] = content[Cterm] = 1.0 (one terminus pair per sequence);
# content[K/R/H/D/E/Y/C] = residue count. The scalar oracle injects IPC 2.0
# pKa overrides (properties._ipc2_isoelectric_point): pos_pKs = {Nterm, K, R, H},
# neg_pKs = {Cterm, D, E, Y} and C only when include_cys and "C" in side_chain.
# --------------------------------------------------------------------------- #


def _charge_raw(sub: Substrate, ph: float, pka_set: PKaSet, include_cys: bool) -> np.ndarray:
    """Vectorized BioPython charge_at_pH over ALL rows — finite for every row,
    invalid rows are NOT masked here. The public wrappers apply NaN via _mask.

    Kept finite-for-all-rows on purpose: the pI bisection (Task 8) reuses this
    over the same count matrix and needs a clean charge value per row.
    """
    c = sub.counts
    # One N-terminus per sequence (count = 1), basic side chains K/R/H.
    pos = 1.0 / (10 ** (ph - pka_set.n_terminus) + 1.0)
    for aa in ("K", "R", "H"):
        pk = pka_set.side_chain[aa]
        pos = pos + c[:, _AA_INDEX[aa]] * (1.0 / (10 ** (ph - pk) + 1.0))

    # One C-terminus per sequence (count = 1), acidic side chains D/E/Y (+C).
    neg = 1.0 / (10 ** (pka_set.c_terminus - ph) + 1.0)
    neg_aas = ["D", "E", "Y"]
    if include_cys and "C" in pka_set.side_chain:
        neg_aas.append("C")
    for aa in neg_aas:
        pk = pka_set.side_chain[aa]
        neg = neg + c[:, _AA_INDEX[aa]] * (1.0 / (10 ** (pk - ph) + 1.0))

    return pos - neg


def charge_at_ph(sub: Substrate, ph: float, pka_set: PKaSet, include_cys: bool = True) -> np.ndarray:
    """Net charge at a given pH, float64, NaN for invalid rows.

    Mirrors properties.charge_at_ph (BioPython IsoelectricPoint.charge_at_pH
    with IPC 2.0 pKa overrides).
    """
    return _mask(_charge_raw(sub, ph, pka_set, include_cys), sub.valid)


def charge_shift(
    sub: Substrate,
    pka_set: PKaSet,
    include_cys: bool = True,
    ph_from: float = 7.4,
    ph_to: float = 6.0,
) -> np.ndarray:
    """ΔCharge = charge(ph_from) - charge(ph_to), float64, NaN for invalid rows.

    Mirrors properties.charge_shift (defaults ph_from=7.4, ph_to=6.0).
    """
    hi = _charge_raw(sub, ph_from, pka_set, include_cys)
    lo = _charge_raw(sub, ph_to, pka_set, include_cys)
    return _mask(hi - lo, sub.valid)


def isoelectric_point(
    sub: Substrate,
    pka_set: PKaSet,
    include_cys: bool = True,
    lo: float = 0.0,
    hi: float = 14.0,
    tol: float = 1e-3,
) -> np.ndarray:
    """Vectorized pI: lockstep fixed-iteration bisection of `_charge_raw` over
    all rows, bit-identical to the scalar `properties.isoelectric_point`.

    Mirrors `properties._bisect_charge_zero(charge_fn, lo, hi, tol)` exactly:

    * SAME bracket [lo, hi] and tol -> the bracket width starts at (hi - lo) and
      halves every iteration regardless of data, so the loop runs the SAME
      data-independent `ceil(log2((hi - lo) / tol))` iterations (== 14 for
      [0, 14], tol 1e-3) for every row;
    * SAME branch test `(f_mid > 0) == (f_lo > 0)` to move `lo` up to `mid`,
      else move `hi` down to `mid`;
    * `_charge_raw` is bit-identical to the scalar charge (charge parity = 0.0).

    The only scalar behaviour not reproduced is the `if f_mid == 0.0: return mid`
    early-out — a midpoint landing on EXACT zero net charge — which is
    astronomically rare and, when it happens, the lockstep loop still converges
    to the same value within `tol`.

    Rows with no zero crossing in [lo, hi] (both endpoints strictly same-sign)
    or invalid rows -> NaN, matching the scalar `None`. Returns float64.
    """
    f = lambda ph: _charge_raw(sub, ph, pka_set, include_cys)  # noqa: E731

    f_lo = f(lo)
    f_hi = f(hi)
    # No zero-crossing when both endpoints share a strict sign (matches scalar).
    crossing = ~(((f_lo > 0) & (f_hi > 0)) | ((f_lo < 0) & (f_hi < 0)))

    n = sub.counts.shape[0]
    lo_a = np.full(n, lo, dtype=np.float64)
    hi_a = np.full(n, hi, dtype=np.float64)
    f_lo_a = np.asarray(f_lo, dtype=np.float64).copy()

    iters = math.ceil(math.log2((hi - lo) / tol))
    for _ in range(iters):
        mid = 0.5 * (lo_a + hi_a)
        f_mid = f(mid)
        go_lo = (f_mid > 0) == (f_lo_a > 0)  # same branch test as scalar
        lo_a = np.where(go_lo, mid, lo_a)
        f_lo_a = np.where(go_lo, f_mid, f_lo_a)
        hi_a = np.where(go_lo, hi_a, mid)

    pi = 0.5 * (lo_a + hi_a)
    pi[~(crossing & sub.valid)] = np.nan
    return pi
