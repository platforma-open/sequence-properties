"""Vectorized compute substrate. A column of sequences -> per-residue count
matrix + length + validity, matching properties.py's scalar cleaning exactly.
Single-threaded by construction (pure numpy elementwise + per-seq counting)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from Bio.Data.IUPACData import protein_weights
from Bio.SeqUtils.ProtParam import ProtParamData

from aa_tables import STANDARD_AAS

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
