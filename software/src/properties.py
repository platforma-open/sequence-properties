"""Per-sequence physicochemical property functions.

All functions accept a single amino-acid string (single-letter codes,
case-insensitive — sequences are upper-cased on entry). Non-standard residues
(B, Z, X, U, J, gaps `-`) are excluded from computation. Stop codon `*` and
zero-length sequences invalidate the whole sequence.

Functions return `float | None`. `None` is the convention for "this property
is NA for this sequence" — the IO layer translates that to an empty TSV cell.

Cysteine handling — three contexts, three rules, per spec:

* `peptide_mode=True`: include Cys as ionizable acid (free thiol assumption).
* `cdr3_mode=True`: include Cys as ionizable acid (CDR3 Cys treated as free
  thiol; intra-CDR3 disulfides not detected).
* `full_chain_mode=True`: exclude Cys from ionizable residue sum (assume
  disulfide-bonded; engineered free Cys not detected).

Ionizable-Cys inclusion is controlled at the call site by selecting the right
pKa-set wrapper (see `charge_at_ph` / `isoelectric_point`).
"""

from __future__ import annotations

import math
from collections.abc import Callable

from aa_tables import (
    AROMATIC_AAS,
    AVG_RESIDUE_MASS,
    EC_DISULFIDE,
    EC_TRP,
    EC_TYR,
    H2O_AVG_MASS,
    KD_SCALE,
    STANDARD_AA_SET,
    STANDARD_AAS,
)
from instability import diwv
from pka_tables import ACIDIC_AAS, BASIC_AAS, PKaSet

# ---------------------------------------------------------------------------
# Sequence cleanup
# ---------------------------------------------------------------------------


def is_invalid_sequence(seq: str | None) -> bool:
    """A sequence is invalid (NA for every property) when empty or contains a
    stop codon `*`. Non-standard residues other than `*` do not invalidate the
    whole sequence — they are filtered out residue-by-residue downstream.
    """
    if seq is None:
        return True
    if not seq:
        return True
    if "*" in seq:
        return True
    return False


def clean_sequence(seq: str) -> str:
    """Upper-case + filter to standard residues only. Returns empty string if
    every residue is non-standard. Caller must check `is_invalid_sequence`
    first to handle the stop-codon case correctly.
    """
    return "".join(c for c in seq.upper() if c in STANDARD_AA_SET)


def effective_length(seq: str) -> int:
    """Count of standard residues. Non-standard / gap characters do not count."""
    return len(clean_sequence(seq))


def aa_counts(seq: str) -> dict[str, int]:
    """Count of each of the 20 standard residues. Non-standard residues are
    silently dropped — they neither contribute to a count nor to any
    denominator."""
    counts = {aa: 0 for aa in STANDARD_AAS}
    for c in seq.upper():
        if c in STANDARD_AA_SET:
            counts[c] += 1
    return counts


def _prepare(seq: str | None) -> str | None:
    """Validate + clean. Returns the cleaned (upper-cased, ambiguity-stripped)
    sequence, or None when the input is invalid OR when nothing standard
    remains after cleaning. Centralises the common prelude every property
    function used to repeat.
    """
    if is_invalid_sequence(seq):
        return None
    cleaned = clean_sequence(seq)
    return cleaned or None


# ---------------------------------------------------------------------------
# Charge and isoelectric point
# ---------------------------------------------------------------------------


def _residue_charge(aa: str, ph: float, pka_set: PKaSet, include_cys: bool) -> float:
    """Charge contribution of a single ionizable side chain at the given pH.

    Acids (D, E, Y, C when included): −1 / (1 + 10^(pKa − pH))
    Bases (H, K, R):                  +1 / (1 + 10^(pH − pKa))

    `include_cys=False` excludes Cys from the ionizable sum. Used for full VH
    / VL chains where Cys is assumed to be in disulfide bonds.
    """
    if aa == "C" and not include_cys:
        return 0.0
    pka = pka_set.get(aa)
    if pka is None:
        return 0.0
    if aa in ACIDIC_AAS:
        return -1.0 / (1.0 + 10.0 ** (pka - ph))
    if aa in BASIC_AAS:
        return 1.0 / (1.0 + 10.0 ** (ph - pka))
    return 0.0


def charge_at_ph(seq: str, ph: float, pka_set: PKaSet, include_cys: bool = True) -> float | None:
    """Net charge of a sequence at a given pH (Henderson-Hasselbalch).

    Returns None when the sequence is invalid or when no standard residues
    remain after non-standard filtering.
    """
    cleaned = _prepare(seq)
    if cleaned is None:
        return None

    # N- and C-terminus: each free terminus contributes one ionizable group.
    # N-terminal amine — basic; C-terminal carboxyl — acidic.
    n_term = 1.0 / (1.0 + 10.0 ** (ph - pka_set.n_terminus))
    c_term = -1.0 / (1.0 + 10.0 ** (pka_set.c_terminus - ph))
    side_chain_total = sum(_residue_charge(c, ph, pka_set, include_cys) for c in cleaned)
    return n_term + c_term + side_chain_total


def _bisect_zero(
    f: Callable[[float], float],
    lo: float,
    hi: float,
    tolerance: float = 0.001,
) -> float | None:
    """Find the pH at which f(pH) = 0 by bisection over [lo, hi]. Returns
    None when both endpoints have the same sign — i.e. no zero crossing in
    range — instead of looping forever or clamping to a boundary. The same-
    sign guard is the load-bearing piece of the spec's NA rule for polybasic /
    polyacidic synthetic sequences.
    """
    f_lo, f_hi = f(lo), f(hi)
    if math.copysign(1.0, f_lo) == math.copysign(1.0, f_hi):
        return None
    while hi - lo > tolerance:
        mid = 0.5 * (lo + hi)
        f_mid = f(mid)
        if f_mid == 0.0:
            return mid
        if math.copysign(1.0, f_mid) == math.copysign(1.0, f_lo):
            lo, f_lo = mid, f_mid
        else:
            hi, f_hi = mid, f_mid
    return 0.5 * (lo + hi)


def isoelectric_point(
    seq: str,
    pka_set: PKaSet,
    include_cys: bool = True,
    tolerance: float = 0.001,
) -> float | None:
    """pI by bisection on `charge_at_ph` over [0, 14]. Returns None when the
    sequence has no zero crossing in range or after non-standard filtering
    yields an empty sequence.
    """
    cleaned = _prepare(seq)
    if cleaned is None:
        return None
    return _bisect_zero(
        lambda ph: charge_at_ph(cleaned, ph, pka_set, include_cys=include_cys) or 0.0,
        0.0,
        14.0,
        tolerance,
    )


# ---------------------------------------------------------------------------
# Bulk-composition properties
# ---------------------------------------------------------------------------


def gravy(seq: str) -> float | None:
    """Kyte-Doolittle GRAVY: mean hydropathy across standard residues. NA when
    no standard residues remain.
    """
    cleaned = _prepare(seq)
    if cleaned is None:
        return None
    return sum(KD_SCALE[c] for c in cleaned) / len(cleaned)


def molecular_weight(seq: str) -> float | None:
    """Average mass (Da) — sum of residue masses + one H₂O. NA when no
    standard residues remain.
    """
    cleaned = _prepare(seq)
    if cleaned is None:
        return None
    return sum(AVG_RESIDUE_MASS[c] for c in cleaned) + H2O_AVG_MASS


def extinction_coefficients(seq: str) -> tuple[float | None, float | None]:
    """Pace et al. extinction coefficients at 280 nm. Returns
    `(ε_oxidized, ε_reduced)`.

    Sequences with no Tyr or Trp emit 0 (not NA) — the caller can interpret
    "0 means A280 quantification not possible" downstream. The tuple is
    (None, None) only when the sequence itself is invalid (empty / stop).
    """
    if is_invalid_sequence(seq):
        return (None, None)
    counts = aa_counts(seq)
    y = counts["Y"]
    w = counts["W"]
    c = counts["C"]
    eps_red = float(y * EC_TYR + w * EC_TRP)
    eps_ox = eps_red + (c // 2) * EC_DISULFIDE
    return (eps_ox, eps_red)


def instability_index(seq: str) -> float | None:
    """Guruprasad instability index. Requires effective length ≥ 10 — emits
    NA otherwise per spec.

    The 10-residue floor uses *effective length* (post non-standard filter),
    not the raw input string. A sequence padded with X to look long is short.
    """
    cleaned = _prepare(seq)
    if cleaned is None or len(cleaned) < 10:
        return None
    n = len(cleaned)
    total = 0.0
    counted = 0
    for i in range(n - 1):
        v = diwv(cleaned[i], cleaned[i + 1])
        if v is None:
            continue  # impossible after _prepare — defensive
        total += v
        counted += 1
    if counted == 0:
        return None
    return (10.0 / n) * total


def aliphatic_index(seq: str) -> float | None:
    """Ikai 1980 aliphatic index. AI = 100 × (X_A + 2.9·X_V + 3.9·(X_I + X_L))
    where X_aa is mole fraction over the cleaned sequence.
    """
    cleaned = _prepare(seq)
    if cleaned is None:
        return None
    n = len(cleaned)
    counts = aa_counts(seq)
    x_a = counts["A"] / n
    x_v = counts["V"] / n
    x_i = counts["I"] / n
    x_l = counts["L"] / n
    return 100.0 * (x_a + 2.9 * x_v + 3.9 * (x_i + x_l))


def aromaticity(seq: str) -> float | None:
    """Aromatic fraction (F + W + Y) / effective_length."""
    cleaned = _prepare(seq)
    if cleaned is None:
        return None
    n = len(cleaned)
    counts = aa_counts(seq)
    arom = sum(counts[a] for a in AROMATIC_AAS)
    return arom / n


def aa_fractions(seq: str) -> dict[str, float] | None:
    """Mole fraction of each of the 20 standard residues. Sums to 1.0 (modulo
    floating-point) when at least one standard residue is present. Returns
    None when the sequence is invalid or has no standard residues.
    """
    if is_invalid_sequence(seq):
        return None
    counts = aa_counts(seq)
    n = sum(counts.values())
    if n == 0:
        return None
    return {aa: counts[aa] / n for aa in STANDARD_AAS}


# ---------------------------------------------------------------------------
# Fv (paired-chain) properties
# ---------------------------------------------------------------------------


def fv_charge(vh: str, vl: str, ph: float, pka_set: PKaSet) -> float | None:
    """Fv net charge at a given pH = charge(VH) + charge(VL). Both chains
    excluded-Cys per the full-chain rule. Returns None if either chain is
    invalid.
    """
    c_vh = charge_at_ph(vh, ph, pka_set, include_cys=False)
    c_vl = charge_at_ph(vl, ph, pka_set, include_cys=False)
    if c_vh is None or c_vl is None:
        return None
    return c_vh + c_vl


def fv_isoelectric_point(
    vh: str,
    vl: str,
    pka_set: PKaSet,
    tolerance: float = 0.001,
) -> float | None:
    """Fv pI = pH where charge(VH, pH) + charge(VL, pH) = 0. Bisection over
    [0, 14] using the per-chain-sum charge function (NOT a concatenated
    string, per spec). Returns None if there is no zero crossing or either
    chain is invalid.
    """
    vh_clean = _prepare(vh)
    vl_clean = _prepare(vl)
    if vh_clean is None or vl_clean is None:
        return None

    def f(ph: float) -> float:
        a = charge_at_ph(vh_clean, ph, pka_set, include_cys=False) or 0.0
        b = charge_at_ph(vl_clean, ph, pka_set, include_cys=False) or 0.0
        return a + b

    return _bisect_zero(f, 0.0, 14.0, tolerance)


def fv_extinction_coefficients(vh: str, vl: str) -> tuple[float | None, float | None]:
    """Fv ε at 280 nm — additive: ε(Fv) = ε(VH) + ε(VL). The per-chain additive
    formula gives the correct disulfide count when each chain has an even
    number of Cys; a chain with an odd count contributes floor(odd/2) bonds.
    """
    if is_invalid_sequence(vh) or is_invalid_sequence(vl):
        return (None, None)
    ox_vh, red_vh = extinction_coefficients(vh)
    ox_vl, red_vl = extinction_coefficients(vl)
    if ox_vh is None or red_vh is None or ox_vl is None or red_vl is None:
        return (None, None)
    return (ox_vh + ox_vl, red_vh + red_vl)


def fv_molecular_weight(vh: str, vl: str) -> float | None:
    """Fv molecular weight = MW(VH) + MW(VL). Each chain contributes its own
    H₂O term, matching the convention that VH and VL are independent
    polypeptides each with free termini.
    """
    a = molecular_weight(vh)
    b = molecular_weight(vl)
    if a is None or b is None:
        return None
    return a + b
