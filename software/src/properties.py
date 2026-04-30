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

Implementation: BioPython's `Bio.SeqUtils.ProtParam.ProteinAnalysis` and
`Bio.SeqUtils.IsoelectricPoint.IsoelectricPoint` provide GRAVY, MW,
aromaticity, instability index, extinction coefficient, charge, and pI per
spec M1 strategy. The IPC 2.0 pKa values are injected onto the IsoelectricPoint
instance after construction (`pos_pKs` / `neg_pKs`), overriding BioPython's
default Bjellqvist values per spec direction.
"""

from __future__ import annotations

from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from aa_tables import STANDARD_AA_SET, STANDARD_AAS
from pka_tables import PKaSet

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

    BioPython's ProtParam/IsoelectricPoint reject any non-standard residue
    with a ValueError, so cleaning to standard-AAs-only is a precondition for
    every BioPython call below.
    """
    if is_invalid_sequence(seq):
        return None
    cleaned = clean_sequence(seq)
    return cleaned or None


# ---------------------------------------------------------------------------
# Charge and isoelectric point
# ---------------------------------------------------------------------------


def _ipc2_isoelectric_point(seq: str, pka_set: PKaSet, include_cys: bool) -> IsoelectricPoint:
    """BioPython IsoelectricPoint with IPC 2.0 pKa overrides on the instance.

    Override `pos_pKs` / `neg_pKs` after construction; BioPython reads these
    dicts at `charge_at_pH` time, so per-instance overrides take effect
    without touching module globals.

    `include_cys=False` omits Cys from `neg_pKs` for the full-chain rule
    (Cys assumed disulfide-bonded).
    """
    ip = IsoelectricPoint(seq)
    ip.pos_pKs = {
        "Nterm": pka_set.n_terminus,
        "K": pka_set.side_chain["K"],
        "R": pka_set.side_chain["R"],
        "H": pka_set.side_chain["H"],
    }
    ip.neg_pKs = {
        "Cterm": pka_set.c_terminus,
        "D": pka_set.side_chain["D"],
        "E": pka_set.side_chain["E"],
        "Y": pka_set.side_chain["Y"],
    }
    if include_cys and "C" in pka_set.side_chain:
        ip.neg_pKs["C"] = pka_set.side_chain["C"]
    return ip


def charge_at_ph(seq: str, ph: float, pka_set: PKaSet, include_cys: bool = True) -> float | None:
    """Net charge of a sequence at a given pH (Henderson-Hasselbalch via
    BioPython's IsoelectricPoint, with IPC 2.0 pKa overrides).

    Returns None when the sequence is invalid or when no standard residues
    remain after non-standard filtering.
    """
    cleaned = _prepare(seq)
    if cleaned is None:
        return None
    return _ipc2_isoelectric_point(cleaned, pka_set, include_cys).charge_at_pH(ph)


def isoelectric_point(
    seq: str,
    pka_set: PKaSet,
    include_cys: bool = True,
) -> float | None:
    """pI as the pH where BioPython's IPC 2.0-overridden charge function = 0.

    BioPython's `IsoelectricPoint.pi()` brackets [4.05, 12.0] internally —
    too narrow for polyacidic / polybasic synthetic sequences whose true pI
    sits outside that range. Bisect [0, 14] locally using BioPython's
    `charge_at_pH`, keeping the formula and pKa BioPython's.

    Returns None when no zero crossing exists in [0, 14], or after
    non-standard filtering yields an empty sequence.
    """
    cleaned = _prepare(seq)
    if cleaned is None:
        return None
    ip = _ipc2_isoelectric_point(cleaned, pka_set, include_cys)
    return _bisect_charge_zero(ip.charge_at_pH)


def _bisect_charge_zero(charge_fn, lo: float = 0.0, hi: float = 14.0, tol: float = 0.001) -> float | None:
    """Bisect for charge_fn(pH) = 0 over [lo, hi]. Returns None when both
    endpoints have the same sign (no zero crossing in range).
    """
    f_lo = charge_fn(lo)
    f_hi = charge_fn(hi)
    if (f_lo > 0 and f_hi > 0) or (f_lo < 0 and f_hi < 0):
        return None
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        f_mid = charge_fn(mid)
        if f_mid == 0.0:
            return mid
        if (f_mid > 0) == (f_lo > 0):
            lo, f_lo = mid, f_mid
        else:
            hi, f_hi = mid, f_mid
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------------------
# Bulk-composition properties (BioPython ProtParam)
# ---------------------------------------------------------------------------


def _protein_analysis(seq: str) -> ProteinAnalysis | None:
    """Cleaned BioPython ProteinAnalysis object, or None when invalid.

    BioPython rejects non-standard residues — clean first.
    """
    cleaned = _prepare(seq)
    if cleaned is None:
        return None
    return ProteinAnalysis(cleaned)


def gravy(seq: str) -> float | None:
    """Kyte-Doolittle GRAVY: mean hydropathy across standard residues."""
    p = _protein_analysis(seq)
    return None if p is None else p.gravy()


def molecular_weight(seq: str) -> float | None:
    """Average mass (Da). BioPython uses the same residue-mass table we did."""
    p = _protein_analysis(seq)
    return None if p is None else p.molecular_weight()


def extinction_coefficients(seq: str) -> tuple[float | None, float | None]:
    """Pace et al. extinction coefficients at 280 nm. Returns
    `(ε_oxidized, ε_reduced)`.

    BioPython's `molar_extinction_coefficient` returns `(reduced, oxidized)` —
    we flip to match the spec column order (oxidized first).

    Sequences with no Tyr or Trp emit (0, 0); invalid sequences emit (None, None).
    """
    p = _protein_analysis(seq)
    if p is None:
        return (None, None)
    reduced, oxidized = p.molar_extinction_coefficient()
    return (float(oxidized), float(reduced))


def instability_index(seq: str) -> float | None:
    """Guruprasad instability index. Spec floor: effective length ≥ 10 (NA
    otherwise) — BioPython has no such floor, so we enforce it here.
    """
    cleaned = _prepare(seq)
    if cleaned is None or len(cleaned) < 10:
        return None
    return ProteinAnalysis(cleaned).instability_index()


def aliphatic_index(seq: str) -> float | None:
    """Ikai 1980 aliphatic index. AI = 100 × (X_A + 2.9·X_V + 3.9·(X_I + X_L))
    where X_aa is mole fraction over the cleaned sequence.

    Custom — BioPython has no aliphatic index helper.
    """
    cleaned = _prepare(seq)
    if cleaned is None:
        return None
    n = len(cleaned)
    counts = aa_counts(cleaned)
    x_a = counts["A"] / n
    x_v = counts["V"] / n
    x_i = counts["I"] / n
    x_l = counts["L"] / n
    return 100.0 * (x_a + 2.9 * x_v + 3.9 * (x_i + x_l))


def aromaticity(seq: str) -> float | None:
    """Aromatic fraction (F + W + Y) / effective_length, via BioPython."""
    p = _protein_analysis(seq)
    return None if p is None else p.aromaticity()


def aa_fractions(seq: str) -> dict[str, float] | None:
    """Mole fraction of each of the 20 standard residues. Sums to 1.0 when at
    least one standard residue is present. Returns None when invalid.

    BioPython's `amino_acids_percent` returns values in PERCENT (sum = 100).
    Convert to fractions to match the PColumn output convention.
    """
    p = _protein_analysis(seq)
    if p is None:
        return None
    pct = p.amino_acids_percent
    return {aa: pct[aa] / 100.0 for aa in STANDARD_AAS}


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
) -> float | None:
    """Fv pI = pH where charge(VH, pH) + charge(VL, pH) = 0.

    Spec: bisect the per-chain charge SUM, not the pI of a concatenated VH+VL
    string — concatenation would drop one terminus pair and add a fake peptide
    bond. Bisection runs locally with BioPython's per-chain charge functions.
    """
    vh_clean = _prepare(vh)
    vl_clean = _prepare(vl)
    if vh_clean is None or vl_clean is None:
        return None

    ip_vh = _ipc2_isoelectric_point(vh_clean, pka_set, include_cys=False)
    ip_vl = _ipc2_isoelectric_point(vl_clean, pka_set, include_cys=False)
    return _bisect_charge_zero(lambda ph: ip_vh.charge_at_pH(ph) + ip_vl.charge_at_pH(ph))


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
