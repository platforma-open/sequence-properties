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
    """Count of standard residues. Non-standard / gap characters do not count.

    Counts in a single pass without allocating the cleaned string — saves
    O(L) intermediate memory on hot stats paths (e.g. per-chain CDR3 medians).
    """
    return sum(1 for c in seq.upper() if c in STANDARD_AA_SET)


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
# Per-sequence cached context
# ---------------------------------------------------------------------------
#
# Building one ProteinAnalysis or one IsoelectricPoint per property read repeats
# the same O(L) AA-count work N times per sequence. SequenceContext computes
# `_prepare` once, then memoises one ProteinAnalysis and one IsoelectricPoint
# per (pka_set, include_cys) — so all property reads on a single sequence share
# the same cached internals. Public top-level functions still build a one-shot
# context, keeping the existing API and tests unchanged.


class SequenceContext:
    """Per-sequence cached state. Build via `SequenceContext.from_seq(seq)` —
    returns None for invalid/empty-after-cleaning input. Methods mirror the
    module's top-level property functions but skip the per-call `_prepare` and
    reuse cached BioPython objects.
    """

    __slots__ = ("cleaned", "_pa", "_ip_cache")

    def __init__(self, cleaned: str) -> None:
        self.cleaned = cleaned
        self._pa: ProteinAnalysis | None = None
        self._ip_cache: dict[tuple[str, bool], IsoelectricPoint] = {}

    @classmethod
    def from_seq(cls, seq: str | None) -> SequenceContext | None:
        cleaned = _prepare(seq)
        return None if cleaned is None else cls(cleaned)

    @property
    def length(self) -> int:
        return len(self.cleaned)

    @property
    def protein_analysis(self) -> ProteinAnalysis:
        if self._pa is None:
            self._pa = ProteinAnalysis(self.cleaned)
        return self._pa

    def isoelectric(self, pka_set: PKaSet, include_cys: bool) -> IsoelectricPoint:
        key = (pka_set.name, include_cys)
        ip = self._ip_cache.get(key)
        if ip is None:
            ip = _ipc2_isoelectric_point(self.cleaned, pka_set, include_cys)
            self._ip_cache[key] = ip
        return ip

    def charge_at_ph(self, ph: float, pka_set: PKaSet, include_cys: bool = True) -> float:
        return self.isoelectric(pka_set, include_cys).charge_at_pH(ph)

    def charge_shift(
        self,
        pka_set: PKaSet,
        include_cys: bool = True,
        ph_from: float = 7.4,
        ph_to: float = 6.0,
    ) -> float:
        ip = self.isoelectric(pka_set, include_cys)
        return ip.charge_at_pH(ph_from) - ip.charge_at_pH(ph_to)

    def isoelectric_point(self, pka_set: PKaSet, include_cys: bool = True) -> float | None:
        return _bisect_charge_zero(self.isoelectric(pka_set, include_cys).charge_at_pH)

    def gravy(self) -> float:
        return self.protein_analysis.gravy()

    def molecular_weight(self) -> float:
        return self.protein_analysis.molecular_weight()

    def aromaticity(self) -> float:
        return self.protein_analysis.aromaticity()

    def extinction_coefficients(self) -> tuple[float, float]:
        reduced, oxidized = self.protein_analysis.molar_extinction_coefficient()
        return (float(oxidized), float(reduced))

    def instability_index(self) -> float | None:
        if self.length < INSTABILITY_MIN_LENGTH:
            return None
        return self.protein_analysis.instability_index()

    def aliphatic_index(self) -> float:
        # Aliphatic index needs only 4 of the 20 counts. `str.count` runs at
        # C speed in one pass per residue; faster than building a 20-key dict
        # via `aa_counts`, and the resulting integers are bit-identical so
        # the divide-and-multiply produces the same float.
        n = self.length
        cleaned = self.cleaned
        x_a = cleaned.count("A") / n
        x_v = cleaned.count("V") / n
        x_i = cleaned.count("I") / n
        x_l = cleaned.count("L") / n
        return 100.0 * (x_a + 2.9 * x_v + 3.9 * (x_i + x_l))

    def aa_fractions(self) -> dict[str, float]:
        pct = self.protein_analysis.amino_acids_percent
        return {aa: pct[aa] / 100.0 for aa in STANDARD_AAS}


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
    ctx = SequenceContext.from_seq(seq)
    return None if ctx is None else ctx.charge_at_ph(ph, pka_set, include_cys)


def charge_shift(
    seq: str,
    pka_set: PKaSet,
    include_cys: bool = True,
    ph_from: float = 7.4,
    ph_to: float = 6.0,
) -> float | None:
    """ΔCharge = charge(seq, ph_from) − charge(seq, ph_to). Henderson-Hasselbalch
    at both pH points using the same pKa set as `charge_at_ph`. Histidine
    dominates the 7.4 → 6.0 window (~−0.46 per His). Returns None when the
    sequence is invalid or no standard residues remain after cleaning.
    """
    ctx = SequenceContext.from_seq(seq)
    return None if ctx is None else ctx.charge_shift(pka_set, include_cys, ph_from, ph_to)


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
    ctx = SequenceContext.from_seq(seq)
    return None if ctx is None else ctx.isoelectric_point(pka_set, include_cys)


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


# Effective-length floor for the Guruprasad instability index (spec R9). Exposed
# so the pipeline can count how many rows fall below the floor without duplicating
# the threshold.
INSTABILITY_MIN_LENGTH = 10


def gravy(seq: str) -> float | None:
    """Kyte-Doolittle GRAVY: mean hydropathy across standard residues."""
    ctx = SequenceContext.from_seq(seq)
    return None if ctx is None else ctx.gravy()


def molecular_weight(seq: str) -> float | None:
    """Average mass (Da). BioPython uses the same residue-mass table we did."""
    ctx = SequenceContext.from_seq(seq)
    return None if ctx is None else ctx.molecular_weight()


def extinction_coefficients(seq: str) -> tuple[float | None, float | None]:
    """Pace et al. extinction coefficients at 280 nm. Returns
    `(ε_oxidized, ε_reduced)`.

    BioPython's `molar_extinction_coefficient` returns `(reduced, oxidized)` —
    we flip to match the spec column order (oxidized first).

    Sequences with no Tyr or Trp emit (0, 0); invalid sequences emit (None, None).
    """
    ctx = SequenceContext.from_seq(seq)
    if ctx is None:
        return (None, None)
    return ctx.extinction_coefficients()


def instability_index(seq: str) -> float | None:
    """Guruprasad instability index. Spec floor: effective length ≥ 10 (NA
    otherwise) — BioPython has no such floor, so we enforce it here.
    """
    ctx = SequenceContext.from_seq(seq)
    return None if ctx is None else ctx.instability_index()


def aliphatic_index(seq: str) -> float | None:
    """Ikai 1980 aliphatic index. AI = 100 × (X_A + 2.9·X_V + 3.9·(X_I + X_L))
    where X_aa is mole fraction over the cleaned sequence.

    Custom — BioPython has no aliphatic index helper.
    """
    ctx = SequenceContext.from_seq(seq)
    return None if ctx is None else ctx.aliphatic_index()


def aromaticity(seq: str) -> float | None:
    """Aromatic fraction (F + W + Y) / effective_length, via BioPython."""
    ctx = SequenceContext.from_seq(seq)
    return None if ctx is None else ctx.aromaticity()


def aa_fractions(seq: str) -> dict[str, float] | None:
    """Mole fraction of each of the 20 standard residues. Sums to 1.0 when at
    least one standard residue is present. Returns None when invalid.

    BioPython's `amino_acids_percent` returns values in PERCENT (sum = 100).
    Convert to fractions to match the PColumn output convention.
    """
    ctx = SequenceContext.from_seq(seq)
    return None if ctx is None else ctx.aa_fractions()


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


def fv_charge_shift(
    vh: str,
    vl: str,
    pka_set: PKaSet,
    ph_from: float = 7.4,
    ph_to: float = 6.0,
) -> float | None:
    """Fv ΔCharge = ΔCharge(VH) + ΔCharge(VL). Equivalent to
    `(charge(VH, ph_from) + charge(VL, ph_from)) − (charge(VH, ph_to) +
    charge(VL, ph_to))`. Both chains excluded-Cys per the full-chain rule.
    Returns None if either chain is invalid.
    """
    s_vh = charge_shift(vh, pka_set, include_cys=False, ph_from=ph_from, ph_to=ph_to)
    s_vl = charge_shift(vl, pka_set, include_cys=False, ph_from=ph_from, ph_to=ph_to)
    if s_vh is None or s_vl is None:
        return None
    return s_vh + s_vl


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
