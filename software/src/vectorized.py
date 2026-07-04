"""Vectorized compute substrate. A column of sequences -> per-residue count
matrix + length + validity, matching properties.py's scalar cleaning exactly.
Single-threaded by construction (pure numpy elementwise + per-seq counting)."""

from __future__ import annotations

import math
import os
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

import numpy as np
from Bio.Data.IUPACData import protein_weights
from Bio.SeqUtils.ProtParam import ProtParamData

from aa_tables import STANDARD_AAS
from pka_tables import PKaSet

_AA_INDEX = {aa: i for i, aa in enumerate(STANDARD_AAS)}

# byte -> AA index lookup (256 entries, int8). Non-standard bytes stay at -1.
# Both upper- and lower-case map to the same index so the table reproduces the
# scalar `.upper()` + standard-AA filter in a single numpy gather. latin-1
# encoding with errors="replace" keeps exactly one byte per Python char (so
# per-char offsets stay aligned) and maps any non-latin-1 char to byte 0x3F
# ('?'), which the table leaves at -1 — i.e. dropped, like any non-standard
# residue. See build_counts for how the table is applied.
_BYTE_TO_AA = np.full(256, -1, dtype=np.int8)
for _i, _aa in enumerate(STANDARD_AAS):
    _BYTE_TO_AA[ord(_aa)] = _i
    _BYTE_TO_AA[ord(_aa.lower())] = _i

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

# Dipeptide instability weights (Guruprasad et al. 1990) as a 20x20 matrix in
# STANDARD_AAS x STANDARD_AAS order: _DIWV[i, j] == DIWV[STANDARD_AAS[i]][...[j]].
# ProteinAnalysis.instability_index() resolves the same ProtParamData.DIWV nested
# dict, so building the matrix from it (rather than retyping) keeps parity exact.
_DIWV = np.array(
    [[ProtParamData.DIWV[a][b] for b in STANDARD_AAS] for a in STANDARD_AAS],
    dtype=np.float64,
)

# Spec R9 floor for the instability index — kept in sync with
# properties.INSTABILITY_MIN_LENGTH (cleaned-sequence length below this -> NA).
_INSTABILITY_MIN_LENGTH = 10

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


def _matvec(counts: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """counts (N,20) · weights (20,) -> (N,).

    Wrapped in np.errstate because numpy's matmul SIMD kernel raises spurious
    overflow / invalid / divide-by-zero RuntimeWarnings on this matvec shape
    even though the result is exact: the flags are read from FP status on
    masked tail SIMD lanes, not from the data. The output is bit-identical to
    the unwrapped `counts @ weights` (verified, maxdiff 0.0) and parity-tested
    to 1e-16; multiply-sum / einsum silence the warning too but change the
    summation order (and the emitted bytes), so they are not used here.
    """
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        return counts @ weights


@dataclass(frozen=True)
class Substrate:
    counts: np.ndarray  # (N, 20) int64, STANDARD_AAS order
    length: np.ndarray  # (N,) int64 effective length
    valid: np.ndarray  # (N,) bool


@dataclass(frozen=True)
class _Cleaned:
    """Vectorized cleaning intermediate, shared by build_counts and
    instability_index. Encodes every candidate sequence (not None/empty/no `*`)
    as one flat AA-index buffer with row boundaries, so per-char Python loops
    are replaced by a single numpy gather over the joined bytes.

    * `cand_rows`   — (C,) int32 original row index of each candidate (into N).
    * `aa_flat`     — (T,) int8 AA index per *standard* residue, concatenated in
                      candidate then residue order (non-standard chars dropped).
    * `row_of_res`  — (T,) int32 candidate-local row id (0..C-1) for each entry of
                      aa_flat — i.e. which candidate that residue belongs to.
    * `lengths`     — (C,) standard-residue count per candidate (== len of its
                      slice of aa_flat).
    """

    n: int
    cand_rows: np.ndarray
    aa_flat: np.ndarray
    row_of_res: np.ndarray
    lengths: np.ndarray


def _clean_vectorized(seqs: list[str | None]) -> _Cleaned:
    """Clean a column of sequences in one numpy pass.

    Candidate = not None, not "", and no stop codon `*` (matches the scalar
    `is_invalid_sequence` short-circuit). Each candidate is latin-1 encoded
    (one byte per char), the bytes are concatenated into a single buffer,
    `np.frombuffer`'d once, and mapped through `_BYTE_TO_AA`; standard residues
    keep their AA index, non-standard chars (mapped to -1) are dropped. The
    candidate-local row id of each *kept* residue is built by counting kept
    residues per candidate (segment-sum of the keep mask) and repeating the
    candidate ids by those counts — so `aa_flat[k]` belongs to candidate
    `row_of_res[k]` without ever materializing a per-CHAR row-id array.
    """
    n = len(seqs)
    cand_rows: list[int] = []
    chunks: list[bytes] = []  # per-candidate latin-1 bytes, joined once below
    char_lengths: list[int] = []  # per-candidate char count (== byte count)
    for i, s in enumerate(seqs):
        if s is None or s == "" or "*" in s:
            continue
        b = s.encode("latin-1", "replace")
        cand_rows.append(i)
        char_lengths.append(len(b))
        chunks.append(b)
    # Single exact-sized join instead of an incremental `bytearray +=` loop:
    # the per-append reallocation/overallocation of a many-hundred-MB buffer was
    # a large RSS transient (and allocator churn over the millions of per-row
    # encode() chunks). `b"".join` allocates the final buffer once.
    buf = b"".join(chunks)
    del chunks

    # int32 indices throughout: cand_rows are row ids into the full N and
    # row_of_res are candidate-local row ids (0..C-1) — both are < 2**31 for any
    # real dataset (<= ~2.1B clones), so int32 is safe and halves these
    # per-residue arrays vs the numpy default int64/intp.
    cand_rows_arr = np.asarray(cand_rows, dtype=np.int32)
    if not buf:
        empty_i = np.empty(0, dtype=np.int32)
        empty_aa = np.empty(0, dtype=np.int8)
        return _Cleaned(
            n=n,
            cand_rows=cand_rows_arr,
            aa_flat=empty_aa,
            row_of_res=empty_i,
            lengths=np.zeros(len(cand_rows_arr), dtype=np.int64),
        )

    char_len_arr = np.asarray(char_lengths, dtype=np.intp)
    # AA index per character of the joined buffer; non-standard chars -> -1.
    # frombuffer is a zero-copy read-only view over `buf`; the gather copies the
    # mapped indices into a fresh array, so `buf` can be dropped right after.
    aa_per_char = _BYTE_TO_AA[np.frombuffer(buf, dtype=np.uint8)]
    del buf
    keep = aa_per_char >= 0
    aa_flat = aa_per_char[keep]
    del aa_per_char
    # Kept-residue count per candidate: segment-sum the keep mask over each
    # candidate's char slice. Every candidate has >=1 char (empty strings are
    # filtered above), so the start offsets have no zero-length segments and
    # reduceat is well defined. This is `lengths` directly and avoids both the
    # per-char row-id array (the dominant transient) and a separate bincount.
    starts = np.empty(len(char_len_arr), dtype=np.intp)
    starts[0] = 0
    np.cumsum(char_len_arr[:-1], out=starts[1:])
    lengths = np.add.reduceat(keep, starts).astype(np.int64)
    del keep
    # Candidate-local row id per KEPT residue: repeat each candidate id by its
    # kept count — produces the (T_kept,) survivor directly, no masking.
    # int32: candidate count C < 2**31 for any real dataset (see above).
    row_of_res = np.repeat(np.arange(len(char_len_arr), dtype=np.int32), lengths)
    return _Cleaned(
        n=n,
        cand_rows=cand_rows_arr,
        aa_flat=aa_flat,
        row_of_res=row_of_res,
        lengths=lengths,
    )


def counts_from_cleaned(cl: _Cleaned) -> Substrate:
    """Build the (N, 20) count substrate from an already-cleaned column.

    Counts are scattered with a single flattened `np.bincount` over the
    (candidate-local-row, AA-index) pairs, then mapped back to the full N rows.
    A candidate is valid iff at least one standard residue remained
    (lengths > 0) — reproducing the scalar `clean_sequence(...) or None` rule.

    Splitting the clean (`_clean_vectorized`) from this derivation lets a caller
    clean a sequence-set ONCE and feed the same `_Cleaned` to both
    `counts_from_cleaned` and `instability_from_cleaned`, avoiding a duplicate
    full-chain clean (a large transient on the full-chain path).
    """
    n = cl.n
    counts = np.zeros((n, 20), dtype=np.int64)
    valid = np.zeros(n, dtype=bool)
    n_cand = len(cl.cand_rows)
    if n_cand:
        if cl.aa_flat.size:
            # 2D count histogram via a single flattened bincount rather than
            # np.add.at's unbuffered scatter (~3x faster, bit-identical for
            # integer counts). Flat index = candidate-row * 20 + AA-index.
            flat = cl.row_of_res.astype(np.int64) * 20 + cl.aa_flat
            cand_counts = np.bincount(flat, minlength=n_cand * 20).reshape(n_cand, 20).astype(np.int64, copy=False)
        else:
            cand_counts = np.zeros((n_cand, 20), dtype=np.int64)
        counts[cl.cand_rows] = cand_counts
        valid[cl.cand_rows] = cl.lengths > 0
    length = counts.sum(axis=1)
    return Substrate(counts=counts, length=length, valid=valid)


def build_counts(seqs: list[str | None]) -> Substrate:
    """Per-residue (N, 20) count matrix + effective length + validity, matching
    properties.py's scalar cleaning exactly.

    Fully vectorized: candidates are cleaned to one flat AA-index buffer
    (`_clean_vectorized`), then `counts_from_cleaned` scatters the counts. Public
    entry point for callers that only need counts; the pipeline's full-chain path
    instead cleans once and shares the `_Cleaned` with `instability_from_cleaned`.
    """
    return counts_from_cleaned(_clean_vectorized(seqs))


# --------------------------------------------------------------------------- #
# Linear properties — pure array ops over the count matrix. Each returns a
# float64 array (or pair of arrays) with NaN for invalid rows, matching the
# scalar properties.py oracle within 1e-6.
# --------------------------------------------------------------------------- #


def gravy(sub: Substrate) -> np.ndarray:
    """Kyte-Doolittle GRAVY: (sum of per-residue hydropathy) / length.

    Mirrors ProteinAnalysis.gravy(), which returns total_gravy / length.
    """
    return _safe_div(_matvec(sub.counts, _KD), sub.length, sub.valid)


def molecular_weight(sub: Substrate) -> np.ndarray:
    """Average mass (Da): sum(residue masses) - (length - 1) * water.

    Mirrors Bio.SeqUtils.molecular_weight(seq, "protein"); water = 18.0153.
    """
    mw = _matvec(sub.counts, _MASS) - (sub.length - 1) * _WATER
    return _mask(mw, sub.valid & (sub.length > 0))


def aromaticity(sub: Substrate) -> np.ndarray:
    """Aromatic mole fraction (F + W + Y) / length."""
    return _safe_div(_matvec(sub.counts, _AROMATIC), sub.length, sub.valid)


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
    red = _matvec(sub.counts, _EXT_RED)
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

    Computes 10**ph ONCE and reuses it for every pKa term via a scalar factor
    (10**(ph-pk) = 10**ph * 10**(-pk); 10**(pk-ph) = 10**pk / 10**ph) instead of a
    separate 10**(ph±pk) per amino acid. `ph` is an (N,) array during the pI
    bisection, so this replaces ~7 array-wide transcendentals with one (~2.4x
    faster on the charge path). The FP drift vs the per-AA form is ~1e-15 — far
    below the 3-dp charge/pI quantization and the 1e-6 parity tolerance, so
    emitted bytes (and the CID) are unchanged.
    """
    c = sub.counts
    ten_ph = 10.0**ph  # single transcendental; ph may be scalar or (N,)
    # One N-terminus per sequence (count = 1), basic side chains K/R/H.
    pos = 1.0 / (ten_ph * (10.0 ** (-pka_set.n_terminus)) + 1.0)
    for aa in ("K", "R", "H"):
        pk = pka_set.side_chain[aa]
        pos = pos + c[:, _AA_INDEX[aa]] * (1.0 / (ten_ph * (10.0 ** (-pk)) + 1.0))

    # One C-terminus per sequence (count = 1), acidic side chains D/E/Y (+C).
    neg = 1.0 / ((10.0**pka_set.c_terminus) / ten_ph + 1.0)
    neg_aas = ["D", "E", "Y"]
    if include_cys and "C" in pka_set.side_chain:
        neg_aas.append("C")
    for aa in neg_aas:
        pk = pka_set.side_chain[aa]
        neg = neg + c[:, _AA_INDEX[aa]] * (1.0 / ((10.0**pk) / ten_ph + 1.0))

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


def _bisect_charge_zero_vec(
    charge_fn,
    valid: np.ndarray,
    lo: float = 0.0,
    hi: float = 14.0,
    tol: float = 1e-3,
) -> np.ndarray:
    """Lockstep fixed-iteration bisection of a row-wise charge function for the
    pH where net charge crosses zero, bit-identical to the scalar
    `properties._bisect_charge_zero(charge_fn, lo, hi, tol)`.

    `charge_fn(ph)` must return a finite-for-all-rows float64 array of net charge
    at the scalar `ph` (e.g. `_charge_raw` for one substrate, or a per-chain SUM
    of `_charge_raw` for the Fv path). `valid` is the per-row validity mask.

    Mirrors the scalar bisection exactly:

    * SAME bracket [lo, hi] and tol -> the bracket width starts at (hi - lo) and
      halves every iteration regardless of data, so the loop runs the SAME
      data-independent `ceil(log2((hi - lo) / tol))` iterations (== 14 for
      [0, 14], tol 1e-3) for every row;
    * SAME branch test `(f_mid > 0) == (f_lo > 0)` to move `lo` up to `mid`,
      else move `hi` down to `mid`.

    The only scalar behaviour not reproduced is the `if f_mid == 0.0: return mid`
    early-out — a midpoint landing on EXACT zero net charge — which is
    astronomically rare and, when it happens, the lockstep loop still converges
    to the same value within `tol`.

    Rows with no zero crossing in [lo, hi] (both endpoints strictly same-sign)
    or invalid rows -> NaN, matching the scalar `None`. Returns float64.
    """
    f_lo = charge_fn(lo)
    f_hi = charge_fn(hi)
    # No zero-crossing when both endpoints share a strict sign (matches scalar).
    crossing = ~(((f_lo > 0) & (f_hi > 0)) | ((f_lo < 0) & (f_hi < 0)))

    n = len(valid)
    lo_a = np.full(n, lo, dtype=np.float64)
    hi_a = np.full(n, hi, dtype=np.float64)
    f_lo_a = np.asarray(f_lo, dtype=np.float64).copy()

    iters = math.ceil(math.log2((hi - lo) / tol))
    for _ in range(iters):
        mid = 0.5 * (lo_a + hi_a)
        f_mid = charge_fn(mid)
        go_lo = (f_mid > 0) == (f_lo_a > 0)  # same branch test as scalar
        lo_a = np.where(go_lo, mid, lo_a)
        f_lo_a = np.where(go_lo, f_mid, f_lo_a)
        hi_a = np.where(go_lo, hi_a, mid)

    pi = 0.5 * (lo_a + hi_a)
    pi[~(crossing & valid)] = np.nan
    return pi


_PARALLEL_MIN_ROWS = 50_000


def _n_workers() -> int:
    """Row-parallel worker count = the allocated core count. The workflow sets
    POLARS_MAX_THREADS to its cpu() request, so this stays in lockstep with the
    quota; defaults to 1 (serial) for local / test runs."""
    try:
        return max(1, int(os.environ.get("POLARS_MAX_THREADS", "1")))
    except ValueError:
        return 1


def _bisect_rows_parallel(make_charge_fn, valid, lo=0.0, hi=14.0, tol=1e-3):
    """Row-parallel wrapper around `_bisect_charge_zero_vec`.

    `make_charge_fn(a, b)` returns a charge function bound to rows [a, b). Rows
    are split into contiguous, data-independent blocks and bisected on separate
    threads (numpy's elementwise pow / where release the GIL), then concatenated
    in row order. Each row's bisection depends only on that row's own counts, so
    the result is bit-identical to the serial path regardless of the block
    count — the byte-stability contract holds independent of worker count. Small
    inputs run serial to avoid thread-dispatch overhead.
    """
    n = len(valid)
    workers = _n_workers()
    if workers <= 1 or n < _PARALLEL_MIN_ROWS:
        return _bisect_charge_zero_vec(make_charge_fn(0, n), valid, lo=lo, hi=hi, tol=tol)
    edges = [(n * i) // workers for i in range(workers + 1)]
    blocks = [(edges[i], edges[i + 1]) for i in range(workers) if edges[i] < edges[i + 1]]

    def run_block(ab):
        a, b = ab
        return _bisect_charge_zero_vec(make_charge_fn(a, b), valid[a:b], lo=lo, hi=hi, tol=tol)

    with ThreadPoolExecutor(max_workers=len(blocks)) as ex:
        parts = list(ex.map(run_block, blocks))
    return np.concatenate(parts)


def isoelectric_point(
    sub: Substrate,
    pka_set: PKaSet,
    include_cys: bool = True,
    lo: float = 0.0,
    hi: float = 14.0,
    tol: float = 1e-3,
) -> np.ndarray:
    """Vectorized pI: lockstep fixed-iteration bisection of `_charge_raw`,
    bit-identical to the scalar `properties.isoelectric_point`.

    Row-parallel across `_n_workers()` threads (see `_bisect_rows_parallel`); the
    per-row bisection is unchanged, so the pI is bit-identical to the scalar pI
    for every valid sequence and independent of the worker count.
    """

    def make(a: int, b: int):
        block = Substrate(counts=sub.counts[a:b], length=sub.length[a:b], valid=sub.valid[a:b])
        return lambda ph: _charge_raw(block, ph, pka_set, include_cys)

    return _bisect_rows_parallel(make, sub.valid, lo=lo, hi=hi, tol=tol)


def fv_isoelectric_point(
    sub_vh: Substrate,
    sub_vl: Substrate,
    pka_set: PKaSet,
    include_cys: bool = False,
    lo: float = 0.0,
    hi: float = 14.0,
    tol: float = 1e-3,
) -> np.ndarray:
    """Vectorized Fv pI: pH where charge(VH, pH) + charge(VL, pH) = 0.

    Mirrors the scalar `properties.fv_isoelectric_point`: bisect the per-chain
    charge SUM, not the pI of a concatenated VH+VL string. The bisected function
    is
    `f(ph) = _charge_raw(vh, ph) + _charge_raw(vl, ph)` — both finite for all
    rows — and the row is valid only where BOTH chains are valid (matches the
    scalar's "None if either chain invalid").

    `sub_vh` and `sub_vl` must be aligned row-for-row (same N, same clone order).
    Returns float64, NaN where either chain is invalid or no zero crossing.
    """
    valid = sub_vh.valid & sub_vl.valid

    def make(a: int, b: int):
        vh = Substrate(counts=sub_vh.counts[a:b], length=sub_vh.length[a:b], valid=sub_vh.valid[a:b])
        vl = Substrate(counts=sub_vl.counts[a:b], length=sub_vl.length[a:b], valid=sub_vl.valid[a:b])
        return lambda ph: _charge_raw(vh, ph, pka_set, include_cys) + _charge_raw(vl, ph, pka_set, include_cys)

    return _bisect_rows_parallel(make, valid, lo=lo, hi=hi, tol=tol)


# --------------------------------------------------------------------------- #
# Instability index — the only ORDER-dependent property. The count matrix is
# insufficient: it needs the consecutive-residue PAIRS of the cleaned sequence.
# We clean each sequence to an index array (matching properties._prepare /
# clean_sequence) and gather the dipeptide weights with numpy — NO BioPython
# ProteinAnalysis object is constructed per row.
#
# Per the Guruprasad et al. 1990 method as implemented by
# Bio.SeqUtils.ProtParam.ProteinAnalysis.instability_index over the CLEANED seq:
#   score = sum(DIWV[seq[i]][seq[i+1]] for i in range(len - 1))
#   index = (10.0 / len(seq)) * score
# The scalar oracle (properties.SequenceContext.instability_index) additionally
# returns None below the spec R9 floor (cleaned length < INSTABILITY_MIN_LENGTH).
# --------------------------------------------------------------------------- #


def instability_from_cleaned(cl: _Cleaned) -> np.ndarray:
    """Guruprasad instability index from an already-cleaned column, float64, NaN
    where the scalar oracle returns None.

    Splitting the clean from this derivation lets the full-chain path reuse the
    same `_Cleaned` for both counts and instability (see `counts_from_cleaned`).
    A row is NaN when the sequence is invalid (None / empty / contains a stop
    codon `*` / nothing standard remains after cleaning) OR when the cleaned
    (standard-AA-only) length is below `_INSTABILITY_MIN_LENGTH` (= 10). Order is
    preserved: each row's value is a pure function of that row's sequence, with
    no dict/set iteration leaking into the output.
    """
    out = np.full(cl.n, np.nan, dtype=np.float64)
    n_cand = len(cl.cand_rows)
    if n_cand == 0:
        return out

    # Consecutive cleaned residues k, k+1 form a dipeptide only when they belong
    # to the SAME candidate (row boundaries in the joined buffer must not leak a
    # pair across two sequences). Mask on row_of_res, gather DIWV weights, then
    # sum per candidate with np.add.at — order-preserving, no per-row Python.
    # aa_flat stays int8 (no intp promotion): _DIWV fancy-indexing accepts an
    # int8 index, and avoiding the int64 copy saves ~8x its size on full chains
    # (the int64 promotion was a multi-hundred-MB transient). The gathered DIWV
    # weights and the per-candidate sum are unchanged, so the result is identical.
    aa = cl.aa_flat
    rows = cl.row_of_res
    if aa.size >= 2:
        same = rows[:-1] == rows[1:]
        left = aa[:-1][same]
        right = aa[1:][same]
        pair_row = rows[:-1][same]
        weights = _DIWV[left, right]
        # Weighted per-candidate sum via bincount instead of np.add.at. Both
        # accumulate weights[i] into bin pair_row[i] in array order, so the float
        # summation order — and the result — is identical, while bincount is
        # ~3x faster. Empty pair_row (no same-row adjacent pair) -> zeros.
        score = np.bincount(pair_row, weights=weights, minlength=n_cand)
    else:
        score = np.zeros(n_cand, dtype=np.float64)

    lengths = cl.lengths
    keep = lengths >= _INSTABILITY_MIN_LENGTH  # excludes length 0 ("nothing standard")
    if keep.any():
        rows_keep = cl.cand_rows[keep]
        out[rows_keep] = (10.0 / lengths[keep]) * score[keep]
    return out


def instability_index(seqs: list[str | None]) -> np.ndarray:
    """Public entry point: clean then derive the instability index. Callers that
    also need the count substrate should clean once with `_clean_vectorized` and
    feed the shared `_Cleaned` to both `counts_from_cleaned` and
    `instability_from_cleaned` instead of calling this and `build_counts`.
    """
    return instability_from_cleaned(_clean_vectorized(seqs))
