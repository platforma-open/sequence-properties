"""Compute properties for every entity in the input TSV given a plan dict.

Two top-level paths:

* `run_peptide` — one row per peptide, computes all 9 scalar properties +
  AA-fraction long-format frame.

* `run_antibody_tcr` — one row per clone. Per-chain CDR3 properties when CDR3
  present; full-chain VDJ properties for chains where all 7 regions were
  observed; Fv columns when both VH and VL have full coverage and receptor
  is IG.

Output column names are the contract from `workflow/src/process.tpl.tengo`.
"""

from __future__ import annotations

import functools
import logging
import os
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from typing import Any

import polars as pl

from aa_tables import STANDARD_AAS
from pka_tables import IPC2_PEPTIDE, IPC2_PROTEIN
from properties import (
    INSTABILITY_MIN_LENGTH,
    SequenceContext,
    _bisect_charge_zero,
    effective_length,
)

log = logging.getLogger(__name__)

PH = 7.0  # All charge values computed at pH 7 (spec default).


# Parallelism. Below this row count the pool's process startup + pickle cost
# outweighs the benefit, so we stay in-process. The threshold also keeps the
# whole existing unit-test suite on the fast sequential path.
# The threshold is compared against ROW count regardless of per-row cost — an
# antibody/TCR clone does several times the work of a peptide (CDR3 × chains +
# full-chain × chains + Fv), so 2000 is intentionally more conservative for
# antibody mode than peptide mode.
_PARALLEL_MIN_ROWS = 2000

# Rows dispatched per pool task. Bounds per-task pickle/IPC overhead without
# starving workers; never overridden today but kept as a tuning knob.
_PARALLEL_CHUNKSIZE = 256


def resolve_workers(workers: int | None) -> int:
    """How many worker processes to use. Explicit arg wins (used by tests and
    by main.py once it reads the platform's CPU allocation). Falls back to the
    PL_COMPUTE_WORKERS env var, then os.cpu_count(). The RESULT never depends on
    this number — only the wall-clock does — so an over- or under-estimate is a
    speed concern, never a correctness one.
    """
    if workers is not None:
        return max(1, int(workers))
    env = os.environ.get("PL_COMPUTE_WORKERS")
    if env and env.isdigit() and int(env) > 0:
        return int(env)
    return max(1, os.cpu_count() or 1)


def _pmap(fn, items: list, workers: int, chunksize: int = _PARALLEL_CHUNKSIZE) -> list:
    """Map fn over items, preserving input order. Sequential below the
    threshold or when workers<=1; otherwise a process pool. ProcessPoolExecutor
    .map() preserves input order, so results align with items by index — the
    property the byte-stable output depends on.
    """
    if workers <= 1 or len(items) < _PARALLEL_MIN_ROWS:
        return [fn(x) for x in items]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        try:
            return list(ex.map(fn, items, chunksize=chunksize))
        except BrokenProcessPool as exc:
            raise RuntimeError(
                f"compute worker pool broke while processing {len(items)} items "
                f"with {workers} workers — a worker process was killed (most likely "
                f"out of memory). Reduce the input size or raise the step's memory."
            ) from exc


# ---------------------------------------------------------------------------
# CID quantization
# ---------------------------------------------------------------------------
#
# Only `charge_*` and `pi_*` outputs depend on a transcendental — `10**x` via
# libm — and only those carry ULP-level variance when the underlying FP path
# changes (libm patch, numpy SIMD reduction strategy, Python → numpy code
# substitution). Every other property is closed-form integer / constant
# arithmetic and bit-exact under IEEE-754 on a single machine.
#
# Rounding to 3 decimals matches the isoelectric_point bisection tolerance of
# 1e-3 — the value's true precision is already 0.0005, so rounding to 0.001
# discards only ULP noise without losing real information. Display format
# (.2f) is even coarser, so users see no change.
#
# The quantization is a *boundary* concern. Internal property functions
# (`charge_at_ph`, `isoelectric_point`, etc.) keep full precision so golden-
# value tests stay sharp. Only the pipeline's emitted DataFrame is rounded.
CID_QUANTIZE_PREFIXES = ("charge_", "chargeShift_", "pi_")
CID_QUANTIZE_DECIMALS = 3


def _quantize_for_cid(df: pl.DataFrame) -> pl.DataFrame:
    cols = [c for c in df.columns if any(c.startswith(p) for p in CID_QUANTIZE_PREFIXES)]
    if not cols:
        return df
    return df.with_columns([pl.col(c).round(CID_QUANTIZE_DECIMALS) for c in cols])


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def run(reads: pl.DataFrame, plan: dict[str, Any], workers: int | None = None) -> dict[str, Any]:
    """Dispatch by mode. `workers` controls parallelism only — output is
    identical for any value (see test_parallel_invariance). Returns a dict with
    `properties`, `aa_fraction`, and `stats` entries (unchanged contract).
    """
    n_workers = resolve_workers(workers)
    mode = plan["mode"]
    if mode == "peptide":
        log.info("Running peptide mode (%d entities, %d workers)", reads.height, n_workers)
        return run_peptide(reads, n_workers)
    log.info(
        "Running antibody/TCR mode (receptor=%s, %d clones, %d workers)",
        plan.get("receptor", "IG"),
        reads.height,
        n_workers,
    )
    return run_antibody_tcr(reads, plan, n_workers)


# ---------------------------------------------------------------------------
# Peptide mode
# ---------------------------------------------------------------------------

PEPTIDE_PROPERTY_COLUMNS = [
    "charge_peptide",
    "chargeShift_peptide",
    "gravy_peptide",
    "mw_peptide",
    "pi_peptide",
    "eox_peptide",
    "ered_peptide",
    "instability_peptide",
    "aliphatic_peptide",
    "aromaticity_peptide",
]


_NA_PEPTIDE_ROW: dict[str, float | None] = dict.fromkeys(PEPTIDE_PROPERTY_COLUMNS)


def _compute_peptide_row(seq: str) -> dict[str, float | None]:
    """All 9 scalar properties for a single peptide. Cys is included as
    ionizable (free thiol assumption) — the IPC 2.0 peptide pKa set is used.

    Uses one `SequenceContext` per sequence so `_prepare`, `ProteinAnalysis`,
    and `IsoelectricPoint(IPC2_PEPTIDE, include_cys=True)` are constructed
    exactly once and shared across all 10 property reads.
    """
    ctx = SequenceContext.from_seq(seq)
    if ctx is None:
        return dict(_NA_PEPTIDE_ROW)
    eox, ered = ctx.extinction_coefficients()
    return {
        "charge_peptide": ctx.charge_at_ph(PH, IPC2_PEPTIDE, include_cys=True),
        "chargeShift_peptide": ctx.charge_shift(IPC2_PEPTIDE, include_cys=True),
        "gravy_peptide": ctx.gravy(),
        "mw_peptide": ctx.molecular_weight(),
        "pi_peptide": ctx.isoelectric_point(IPC2_PEPTIDE, include_cys=True),
        "eox_peptide": eox,
        "ered_peptide": ered,
        "instability_peptide": ctx.instability_index(),
        "aliphatic_peptide": ctx.aliphatic_index(),
        "aromaticity_peptide": ctx.aromaticity(),
    }


def _compute_peptide_row_from_ctx(ctx: SequenceContext) -> dict[str, float | None]:
    """Variant that takes a pre-built context — used when `run_peptide` already
    constructed one to share with the AA-fraction pass.
    """
    eox, ered = ctx.extinction_coefficients()
    return {
        "charge_peptide": ctx.charge_at_ph(PH, IPC2_PEPTIDE, include_cys=True),
        "chargeShift_peptide": ctx.charge_shift(IPC2_PEPTIDE, include_cys=True),
        "gravy_peptide": ctx.gravy(),
        "mw_peptide": ctx.molecular_weight(),
        "pi_peptide": ctx.isoelectric_point(IPC2_PEPTIDE, include_cys=True),
        "eox_peptide": eox,
        "ered_peptide": ered,
        "instability_peptide": ctx.instability_index(),
        "aliphatic_peptide": ctx.aliphatic_index(),
        "aromaticity_peptide": ctx.aromaticity(),
    }


def _peptide_worker(seq: str) -> tuple[dict[str, float | None], list[float | None] | None]:
    """Picklable per-peptide unit: (properties row, 20 AA fractions in
    STANDARD_AAS order | None). One SequenceContext per sequence, shared across
    all 11 reads — same sharing the old inline loop relied on.
    """
    ctx = SequenceContext.from_seq(seq)
    if ctx is None:
        return (dict(_NA_PEPTIDE_ROW), None)
    props = _compute_peptide_row_from_ctx(ctx)
    fr = ctx.aa_fractions()
    return (props, [fr[aa] for aa in STANDARD_AAS])


def run_peptide(reads: pl.DataFrame, workers: int = 1) -> dict[str, Any]:
    """Compute peptide-mode outputs. Per-sequence work runs through `_pmap`
    (sequential or pooled); results are reassembled in input order so the
    serialized output is byte-identical regardless of worker count.
    """
    keys = reads["entity_key"].to_list()
    seqs = reads["sequence"].to_list()
    n = len(seqs)

    log.info("Computing peptide properties + AA fractions (%d sequences)", n)
    results = _pmap(_peptide_worker, seqs, workers)

    prop_cols: dict[str, list[Any]] = {"entity_key": list(keys), **{c: [] for c in PEPTIDE_PROPERTY_COLUMNS}}
    aa_entity: list[str] = []
    aa_amino: list[str] = []
    aa_value: list[float | None] = []
    for k, (props, fractions) in zip(keys, results):
        for c in PEPTIDE_PROPERTY_COLUMNS:
            prop_cols[c].append(props[c])
        if fractions is None:
            for aa in STANDARD_AAS:
                aa_entity.append(k)
                aa_amino.append(aa)
                aa_value.append(None)
        else:
            for aa, val in zip(STANDARD_AAS, fractions):
                aa_entity.append(k)
                aa_amino.append(aa)
                aa_value.append(val)
    properties = pl.DataFrame(
        prop_cols,
        schema={"entity_key": pl.Utf8, **{c: pl.Float64 for c in PEPTIDE_PROPERTY_COLUMNS}},
    )
    aa_fraction = pl.DataFrame(
        {"entity_key": aa_entity, "aminoAcid": aa_amino, "value": aa_value},
        schema={"entity_key": pl.Utf8, "aminoAcid": pl.Utf8, "value": pl.Float64},
    )

    has_below_floor = any(0 < effective_length(s) < INSTABILITY_MIN_LENGTH for s in seqs if s)
    stats = {
        "medianCdr3Length": {},
        "hasPeptideBelowInstabilityFloor": has_below_floor,
    }

    return {
        "properties": _quantize_for_cid(properties),
        "aa_fraction": aa_fraction,
        "stats": stats,
    }


# ---------------------------------------------------------------------------
# Antibody / TCR mode
# ---------------------------------------------------------------------------

CDR3_PROPS = ("charge", "chargeShift", "gravy")
FULL_CHAIN_PROPS = (
    "charge",
    "pi",
    "gravy",
    "mw",
    "eox",
    "ered",
    "instability",
    "aliphatic",
    "aromaticity",
)
FV_PROPS = ("charge", "chargeShift", "pi", "eox", "ered", "mw")

REQUIRED_FEATURES = ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")


_NA_CDR3_ROW: dict[str, float | None] = dict.fromkeys(CDR3_PROPS)
_NA_FULL_CHAIN_ROW: dict[str, float | None] = dict.fromkeys(FULL_CHAIN_PROPS)
_NA_FV_ROW: dict[str, float | None] = dict.fromkeys(FV_PROPS)


def _compute_cdr3_row(cdr3: str) -> dict[str, float | None]:
    """CDR3 charge, ΔCharge, and GRAVY. CDR3 uses the IPC 2.0 peptide pKa set
    with Cys included as ionizable (per spec — CDR3 Cys treated as free thiol).
    """
    ctx = SequenceContext.from_seq(cdr3)
    if ctx is None:
        return dict(_NA_CDR3_ROW)
    return {
        "charge": ctx.charge_at_ph(PH, IPC2_PEPTIDE, include_cys=True),
        "chargeShift": ctx.charge_shift(IPC2_PEPTIDE, include_cys=True),
        "gravy": ctx.gravy(),
    }


def _compute_full_chain_row(chain_seq: str) -> dict[str, float | None]:
    """Full-chain (VH / VL etc.) — protein pKa set, Cys excluded from
    ionisation (assumed disulfide-bonded). One context, one ProteinAnalysis,
    one IsoelectricPoint shared across all 9 reads.
    """
    return _compute_full_chain_row_from_ctx(SequenceContext.from_seq(chain_seq))


def _compute_full_chain_row_from_ctx(ctx: SequenceContext | None) -> dict[str, float | None]:
    if ctx is None:
        return dict(_NA_FULL_CHAIN_ROW)
    eox, ered = ctx.extinction_coefficients()
    return {
        "charge": ctx.charge_at_ph(PH, IPC2_PROTEIN, include_cys=False),
        "pi": ctx.isoelectric_point(IPC2_PROTEIN, include_cys=False),
        "gravy": ctx.gravy(),
        "mw": ctx.molecular_weight(),
        "eox": eox,
        "ered": ered,
        "instability": ctx.instability_index(),
        "aliphatic": ctx.aliphatic_index(),
        "aromaticity": ctx.aromaticity(),
    }


def _compute_fv_row(vh: str, vl: str) -> dict[str, float | None]:
    """Fv columns — IPC 2.0 protein set, Cys-excluded. pI uses the per-chain
    sum of charge functions (NOT a concatenated string), per spec. Fv
    ΔCharge = ΔCharge(VH) + ΔCharge(VL).

    Builds one context per chain so the chain-level full-chain pass and the
    Fv pass share their `IsoelectricPoint(IPC2_PROTEIN, include_cys=False)` —
    the same IP serves both `charge_at_ph(7.0)` and the bisection here.
    """
    return _compute_fv_row_from_ctx(SequenceContext.from_seq(vh), SequenceContext.from_seq(vl))


def _compute_fv_row_from_ctx(
    vh_ctx: SequenceContext | None,
    vl_ctx: SequenceContext | None,
) -> dict[str, float | None]:
    if vh_ctx is None or vl_ctx is None:
        return dict(_NA_FV_ROW)
    ox_vh, red_vh = vh_ctx.extinction_coefficients()
    ox_vl, red_vl = vl_ctx.extinction_coefficients()
    fn_vh = vh_ctx.isoelectric(IPC2_PROTEIN, include_cys=False).charge_at_pH
    fn_vl = vl_ctx.isoelectric(IPC2_PROTEIN, include_cys=False).charge_at_pH
    return {
        "charge": fn_vh(PH) + fn_vl(PH),
        "chargeShift": (fn_vh(7.4) - fn_vh(6.0)) + (fn_vl(7.4) - fn_vl(6.0)),
        "pi": _bisect_charge_zero(lambda ph: fn_vh(ph) + fn_vl(ph)),
        "eox": ox_vh + ox_vl,
        "ered": red_vh + red_vl,
        "mw": vh_ctx.molecular_weight() + vl_ctx.molecular_weight(),
    }


def _reconstruct_chain(row: dict[str, str], chain: str) -> str | None:
    """Concatenate FR1+CDR1+FR2+CDR2+FR3+CDR3+FR4. Returns None if any
    region is missing (empty string in input).
    """
    parts = []
    for feat in REQUIRED_FEATURES:
        col = f"{chain}_{feat}"
        v = row.get(col, "")
        if not v:
            return None
        parts.append(v)
    return "".join(parts)


def _planned_output_columns(plan: dict[str, Any]) -> list[str]:
    """Output column order — matches process.tpl.tengo's xsv import expectations.

    This and `_compute_row_for` are sibling functions: the column list and
    per-row population both walk (chains × CDR3_PROPS), (fullChains × FULL_CHAIN_PROPS),
    and conditionally FV_PROPS. Property name tuples (CDR3_PROPS, FULL_CHAIN_PROPS,
    FV_PROPS) are the single source of truth — both functions consume them.
    """
    cols: list[str] = []
    for ch in plan.get("chains", []):
        cols.extend(f"{p}_{ch}_CDR3" for p in CDR3_PROPS)
    for ch in plan.get("fullChains", []):
        cols.extend(f"{p}_{ch}_VDJRegion" for p in FULL_CHAIN_PROPS)
    if plan.get("hasFv"):
        cols.extend(f"{p}_Fv" for p in FV_PROPS)
    return cols


def _compute_row_for(record: dict[str, Any], plan: dict[str, Any]) -> dict[str, Any]:
    """Build the output row for one input record. The set of populated
    columns matches `_planned_output_columns(plan)` exactly — both are
    driven by the same plan keys and shared property tuples.

    For each chain, the reconstructed full-chain context is held and reused
    by the Fv pass below, so VH/VL `IsoelectricPoint(IPC2_PROTEIN, False)` is
    constructed once per clone instead of twice.
    """
    out: dict[str, Any] = {"entity_key": record["entity_key"]}

    # CDR3 per chain — empty cell ⇒ NA for this clone, not for the column.
    for ch in plan.get("chains", []):
        cdr3_props = _compute_cdr3_row(record.get(f"{ch}_CDR3") or "")
        for p in CDR3_PROPS:
            out[f"{p}_{ch}_CDR3"] = cdr3_props[p]

    # Full chain — reconstruct then compute. NA per-clone if any of the
    # seven regions is empty for that clone. Cache contexts for Fv reuse.
    chain_ctx: dict[str, SequenceContext | None] = {}
    for ch in plan.get("fullChains", []):
        reconstructed = _reconstruct_chain(record, ch)
        ctx = SequenceContext.from_seq(reconstructed) if reconstructed is not None else None
        chain_ctx[ch] = ctx
        full_props = _compute_full_chain_row_from_ctx(ctx)
        for p in FULL_CHAIN_PROPS:
            out[f"{p}_{ch}_VDJRegion"] = full_props[p]

    # Fv — only when both VH and VL fully reconstructed for this clone.
    # Reuses the per-chain contexts from the full-chain pass above so the
    # IPC2_PROTEIN/include_cys=False IsoelectricPoint is shared.
    if plan.get("hasFv"):
        fv = _compute_fv_row_from_ctx(chain_ctx.get("A"), chain_ctx.get("B"))
        for p in FV_PROPS:
            out[f"{p}_Fv"] = fv[p]

    return out


def _median_cdr3_length_by_chain(reads: pl.DataFrame, chains: list[str]) -> dict[str, float]:
    """Median effective length of CDR3 sequences per chain.

    Only chains with at least one non-empty CDR3 in the dataset appear in the
    result. Effective length excludes ambiguity codes — matches the convention
    used by all property functions.
    """
    out: dict[str, float] = {}
    for ch in chains:
        col = f"{ch}_CDR3"
        if col not in reads.columns:
            continue
        lengths = [effective_length(s) for s in reads[col].to_list() if s]
        if not lengths:
            continue
        lengths.sort()
        n = len(lengths)
        mid = n // 2
        out[ch] = float(lengths[mid]) if n % 2 == 1 else 0.5 * (lengths[mid - 1] + lengths[mid])
    return out


def _antibody_worker(record: dict, plan: dict[str, Any]) -> dict[str, Any]:
    """Picklable per-clone unit. `plan` is bound per-call via functools.partial;
    it is a plain dict and pickles cleanly. Delegates to the existing
    _compute_row_for so the computation is unchanged.
    """
    return _compute_row_for(record, plan)


def run_antibody_tcr(reads: pl.DataFrame, plan: dict[str, Any], workers: int = 1) -> dict[str, Any]:
    chains = plan.get("chains", [])
    full_chains = plan.get("fullChains", [])
    n = reads.height
    if chains:
        log.info("Computing CDR3 properties for chains %s (%d clones)", list(chains), n)
    if full_chains:
        log.info("Reconstructing full chains %s and computing full-chain properties", list(full_chains))
    if plan.get("hasFv"):
        log.info("Computing Fv properties (paired VH+VL)")

    out_cols = _planned_output_columns(plan)
    records = list(reads.iter_rows(named=True))
    worker = functools.partial(_antibody_worker, plan=plan)
    rows = _pmap(worker, records, workers)

    columns: dict[str, list[Any]] = {"entity_key": [], **{c: [] for c in out_cols}}
    for row in rows:
        columns["entity_key"].append(row["entity_key"])
        for c in out_cols:
            columns[c].append(row[c])
    schema = {"entity_key": pl.Utf8, **{c: pl.Float64 for c in out_cols}}
    properties = pl.DataFrame(columns, schema=schema)
    aa_fraction = pl.DataFrame(schema={"entity_key": pl.Utf8, "aminoAcid": pl.Utf8, "value": pl.Float64})
    stats = {"medianCdr3Length": _median_cdr3_length_by_chain(reads, chains)}
    return {
        "properties": _quantize_for_cid(properties),
        "aa_fraction": aa_fraction,
        "stats": stats,
    }
