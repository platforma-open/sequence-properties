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

import logging
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


def run(reads: pl.DataFrame, plan: dict[str, Any]) -> dict[str, Any]:
    """Dispatch by mode. Returns a dict with three entries:

    - `properties` (DataFrame): one row per entity, columns per the plan.
    - `aa_fraction` (DataFrame): long-format (entity_key, aminoAcid, value).
      Empty body when mode is not peptide.
    - `stats` (dict): dataset-level stats consumed by the workflow info layer
      (e.g. R11c VHH detection — median CDR-H3 length per chain;
      R9 — peptide count below the Instability Index length floor).
    """
    mode = plan["mode"]
    if mode == "peptide":
        log.info("Running peptide mode (%d entities)", reads.height)
        return run_peptide(reads)
    log.info(
        "Running antibody/TCR mode (receptor=%s, %d clones)",
        plan.get("receptor", "IG"),
        reads.height,
    )
    return run_antibody_tcr(reads, plan)


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


def run_peptide(reads: pl.DataFrame) -> dict[str, Any]:
    """Compute peptide-mode outputs.

    Builds one `SequenceContext` per sequence and reuses it for both the
    properties row and the AA-fraction rows — so each sequence is `_prepare`d
    once and BioPython objects are constructed once across all 11 reads.
    Accumulates columnar arrays directly into a dict-of-lists (one allocation
    per column, vs. one dict per row) and constructs the DataFrame from those.
    """
    keys = reads["entity_key"].to_list()
    seqs = reads["sequence"].to_list()
    n = len(seqs)

    log.info("Computing peptide properties + AA fractions (%d sequences)", n)
    prop_cols: dict[str, list[Any]] = {"entity_key": [], **{c: [] for c in PEPTIDE_PROPERTY_COLUMNS}}
    aa_entity: list[str] = []
    aa_amino: list[str] = []
    aa_value: list[float | None] = []
    for k, s in zip(keys, seqs):
        prop_cols["entity_key"].append(k)
        ctx = SequenceContext.from_seq(s)
        if ctx is None:
            for c in PEPTIDE_PROPERTY_COLUMNS:
                prop_cols[c].append(None)
            # Emit one row per std AA with NA value, so the 2-axis PColumn
            # keeps a uniform shape across entities.
            for aa in STANDARD_AAS:
                aa_entity.append(k)
                aa_amino.append(aa)
                aa_value.append(None)
        else:
            row = _compute_peptide_row_from_ctx(ctx)
            for c in PEPTIDE_PROPERTY_COLUMNS:
                prop_cols[c].append(row[c])
            fractions = ctx.aa_fractions()
            for aa in STANDARD_AAS:
                aa_entity.append(k)
                aa_amino.append(aa)
                aa_value.append(fractions[aa])
    properties = pl.DataFrame(
        prop_cols,
        schema={"entity_key": pl.Utf8, **{c: pl.Float64 for c in PEPTIDE_PROPERTY_COLUMNS}},
    )
    aa_fraction = pl.DataFrame(
        {"entity_key": aa_entity, "aminoAcid": aa_amino, "value": aa_value},
        schema={"entity_key": pl.Utf8, "aminoAcid": pl.Utf8, "value": pl.Float64},
    )

    # R9 — flag whether any *real* peptide falls below the Instability Index
    # floor. `if s` filters None / "" so the banner does not fire on empty
    # cells (which are no peptide, not a short peptide); `0 < effective_length`
    # filters sequences that clean to empty (e.g. all-non-standard residues).
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


def run_antibody_tcr(reads: pl.DataFrame, plan: dict[str, Any]) -> dict[str, Any]:
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
    columns: dict[str, list[Any]] = {"entity_key": [], **{c: [] for c in out_cols}}
    for record in reads.iter_rows(named=True):
        row = _compute_row_for(record, plan)
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
