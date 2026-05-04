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
    aa_fractions,
    aliphatic_index,
    aromaticity,
    charge_at_ph,
    effective_length,
    extinction_coefficients,
    fv_charge,
    fv_extinction_coefficients,
    fv_isoelectric_point,
    fv_molecular_weight,
    gravy,
    instability_index,
    isoelectric_point,
    molecular_weight,
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
CID_QUANTIZE_PREFIXES = ("charge_", "pi_")
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
      (e.g. R11c VHH detection — median CDR-H3 length per chain).
    """
    mode = plan["mode"]
    if mode == "peptide":
        log.info("Running peptide mode (%d entities)", reads.height)
        out = run_peptide(reads)
        out["stats"] = {"medianCdr3Length": {}}
        return out
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
    "gravy_peptide",
    "mw_peptide",
    "pi_peptide",
    "eox_peptide",
    "ered_peptide",
    "instability_peptide",
    "aliphatic_peptide",
    "aromaticity_peptide",
]


def _compute_peptide_row(seq: str) -> dict[str, float | None]:
    """All 9 scalar properties for a single peptide. Cys is included as
    ionizable (free thiol assumption) — the IPC 2.0 peptide pKa set is used.
    """
    eox, ered = extinction_coefficients(seq)
    return {
        "charge_peptide": charge_at_ph(seq, PH, IPC2_PEPTIDE, include_cys=True),
        "gravy_peptide": gravy(seq),
        "mw_peptide": molecular_weight(seq),
        "pi_peptide": isoelectric_point(seq, IPC2_PEPTIDE, include_cys=True),
        "eox_peptide": eox,
        "ered_peptide": ered,
        "instability_peptide": instability_index(seq),
        "aliphatic_peptide": aliphatic_index(seq),
        "aromaticity_peptide": aromaticity(seq),
    }


def run_peptide(reads: pl.DataFrame) -> dict[str, pl.DataFrame]:
    """Compute peptide-mode outputs."""
    keys = reads["entity_key"].to_list()
    seqs = reads["sequence"].to_list()

    log.info("Computing peptide scalar properties (%d sequences)", len(seqs))
    rows = [{"entity_key": k, **_compute_peptide_row(s)} for k, s in zip(keys, seqs)]
    properties = pl.DataFrame(rows, schema={"entity_key": pl.Utf8, **{c: pl.Float64 for c in PEPTIDE_PROPERTY_COLUMNS}})

    log.info("Computing AA fractions (%d sequences)", len(seqs))
    aa_rows: list[dict[str, Any]] = []
    for k, s in zip(keys, seqs):
        fractions = aa_fractions(s)
        if fractions is None:
            # Emit one row per std AA with NA value, so the 2-axis PColumn
            # keeps a uniform shape across entities.
            for aa in STANDARD_AAS:
                aa_rows.append({"entity_key": k, "aminoAcid": aa, "value": None})
        else:
            for aa in STANDARD_AAS:
                aa_rows.append({"entity_key": k, "aminoAcid": aa, "value": fractions[aa]})
    aa_fraction = pl.DataFrame(
        aa_rows,
        schema={"entity_key": pl.Utf8, "aminoAcid": pl.Utf8, "value": pl.Float64},
    )

    return {"properties": _quantize_for_cid(properties), "aa_fraction": aa_fraction}


# ---------------------------------------------------------------------------
# Antibody / TCR mode
# ---------------------------------------------------------------------------

CDR3_PROPS = ("charge", "gravy")
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
FV_PROPS = ("charge", "pi", "eox", "ered", "mw")

REQUIRED_FEATURES = ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")


def _compute_cdr3_row(cdr3: str) -> dict[str, float | None]:
    """CDR3 charge + GRAVY. CDR3 uses the IPC 2.0 peptide pKa set with Cys
    included as ionizable (per spec — CDR3 Cys treated as free thiol).
    """
    return {
        "charge": charge_at_ph(cdr3, PH, IPC2_PEPTIDE, include_cys=True),
        "gravy": gravy(cdr3),
    }


def _compute_full_chain_row(chain_seq: str) -> dict[str, float | None]:
    """Full-chain (VH / VL etc.) — protein pKa set, Cys excluded from
    ionisation (assumed disulfide-bonded).
    """
    eox, ered = extinction_coefficients(chain_seq)
    return {
        "charge": charge_at_ph(chain_seq, PH, IPC2_PROTEIN, include_cys=False),
        "pi": isoelectric_point(chain_seq, IPC2_PROTEIN, include_cys=False),
        "gravy": gravy(chain_seq),
        "mw": molecular_weight(chain_seq),
        "eox": eox,
        "ered": ered,
        "instability": instability_index(chain_seq),
        "aliphatic": aliphatic_index(chain_seq),
        "aromaticity": aromaticity(chain_seq),
    }


def _compute_fv_row(vh: str, vl: str) -> dict[str, float | None]:
    """Fv columns — IPC 2.0 protein set, Cys-excluded. pI uses the per-chain
    sum of charge functions (NOT a concatenated string), per spec.
    """
    eox, ered = fv_extinction_coefficients(vh, vl)
    return {
        "charge": fv_charge(vh, vl, PH, IPC2_PROTEIN),
        "pi": fv_isoelectric_point(vh, vl, IPC2_PROTEIN),
        "eox": eox,
        "ered": ered,
        "mw": fv_molecular_weight(vh, vl),
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
    """
    out: dict[str, Any] = {"entity_key": record["entity_key"]}

    # CDR3 per chain — empty cell ⇒ NA for this clone, not for the column.
    for ch in plan.get("chains", []):
        cdr3 = record.get(f"{ch}_CDR3") or ""
        cdr3_props = _compute_cdr3_row(cdr3) if cdr3 else dict.fromkeys(CDR3_PROPS)
        for p in CDR3_PROPS:
            out[f"{p}_{ch}_CDR3"] = cdr3_props[p]

    # Full chain — reconstruct then compute. NA per-clone if any of the
    # seven regions is empty for that clone.
    reconstructed: dict[str, str | None] = {}
    for ch in plan.get("fullChains", []):
        reconstructed[ch] = _reconstruct_chain(record, ch)
        full_props = (
            _compute_full_chain_row(reconstructed[ch])
            if reconstructed[ch] is not None
            else dict.fromkeys(FULL_CHAIN_PROPS)
        )
        for p in FULL_CHAIN_PROPS:
            out[f"{p}_{ch}_VDJRegion"] = full_props[p]

    # Fv — only when both VH and VL fully reconstructed for this clone.
    if plan.get("hasFv"):
        vh = reconstructed.get("A")
        vl = reconstructed.get("B")
        fv = _compute_fv_row(vh, vl) if vh and vl else dict.fromkeys(FV_PROPS)
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
    rows = [_compute_row_for(record, plan) for record in reads.iter_rows(named=True)]
    out_cols = _planned_output_columns(plan)
    schema = {"entity_key": pl.Utf8, **{c: pl.Float64 for c in out_cols}}
    properties = pl.DataFrame(rows, schema=schema)
    aa_fraction = pl.DataFrame(schema={"entity_key": pl.Utf8, "aminoAcid": pl.Utf8, "value": pl.Float64})
    stats = {"medianCdr3Length": _median_cdr3_length_by_chain(reads, chains)}
    return {
        "properties": _quantize_for_cid(properties),
        "aa_fraction": aa_fraction,
        "stats": stats,
    }
