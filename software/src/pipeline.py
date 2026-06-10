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

import numpy as np
import polars as pl

import vectorized as vec
from aa_tables import STANDARD_AAS
from pka_tables import IPC2_PEPTIDE, IPC2_PROTEIN
from properties import (
    INSTABILITY_MIN_LENGTH,
    effective_length,
)

log = logging.getLogger(__name__)

PH = 7.0  # All charge values computed at pH 7 (spec default).


# ---------------------------------------------------------------------------
# CID quantization
# ---------------------------------------------------------------------------
#
# Determinism contract: "quantized-equal". Every emitted float column is
# rounded at the output boundary to one digit below its display precision (and
# at least 4 dp). This serves two ends: users see no change (rounding stays
# below the `.Nf` display format), and the upcoming BioPython → numpy code
# substitution's floating-point drift is absorbed, so the workflow's content-
# addressable id stays byte-stable across a same-machine re-run.
#
# `charge_*`, `chargeShift_*`, and `pi_*` are the only outputs that depend on a
# transcendental (`10**x` via libm), so they carry true ULP-level FP variance
# today. They round to 3 dp, matching the isoelectric_point bisection tolerance
# of 1e-3 (the value's real precision is ~0.0005). The other families are
# closed-form arithmetic — bit-exact under IEEE-754 on one machine now, but
# rounded anyway so the contract holds uniformly once the math is vectorized.
#
# `eox_*`/`ered_*` are integer-valued; rounding to 0 dp is exact and never
# perturbs them.
#
# The quantization is a *boundary* concern. Internal property functions
# (`charge_at_ph`, `isoelectric_point`, etc.) keep full precision so golden-
# value tests stay sharp. Only the pipeline's emitted DataFrame is rounded.
#
# `CID_QUANTIZE_PREFIXES` / `CID_QUANTIZE_DECIMALS` remain as documentation of
# the charge/chargeShift/pi family's 3-dp contract (other tests assert them).
CID_QUANTIZE_PREFIXES = ("charge_", "chargeShift_", "pi_")
CID_QUANTIZE_DECIMALS = 3

# Per-column-family rounding at the output boundary. Keys match the TSV column
# *prefixes* the pipeline emits; the first matching prefix wins. Chosen below
# each property's display format (columns.lib.tengo) with headroom over the
# BioPython->numpy FP drift the upcoming vectorization introduces.
CID_QUANTIZE_DECIMALS_BY_PREFIX: tuple[tuple[str, int], ...] = (
    ("charge_", 3),
    ("chargeShift_", 3),
    ("pi_", 3),  # .2f display, unchanged
    ("instability_", 4),  # .2f display
    ("mw_", 4),
    ("aliphatic_", 4),  # .1f display
    ("gravy_", 5),
    ("aromaticity_", 5),  # .3f display
    ("eox_", 0),
    ("ered_", 0),  # .0f, integer-valued
)


def _decimals_for(col: str) -> int | None:
    for prefix, dp in CID_QUANTIZE_DECIMALS_BY_PREFIX:
        if col.startswith(prefix):
            return dp
    return None


def _quantize_for_cid(df: pl.DataFrame) -> pl.DataFrame:
    exprs = []
    for c in df.columns:
        dp = _decimals_for(c)
        if dp is not None and df.schema[c] == pl.Float64:
            # Round, then canonicalize signed zero to `+0.0`. A property whose
            # true value is ~0 (e.g. GRAVY of a charge-balanced chain) lands on a
            # sub-ULP residual whose SIGN depends on summation order — the scalar
            # BioPython path sums in residue order, the vectorized `counts @ KD`
            # in AA-index order, so the same numeric zero rounds to `-0.0` on one
            # and `+0.0` on the other. Both are numerically equal, but the TSV
            # writer emits different bytes (`-0.0` vs `0.0`). `-0.0 == 0.0` is
            # True in polars, so the `when` maps BOTH signed zeros to `+0.0`
            # (note: `col + 0.0` does NOT canonicalize in polars — it preserves
            # the negative-zero bit). This makes the content-addressable id
            # insensitive to FP-residual-sign drift — the same determinism
            # guarantee the rounding already gives the other digits. round(null)
            # is null and `null == 0.0` is null, so NA cells stay null/empty.
            rounded = pl.col(c).round(dp)
            exprs.append(pl.when(rounded == 0.0).then(pl.lit(0.0)).otherwise(rounded).alias(c))
    return df.with_columns(exprs) if exprs else df


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


def _na_to_null(values: np.ndarray) -> pl.Series:
    """Wrap a numpy float64 array as a Float64 polars Series with NaN -> null.

    THE #1 vectorization trap: the vectorized functions emit `np.nan` for NA
    rows, but the golden's NA cells are EMPTY. A NaN survives to TSV as the
    literal "NaN", not an empty cell. `fill_nan(None)` converts NaN -> null so
    the IO layer writes an empty cell, matching the scalar `None` -> empty path.
    """
    return pl.Series(values=values, dtype=pl.Float64).fill_nan(None)


def run_peptide(reads: pl.DataFrame) -> dict[str, Any]:
    """Compute peptide-mode outputs on the vectorized engine.

    Builds the per-residue count substrate once over the `sequence` column, then
    derives all 10 peptide properties + the AA-fraction matrix with array ops.
    charge / chargeShift / pi use the IPC 2.0 peptide pKa set with Cys included
    as ionizable (free thiol); the rest are the linear / instability funcs.

    NaN -> null is applied to every emitted float column (and the AA-fraction
    `value`) so NA rows render as empty cells, byte-matching the scalar path.
    """
    keys = reads["entity_key"].to_list()
    seqs = reads["sequence"].to_list()
    n = len(seqs)

    log.info("Computing peptide properties + AA fractions (%d sequences)", n)
    sub = vec.build_counts(seqs)

    eox, ered = vec.extinction(sub)
    prop_arrays: dict[str, np.ndarray] = {
        "charge_peptide": vec.charge_at_ph(sub, PH, IPC2_PEPTIDE, include_cys=True),
        "chargeShift_peptide": vec.charge_shift(sub, IPC2_PEPTIDE, include_cys=True),
        "gravy_peptide": vec.gravy(sub),
        "mw_peptide": vec.molecular_weight(sub),
        "pi_peptide": vec.isoelectric_point(sub, IPC2_PEPTIDE, include_cys=True),
        "eox_peptide": eox,
        "ered_peptide": ered,
        "instability_peptide": vec.instability_index(seqs),
        "aliphatic_peptide": vec.aliphatic_index(sub),
        "aromaticity_peptide": vec.aromaticity(sub),
    }
    properties = pl.DataFrame(
        {
            "entity_key": pl.Series(values=keys, dtype=pl.Utf8),
            **{c: _na_to_null(prop_arrays[c]) for c in PEPTIDE_PROPERTY_COLUMNS},
        }
    )

    # AA-fraction long frame from the (N, 20) mole-fraction matrix: 20 rows per
    # entity (one per STANDARD_AAS), value NaN -> null for invalid entities so
    # the 2-axis PColumn keeps a uniform shape across entities.
    fractions = vec.aa_fractions(sub)  # (N, 20), NaN for invalid rows
    aa_entity = [k for k in keys for _ in STANDARD_AAS]
    aa_amino = list(STANDARD_AAS) * n
    aa_value = fractions.reshape(-1)  # row-major: entity 0's 20 AAs, then entity 1's, ...
    aa_fraction = pl.DataFrame(
        {
            "entity_key": pl.Series(values=aa_entity, dtype=pl.Utf8),
            "aminoAcid": pl.Series(values=aa_amino, dtype=pl.Utf8),
            "value": _na_to_null(aa_value),
        }
    )
    # Quantize at the output boundary like the property columns. aaFraction
    # displays as .3f; round to 5 dp (one digit below display, headroom for the
    # vectorization's FP drift). round(null) == null, so NA cells stay empty.
    aa_fraction = aa_fraction.with_columns(pl.col("value").round(5))

    # R9 — flag whether any *real* peptide falls below the Instability Index
    # floor. `if s` filters None / "" so the banner does not fire on empty
    # cells (no peptide, not a short peptide); `0 < effective_length` filters
    # sequences that clean to empty (all-non-standard residues). NOTE: this is
    # the scalar `effective_length` count, NOT `sub.length`: a stop-codon cell
    # like "ACDE*FGHI" is substrate-INVALID (length 0) but still has 8 standard
    # residues, which the scalar oracle counts toward the floor. Using
    # effective_length keeps this stat bit-identical to the scalar path.
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

    Walks (chains × CDR3_PROPS), (fullChains × FULL_CHAIN_PROPS), and
    conditionally FV_PROPS. The property name tuples (CDR3_PROPS,
    FULL_CHAIN_PROPS, FV_PROPS) are the single source of truth — `run_antibody_tcr`
    populates exactly these columns in this order.
    """
    cols: list[str] = []
    for ch in plan.get("chains", []):
        cols.extend(f"{p}_{ch}_CDR3" for p in CDR3_PROPS)
    for ch in plan.get("fullChains", []):
        cols.extend(f"{p}_{ch}_VDJRegion" for p in FULL_CHAIN_PROPS)
    if plan.get("hasFv"):
        cols.extend(f"{p}_Fv" for p in FV_PROPS)
    return cols


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


def _column_or_empty(reads: pl.DataFrame, col: str) -> list[str | None]:
    """Return `col` as a python list, or a list of empty strings (one per row)
    when the column is absent from the reads. Mirrors the old per-row
    `record.get(col, "")` / `record.get(col) or ""` defence against a plan that
    names a chain/region the input TSV does not carry.
    """
    if col not in reads.columns:
        return [""] * reads.height
    return reads[col].to_list()


def _reconstruct_chain_column(reads: pl.DataFrame, chain: str) -> list[str | None]:
    """Reconstruct the full chain for every clone: FR1+CDR1+...+FR4, None if any
    region is missing (empty cell OR absent column).

    Pulls each region column once, then delegates the per-clone rule to
    `_reconstruct_chain` so the reconstruction logic has a single source of
    truth. A region column absent from the reads is materialised as empty
    strings, so the missing-region -> None rule fires the same way.
    """
    cols = {feat: _column_or_empty(reads, f"{chain}_{feat}") for feat in REQUIRED_FEATURES}
    out: list[str | None] = []
    for i in range(reads.height):
        row = {f"{chain}_{feat}": cols[feat][i] for feat in REQUIRED_FEATURES}
        out.append(_reconstruct_chain(row, chain))
    return out


def run_antibody_tcr(reads: pl.DataFrame, plan: dict[str, Any]) -> dict[str, Any]:
    """Compute antibody/TCR-mode outputs on the vectorized engine.

    * CDR3 (per `plan["chains"]`): build_counts on the `{chain}_CDR3` column ->
      charge / chargeShift / gravy with the IPC 2.0 peptide pKa set, Cys
      included (CDR3 Cys treated as free thiol). NaN where the CDR3 cell is
      empty / invalid.
    * Full chain (per `plan["fullChains"]`): reconstruct FR1..FR4 -> build_counts
      -> charge / pi (IPC 2.0 protein, Cys excluded) + gravy / mw / eox / ered /
      instability / aliphatic / aromaticity. NaN where ANY region is missing.
    * Fv (when `plan["hasFv"]`): per-chain SUMS over VH=A, VL=B — charge,
      chargeShift, eox, ered, mw additive; pi = bisection of the SUMMED charge
      function. NaN where EITHER chain is not fully reconstructed.

    Every emitted float column has NaN -> null applied so NA rows render as
    empty cells, byte-matching the scalar path. Column order is
    `_planned_output_columns(plan)`.
    """
    chains = plan.get("chains", [])
    full_chains = plan.get("fullChains", [])
    n = reads.height
    if chains:
        log.info("Computing CDR3 properties for chains %s (%d clones)", list(chains), n)
    if full_chains:
        log.info("Reconstructing full chains %s and computing full-chain properties", list(full_chains))
    if plan.get("hasFv"):
        log.info("Computing Fv properties (paired VH+VL)")

    series: dict[str, pl.Series] = {"entity_key": pl.Series(values=reads["entity_key"].to_list(), dtype=pl.Utf8)}

    # CDR3 per chain — IPC2_PEPTIDE, Cys included. Empty cell -> NA for this
    # clone (build_counts marks it invalid -> the funcs emit NaN -> null).
    for ch in chains:
        sub = vec.build_counts(_column_or_empty(reads, f"{ch}_CDR3"))
        cdr3_arrays = {
            "charge": vec.charge_at_ph(sub, PH, IPC2_PEPTIDE, include_cys=True),
            "chargeShift": vec.charge_shift(sub, IPC2_PEPTIDE, include_cys=True),
            "gravy": vec.gravy(sub),
        }
        for p in CDR3_PROPS:
            series[f"{p}_{ch}_CDR3"] = _na_to_null(cdr3_arrays[p])

    # Full chain — reconstruct then compute. IPC2_PROTEIN, Cys excluded. NA per
    # clone where any region is missing (reconstruction None -> invalid row).
    # Cache the per-chain substrates for Fv reuse below.
    chain_subs: dict[str, vec.Substrate] = {}
    chain_seqs: dict[str, list[str | None]] = {}
    for ch in full_chains:
        seqs = _reconstruct_chain_column(reads, ch)
        sub = vec.build_counts(seqs)
        chain_subs[ch] = sub
        chain_seqs[ch] = seqs
        eox, ered = vec.extinction(sub)
        full_arrays = {
            "charge": vec.charge_at_ph(sub, PH, IPC2_PROTEIN, include_cys=False),
            "pi": vec.isoelectric_point(sub, IPC2_PROTEIN, include_cys=False),
            "gravy": vec.gravy(sub),
            "mw": vec.molecular_weight(sub),
            "eox": eox,
            "ered": ered,
            "instability": vec.instability_index(seqs),
            "aliphatic": vec.aliphatic_index(sub),
            "aromaticity": vec.aromaticity(sub),
        }
        for p in FULL_CHAIN_PROPS:
            series[f"{p}_{ch}_VDJRegion"] = _na_to_null(full_arrays[p])

    # Fv — per-chain sums over VH=A, VL=B. charge/chargeShift/eox/ered/mw are
    # additive (NaN propagates if either chain invalid); pi bisects the SUMMED
    # charge function. Reuses the cached full-chain substrates so the chains are
    # reconstructed/counted once. Mirrors scalar `_compute_fv_row_from_ctx`.
    if plan.get("hasFv"):
        sub_vh = chain_subs["A"]
        sub_vl = chain_subs["B"]
        eox_vh, ered_vh = vec.extinction(sub_vh)
        eox_vl, ered_vl = vec.extinction(sub_vl)
        fv_arrays = {
            "charge": vec.charge_at_ph(sub_vh, PH, IPC2_PROTEIN, include_cys=False)
            + vec.charge_at_ph(sub_vl, PH, IPC2_PROTEIN, include_cys=False),
            "chargeShift": vec.charge_shift(sub_vh, IPC2_PROTEIN, include_cys=False)
            + vec.charge_shift(sub_vl, IPC2_PROTEIN, include_cys=False),
            "pi": vec.fv_isoelectric_point(sub_vh, sub_vl, IPC2_PROTEIN, include_cys=False),
            "eox": eox_vh + eox_vl,
            "ered": ered_vh + ered_vl,
            "mw": vec.molecular_weight(sub_vh) + vec.molecular_weight(sub_vl),
        }
        for p in FV_PROPS:
            series[f"{p}_Fv"] = _na_to_null(fv_arrays[p])

    out_cols = _planned_output_columns(plan)
    properties = pl.DataFrame({"entity_key": series["entity_key"], **{c: series[c] for c in out_cols}})
    aa_fraction = pl.DataFrame(schema={"entity_key": pl.Utf8, "aminoAcid": pl.Utf8, "value": pl.Float64})
    stats = {"medianCdr3Length": _median_cdr3_length_by_chain(reads, chains)}
    return {
        "properties": _quantize_for_cid(properties),
        "aa_fraction": aa_fraction,
        "stats": stats,
    }
