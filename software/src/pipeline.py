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
from collections import Counter
from collections.abc import Callable, Iterable, Iterator
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


def _round_and_canonicalize_zero(col: str, dp: int) -> pl.Expr:
    """Round `col` to `dp` decimals, then canonicalize signed zero to `+0.0`.

    This is the single definition of "how an emitted float is made CID-stable",
    shared by the property frame (`_quantize_for_cid`) and the aa_fraction
    `value` column. A property whose true value is ~0 (e.g. GRAVY of a
    charge-balanced chain) lands on a sub-ULP residual whose SIGN depends on
    summation order — the scalar BioPython path sums in residue order, the
    vectorized `counts @ KD` in AA-index order, so the same numeric zero rounds
    to `-0.0` on one and `+0.0` on the other. Both are numerically equal, but the
    TSV writer emits different bytes (`-0.0` vs `0.0`). `-0.0 == 0.0` is True in
    polars, so the `when` maps BOTH signed zeros to `+0.0` (note: `col + 0.0`
    does NOT canonicalize in polars — it preserves the negative-zero bit). This
    makes the content-addressable id insensitive to FP-residual-sign drift — the
    same determinism guarantee the rounding already gives the other digits.
    round(null) is null and `null == 0.0` is null, so NA cells stay null/empty.
    """
    rounded = pl.col(col).round(dp)
    return pl.when(rounded == 0.0).then(pl.lit(0.0)).otherwise(rounded).alias(col)


def _quantize_for_cid(df: pl.DataFrame) -> pl.DataFrame:
    exprs = [
        _round_and_canonicalize_zero(c, dp)
        for c in df.columns
        if (dp := _decimals_for(c)) is not None and df.schema[c] == pl.Float64
    ]
    return df.with_columns(exprs) if exprs else df


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


# Rows per compute chunk. Every property is a per-row function, so peak memory
# scales with the chunk, not the dataset: the per-residue buffers (cleaned AA
# indices, row ids, count matrices) and both output frames exist for one chunk
# at a time. The input reader bounds a chunk's sequence bytes
# (`io_layer.READ_BLOCK_BYTES`).
DEFAULT_CHUNK_ROWS = 100_000

AA_FRACTION_SCHEMA = {"entity_key": pl.Utf8, "aminoAcid": pl.Utf8, "value": pl.Float64}

# Receives each chunk's (properties, aa_fraction) frames, in write order.
Emit = Callable[[pl.DataFrame, pl.DataFrame], None]


# The whole input as one frame, or its consecutive row batches.
Reads = pl.DataFrame | Iterable[pl.DataFrame]


def _row_chunks(reads: Reads, chunk_rows: int) -> Iterator[pl.DataFrame]:
    """Yield the input rows in input order, at most `chunk_rows` rows per frame.
    A frame with zero rows yields itself, so callers still build correctly
    typed (empty) outputs."""
    if chunk_rows < 1:
        raise ValueError(f"chunk_rows must be >= 1, got {chunk_rows}")
    for frame in [reads] if isinstance(reads, pl.DataFrame) else reads:
        if frame.height == 0:
            yield frame
        for offset in range(0, frame.height, chunk_rows):
            yield frame.slice(offset, chunk_rows)


def run_streaming(
    reads: Reads, plan: dict[str, Any], emit: Emit, chunk_rows: int = DEFAULT_CHUNK_ROWS
) -> dict[str, Any]:
    """Dispatch by mode. Rows are computed `chunk_rows` at a time in input
    order; each chunk's frames go to `emit`:

    - properties: one row per entity, columns per the plan.
    - aa_fraction: long-format (entity_key, aminoAcid, value), entities in
      input order, 20 rows each in aminoAcid order. Zero rows when mode is not
      peptide.

    Returns the dataset-level stats consumed by the workflow info layer (e.g.
    R11c VHH detection — median CDR3 length per chain; R9 — peptide count below
    the Instability Index length floor). Outputs do not depend on `chunk_rows`.
    """
    mode = plan["mode"]
    if mode == "peptide":
        log.info("Running peptide mode (chunk %d rows)", chunk_rows)
        return run_peptide(reads, emit, chunk_rows)
    log.info("Running antibody/TCR mode (receptor=%s, chunk %d rows)", plan.get("receptor", "IG"), chunk_rows)
    return run_antibody_tcr(reads, plan, emit, chunk_rows)


def run(reads: Reads, plan: dict[str, Any], chunk_rows: int = DEFAULT_CHUNK_ROWS) -> dict[str, Any]:
    """`run_streaming` collected into whole frames: a dict with `properties`
    and `aa_fraction` (DataFrames in write order) and `stats`."""
    properties: list[pl.DataFrame] = []
    aa_fraction: list[pl.DataFrame] = []

    def collect(props: pl.DataFrame, aa: pl.DataFrame) -> None:
        properties.append(props)
        aa_fraction.append(aa)

    stats = run_streaming(reads, plan, collect, chunk_rows)
    return {"properties": pl.concat(properties), "aa_fraction": pl.concat(aa_fraction), "stats": stats}


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


def run_peptide(reads: Reads, emit: Emit, chunk_rows: int = DEFAULT_CHUNK_ROWS) -> dict[str, Any]:
    """Compute peptide-mode outputs on the vectorized engine, `chunk_rows` rows
    at a time (see `_peptide_chunk`), and return the stats. The R9 below-floor
    flag is OR-ed across chunks.
    """
    log.info("Computing peptide properties + AA fractions")
    has_below_floor = False
    rows = 0
    for chunk in _row_chunks(reads, chunk_rows):
        props, aa_fraction, below = _peptide_chunk(chunk)
        emit(props, aa_fraction)
        has_below_floor = has_below_floor or below
        rows += chunk.height
    log.info("Computed %d sequences", rows)
    return {
        "medianCdr3Length": {},
        "hasPeptideBelowInstabilityFloor": has_below_floor,
    }


def _aa_fraction_long(keys: pl.Series, sub: vec.Substrate) -> pl.DataFrame:
    """Long AA-fraction frame for one chunk, entities in `keys` order: 20 rows
    per entity, one per STANDARD_AAS, value NaN -> null for invalid entities so
    the 2-axis PColumn keeps a uniform shape across entities. Each entity's
    rows are in ascending aminoAcid order (STANDARD_AAS is ascending).
    """
    n = len(keys)
    fractions = vec.aa_fractions(sub)  # (n, 20), NaN for invalid rows
    frame = pl.DataFrame(
        {
            "entity_key": keys.gather(np.repeat(np.arange(n), len(STANDARD_AAS))),
            "aminoAcid": pl.Series(values=list(STANDARD_AAS) * n, dtype=pl.Utf8),
            "value": _na_to_null(fractions.reshape(-1)),  # row-major: entity 0's 20 AAs, ...
        },
        schema=AA_FRACTION_SCHEMA,
    )
    # Quantized at the output boundary through the shared CID-stability
    # helper. aaFraction displays as .3f; round to 5 dp. Fractions are >= 0,
    # so the signed-zero canonicalization is a no-op here.
    return frame.with_columns(_round_and_canonicalize_zero("value", 5))


def _peptide_chunk(reads: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame, bool]:
    """Peptide-mode outputs for one row chunk: (quantized properties, long
    AA-fraction frame, whether any peptide is below the instability floor).

    Builds the per-residue count substrate once over the `sequence` column, then
    derives all 10 peptide properties with array ops. charge / chargeShift / pi
    use the IPC 2.0 peptide pKa set with Cys included as ionizable (free thiol);
    the rest are the linear / instability funcs.

    NaN -> null is applied to every emitted float column so NA rows render as
    empty cells, byte-matching the scalar path.
    """
    key_series = reads["entity_key"].cast(pl.Utf8)
    seqs = reads["sequence"].to_list()

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
            "entity_key": key_series,
            **{c: _na_to_null(prop_arrays[c]) for c in PEPTIDE_PROPERTY_COLUMNS},
        }
    )

    # R9 — flag whether any *real* peptide falls below the Instability Index
    # floor. `if s` filters None / "" so the banner does not fire on empty
    # cells (no peptide, not a short peptide); `0 < effective_length` filters
    # sequences that clean to empty (all-non-standard residues). NOTE: this is
    # the scalar `effective_length` count, NOT `sub.length`: a stop-codon cell
    # like "ACDE*FGHI" is substrate-INVALID (length 0) but still has 8 standard
    # residues, which the scalar oracle counts toward the floor. Using
    # effective_length keeps this stat bit-identical to the scalar path.
    has_below_floor = any(0 < effective_length(s) < INSTABILITY_MIN_LENGTH for s in seqs if s)
    return _quantize_for_cid(properties), _aa_fraction_long(key_series, sub), has_below_floor


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


def _count_cdr3_lengths(counts: dict[str, Counter[int]], reads: pl.DataFrame, chains: list[str]) -> None:
    """Add this chunk's non-empty CDR3 effective lengths to `counts`, per
    chain. Effective length excludes ambiguity codes — matches the convention
    used by all property functions. A chain with no `{chain}_CDR3` column or no
    non-empty CDR3 gets no entry.
    """
    for ch in chains:
        col = f"{ch}_CDR3"
        if col not in reads.columns:
            continue
        lengths = [effective_length(s) for s in reads[col].to_list() if s]
        if lengths:
            counts.setdefault(ch, Counter()).update(lengths)


def _median_cdr3_length_by_chain(counts: dict[str, Counter[int]]) -> dict[str, float]:
    """Median CDR3 effective length per chain, from `_count_cdr3_lengths`
    counts. An even count averages the two middle lengths.
    """
    out: dict[str, float] = {}
    for ch, by_length in counts.items():
        n = by_length.total()
        # 0-based ranks of the middle value(s) in the sorted lengths.
        ranks = (n // 2,) if n % 2 == 1 else (n // 2 - 1, n // 2)
        middle: list[int] = []
        seen = 0
        for length in sorted(by_length):
            seen += by_length[length]
            while len(middle) < len(ranks) and ranks[len(middle)] < seen:
                middle.append(length)
        out[ch] = float(middle[0]) if n % 2 == 1 else 0.5 * (middle[0] + middle[1])
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

    Vectorized with polars: each region is normalised so that both null and the
    empty string map to null (the `if not v` "missing region" rule), then
    `concat_str(ignore_nulls=False)` joins the seven regions — yielding null
    whenever ANY region is null/empty (any-missing -> None) without a per-row
    Python loop. A region column absent from the reads is treated as all-empty,
    so its rows reconstruct to None, same as the per-row path.
    """
    parts = []
    for feat in REQUIRED_FEATURES:
        col = f"{chain}_{feat}"
        if col not in reads.columns:
            # Absent column -> every row missing this region -> chain is None.
            return [None] * reads.height
        # Empty string counts as missing (== None for the any-missing rule).
        parts.append(pl.when(pl.col(col) == "").then(None).otherwise(pl.col(col)))
    reconstructed = reads.select(pl.concat_str(parts, ignore_nulls=False).alias("_chain"))["_chain"]
    return reconstructed.to_list()


def run_antibody_tcr(
    reads: Reads, plan: dict[str, Any], emit: Emit, chunk_rows: int = DEFAULT_CHUNK_ROWS
) -> dict[str, Any]:
    """Compute antibody/TCR-mode outputs on the vectorized engine, `chunk_rows`
    rows at a time (see `_antibody_tcr_chunk`), and return the stats. The median
    CDR3 length stat comes from per-chain length counts added up over the
    chunks.
    """
    chains = plan.get("chains", [])
    full_chains = plan.get("fullChains", [])
    if chains:
        log.info("Computing CDR3 properties for chains %s", list(chains))
    if full_chains:
        log.info("Reconstructing full chains %s and computing full-chain properties", list(full_chains))
    if plan.get("hasFv"):
        log.info("Computing Fv properties (paired VH+VL)")

    no_aa_fraction = pl.DataFrame(schema=AA_FRACTION_SCHEMA)
    cdr3_lengths: dict[str, Counter[int]] = {}
    rows = 0
    for chunk in _row_chunks(reads, chunk_rows):
        emit(_antibody_tcr_chunk(chunk, plan), no_aa_fraction)
        _count_cdr3_lengths(cdr3_lengths, chunk, chains)
        rows += chunk.height
    log.info("Computed %d clones", rows)
    return {"medianCdr3Length": _median_cdr3_length_by_chain(cdr3_lengths)}


def _antibody_tcr_chunk(reads: pl.DataFrame, plan: dict[str, Any]) -> pl.DataFrame:
    """Quantized antibody/TCR property frame for one row chunk.

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

    series: dict[str, pl.Series] = {"entity_key": reads["entity_key"].cast(pl.Utf8)}

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
    # Cache the per-chain substrates for Fv reuse below. Each chain is cleaned
    # ONCE: the shared `_Cleaned` feeds both the count substrate and the
    # instability index, so the full chains are not cleaned twice (the clean's
    # flat AA buffer is the dominant full-chain transient). The reconstructed
    # `seqs` list and the `_Cleaned` are dropped at the end of each iteration —
    # Fv reuses only the substrates (chain_subs), never the raw sequences.
    chain_subs: dict[str, vec.Substrate] = {}
    for ch in full_chains:
        seqs = _reconstruct_chain_column(reads, ch)
        cleaned = vec._clean_vectorized(seqs)
        del seqs  # only the cleaned form is needed past this point
        sub = vec.counts_from_cleaned(cleaned)
        chain_subs[ch] = sub
        eox, ered = vec.extinction(sub)
        full_arrays = {
            "charge": vec.charge_at_ph(sub, PH, IPC2_PROTEIN, include_cys=False),
            "pi": vec.isoelectric_point(sub, IPC2_PROTEIN, include_cys=False),
            "gravy": vec.gravy(sub),
            "mw": vec.molecular_weight(sub),
            "eox": eox,
            "ered": ered,
            "instability": vec.instability_from_cleaned(cleaned),
            "aliphatic": vec.aliphatic_index(sub),
            "aromaticity": vec.aromaticity(sub),
        }
        for p in FULL_CHAIN_PROPS:
            series[f"{p}_{ch}_VDJRegion"] = _na_to_null(full_arrays[p])

    # Fv — per-chain sums over VH=A, VL=B. charge/chargeShift/eox/ered/mw are
    # additive (NaN propagates if either chain invalid); pi bisects the SUMMED
    # charge function. Reuses the cached full-chain substrates so the chains are
    # reconstructed/counted once. Mirrors the scalar Fv anchors:
    # `properties.fv_isoelectric_point` (pi) and the additive `properties.fv_charge`
    # / `properties.fv_molecular_weight` / `properties.fv_extinction_coefficients`.
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
        del eox_vh, ered_vh, eox_vl, ered_vl, sub_vh, sub_vl, fv_arrays

    # The full-chain substrates (each an (N, 20) count matrix + length + valid)
    # are the last large live arrays besides the emitted Series. Nothing past
    # this point reads them, so drop them before assembling the output frame.
    chain_subs.clear()

    out_cols = _planned_output_columns(plan)
    properties = pl.DataFrame({"entity_key": series["entity_key"], **{c: series[c] for c in out_cols}})
    return _quantize_for_cid(properties)
