---
'@platforma-open/milaboratories.sequence-properties.software': patch
'@platforma-open/milaboratories.sequence-properties.workflow': patch
---

Parallelize the compute step and size its resources to the input, keeping output byte-identical.

Memory now scales with the input: `memFormula(clamp(max(lineCount * perRow, size * 64), 2GiB, 64GiB))` replaces the flat `mem("12GiB")`. The two terms cover two cost sources — a row-scaled term for per-clone structures (count matrix, output columns, aa_fraction; measured ~3 KiB/row peptide, ~4.3 KiB/row antibody), and a residue-scaled term (`size` ≈ total residues) for the O(residues) cleaning transients. The size term keeps long-sequence/low-row inputs (amplicon variants) off the 2 GiB floor, where a `lineCount`-only formula would OOM them. Large datasets get the RAM they need; small ones stop over-reserving. A `fallback: "12GiB"` preserves the old behavior on backends that cannot evaluate formulas.

CPU scales with the input too: both modes request up to 4 cores (`cpuFormula`), matching their measured thread speedups (peptide ~1.7x, antibody ~1.5x). `POLARS_MAX_THREADS` — which sizes both the polars pool and the numpy pI-bisection workers — takes the cores the backend actually grants, via the `{system.cpu}` command expression, so a sub-cap allocation never oversubscribes its quota. A static `env` covers backends without command expressions.

Determinism holds regardless of core count. Polars parallelizes only order-stable work (CSV read, chain reconstruction, aa_fraction reshape, the unique-key sort, CSV write), and `main.py` pins the BLAS/OpenMP intra-op threads to 1 so the one order-sensitive step — the `counts @ weights` reduction — never splits across threads. The emitted bytes, and the resource CID, stay identical.

Compute-engine optimizations, all output-preserving:

- `_charge_raw` computes `10**ph` once and reuses it via scalar factors instead of a per-amino-acid `10**(ph±pk)` — ~2.4x faster on the charge/pI path (FP drift ~1e-15, absorbed by the 3-dp quantization).
- Peptide mode cleans its column once, sharing the intermediate between the count substrate and the instability index instead of cleaning twice.
- The count and instability scatters use `np.bincount` instead of the slower unbuffered `np.add.at` — bit-identical, ~3x faster.
- The per-chain median-CDR3-length stat uses a vectorized polars `count_matches` instead of one Python `effective_length` call per row.
- Peptide `aa_fraction` uses a polars unpivot instead of a 20×N Python-list explosion, cutting peak memory on the peptide path.
- The pI bisection (`isoelectric_point` / `fv_isoelectric_point`) runs row-parallel across threads — contiguous, data-independent blocks, each depending only on its own rows, so the result is bit-identical regardless of worker count.

Together these cut peptide compute ~2.3x and make antibody scale from ~1.1x to ~1.5x on 4 cores (10.0s → 4.6s at 1M rows).

Output stays within the quantized-equal contract, verified byte-for-byte between 1 and 4 threads by an extended determinism test and against the pre-change implementation (peptide and antibody, including edge cases).
