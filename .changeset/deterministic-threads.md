---
'@platforma-open/milaboratories.sequence-properties.software': patch
'@platforma-open/milaboratories.sequence-properties.workflow': patch
---

Parallelize the compute step and scale its resource requests to the input size, while keeping output byte-identical.

The workflow now sizes memory and CPU from the input via `exec` resource formulas instead of a flat `mem("12GiB")`/`cpu(1)`. Memory scales with clone count — `memFormula(clamp(lineCount * perRow, 2GiB, 64GiB))` — using per-row constants calibrated from measured peak RSS (~3 KiB/row peptide, ~4.3 KiB/row antibody), so large datasets get the RAM they need instead of being OOM-killed while small ones stop over-reserving. Peptide/amplicon mode additionally scales CPU up to 4 cores (`cpuFormula`), matching its measured ~1.7x thread speedup; antibody/TCR stays single-core because its wall time is dominated by the single-threaded numpy path (~1.1x even on 8 cores). A `fallback: "12GiB"` preserves the previous behavior on backends that cannot evaluate formulas.

Determinism: polars parallelizes only order-stable work (CSV read, chain reconstruction, aa_fraction reshape, the unique-key sort, CSV write), and `main.py` pins the BLAS/OpenMP intra-op threads to 1 so the only order-sensitive step — the `counts @ weights` numpy reduction — never splits across threads. The emitted bytes (and the resource CID) are therefore independent of the core count.

Peptide-mode `aa_fraction` is rebuilt with a polars unpivot instead of a 20×N Python-list explosion, cutting peak memory on the peptide path.

Three single-threaded compute optimizations, all output-preserving: (1) `_charge_raw` computes `10**ph` once and reuses it via scalar factors instead of a per-amino-acid `10**(ph±pk)` (~2.4x faster on the charge/pI path; FP drift ~1e-15, absorbed by the 3-dp quantization); (2) peptide mode cleans its column once and shares the intermediate between the count substrate and the instability index, instead of cleaning twice; (3) the count and instability scatters use `np.bincount` instead of the slower unbuffered `np.add.at` (bit-identical, ~3x faster). Together these cut peptide compute ~2.3x and antibody wall time ~25% at 1M rows.

Output is unchanged within the quantized-equal contract, verified byte-for-byte between 1 and 4 threads by an extended determinism test and against the pre-change implementation (peptide and antibody, including edge cases).
