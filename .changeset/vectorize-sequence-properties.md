---
'@platforma-open/milaboratories.sequence-properties.software': patch
'@platforma-open/milaboratories.sequence-properties.workflow': patch
---

Vectorize the property computation — single-threaded numpy/polars array math replaces per-row BioPython. ~5x faster at 50k clones; output unchanged within the quantized-equal contract (signed-zero canonicalized).

Reduce the vectorized engine's peak memory (int8/int32 indices, share the clean intermediate between counts and instability, drop dead retention) — ~4.85 GB/1M clones vs ~8.5 GB before. Raise the compute-properties step's `mem()` from 4GiB to 16GiB to cover ~2M clones at the new footprint with headroom.
