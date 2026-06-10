---
'@platforma-open/milaboratories.sequence-properties.software': patch
---

Vectorize the property computation — single-threaded numpy/polars array math replaces per-row BioPython. ~5x faster at 50k clones; output unchanged within the quantized-equal contract (signed-zero canonicalized).
