---
'@platforma-open/milaboratories.sequence-properties.software': patch
---

Performance: cache per-sequence `_prepare`, `ProteinAnalysis`, and `IsoelectricPoint` via a `SequenceContext` so each sequence does the BioPython setup work once instead of per-property. Pipeline reuses the full-chain context for the Fv pass (one `IsoelectricPoint(IPC2_PROTEIN, include_cys=False)` shared between `charge_at_pH(7.0)` and pI bisection per chain). DataFrames are now built columnarly (dict-of-lists) instead of via list-of-dicts. Output is byte-identical to pre-refactor on the corpus tests; ~1.7× faster on per-property micro-bench, end-to-end ~40k peptides/s and ~10k antibody-clones/s with full-chain + Fv. CDR3 chain-mode byte-stability tests added.
