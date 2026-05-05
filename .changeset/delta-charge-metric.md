---
'@platforma-open/milaboratories.sequence-properties.software': minor
'@platforma-open/milaboratories.sequence-properties.workflow': minor
'@platforma-open/milaboratories.sequence-properties': minor
---

Add ΔCharge (pH 7.4 → 6.0) metric — `pl7.app/chargeShift` — emitted at peptide, CDR3 (per chain), and Fv scopes. Captures pH-switching capacity (FcRn recycling, endosomal release); negative values mean the molecule gains positive charge on acidification, the productive direction for histidine-driven pH switching. Histidine dominates the metric (~−0.46 per His; pKa ~6.0 sits in the window). Domain carries the pH endpoints (`pl7.app/pH/from`, `pl7.app/pH/to`) so additional pH pairs can land later without breaking the v1 column identity. Default-visible alongside the static charge column at each scope; not marked `isScore` (interpretive, not a Lead Selection ranking criterion).

Performance: cache per-sequence `_prepare`, `ProteinAnalysis`, and `IsoelectricPoint` via a `SequenceContext` so each sequence does the BioPython setup work once instead of per-property. Pipeline reuses the full-chain context for the Fv pass (one `IsoelectricPoint(IPC2_PROTEIN, include_cys=False)` shared between `charge_at_pH(7.0)` and pI bisection per chain). DataFrames are built columnarly (dict-of-lists) instead of via list-of-dicts. Output is byte-identical to pre-refactor on the corpus tests; ~1.7× faster on per-property micro-bench, end-to-end ~40k peptides/s and ~10k antibody-clones/s with full-chain + Fv. CDR3 chain-mode byte-stability tests added.
