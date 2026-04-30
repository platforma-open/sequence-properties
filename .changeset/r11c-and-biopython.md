---
'@platforma-open/MiLaboratories.sequence-properties.software': minor
'@platforma-open/MiLaboratories.sequence-properties.workflow': minor
'@platforma-open/MiLaboratories.sequence-properties': minor
---

Resolve missing items from the comprehensive spec review and migrate property computation to BioPython.

- R11c VHH/single-domain antibody detection: Python computes median CDR-H3 length per chain; the workflow emits the spec info message when receptor is IG, only chain A has CDR3, and the median is ≥16 aa.
- Adds missing description annotations on `ered_peptide`, `aliphatic_peptide`, and `ered_Fv`. Differentiates CDR-L3 description from CDR-H3 per `pcolumn-spec.md`.
- Migrates `properties.py` to wrap `Bio.SeqUtils.ProtParam` and `IsoelectricPoint` with IPC 2.0 pKa overrides on the instance. Wrapper layer documented in `docs/reviews/biopython-tradeoffs.md`. Deletes `instability.py` (BioPython has the same DIWV table).
- Architectural shift: info-blob assembly moves to `process.tpl.tengo` since it depends on Python's stats output. Adds a new `stats.json` output to the Python step.
- MW values shift by ~0.07–0.13 Da on antibody chains (BioPython's mass table), below the `.1f` display precision; charge / pI numerics unchanged within the existing 3-dp CID quantization.
