---
'@platforma-open/MiLaboratories.sequence-properties.workflow': minor
'@platforma-open/MiLaboratories.sequence-properties.software': minor
'@platforma-open/MiLaboratories.sequence-properties.model': patch
'@platforma-open/MiLaboratories.sequence-properties': minor
---

Address punchlist follow-ups from the comprehensive spec review.

- **IPC 2.0 pKa correction (CRITICAL).** Verified `pka_tables.py` against the official IPC 2.0 site (`http://ipc2-isoelectric-point.org/theory.html`). Previous values were the IPC 1.0 peptide / protein scales (Kozlowski 2016) with labels swapped between the two contexts — every pI and charge result emitted to date was miscalibrated against the spec's intended IPC 2.0 reference. Numerical impact (VH/VL/Fv): pI shifts of 0.6–1.7 units, VH charge sign flip. CIDs for `charge_*` / `pi_*` columns will all change on next run.
- **`blockId` mutation fix.** `process.tpl.tengo` previously mutated `col.spec.domain` in place when stamping `pl7.app/blockId` on export columns; because the column references were shared with `outputSpecs.columns`, the stamp leaked to `propertiesPf`. Export columns are now built via fresh dict construction.
- **R11b silent fallthrough fix.** A chain with 1-6 of 7 regions but no CDR3 fell through every branch in `main.tpl.tengo` chain-collection loop; the user saw no per-chain output and no explanation. New `else` branch emits a "CDR3 absent" info message.
- **`defaultBlockLabel` wiring.** Model `.title()` now reads `ctx.data.defaultBlockLabel ?? "Sequence Properties"`, reflecting the selected input dataset in the block title. Field removed from `BlockArgs` (UI-only state).
- **`spec-deviations.md`** documents SD-004…SD-007 covering peptide_seq column name, omitted receptor_type column, last-axis selection, and `pl7.app/vdj/clonotypeKey` alias.
