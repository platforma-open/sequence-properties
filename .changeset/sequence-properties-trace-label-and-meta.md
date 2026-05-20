---
'@platforma-open/milaboratories.sequence-properties.workflow': minor
'@platforma-open/milaboratories.sequence-properties.model': minor
'@platforma-open/milaboratories.sequence-properties.ui': minor
'@platforma-open/milaboratories.sequence-properties': minor
'@platforma-open/milaboratories.sequence-properties.test': patch
---

Per-instance trace label, broader plot pickers, locked-in test coverage.

- **Trace label is per-instance.** The workflow's `pl7.app/trace.label` resolves to `customBlockLabel || defaultBlockLabel || "Sequence Properties"` (centralised in `model/src/label.ts`). Two sequence-properties blocks on the same dataset show distinguishable entries in Lead Selection and other downstream pickers once the user customises the `PlBlockPage` subtitle. Same pattern as clonotype-clustering and titeseq-analysis PR #13.
- **Scatter and Histogram metadata pickers accept own-block columns.** Filter, Grouping/Color, Highlight, Size, Tab, Tooltip, Label, and Additional-curves now treat every column in the property pframe as a candidate — own scalars and upstream metadata alike. Users can color the Property Relationships scatter by Aromaticity while plotting Charge vs Hydrophobicity. X/Y axis defaults unchanged.
- **Migration backfill.** A new `Ver_2026_05_18` step fills the new label fields onto projects tagged at the deployed `Ver_2026_05_05`, preserving any interim-deployed value via `?? ""`. Without the split, already-V2 projects would skip the migration and the workflow would receive `args.customBlockLabel === undefined`.
- **Test coverage.** Model vitest locks the resolution chain (6 cases) and the migration backfill (4 cases). The two-instance dedup integration test now runs in CI alongside a subprocess-based Python byte-compare test. `build.yaml` enables `test: true` so block-level tests exercise on every PR.
