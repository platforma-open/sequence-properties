---
'@platforma-open/milaboratories.sequence-properties.workflow': patch
'@platforma-open/milaboratories.sequence-properties.model': patch
'@platforma-open/milaboratories.sequence-properties.ui': patch
'@platforma-open/milaboratories.sequence-properties': patch
---

Refactor to fix CID conflicts between block instances with identical inputs and to align outputs with exports.

- Workflow: split per-block work out of the deferred render template so it dedups across block instances. The process template is renamed `info.tpl.tengo` and now only assembles the info blob from Python stats; xsv.importFile and pFrame building moved into `main.tpl.tengo` where `blockId`/trace are stamped at the spec layer.
- Unified pFrame: the same pFrame is now published as both the block's UI output and the result-pool export — one canonical resource, two consumers. All calculated properties (charge, GRAVY, MW, pI, extinction coefficients, instability, aliphatic, aromaticity, ΔCharge) are exported, not just those flagged `isScore`.
- Columns library: drop bespoke `pl7.app/isOutput` annotation in favor of `pl7.app/trace` for own-block identification in the UI. Public surface simplified to `buildColumns`, `aaFractionColumn`, `cloneSpec`.
- Table: downgrade to PlAgDataTableV2 with the block's own columns only — fixed, predictable column list — until V3 default-visibility rules stabilize for the multi-source layout this block needs.
- Clone-id / variant-key axis now visible by default in the table — it's the join key the user reads against the property values.
- Plot pages: `dataColumnPredicate` filters by trace instead of `isOutput`.
- Description: broadened to reflect the block's general utility beyond Lead Selection ranking.
