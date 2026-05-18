---
'@platforma-open/milaboratories.sequence-properties.workflow': minor
'@platforma-open/milaboratories.sequence-properties.model': minor
'@platforma-open/milaboratories.sequence-properties.ui': minor
'@platforma-open/milaboratories.sequence-properties': minor
---

Per-instance trace label. The workflow resolves `pl7.app/trace.label` to `customBlockLabel || defaultBlockLabel` instead of the hardcoded "Sequence Properties". Two sequence-properties blocks on the same dataset now show distinguishable entries in Lead Selection and other downstream pickers when the user customizes the `PlBlockPage` subtitle. Mirrors titeseq-analysis PR #13.
