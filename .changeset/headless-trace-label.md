---
"@platforma-open/milaboratories.sequence-properties.model": patch
"@platforma-open/milaboratories.sequence-properties.ui": patch
"@platforma-open/milaboratories.sequence-properties.workflow": patch
---

Compute the block label without the UI

The trace label is now the custom block label or "Sequence Properties", computed by the
model, so opening the block no longer changes its arguments and a block that was never
opened gets the same label. The input dataset (or subset) name moves to the subtitle only,
written when the user picks the input; it was already part of the inherited upstream trace.
Existing blocks show as stale once after the update because the arguments change.
