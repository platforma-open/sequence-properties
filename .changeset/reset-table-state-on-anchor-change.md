---
'@platforma-open/milaboratories.sequence-properties.ui': patch
---

Reset the properties-table state (column visibility, sort, filters) when the user picks a different input dataset. Persisted state from the previous anchor referenced columns that no longer exist in the new view, which made every column appear default-visible after a switch. Implemented via a `setInput` callback wired to the dropdown's `@update:model-value` (matches the convention in `clonotype-clustering` and avoids the false-positive resets that a `watch` on `inputAnchor` would trigger when the SDK replaces the entire `data` object on server patches).
