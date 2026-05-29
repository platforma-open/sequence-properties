---
'@platforma-open/milaboratories.sequence-properties.ui': minor
'@platforma-open/milaboratories.sequence-properties': minor
---

Closeable info messages on the Main tab — session-only variant. The
advisory alerts emitted by the workflow now show a close button.
Dismissals live in a module-scope UI ref: they persist across in-block
navigation and reset on project close, block reload, or app restart. No
`BlockData` change; no migration; the block model is unchanged.

This is the minimal-surface alternative to a persisted-dismissal
approach. Pick this when you want "acknowledged for this session, fresh
on reopen" semantics and don't want to extend the block schema.
