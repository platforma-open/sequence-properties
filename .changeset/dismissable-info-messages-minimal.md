---
'@platforma-open/milaboratories.sequence-properties.ui': minor
'@platforma-open/milaboratories.sequence-properties': minor
---

Info-message alerts on the Main tab now have a close button. Dismissals
follow the firing: hidden while the message is in the workflow output,
surface fresh if the workflow stops emitting it and later re-emits.
State lives in a module-scope UI ref — persists across in-block
navigation, resets on project close, block reload, or app restart.
