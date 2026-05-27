---
'@platforma-open/milaboratories.sequence-properties.model': minor
'@platforma-open/milaboratories.sequence-properties.ui': minor
'@platforma-open/milaboratories.sequence-properties': minor
---

Closeable info messages on the Main tab. The advisory alerts emitted by
the workflow (VHH detection, partial-region inputs, peptide-instability
floor, etc.) now show a close button. Dismissals persist in
`BlockData.dismissedInfoMessages` — server-side, across project reopens
and clients. The Settings modal includes a "Reset dismissed info
messages" action; a small "N hidden — show all" footer surfaces the
dismissed state in-context.

Model schema: new `Ver_2026_05_27` migration step backfills
`dismissedInfoMessages: []` on existing projects.
