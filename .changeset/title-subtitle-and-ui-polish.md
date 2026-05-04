---
'@platforma-open/MiLaboratories.sequence-properties.model': patch
'@platforma-open/MiLaboratories.sequence-properties.ui': patch
'@platforma-open/MiLaboratories.sequence-properties.workflow': patch
'@platforma-open/MiLaboratories.sequence-properties.software': patch
'@platforma-open/MiLaboratories.sequence-properties': patch
---

Match block title/subtitle convention; tighten user-facing text.

- Block title is now the static "Sequence Properties"; the selected input
  dataset becomes the subtitle. Matches the convention in cdr3-spectratype,
  repertoire-diversity, and clonotype-browser.
- Main page table gains `show-export-button` (CSV export) and a not-ready
  prompt before configuration; input dropdown gets a tooltip describing
  supported inputs and modality auto-detection.
- Removed PColumn-domain namespace tokens from user-facing messages:
  the receptor-not-detected info message and peptide-mode panic no longer
  mention `pl7.app/vdj/receptor` or `pl7.app/sequence` domain shape.
- Python progress log says "full-chain properties" instead of the internal
  `VDJRegion` domain value.
