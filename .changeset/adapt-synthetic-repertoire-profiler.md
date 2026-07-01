---
'@platforma-open/milaboratories.sequence-properties.model': patch
'@platforma-open/milaboratories.sequence-properties.workflow': patch
'@platforma-open/milaboratories.sequence-properties.ui': patch
---

Support `synthetic-repertoire-profiler` (amplicon) variant datasets:

- `detectMode` now recognizes the profiler's `pl7.app/variantKey` axis (axis domain `pl7.app/repertoire/extractionRunId`) as a new `"amplicon"` mode, instead of falling through to the "no recognized sequence key axis" panic.
- Amplicon runs the same whole-sequence physicochemical computation as peptide mode (the Python engine runs it under its `peptide` path), reading the whole-variant amino-acid sequence (`pl7.app/feature: "amplicon-sequence"`).
- Output property columns (and the AA-fraction column) are labeled with the `amplicon-sequence` feature instead of `peptide`, so they attach to the correct entity.
- The UI's default scatter/histogram axes are amplicon-aware (charge / hydrophobicity on the `amplicon-sequence` feature).

Per-region properties for amplicon (using the profiler's region subsequences) are out of scope here — whole-sequence descriptors only.
