---
"@platforma-open/milaboratories.sequence-properties.workflow": patch
"@platforma-open/milaboratories.sequence-properties": patch
---

Score centroid datasets from Sequence Clustering when the clustered dataset came from Amplicon Profiling. In peptide and amplicon modes the block now reads the amino-acid sequence column marked `pl7.app/isMainSequence`, and falls back to the `peptide` / `amplicon-sequence` feature columns for data without the marker. Before, a centroid dataset clustered on a region failed with "amplicon mode detected but no whole-variant amino-acid sequence column".

Update the SDK to 1.88 (model 1.88.0, ui-vue 1.88.3, workflow-tengo 6.22.0) and the toolchain to Node 24.22.1 and pnpm 12.
