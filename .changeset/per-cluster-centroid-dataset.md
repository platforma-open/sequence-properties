---
"@platforma-open/milaboratories.sequence-properties.workflow": patch
"@platforma-open/milaboratories.sequence-properties.model": patch
---

Add support for per-cluster centroid datasets (clonotype-clustering). The model accepts a 2-axis `[sampleId, centroidId]` anchor, and the workflow processes the synthetic per-cluster consensus sequence in peptide mode (the upstream "Export consensus sequences as a dataset" feature is peptide-only, so the centroid sequence column is a peptide amino-acid column).
