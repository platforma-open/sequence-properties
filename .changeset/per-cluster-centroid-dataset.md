---
"@platforma-open/milaboratories.sequence-properties.workflow": patch
"@platforma-open/milaboratories.sequence-properties.model": patch
---

Add support for per-cluster centroid datasets (clonotype-clustering). The model accepts a 2-axis `[sampleId, centroidId]` anchor, and the workflow treats the synthetic per-cluster centroid sequence as legacy-bulk VDJ for sequence pulling.
