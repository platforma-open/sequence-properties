---
'@platforma-open/milaboratories.sequence-properties.software': patch
'@platforma-open/milaboratories.sequence-properties.workflow': patch
---

Parallelize the property computation across CPU cores and raise the compute step's CPU/memory request, so large datasets compute on the first run without the previous single-core bottleneck. Output is byte-identical to the prior version (worker count never affects results).
