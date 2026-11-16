---
'@platforma-open/milaboratories.sequence-properties.workflow': patch
'@platforma-open/milaboratories.sequence-properties': patch
---

Size the compute-properties step's memory request from the input's clone count instead of a fixed 48 GiB, so multi-million-clone inputs get the memory they need and small inputs request less. The request has a 12 GiB floor; backends that cannot evaluate resource formulas keep the static 48 GiB.
