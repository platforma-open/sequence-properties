---
'@platforma-open/milaboratories.sequence-properties.workflow': patch
'@platforma-open/milaboratories.sequence-properties': patch
---

Size the compute-properties step's memory request from the input's clone count instead of a fixed 12 GiB, so multi-million-clone inputs get the memory they need. The request keeps 12 GiB as its floor.
