---
"@platforma-open/milaboratories.sequence-properties.software": patch
"@platforma-open/milaboratories.sequence-properties.workflow": patch
"@platforma-open/milaboratories.sequence-properties": patch
---

Read, compute and write the table in bounded batches when the input is in entity key order (the order the workflow provides), and sort it in memory only when it is not, so peak memory no longer grows with the number of clones or the sequence length: under 1 GiB for 100M paired antibodies and for 200M peptides. The compute step requests a fixed 12 GiB, and the table export and imports use the SDK's input-volume sizing.
