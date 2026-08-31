---
'@platforma-open/milaboratories.sequence-properties.workflow': minor
'@platforma-open/milaboratories.sequence-properties.ui': minor
'@platforma-open/milaboratories.sequence-properties.software': minor
'@platforma-open/milaboratories.sequence-properties': minor
---

Fix inverted TCR chain labels. Chain letter A is the D-recombining chain, so for TCRAB it is beta and for TCRGD it is delta — CDR3 and full-chain column labels, and the partial-coverage info messages, named the other chain. PColumn names and chain domains are unchanged; only labels move.

Per-chain CDR3 descriptions no longer name antibody loops (CDR-H3 / CDR-L3) on TCR input.

Adds a Tengo unit test covering the chain-letter to label mapping for all three receptors, and wires `pl-tengo test` into the workflow package so it runs.

Orders paired-chain columns by the receptor's spoken naming — alpha before beta, gamma before delta — instead of by chain slot. IG keeps heavy before light. The default scatter and histogram source stays on chain A.
