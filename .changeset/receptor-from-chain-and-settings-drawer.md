---
'@platforma-open/milaboratories.sequence-properties.workflow': patch
'@platforma-open/milaboratories.sequence-properties.model': patch
'@platforma-open/milaboratories.sequence-properties.ui': patch
'@platforma-open/milaboratories.sequence-properties.software': patch
---

**Workflow (SD-008):** Derive receptor from `pl7.app/vdj/chain` when `pl7.app/vdj/receptor` is absent. Bulk MiXCR axes carry the chain key (`IGHeavy`, `IGLight`, `TCRAlpha`, `TCRBeta`, `TCRGamma`, `TCRDelta`) without the receptor key, which previously fired R13b's "Receptor type not detected" warning on every bulk run. Detection precedence now: axis receptor → axis chain → per-column receptor → per-column chain → IG default + warning. Inputs that carry receptor explicitly are unaffected.

**Model:** Default-visible single upstream amino-acid sequence column matching the analysed coverage tier — peptide (`pl7.app/sequence`, feature=peptide) for peptide mode; full-chain VDJRegion (`pl7.app/vdj/sequence`, feature=VDJRegion or VDJRegionInFrame, chain A) for `full_chain` tier; CDR3 (`pl7.app/vdj/sequence`, feature=CDR3, chain A) for `cdr3_only` and `partial` tiers. Chain B (light / beta / delta) and secondary alleles stay available via the column picker.

**UI:** Move the "Input dataset" picker into a Settings slide-out drawer (matches the convention used by clonotype-clustering and other sibling blocks). Drawer auto-opens on first load when no input is selected, auto-closes when the workflow starts running.

**Software:** Switch the Python runenv from `runenv-python-3:3.12.10` to `runenv-python-3:3.12.10-scientific-slim`. The scientific-slim image now ships biopython, so the per-block dependency install is faster and matches the convention used by sibling python blocks (e.g. `titeseq-analysis`).
