---
'@platforma-open/milaboratories.sequence-properties': minor
'@platforma-open/milaboratories.sequence-properties.workflow': minor
---

Read the input's declared `pl7.app/modality` instead of guessing from the run-id key

`synthetic-repertoire-profiler` runs one pipeline over both antibody/TCR parents and designed
libraries, and everything it emits sits on the modality-neutral `pl7.app/variantKey` axis. It
keeps the same `pl7.app/repertoire/extractionRunId` key in both cases, so a VDJ run used to be
read as amplicon and scanned as one whole sequence — its FR1–FR4 regions, CDR3 included, were
never used.

The profiler now declares which kind of repertoire it produced in the entity-axis domain, and
this block reads that declaration: `vdj` takes the per-region antibody/TCR path, `amplicon` the
whole-sequence path. A VDJ run therefore gets CDR3 charge and hydrophobicity, and full-chain
properties, the same as MiXCR input does.

Datasets with no declaration keep the previous behaviour exactly, so projects made before the
declaration landed are unaffected.

Two smaller things follow from it. A declared VDJ input carries no receptor — the profiler has
none to emit, since germline auto-detect builds its reference from the user's own parent
sequences — so the "receptor not detected" notice no longer tells those users to pick a MiXCR
preset they never used; it says the labels defaulted to the antibody convention and that only
labels are affected. And modality detection moved into its own workflow module so it could be
unit-tested; `pnpm test` on the workflow package now runs those tests.
