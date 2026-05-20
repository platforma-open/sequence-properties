---
'@platforma-open/milaboratories.sequence-properties.model': patch
'@platforma-open/milaboratories.sequence-properties.test': patch
---

Lock label-resolution behaviour with model-package vitest. `model/src/label.ts` centralises `resolveSubtitle` and `resolveTraceLabel`; `model/src/label.test.ts` and `model/src/dataModel.test.ts` cover the resolution chain and the `Ver_2026_05_18` migration backfill respectively. Stages SC IG fastq fixtures so the two-instance dedup integration test runs in CI instead of self-skipping.
