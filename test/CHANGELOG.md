# @platforma-open/MiLaboratories.sequence-properties.test

## 1.1.8

### Patch Changes

- aa281ac: Per-instance trace label, broader plot pickers, locked-in test coverage.

  - **Trace label is per-instance.** The workflow's `pl7.app/trace.label` resolves to `customBlockLabel || defaultBlockLabel || "Sequence Properties"` (centralised in `model/src/label.ts`). Two sequence-properties blocks on the same dataset show distinguishable entries in Lead Selection and other downstream pickers once the user customises the `PlBlockPage` subtitle. Same pattern as clonotype-clustering and titeseq-analysis PR #13.
  - **Scatter and Histogram metadata pickers accept own-block columns.** Filter, Grouping/Color, Highlight, Size, Tab, Tooltip, Label, and Additional-curves now treat every column in the property pframe as a candidate — own scalars and upstream metadata alike. Users can color the Property Relationships scatter by Aromaticity while plotting Charge vs Hydrophobicity. X/Y axis defaults unchanged.
  - **Migration backfill.** A new `Ver_2026_05_18` step fills the new label fields onto projects tagged at the deployed `Ver_2026_05_05`, preserving any interim-deployed value via `?? ""`. Without the split, already-V2 projects would skip the migration and the workflow would receive `args.customBlockLabel === undefined`.
  - **Test coverage.** Model vitest locks the resolution chain (6 cases) and the migration backfill (4 cases). A subprocess-based Python byte-compare test guards Python output determinism. `build.yaml` enables `test: true` so block-level tests exercise on every PR.

- Updated dependencies [aa281ac]
  - @platforma-open/milaboratories.sequence-properties.model@1.3.0

## 1.1.7

### Patch Changes

- Updated dependencies [8646592]
  - @platforma-open/milaboratories.sequence-properties.model@1.2.4

## 1.1.6

### Patch Changes

- Updated dependencies [7cb0850]
  - @platforma-open/milaboratories.sequence-properties.model@1.2.3

## 1.1.5

### Patch Changes

- Updated dependencies [73b4e4d]
  - @platforma-open/milaboratories.sequence-properties.model@1.2.2

## 1.1.4

### Patch Changes

- Updated dependencies [0c6ed9c]
  - @platforma-open/milaboratories.sequence-properties.model@1.2.1

## 1.1.3

### Patch Changes

- Updated dependencies [c6e975d]
  - @platforma-open/milaboratories.sequence-properties.model@1.2.0

## 1.1.2

### Patch Changes

- Updated dependencies [9d628d9]
  - @platforma-open/milaboratories.sequence-properties.model@1.1.2

## 1.1.1

### Patch Changes

- Updated dependencies [bb07f98]
  - @platforma-open/milaboratories.sequence-properties.model@1.1.1

## 1.1.0

### Minor Changes

- 1059d80: Initial release of the Sequence Properties block.

  Computes physico-chemical properties (charge, pI, GRAVY, MW, extinction
  coefficients, instability and aliphatic indices, aromaticity, AA composition)
  for peptide and antibody/TCR sequence inputs. The block auto-detects modality
  from the input axes and degrades gracefully on partial coverage: CDR3
  properties when CDR3 is present, full-chain VH/VL when all seven IMGT regions
  are exported, and Fv-level properties when both chains reconstruct. An R11c
  heuristic flags likely VHH/single-domain inputs.

  Property math uses BioPython ProtParam + IsoelectricPoint with IPC 2.0 pKa
  overrides — peptide set for peptide and CDR3 inputs, protein set for full
  VH/VL. Charge and pI round to 3 decimals at the output boundary; combined
  with sorted Tengo iteration, canonical-JSON resources, and sorted TSV writes,
  output bytes hash identically across runs so the block joins the dedup path.

  M3 validation is locked down by `tests/unit/test_m3_validation.py` (38 cases:
  ≥5 VH pI, ≥2 VL pI, Fv on ≥2 paired chains, ≥10 CDR-H3 charge, ≥3 CDR-L3
  charge, ≥3 VH aliphatic) against pinned IPC 2.0 webserver values and an
  independent Henderson-Hasselbalch reference.

  Block title is the static "Sequence Properties"; the selected input dataset
  appears as the subtitle.

### Patch Changes

- Updated dependencies [1059d80]
  - @platforma-open/MiLaboratories.sequence-properties.model@1.1.0
