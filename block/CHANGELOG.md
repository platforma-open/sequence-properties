## 1.0.0

## 1.1.2

### Patch Changes

- Updated dependencies [9d628d9]
  - @platforma-open/milaboratories.sequence-properties.workflow@1.1.2
  - @platforma-open/milaboratories.sequence-properties.model@1.1.2
  - @platforma-open/milaboratories.sequence-properties.ui@1.1.2

## 1.1.1

### Patch Changes

- bb07f98: Rename all package scopes from `MiLaboratories.sequence-properties` to `milaboratories.sequence-properties`. npm registry rejects new package names with uppercase letters, which blocked the first publish. Lowercase form aligns with the existing `@platforma-open/milaboratories.*` convention used by sibling blocks. Also corrects the GitHub URL in the block manifest to point at the actual repo (`platforma-open/sequence-properties`).
- Updated dependencies [bb07f98]
  - @platforma-open/milaboratories.sequence-properties.workflow@1.1.1
  - @platforma-open/milaboratories.sequence-properties.model@1.1.1
  - @platforma-open/milaboratories.sequence-properties.ui@1.1.1

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
  - @platforma-open/MiLaboratories.sequence-properties.ui@1.1.0
  - @platforma-open/MiLaboratories.sequence-properties.workflow@1.1.0

Initial scaffold of `sequence-properties` block: model (BlockModelV3), Tengo workflow with modality detection and region collection, UI (PlBlockPage + PlAgDataTableV2), Python stub for property computation.
