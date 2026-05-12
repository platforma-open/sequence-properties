## 1.0.0

## 1.3.2

### Patch Changes

- Updated dependencies [73b4e4d]
  - @platforma-open/milaboratories.sequence-properties.workflow@1.2.1
  - @platforma-open/milaboratories.sequence-properties.model@1.2.2
  - @platforma-open/milaboratories.sequence-properties.ui@1.2.2

## 1.3.1

### Patch Changes

- Updated dependencies [0c6ed9c]
  - @platforma-open/milaboratories.sequence-properties.model@1.2.1
  - @platforma-open/milaboratories.sequence-properties.ui@1.2.1

## 1.3.0

### Minor Changes

- c6e975d: Add Scatterplot and Histogram tabs to the Sequence Properties block. Both
  panels read the existing `propertiesPf` p-frame and pick modality-aware
  defaults: peptide charge / hydrophobicity in peptide mode, chain "A" CDR3
  charge / hydrophobicity in antibody/TCR mode. Axis pickers list every
  numeric scalar PColumn emitted by the run, excluding the 2-axis AA fraction
  column. R21 / R21a reference line at GRAVY = 0 deferred — see
  docs/spec-deviations.md SD-009.

### Patch Changes

- Updated dependencies [c6e975d]
  - @platforma-open/milaboratories.sequence-properties.model@1.2.0
  - @platforma-open/milaboratories.sequence-properties.ui@1.2.0

## 1.2.0

### Minor Changes

- a3eaf4b: Add ΔCharge (pH 7.4 → 6.0) metric — `pl7.app/chargeShift` — emitted at peptide, CDR3 (per chain), and Fv scopes. Captures pH-switching capacity (FcRn recycling, endosomal release); negative values mean the molecule gains positive charge on acidification, the productive direction for histidine-driven pH switching. Histidine dominates the metric (~−0.46 per His; pKa ~6.0 sits in the window). Domain carries the pH endpoints (`pl7.app/pH/from`, `pl7.app/pH/to`) so additional pH pairs can land later without breaking the v1 column identity. Default-visible alongside the static charge column at each scope; not marked `isScore` (interpretive, not a Lead Selection ranking criterion).

  Performance: cache per-sequence `_prepare`, `ProteinAnalysis`, and `IsoelectricPoint` via a `SequenceContext` so each sequence does the BioPython setup work once instead of per-property. Pipeline reuses the full-chain context for the Fv pass (one `IsoelectricPoint(IPC2_PROTEIN, include_cys=False)` shared between `charge_at_pH(7.0)` and pI bisection per chain). DataFrames are built columnarly (dict-of-lists) instead of via list-of-dicts. Output is byte-identical to pre-refactor on the corpus tests; ~1.7× faster on per-property micro-bench, end-to-end ~40k peptides/s and ~10k antibody-clones/s with full-chain + Fv. CDR3 chain-mode byte-stability tests added.

### Patch Changes

- Updated dependencies [a3eaf4b]
  - @platforma-open/milaboratories.sequence-properties.workflow@1.2.0

## 1.1.3

### Patch Changes

- Updated dependencies [74fbd0d]
  - @platforma-open/milaboratories.sequence-properties.ui@1.1.3

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
