# @platforma-open/MiLaboratories.sequence-properties.model

## 1.2.1

### Patch Changes

- 0c6ed9c: Minor fixes

## 1.2.0

### Minor Changes

- c6e975d: Add Scatterplot and Histogram tabs to the Sequence Properties block. Both
  panels read the existing `propertiesPf` p-frame and pick modality-aware
  defaults: peptide charge / hydrophobicity in peptide mode, chain "A" CDR3
  charge / hydrophobicity in antibody/TCR mode. Axis pickers list every
  numeric scalar PColumn emitted by the run, excluding the 2-axis AA fraction
  column. R21 / R21a reference line at GRAVY = 0 deferred — see
  docs/spec-deviations.md SD-009.

## 1.1.2

### Patch Changes

- 9d628d9: **Workflow (SD-008):** Derive receptor from `pl7.app/vdj/chain` when `pl7.app/vdj/receptor` is absent. Bulk MiXCR axes carry the chain key (`IGHeavy`, `IGLight`, `TCRAlpha`, `TCRBeta`, `TCRGamma`, `TCRDelta`) without the receptor key, which previously fired R13b's "Receptor type not detected" warning on every bulk run. Detection precedence now: axis receptor → axis chain → per-column receptor → per-column chain → IG default + warning. Inputs that carry receptor explicitly are unaffected.

  **Model:** Default-visible single upstream amino-acid sequence column matching the analysed coverage tier — peptide (`pl7.app/sequence`, feature=peptide) for peptide mode; full-chain VDJRegion (`pl7.app/vdj/sequence`, feature=VDJRegion or VDJRegionInFrame, chain A) for `full_chain` tier; CDR3 (`pl7.app/vdj/sequence`, feature=CDR3, chain A) for `cdr3_only` and `partial` tiers. Chain B (light / beta / delta) and secondary alleles stay available via the column picker.

  **UI:** Move the "Input dataset" picker into a Settings slide-out drawer (matches the convention used by clonotype-clustering and other sibling blocks). Drawer auto-opens on first load when no input is selected, auto-closes when the workflow starts running.

  **Software:** Switch the Python runenv from `runenv-python-3:3.12.10` to `runenv-python-3:3.12.10-scientific-slim`. The scientific-slim image now ships biopython, so the per-block dependency install is faster and matches the convention used by sibling python blocks (e.g. `titeseq-analysis`).

## 1.1.1

### Patch Changes

- bb07f98: Rename all package scopes from `MiLaboratories.sequence-properties` to `milaboratories.sequence-properties`. npm registry rejects new package names with uppercase letters, which blocked the first publish. Lowercase form aligns with the existing `@platforma-open/milaboratories.*` convention used by sibling blocks. Also corrects the GitHub URL in the block manifest to point at the actual repo (`platforma-open/sequence-properties`).

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
