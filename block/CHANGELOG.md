## 1.0.0

## 1.4.0

### Minor Changes

- aa281ac: Per-instance trace label, broader plot pickers, locked-in test coverage.

  - **Trace label is per-instance.** The workflow's `pl7.app/trace.label` resolves to `customBlockLabel || defaultBlockLabel || "Sequence Properties"` (centralised in `model/src/label.ts`). Two sequence-properties blocks on the same dataset show distinguishable entries in Lead Selection and other downstream pickers once the user customises the `PlBlockPage` subtitle. Same pattern as clonotype-clustering and titeseq-analysis PR #13.
  - **Scatter and Histogram metadata pickers accept own-block columns.** Filter, Grouping/Color, Highlight, Size, Tab, Tooltip, Label, and Additional-curves now treat every column in the property pframe as a candidate — own scalars and upstream metadata alike. Users can color the Property Relationships scatter by Aromaticity while plotting Charge vs Hydrophobicity. X/Y axis defaults unchanged.
  - **Migration backfill.** A new `Ver_2026_05_18` step fills the new label fields onto projects tagged at the deployed `Ver_2026_05_05`, preserving any interim-deployed value via `?? ""`. Without the split, already-V2 projects would skip the migration and the workflow would receive `args.customBlockLabel === undefined`.
  - **Test coverage.** Model vitest locks the resolution chain (6 cases) and the migration backfill (4 cases). A subprocess-based Python byte-compare test guards Python output determinism. `build.yaml` enables `test: true` so block-level tests exercise on every PR.

### Patch Changes

- Updated dependencies [aa281ac]
  - @platforma-open/milaboratories.sequence-properties.workflow@1.3.0
  - @platforma-open/milaboratories.sequence-properties.model@1.3.0
  - @platforma-open/milaboratories.sequence-properties.ui@1.3.0

## 1.3.4

### Patch Changes

- Updated dependencies [8646592]
  - @platforma-open/milaboratories.sequence-properties.ui@1.2.4
  - @platforma-open/milaboratories.sequence-properties.model@1.2.4

## 1.3.3

### Patch Changes

- 7cb0850: Refactor to fix CID conflicts between block instances with identical inputs and to align outputs with exports.

  - Workflow: split per-block work out of the deferred render template so it dedups across block instances. The process template is renamed `info.tpl.tengo` and now only assembles the info blob from Python stats; xsv.importFile and pFrame building moved into `main.tpl.tengo` where `blockId`/trace are stamped at the spec layer.
  - Unified pFrame: the same pFrame is now published as both the block's UI output and the result-pool export — one canonical resource, two consumers. All calculated properties (charge, GRAVY, MW, pI, extinction coefficients, instability, aliphatic, aromaticity, ΔCharge) are exported, not just those flagged `isScore`.
  - Columns library: drop bespoke `pl7.app/isOutput` annotation in favor of `pl7.app/trace` for own-block identification in the UI. Public surface simplified to `buildColumns`, `aaFractionColumn`, `cloneSpec`.
  - Table: downgrade to PlAgDataTableV2 with the block's own columns only — fixed, predictable column list — until V3 default-visibility rules stabilize for the multi-source layout this block needs.
  - Clone-id / variant-key axis now visible by default in the table — it's the join key the user reads against the property values.
  - Plot pages: `dataColumnPredicate` filters by trace instead of `isOutput`.
  - Description: broadened to reflect the block's general utility beyond Lead Selection ranking.

- Updated dependencies [7cb0850]
  - @platforma-open/milaboratories.sequence-properties.workflow@1.2.2
  - @platforma-open/milaboratories.sequence-properties.model@1.2.3
  - @platforma-open/milaboratories.sequence-properties.ui@1.2.3

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
