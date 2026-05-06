# Visualizations — implementation plan (PR #118 in docs/text)

Spec source: `docs/text/work/projects/sequence-properties/README.md` §Visualizations + R18–R21a.

## Summary of the spec change

Add two graph-maker panels reading the existing `propertiesPf` p-frame:
- **Scatterplot** — defaults adapt to detected modality (peptide charge vs hydrophobicity, or chain-A CDR3 charge vs hydrophobicity for IG/TCRAB/TCRGD).
- **Histogram** — default metric also modality-adaptive (peptide hydrophobicity / chain-A CDR3 hydrophobicity).
- **Axis pickers** enumerate every numeric scalar PColumn emitted by the run; exclude the 2-axis `pl7.app/aaFraction` column (R7) and discrete axes.
- **Fallback** when defaults absent: first one/two numeric scalars in workflow emission order (peptide → CDR3 → full-chain → Fv).
- **Reference line** at 0 on a hydrophobicity scatter axis: inject `significantLines: [0]` block-locally.
- **Histogram reference line**: defer — `composeHistogramSettings.ts` does not consume `significantLines` today (verified). No fallback display logic; line just doesn't render.
- **Phase 2 (out of scope):** cross-component selection, lasso → table, region shading.

## Workspace patterns to align with

- Reference pages: `blocks/clonotype-enrichment/ui/src/pages/ScatterPage.vue`, `blocks/clonotype-clustering/ui/src/pages/HistogramPage.vue`, `blocks/titeseq-analysis/ui/src/pages/KDDistributionPage.vue`.
- Graph state shape stored in `BlockData`: `GraphMakerState` from `@milaboratories/graph-maker`.
- PFrame for graph-maker: `createPFrameForGraphs(ctx, pCols)` exposed via `outputWithStatus(..., (ctx): PFrameHandle | undefined => …)`.
- Default-options use `PredefinedGraphOption<'scatterplot'|'histogram'>` with `selectedSource: PColumnSpec` looked up by `name`+`domain`.
- `dataColumnPredicate` filters columns shown in axis pickers.
- **No workspace block today renders a graph-maker panel and `PlAgDataTableV2` on the same page** — the convention is one section per panel. Spec says "alongside the properties table on the block's main page"; need to confirm layout (see Q1).

## Plan

### 1. Workflow

**No changes.** All required PColumns already emitted by `process.tpl.tengo` (`pl7.app/charge`, `pl7.app/hydrophobicity`, `pl7.app/chargeShift`, AA composition fractions, `pl7.app/aaFraction`). Workflow emission order in `process.tpl.tengo` is already peptide → CDR3 → full-chain → Fv (lines ~118–354), satisfying R19a/R20a fallback ordering.

### 2. Model (`model/src/`)

#### `types.ts`
Add to `BlockData`:
```ts
import type { GraphMakerState } from '@milaboratories/graph-maker';

graphStateScatter: GraphMakerState;
graphStateHistogram: GraphMakerState;
```

#### `dataModel.ts`
Bump data-model version (e.g. `"Ver_2026_05_05"`) with a migration that adds the two new graph-state fields initialised to a minimal default `GraphMakerState`. Keep `"Ver_2026_04_28"` registered.

#### `index.ts`
Add three new outputs:

```ts
.outputWithStatus("propertiesPfHandle", (ctx): PFrameHandle | undefined => {
  const pCols = ctx.outputs?.resolve("propertiesPf")?.getPColumns();
  if (pCols === undefined) return undefined;
  return createPFrameForGraphs(ctx, pCols);
})

.output("propertiesPfCols", (ctx) => {
  const pCols = ctx.outputs?.resolve("propertiesPf")?.getPColumns();
  if (pCols === undefined) return undefined;
  return pCols.map(c => ({ columnId: c.id, spec: c.spec }) satisfies PColumnIdAndSpec);
})
```

The UI uses `propertiesPfCols` for default-axis lookup and column-predicate filtering, and `propertiesPfHandle` to feed the chart.

If layout is **separate sections** (Q1), extend `.sections(...)` with two more entries:
```ts
.sections(() => [
  { type: 'link' as const, href: '/' as const, label: 'Properties' },
  { type: 'link' as const, href: '/scatter' as const, label: 'Scatterplot' },
  { type: 'link' as const, href: '/histogram' as const, label: 'Histogram' },
])
```

### 3. UI (`ui/src/`)

#### Numeric-scalar predicate (R18a)
```ts
const NUMERIC = new Set(['Int', 'Long', 'Float', 'Double']);
const dataColumnPredicate = (spec: PColumnSpec) =>
  NUMERIC.has(spec.valueType)
  && spec.axesSpec.length === 1            // excludes 2-axis aaFraction (R7)
  && spec.name !== 'pl7.app/aaFraction';   // belt-and-braces
```

#### Default lookup (R19, R20)
Read modality from `app.model.outputs.info` (existing `WorkflowInfo` carries `mode` + `receptor` + `coverageTier`). Map:
- `mode === 'peptide'` → look up by `name === 'pl7.app/charge'/'hydrophobicity'` with `domain['pl7.app/feature'] === 'peptide'`.
- antibody/TCR (any other mode) → look up by `name === 'pl7.app/charge'/'hydrophobicity'` with `domain['pl7.app/feature'] === 'CDR3'` and `domain['pl7.app/vdj/scClonotypeChain'] === 'A'`.

The label naming (`CDR-H3 / CDR-α3 / CDR-γ3`) is already encoded in column annotations by the workflow (R13a) — graph-maker reads label from spec annotations, no UI-side label work required.

#### Fallback (R19a, R20a)
When the modality default is not present in `propertiesPfCols`, take the first one (histogram) or two (scatter) `PColumnIdAndSpec` entries that pass `dataColumnPredicate`, preserving the workflow's emission order.

When fewer than required, render via `statusText.noPframe.title = 'Select X and Y axes to plot'` (scatter) / `'Select a metric to plot'` (histogram).

#### `defaultOptions` shape
```ts
const scatterDefaults = computed((): PredefinedGraphOption<'scatterplot'>[] | null => {
  const cols = app.model.outputs.propertiesPfCols;
  const info = app.model.outputs.info;
  if (!cols || !info) return null;

  const xSpec = pickDefaultX(cols, info) ?? cols.filter(c => isScalar(c.spec))[0]?.spec;
  const ySpec = pickDefaultY(cols, info) ?? cols.filter(c => isScalar(c.spec))[1]?.spec;
  if (!xSpec || !ySpec) return null;

  return [
    { inputName: 'x', selectedSource: xSpec },
    { inputName: 'y', selectedSource: ySpec },
  ];
});
```

#### Significant lines at GRAVY = 0 (R21)
The block writes `axesSettings.axisX.significantLines = [0]` (or axisY) on the persisted `graphStateScatter` whenever the currently-selected source on that axis is a `pl7.app/hydrophobicity` column; clears it otherwise. Implemented as a `watchEffect` over `app.model.data.graphStateScatter.optionsState.components.x/y.selectorStates[0].selectedSource`. This is `data → data`, not `output → data`, so it's outside the canonical hairpin shape (per `harnesses/block-dev/hairpin.md`) — flagged as a deviation candidate; document in `docs/spec-deviations.md` if it stays. Multi-client races converge because the written value is deterministic from the selected source.

#### Histogram reference line (R21a)
Verified at planning time: `composeHistogramSettings.ts` does not consume `significantLines`. **Skip the injection on the histogram path.** No runtime error; no fallback display. Add a comment pointing at `composeHistogramSettings.ts` so the next reviewer understands why the symmetry is broken.

#### Pages
Two new files: `ui/src/pages/ScatterPage.vue`, `ui/src/pages/HistogramPage.vue`. Both follow the pattern in `clonotype-enrichment/ScatterPage.vue` / `clonotype-clustering/HistogramPage.vue`.

`ui/src/app.ts` adds the routes:
```ts
routes: {
  '/': () => MainPage,
  '/scatter': () => ScatterPage,
  '/histogram': () => HistogramPage,
}
```

`MainPage.vue` is unchanged in this layout.

### 4. Tests

Workflow tests do not change (no new outputs that need backend testing). Optional UI smoke test deferred — the block's existing `test/src/wf.test.ts` covers PColumn shape, which is what these panels rely on.

### 5. Changeset

`.changeset/<name>.md`:
```
---
'@platforma-open/milaboratories.sequence-properties.model': minor
'@platforma-open/milaboratories.sequence-properties.ui': minor
'@platforma-open/milaboratories.sequence-properties': minor
---

Add scatterplot and histogram graph-maker panels with modality-aware defaults.
```

(Root package included because the change is `minor`, not `patch`.)

### 6. Build / verify

`pnpm run build:dev` → reload via pl MCP server (`update_block` on the existing dev project) → verify scatter + histogram render with peptide and antibody fixtures, axis dropdowns enumerate scalars, AA fraction column absent from menu, hydrophobicity = 0 line shows on scatter when selected.

## Decisions (operator confirmed 2026-05-05)

1. **Layout — separate sections.** Three section links: `Properties` (existing main page with the table), `Scatterplot`, `Histogram`. Two new pages under `ui/src/pages/`.
2. **R18a numeric-scalar filter — `valueType ∈ {Int, Long, Float, Double} && axesSpec.length === 1 && name !== 'pl7.app/aaFraction'`.** The 2-axis check excludes the AA fraction column; the explicit name check is belt-and-braces. Sufficient for the current PColumn set.
3. **R21 reference line — deferred.** Graph-maker has no path to inject `significantLines` on a data-column axis today (verified in `composeScatterplotSettings.ts`, `getAxesDataFromForms.ts`, `composeHistogramSettings.ts`). Logged as `SD-009` in `docs/spec-deviations.md`. Block ships the two panels with no reference line on either chart. Pick up R21/R21a when the platform-side threshold extension named in the spec lands.
4. **Data-model migration version — `Ver_2026_05_05`.** Chained migration via `DataModelBuilder.add("Ver_2026_05_05", prev => ({...prev, graphStateScatter: <default>, graphStateHistogram: <default>}))`. Existing `Ver_2026_04_28` stays loadable.

## Out of plan (until operator unblocks)

- R21 / R21a reference line implementation (covered by SD-009).
- Cross-component selection (table ↔ scatter / histogram) — Phase 2 per spec.
- Lasso → table selection — Phase 2 per spec.
- Region/quadrant shading — Phase 2 per spec.

I'll wait for the green light before starting code changes.
