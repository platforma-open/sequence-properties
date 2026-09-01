# Sequence Properties — Spec Deviations

Implementation choices that diverge from the spec at
`docs/text/work/projects/sequence-properties/`. Spec text stays authoritative;
this file records where the code chose differently and why.

---

## SD-001: Skip secondary alleles in single-cell paired data

**Status:** applied
**Date:** 2026-04-29
**Affected file:** `workflow/src/main.tpl.tengo` (chain-collection loop)

### Symptom

The workflow panicked before the Python step ran:

```
header 'A_CDR3' is not unique
@platforma-sdk/workflow-tengo:pframes.xsv-builder:268
```

Two columns mapped to the TSV header `A_CDR3`. Both carried
`pl7.app/vdj/feature: "CDR3"` and `pl7.app/vdj/scClonotypeChain: "A"`. They
differed only on `pl7.app/vdj/scClonotypeChain/index`: one was `"primary"`,
the other `"secondary"`.

### Root cause

MiXCR's single-cell output emits primary and secondary alleles per chain when
a cell expresses two of the same receptor. The original chain-collection loop
built each TSV header from `chain + "_" + feat` and ignored
`pl7.app/vdj/scClonotypeChain/index`, so both alleles claimed the same header.

The spec talks about `pl7.app/vdj/scClonotypeChain` (A/B) as the only chain
disambiguator and never mentions the `/index` domain key. Single-cell paired
alleles are a real-world data shape the spec did not anticipate.

### Trigger

- Dataset: TinySC (Partial), 10x Genomics single-cell IG.
- Modality detected: `antibody_tcr_legacy_sc`.
- MiXCR clonotyping run ID: `7b9a7759-5914-4282-bc00-777c50a983f1`.

### Impact (before fix)

The Python step never ran. Properties were not computed. Bulk MiXCR data and
peptide data were unaffected — the `/index` key is single-cell-specific.

### Options considered

**A. Filter to primary alleles only. [chosen]**
Skip columns where `pl7.app/vdj/scClonotypeChain/index` is present and not
`"primary"`. Two-line change at the chain-collection loop. Matches the spec's
single-allele-per-chain mental model. The primary allele is the canonical
input for developability scoring; secondary alleles are usually expression
artifacts.

**B. Emit per-allele output columns.**
Disambiguate the header (`A_CDR3_primary` vs `A_CDR3_secondary`) and double
the output schema so every property carries primary/secondary variants. Larger
scope, schema change, downstream consumers (Lead Selection) re-keyed, spec
sign-off required.

**C. Surface the collision as a block-level info message and skip the chain.**
Detect duplicate primary+secondary columns, emit "single-cell paired alleles
detected — chain X skipped" via `infoMessages`, and proceed without that
chain. Conservative — produces a partial result instead of a hard failure but
hides data the user supplied.

### Decision

**A.** Smallest change, no schema impact, no downstream re-keying. Revisit if
a customer asks for secondary-allele properties — at that point write a
separate spec for B.

### Implementation

```tengo
// workflow/src/main.tpl.tengo, inside the vdjCols for-loop:
idx := d["pl7.app/vdj/scClonotypeChain/index"]
if idx != undefined && idx != "primary" { continue }
```

### References

- Spec sections touched: `README.md` Requirements R4, R11, R13b.
- Domain key source: `mixcr-clonotyping/workflow/src/process.tpl.tengo` (where
  `pl7.app/vdj/scClonotypeChain/index` is set on single-cell paired output).
- Pattern block also lacks this filter, so the same bug likely affects
  `antibody-sequence-liabilities` on paired single-cell data.

---

## SD-002: Treat `FR4InFrame` as `FR4`

**Status:** applied
**Date:** 2026-04-29
**Affected file:** `workflow/src/main.tpl.tengo` (chain-collection loop)

### Symptom

Every chain reported `6 of 7 required regions found`. Full-chain reconstruction
never ran. Output stayed `coverageTier: cdr3_only` even on full-VDJ MiXCR
presets where every region is exported.

### Root cause

MiXCR exports the FR4 region with `pl7.app/vdj/feature: "FR4InFrame"` (the
in-frame-filtered translation), not `"FR4"`. Spec R4 lists the seven required
regions with the literal name `"FR4"`, and our loop checks
`contains(REQUIRED_FEATURES, feat)` against that literal. Every `FR4InFrame`
column is therefore filtered out, leaving 6 of 7 regions per chain.

Confirmed via `mcp__pl__list_columns` on the MiXCR Clonotyping output PFrame:
all six MiXCR chain slots (Heavy / Light / Alpha / Beta / Gamma / Delta) emit
`<chain> FR4InFrame aa Primary` columns alongside the other six regions.

### Trigger

Any MiXCR preset that emits all VDJ regions. Confirmed on
`10x-sc-xcr-vdj` (5' single-cell). Spec did not anticipate the
`InFrame` suffix on FR4.

### Impact (before fix)

Full-chain properties (charge, pI, hydrophobicity, MW, extinction coefficients,
instability index, aliphatic index, aromaticity) and Fv properties were never
emitted on real MiXCR data. Block silently degraded to CDR3-only output even
when full-chain coverage was available.

### Options considered

**A. Normalise `FR4InFrame` → `FR4` in the chain-collection loop. [chosen]**
One line at the feat-extraction site. Keeps the rest of the workflow, the
seqTable header naming, and the Python 7-region concatenation unchanged.
Treats the in-frame variant as the canonical FR4.

**B. Use `VDJRegionInFrame` directly as the reconstructed full chain.**
MiXCR pre-assembles the full variable region into `VDJRegionInFrame`. Passing
that single column to Python avoids the 7-region concat in both Tengo and
Python. Larger refactor: drops `_reconstruct_chain` in `pipeline.py`,
removes 7 region columns from the seqTable, simplifies the TSV contract.
Better long-term shape; defer to a focused refactor.

**C. Add `FR4InFrame` as a distinct required feature.**
Would force Python to recognise `FR4InFrame` headers and rebuild around them.
More churn for no benefit — the data is the same, only the label differs.

### Decision

**A.** One-line normalisation, no schema impact, restores full-chain mode on
real data. Revisit B as a clean-up when there is appetite for refactoring the
Tengo↔Python TSV contract.

### Implementation

```tengo
// workflow/src/main.tpl.tengo, inside the vdjCols for-loop, after feat extraction:
if feat == "FR4InFrame" { feat = "FR4" }
```

### References

- Spec sections touched: `README.md` Requirement R4 (region list), `pcolumn-spec.md`
  full-chain section (which assumes `FR4` is present).
- MiXCR feature names verified via `mcp__pl__list_columns` on PFrame
  `cfeb3c1d5363b47b0e9305ec7dbdda18c9c34d211da478ab87c5f6b8ed9759f4` (MiXCR
  Clonotyping output, TinyTrees project).

---

## SD-003: Read receptor from clonotypeKey axis domain

**Status:** applied
**Date:** 2026-04-29
**Affected file:** `workflow/src/main.tpl.tengo` (receptor initialisation)

### Symptom

Block emitted the info message:

> Receptor type not found on input columns (pl7.app/vdj/receptor); defaulting
> to antibody (IG) labels.

…on a TCR Alpha/Beta dataset. Output labels then used antibody conventions
(VH / VL, CDR-H3 / CDR-L3) on TCR data.

### Root cause

MiXCR puts `pl7.app/vdj/receptor` on the **clonotypeKey axis domain** (the
input anchor's secondary axis), not on the per-region sequence column domains.
The original loop checked only `d["pl7.app/vdj/receptor"]` on each column's
`spec.domain`, so the lookup always returned `undefined` and the workflow
defaulted to IG.

The receptor IS on the input anchor itself — confirmed via the
`inputSpec` output:

```
axesSpec[1].domain = {
  "pl7.app/vdj/clonotypingRunId": "...",
  "pl7.app/vdj/receptor": "TCRAB",     ← canonical location
  ...
}
```

### Trigger

Any TCR (or any) MiXCR Clonotyping output. Confirmed on TinyTrees /
`10x-sc-xcr-vdj` running TCRAB.

### Impact (before fix)

User-facing labels were wrong on TCR data: antibody nomenclature (VH, VL,
CDR-H3, CDR-L3) instead of TCR (Vα, Vβ, CDR-α3, CDR-β3). Properties were
computed correctly; only the labels misled.

### Options considered

**A. Read receptor from `keyAxisSpec.domain` before the loop. [chosen]**
Matches MiXCR's actual data shape. Keeps the per-column lookup as a fallback
for any future producer that emits receptor on the column domain instead.

**B. Drop the per-column lookup entirely.**
Cleaner if axis-domain is the only contract. Risks breaking compatibility
with any non-MiXCR producer that emits receptor per-column. Defer.

**C. Require an explicit user override in block args.**
Manual workaround. Bad UX — receptor is determinable from data.

### Decision

**A.** Two-source receptor resolution: axis domain first, per-column second.
Removes a false warning on every MiXCR run and produces correct TCR labels
without spec or schema change.

### Implementation

```tengo
// workflow/src/main.tpl.tengo, before the chain-collection loop:
if keyAxisSpec.domain != undefined {
    axisR := keyAxisSpec.domain["pl7.app/vdj/receptor"]
    if axisR == "IG" || axisR == "TCRAB" || axisR == "TCRGD" {
        receptor = axisR
        receptorSeen = true
    }
}
```

The existing per-column check inside the loop is kept as a fallback.

### References

- Spec sections touched: `README.md` Requirement R13b (receptor detection).
- Receptor-to-label mapping unchanged: `process.tpl.tengo` `labelFragments()`.
- Verified via `mcp__pl__get_block_state` `inputSpec` output on the
  sequence-properties block in TinyTrees project (NG:0x388003).

---

## SD-005: TSV Antibody Schema Omits `receptor_type` Column

**Status:** applied
**Date:** 2026-04-29
**Affected file:** `workflow/src/main.tpl.tengo`, `software/src/pipeline.py`

### Symptom

Spec contract drift: the spec L460 lists `receptor_type` as a TSV column
("literal string `IG`, `TCRAB`, or `TCRGD` — taken from the
`pl7.app/vdj/receptor` domain annotation"). The implementation does not
emit `receptor_type` in the TSV; receptor flows via `plan.json` instead.

### Root cause

`receptor_type` is constant per dataset (every clone shares the same
receptor). Encoding it as a per-row TSV column would duplicate the value
N times — once per clonotype. `plan.json` is the natural carrier for
per-run scalars.

### Trigger

Every antibody/TCR run.

### Impact

None functionally — Python reads receptor from `plan.json` (`run_antibody_tcr`).
The Tengo side passes it through `plan` already.

### Options considered

**A. Carry receptor in `plan.json`. [chosen]**
Matches the per-run shape. Smaller TSV, single source of truth. Already
implemented.

**B. Emit `receptor_type` per row.**
Spec-literal but redundant — N copies of one value.

**C. Carry receptor in both.**
Two sources of truth, drift risk.

### Decision

**A.** plan.json is the right shape for per-run constants.

### References

- Spec section: `README.md` L460 (Antibody mode TSV columns).
- Plan schema: `main.tpl.tengo` (`plan := { mode, receptor, chains, ... }`).

---

## SD-007: Accept `pl7.app/vdj/clonotypeKey` Alongside `pl7.app/vdj/cloneId`

**Status:** applied
**Date:** 2026-04-29
**Affected file:** `model/src/index.ts` (`inputAnchorSpecs`),
`workflow/src/main.tpl.tengo` (`detectMode`)

### Symptom

Spec README enumerates the legacy MiXCR bulk anchor as
`{ axes: [..., { name: "pl7.app/vdj/cloneId" }] }`. The implementation
accepts `pl7.app/vdj/clonotypeKey` as an additional anchor and treats it
as the same modality.

### Root cause

Current MiXCR output emits `pl7.app/vdj/clonotypeKey` for what was
historically `pl7.app/vdj/cloneId`. Accepting both lets the block work on
both archived and current MiXCR clonotyping outputs without a forced
migration on the data side.

### Trigger

Any MiXCR clonotyping output produced after the `cloneId →
clonotypeKey` rename.

### Impact

None — both axis names route to `antibody_tcr_legacy_bulk` modality,
identical downstream handling.

### Options considered

**A. Accept both. [chosen]**
Forward and backward compatible. Two `inputAnchorSpecs` entries; one
extra branch in `detectMode`.

**B. Drop `cloneId`, accept only `clonotypeKey`.**
Breaks ingest of archived MiXCR runs. Avoid until a stated migration
window.

**C. Drop `clonotypeKey`, accept only `cloneId`.**
Breaks every current MiXCR run. Not viable.

### Decision

**A.** Compatibility wins for negligible code complexity.

### References

- Anchor specs: `model/src/index.ts:inputAnchorSpecs`.
- Modality detection: `main.tpl.tengo:detectMode`.

---

## SD-008: Derive Receptor From `pl7.app/vdj/chain` When Receptor Annotation Is Absent

**Status:** applied
**Date:** 2026-05-04
**Affected file:** `workflow/src/main.tpl.tengo` (`chainToReceptor`, axis-domain detection block, per-column fallback loop)

### Symptom

The R13b warning fired on every bulk MiXCR run — `IGHeavy`, `IGLight`,
`TCRAlpha`, `TCRBeta`, `TCRGamma`, and `TCRDelta` anchors all surfaced
"Receptor type not detected on the input dataset; defaulting to antibody
labels." even when the underlying chain identity was unambiguous.

### Root cause

Bulk MiXCR's `clonotypes.byCloneKeyBySample/<chain>/umi-count` columns
expose `pl7.app/vdj/chain` on the `clonotypeKey` axis domain (e.g.
`"IGHeavy"`, `"IGLight"`, `"TCRAlpha"`) but do not stamp
`pl7.app/vdj/receptor`. SD-003 fixed receptor detection for single-cell
runs by reading the receptor key on the axis domain, but bulk runs lack
that key entirely, so detection fell through to the IG default and the
R13b warning fired regardless of whether the chain was IG or TCR.

### Trigger

- Any bulk MiXCR run.
- Project: TinyTrees (`NG:0x388003`), bulk MiXCR block
  `da88a5ef-37e6-4a25-9c9d-fcd4713dc4ee` — observed 2026-05-04 with
  IG Heavy and IG Light anchors.

### Impact (before fix)

- TCR bulk runs misreported as antibody (chain labels "heavy"/"light"
  instead of "alpha"/"beta") because `receptor` defaulted to `IG`.
- The R13b warning showed in every bulk-mode block, including
  unambiguously-IG runs, training users to ignore it.
- TCRGD-specific labelling and the γδ message never fired on bulk γδ
  TCR data.

### Options considered

**A. Derive receptor from chain when receptor key is absent. [chosen]**
The MiXCR chain enum maps unambiguously to a receptor:
`IGHeavy`/`IGLight`/`IGKappa`/`IGLambda` → `IG`,
`TCRAlpha`/`TCRBeta` → `TCRAB`,
`TCRGamma`/`TCRDelta` → `TCRGD`. Adds a small helper plus a fallback
inside the existing axis-domain and per-column receptor blocks.

**B. Require MiXCR to emit `pl7.app/vdj/receptor` on bulk axes.**
Correct long-term; out of seqprops's scope and blocks every existing
bulk MiXCR output.

**C. Suppress the warning on bulk mode.**
Hides the symptom but leaves the receptor wrong (TCR misreported as IG),
breaking γδ labelling and γδ heads-up messages.

### Decision

**A.** Receptor detection precedence is now:
1. Axis-domain `pl7.app/vdj/receptor` (SD-003).
2. Axis-domain `pl7.app/vdj/chain` → derived receptor.
3. Per-column `pl7.app/vdj/receptor` (legacy column-domain check).
4. Per-column `pl7.app/vdj/chain` → derived receptor.
5. Default `IG` + R13b warning when nothing matches.

Behaviour preserved on inputs that DO carry receptor — the explicit
key still takes precedence over the derived one.

### Implementation

`chainToReceptor` helper in `workflow/src/main.tpl.tengo`. Two fallback
inserts: the axis-domain block (around the SD-003 site) and the
per-column loop. `receptorSeen` set when the derivation succeeds, so
the R13b warning only fires when neither receptor nor a recognised
chain is present.

### References

- Spec sections touched: `README.md` Requirement R13b (receptor detection).
- MiXCR chain enum verified via `mcp__pl__query_table` on the bulk QC
  pt (`reports/bulk/clonotypesByChain/{IGHeavy,IGLight,TCRAlpha,TCRBeta,TCRGamma,TCRDelta}`).
- Predecessor: SD-003 (receptor on axis domain for single-cell).

---

## SD-009: Defer R21 Reference Line At GRAVY = 0

**Status:** applied
**Date:** 2026-05-05
**Affected file:** `ui/src/pages/ScatterPage.vue` (line not rendered),
`ui/src/pages/HistogramPage.vue` (line not rendered)

### Symptom

Spec R21 calls for `significantLines: [0]` on the scatterplot axis whenever a
hydrophobicity column is plotted, marking the hydrophobic / hydrophilic divide.
Spec R21a calls for the same on the histogram metric axis when the hook
exists. The implementation ships scatter and histogram panels without the
reference line on either chart.

### Root cause

Graph-maker has no path to inject `significantLines` on a data-column axis
today. Verified at `core/visualizations/packages/graph-maker/src/`:

- `composeScatterplotSettings.ts:applyChartInfoFromAnnotations` only reads
  `Annotation.Graph.Thresholds` from the **grouping** column's spec
  (lines ~82–113), not from the X or Y selected source.
- `getAxesDataFromForms.ts:getAxesDataFromFormsScatterplot` propagates
  `axesFormsData.axisX.significantLinesStyle` to the rendered axis but does
  not carry a `significantLines: number[]` array — that field does not exist
  on `AxesState.axisX/axisY` (`constantsCommon.ts` `AxesState`).
- `composeHistogramSettings.ts` does not consume `significantLines` at all
  (the histogram path has no thresholds wiring).

Spec R21 explicitly notes this: "graph-maker today does not read thresholds
from the X/Y data column directly. A platform-side extension to read
thresholds from data columns is tracked separately and is not part of this
block's spec."

### Trigger

Every scatter + histogram render. The reference line never appears regardless
of which column is selected on which axis.

### Impact

Visual cue at hydrophobic / hydrophilic divide is missing. Properties are
computed and plotted correctly; only the divide marker is absent. Users can
still read the value at zero off the axis ticks.

### Options considered

**A. Defer R21 and R21a entirely. [chosen]**
Skip the line on both panels. No graph-maker changes, no data-model
annotations. Land the panels now; pick up the line when graph-maker grows
the affordance.

**B. Extend `composeScatterplotSettings` to read thresholds from X/Y data
columns and annotate hydrophobicity columns with `Annotation.Graph.Thresholds
= [{value: 0}]` in this block's workflow.**
Reusable by other blocks. Spec carves this out as "platform-side extension...
tracked separately, not in this block's spec" — doing it here expands scope
into `core/visualizations` and needs visualizations-team sign-off. Defer.

**C. Add a block-local `axisInjections` prop to `GraphMaker` so blocks can
pass `significantLines` directly without column annotation.**
Matches spec wording ("block-local scope, not a PColumn annotation") most
literally. New graph-maker API surface, design review on prop shape, same
`core/visualizations` touch as B. Defer.

### Decision

**A.** The spec already permits R21a to defer; extending the same posture to
R21 keeps this block's scope contained. Revisit when the platform-side
extension named in the spec lands, or on explicit ask to scope graph-maker
work into this block.

### Implementation

No code injects `significantLines`. Pages call `GraphMaker` with default
options that select hydrophobicity columns when modality dictates; the
chart renders without the reference line.

### References

- Spec sections touched: `README.md` Requirements R21, R21a; Visualizations
  §Reference line at GRAVY = 0.
- Graph-maker render path verified in
  `core/visualizations/packages/graph-maker/src/utils/createChartSettingsForRender/composeScatterplotSettings.ts`,
  `composeHistogramSettings.ts`, and `getAxesDataFromForms.ts`.
- `AxesState` shape: `core/visualizations/packages/graph-maker/src/constantsCommon.ts`.

---

## SD-010: Chain Letter "A" Is The D-Recombining Chain, Not Alpha/Gamma

**Status:** applied
**Date:** 2026-08-31
**Affected file:** `workflow/src/columns.lib.tengo` (`labelFragments`),
`workflow/src/main.tpl.tengo` (`chainLabel`)

### Symptom

On paired αβ TCR input the CDR3 and full-chain labels named the wrong chain:
the column labelled `CDR-α3 Net Charge (pH 7)` carried the beta chain's
values, and `CDR-β3` carried alpha's. Same inversion on γδ (`CDR-γ3` showed
delta), and in the R11b partial-coverage info messages ("…found for alpha
chain" on beta data).

### Root cause

The spec states the opposite mapping. `README.md` R13 (L210-212) and
`pcolumn-spec.md` (L204-218) both assert that for `TCRAB` chain `"A"` is TRA
(alpha) and for `TCRGD` chain `"A"` is TRG (gamma), with an explicit note:
*"for αβ TCR, chain 'A' is alpha, not beta."* The implementation followed the
spec.

The producer does the reverse. MiXCR orders a receptor's chain slots by
diversity — the D-recombining chain first — and the slot index is what becomes
the `pl7.app/vdj/scClonotypeChain` letter:

```tengo
// mixcr-clonotyping/workflow/src/process.tpl.tengo:39-44
// Chain with higher diversity go first
"IG":    { chains: ["IGHeavy", "IGLight"] },
"TCRAB": { chains: ["TCRBeta", "TCRAlpha"] },
"TCRGD": { chains: ["TCRDelta", "TCRGamma"] }
```

So `"A"` is heavy / beta / delta. The spec's cited source
(`antibody-tcr-lead-selection/workflow/src/utils.lib.tengo`) is itself
inverted for TCR and is the origin of the error. R13 carried a pre-M1 gate
requiring this mapping be verified against real MiXCR output; this is that
verification.

### Trigger

Any paired TCR dataset — αβ or γδ. Confirmed on single-cell TCRAB.

### Impact (before fix)

Labels only. PColumn `name` values and the `pl7.app/vdj/scClonotypeChain`
domain (`"A"`/`"B"`) were correct throughout, so downstream blocks selecting
by spec were unaffected. A user reading the table, or picking a default
scatter axis by label, saw one chain's values attributed to the other.

### Options considered

**A. Follow the producer; correct the labels. [chosen]**
Two label maps change. Names and domains untouched, so no schema impact and
no downstream re-keying.

**B. Follow the spec and leave the labels as they were.**
Requires MiXCR to renumber its chain slots — out of scope, and would break
every block already reading the A/B convention.

**C. Emit the concrete chain name instead of a slot-derived label.**
Removes the ambiguity at the source, but changes every TCR column label and
needs the concrete chain on paired input, which single-cell does not carry
per column.

### Decision

**A.** Three shipped components independently encode A = D-recombining chain:
the producer (`mixcr-clonotyping` above), `import-vdj-data`
(`bare-set-specs.lib.tengo:64-71`, `PAIRED_CHAIN_DOMAIN`), and
`antibody-sequence-liabilities` (`main.tpl.tengo:25-29`, with the same
rationale in a comment). The SDK naming-conventions guide agrees. The spec and
`antibody-tcr-lead-selection` are the outliers.

### Implementation

`labelFragments(receptor, chain)` in `columns.lib.tengo` and `chainLabel(ch)`
in `main.tpl.tengo`. Covered by `workflow/src/columns.test.tengo`, which pins
the slot→label mapping for all three receptors and asserts the chain domain
stays independent of the label.

### References

- Spec sections requiring correction: `README.md` R13 (L210-212), the R13a
  label table (L218), L451, the γδ edge-case row (L587); `pcolumn-spec.md`
  L204-218.
- Producer: `blocks/mixcr-clonotyping/workflow/src/process.tpl.tengo:39-44`.
- Public docs: `docs/docs.platforma.bio/docs/30-sdk/100-vdj-guides/60-naming-conventions.md`.
- Still inverted, tracked separately:
  `antibody-tcr-lead-selection/workflow/src/utils.lib.tengo:702-705`.

---

## SD-011: Table Order Follows Spoken Chain Naming, Not Slot Order

**Status:** applied
**Date:** 2026-08-31
**Affected file:** `workflow/src/columns.lib.tengo` (`displaysFirst`,
`buildCdr3Columns`, `buildFullChainColumns`)

### Root cause

`pl7.app/table/orderPriority` was keyed to the chain *slot* — `"A"` always took
the higher band — and the spec fixes those numbers slot-wise from the IG
perspective (`pcolumn-spec.md:295`: `"67000"  // 66000 for light chain`). For IG
that is invisible, since heavy is both the first-named and the slot-`A` chain.
For TCR the two rules diverge, because slot `A` is the D-recombining chain
(SD-010) — beta, not alpha. The table therefore led with `CDR-β3`, leaking
MiXCR's diversity-first slot assignment into the column order exactly as the
labels did before SD-010.

### Options considered

**A. Order by the receptor's spoken naming. [chosen]** Assign the existing
bands by a (receptor, chain) display rank instead of the slot letter. IG keeps
heavy-then-light; no name, domain or value changes.
**B. Leave ordering on the slot.** Matches the spec's literal numbers, but
keeps a producer artifact in front of the scientist.
**C. Order by slot and rename labels to match.** Rejected — re-introduces the
SD-010 bug.

### Decision

**A.** Slot identity belongs to the producer; the label and the reading order
are the user-facing surface. SD-010 established that for labels, and column
position is the same surface reached a different way. Nothing cross-block is
broken: `import-vdj-data` keys `orderPriority` per region, identical for both
chains (`bare-set-specs.lib.tengo:331`), so its chain order is incidental.

### Deliberately not changed: the default plot axis

`ui/src/utils/scalarColumns.ts` still selects `scClonotypeChain: "A"` for the
default scatter and histogram source, per R19/R20. Table order is a reading
convention; the default plotted chain is an analytical one, and beta/delta carry
the greater CDR3 diversity. The chart labels its own axes, so the divergence is
unambiguous. The R19a/R20a fallback reads *emission* order — unchanged here — so
it also still lands on chain `A`.

### Implementation

`displaysFirst(receptor, chain)` in `columns.lib.tengo`, pinned per receptor
across both bands by `Test_buildColumns_tableOrderFollowsSpokenNaming` in
`workflow/src/columns.test.tengo`.

### References

- Spec requiring correction: `pcolumn-spec.md:295` and the slot-keyed
  orderPriority values through L240-L420.
- Predecessor: SD-010. Default-axis spec: `README.md` R19, R20, R19a, R20a.

---

## SD-012: Derive The Chain Slot From The Locus On Bulk Input

**Status:** applied
**Date:** 2026-09-01
**Affected file:** `workflow/src/columns.lib.tengo` (`chainToSlot`),
`workflow/src/main.tpl.tengo` (per-column chain resolution)

### Root cause

Spec R13 reads chain identity from `pl7.app/vdj/scClonotypeChain`. Bulk MiXCR
does not emit that key at all: it names the locus on the key axis via
`pl7.app/vdj/chain` (`"IGHeavy"`, `"TCRAlpha"`, ...), one locus per dataset.
The code filled the gap by assuming slot `"A"` for every bulk column.

While slot `A` was believed to be alpha (pre-SD-010), that assumption happened
to label bulk alpha datasets correctly and bulk beta ones wrongly. SD-010
corrected the slot semantics, which flipped the victims rather than removing
them: `TCRAlpha`, `TCRGamma`, `IGLight`, `IGKappa` and `IGLambda` inputs were
labelled as their paired partner, and the R11b coverage messages named the
wrong chain with them.

SD-008 already derives the *receptor* from the same `pl7.app/vdj/chain` key.
The locus determines the slot just as unambiguously, so the information needed
was present and discarded.

### Options considered

**A. Derive the slot from the locus. [chosen]** `chainToSlot` maps the eight
MiXCR loci onto the slot MiXCR itself seats them in, per its diversity-first
`receptorInfos` ordering. Unknown or absent loci keep the previous `"A"`
default, so non-MiXCR producers are unaffected.
**B. Emit no chain domain on bulk input.** Truer to the data, since a bulk set
has no pairing, but changes the column shape and breaks consumers keying on it.
**C. Keep assuming `"A"`.** Leaves half of all bulk receptors mislabelled.

### Decision

**A.** The derivation already existed for the receptor (SD-008); extending it
to the slot uses the same key and the same MiXCR rule.

### Emitted domain changes on affected datasets

This is not label-only. On bulk `TCRAlpha`, `TCRGamma`, `IGLight`, `IGKappa`
and `IGLambda` input the emitted `pl7.app/vdj/scClonotypeChain` domain moves
from `"A"` to `"B"`, which changes PColumn identity for those datasets. The new
value is the correct one, and a consumer that matched the old `"A"` was
matching a chain that was never there. Bulk `IGHeavy`, `TCRBeta` and `TCRDelta`
are unchanged, as is all paired single-cell input.

### Implementation

`chainToSlot(chain)` in `columns.lib.tengo`, applied in `main.tpl.tengo` where
`scClonotypeChain` is absent, preferring a per-column `pl7.app/vdj/chain` and
falling back to the key axis. Covered by
`Test_chainToSlot_locusSeatsTheDiverseChainInA`.

### References

- Predecessors: SD-008 (receptor from the same key), SD-010 (slot semantics).
- Slot ordering: `blocks/mixcr-clonotyping/workflow/src/process.tpl.tengo:39-44`.
