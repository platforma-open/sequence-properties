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
