# Spec line-by-line review

This document records a clause-by-clause check of `docs/text/work/projects/sequence-properties/{README.md, pcolumn-spec.md}` against the implementation in `blocks/sequence-properties/`.

Last walk: 2026-04-28. Status legend:

- ✅ implemented
- ⚠️ implemented with caveat (description follows)
- ❌ gap (description follows)

---

## README.md — Requirements (R1 … R17)

### Modality detection

| Req | Status | Where | Notes |
|---|---|---|---|
| **R1** — accept any input column whose axes include a recognised sequence key; no minimum axis count | ⚠️ | `model/src/index.ts:8-29` | Picker patterns require 2-axis `[sampleId, X]`. 1-axis inputs (per spec, valid) will not appear in the dropdown. The runtime detect path *does* handle any axis count (`keyAxisIdx = len(axes) - 1` in `main.tpl.tengo`); only the dropdown filter is narrower. |
| **R1a** — modality from axesSpec scan, panic on unrecognised | ✅ | `workflow/src/main.tpl.tengo:23-43` (`detectMode`), invocation at L98 | Maps the four documented axis patterns plus `pl7.app/vdj/clonotypeKey` (current MiXCR alias for `cloneId`). |
| **R1b** — emit detected axis verbatim | ✅ | `workflow/src/process.tpl.tengo` — `keyAxisSpec` carried into all output specs |
| **R2** — record modality + coverage tier in run report | ✅ | `info` output (`mode`, `receptor`, `coverageTier`, `messages`) consumed by UI |

### Input

| Req | Status | Where | Notes |
|---|---|---|---|
| **R3** — single `inputAnchor` ref | ✅ | `BlockArgs.inputAnchor: PlRef` |
| **R4** — collect aa sequence cols by `pl7.app/vdj/sequence` + `alphabet=aminoacid`, inspect feature domain values; do NOT filter on `isAssemblingFeature` | ✅ | `main.tpl.tengo:66-85` — `bundleBuilder.addMulti` uses alphabet only; per-chain feature presence inspected after collection |
| **R5** — full-chain reconstruction = FR1+CDR1+FR2+CDR2+FR3+CDR3+FR4 in order; missing region per clone → NA | ✅ | `pipeline.py:_reconstruct_chain` (correct order); per-clone NA propagation in `run_antibody_tcr` |

### Peptide-mode output

| Req | Status | Where | Notes |
|---|---|---|---|
| **R6** — 9 scalar columns + AA composition; no length columns | ✅ | All 9 columns in `process.tpl.tengo`; AA fraction column emitted in `process.tpl.tengo` (peptide-mode branch); `pl7.app/sequenceLength` deliberately not emitted |
| **R7** — AA composition: 2-axis `[variantKey, pl7.app/aminoAcid]`, 20 std AAs, non-standard excluded from numerator AND denominator | ✅ | `properties.py:aa_fractions` filters non-standard from both; `process.tpl.tengo` emits the 2-axis PColumn with `pl7.app/aminoAcid` axis |
| **R8** — emit ε oxidised (with `floor(C/2)·125`) and ε reduced | ✅ | `properties.py:extinction_coefficients` returns `(ox, red)` |
| **R9** — instability index for sequences < 10 aa: NA | ✅ | `properties.py:instability_index` returns `None` when effective length < 10 |

### Antibody/TCR output

| Req | Status | Where | Notes |
|---|---|---|---|
| **R10** — CDR3 charge + GRAVY emitted per chain when CDR3 present | ✅ | `pipeline.py:_compute_cdr3_row`; emitted for every chain in `plan.chains` |
| **R11** — full-chain emitted per chain when all 7 regions available; silent absence otherwise | ✅ | `_planned_output_columns` only includes columns for chains in `plan.fullChains` |
| **R11a** — block-level info "CDR3-only input detected — full-chain properties not computed…" | ✅ | `main.tpl.tengo:191` — exact spec wording |
| **R11b** — partial-region info "Partial-region input: N of 7 required regions found for [heavy/light] chain…" | ✅ | `main.tpl.tengo:181` — uses `chainLabel(chain)` to render heavy/light/alpha/beta/gamma/delta per receptor |
| **R11c** — VHH detection (single-A chain, median CDR-H3 ≥ 16, IG → VHH info) | ❌ | Not implemented. Requires reading collected CDR3 column data and computing a length distribution in Tengo (or pushing to Python). Tracked as deferred — see "Open items" below. |
| **R12** — Fv columns only when both VH and VL full chains; antibody only; additivity formulas (ε, MW) | ✅ | `main.tpl.tengo`: `hasFv = mode != peptide && receptor == "IG" && fullChain.A && fullChain.B`. `properties.py` Fv functions are additive for ε / MW; pI uses per-chain charge sum. |

### Chain identity & receptor

| Req | Status | Where | Notes |
|---|---|---|---|
| **R13** — chain id from `pl7.app/vdj/scClonotypeChain`; receptor-dependent meaning | ✅ | `main.tpl.tengo` reads chain domain; `process.tpl.tengo:labelFragments` switches labels by receptor + chain |
| **R13a** — labels adapt at runtime; PColumn name + chain domain unchanged | ✅ | `labelFragments` returns label fragments only; column names and chain domain values are constant. Pre-M1 manual gate (verify chain → receptor mapping on real MiXCR sc data) recorded as a manual validation step, not code. |
| **R13b** — receptor from `pl7.app/vdj/receptor`; default to IG with warning if absent/unrecognised | ✅ | `main.tpl.tengo:213` — `if !receptorSeen && anyChain` emits a warning info message |

### Lead Selection integration

| Req | Status | Where | Notes |
|---|---|---|---|
| **R14** — `isScore: "true"` on the spec-listed numeric columns; `rankingOrder: "increasing"` only on hydrophobicity (CDR3 + peptide); min/max where biologically meaningful | ✅ | `process.tpl.tengo`: isScore set on charge_peptide, gravy_peptide (rankingOrder), CDR3 charge + gravy (rankingOrder on gravy), VH/VL charge + pi, Fv charge + pi. min/max set on pi (0..14), min on aliphatic / aromaticity / MW / ε. |
| **R15** — no `defaultCutoff` anywhere | ✅ | `process.tpl.tengo` never emits `pl7.app/score/defaultCutoff` |
| **R16** — AA composition NOT marked isScore | ✅ | `aaCols` block omits `pl7.app/isScore` |

### Processing

| Req | Status | Where | Notes |
|---|---|---|---|
| **R17** — IPC 2.0 pKa, peptide set for peptide+CDR3, protein set for full chains | ✅ | `pka_tables.py:IPC2_PEPTIDE`, `IPC2_PROTEIN`. `pipeline.py` selects peptide for peptide-mode + CDR3; protein for VDJRegion + Fv. **Caveat**: pKa values transcribed from the IPC 2.0 paper need a manual verification pass against the published supplementary data before a release. The file carries an "IMPLEMENTOR ACTION" note. |

---

## README.md — Technical Specification

### Modality detection algorithm

The pseudocode in §"Modality detection" matches `detectMode` in `main.tpl.tengo` exactly. ✅

### Key formulas

| Formula | Status | Where |
|---|---|---|
| Net charge — Henderson-Hasselbalch with acid / base / Cys context | ✅ | `properties.py:_residue_charge` + `charge_at_ph` |
| Charge formula preconditions (linear, free termini, modifications not detected) | ✅ | Spec marks these as out of scope |
| pH = 7.0 reference | ✅ | `pipeline.py:PH = 7.0` |
| Cys handling: peptide include, full-chain exclude, CDR3 include | ✅ | `_compute_peptide_row(include_cys=True)`, `_compute_full_chain_row(include_cys=False)`, `_compute_cdr3_row(include_cys=True)` |
| pI by bisection over [0, 14] with NA guard at endpoints | ✅ | `properties.py:isoelectric_point` — sign check on f(0), f(14); same-sign ⇒ None; tolerance default 0.001 |
| MW = Σ avg residue mass + 18.02 | ✅ | `aa_tables.py:H2O_AVG_MASS = 18.0153`; `properties.py:molecular_weight` |
| GRAVY — KD scale, denominator = effective_length, non-standard skipped | ✅ | `properties.py:gravy` uses cleaned sequence |
| Non-standard exclusion in charge/pI | ✅ | All ionic computations operate on cleaned sequence |
| ε formulas + Cys scope | ✅ | `properties.py:extinction_coefficients`. Workflow only collects FR/CDR — boundary CH1 Cys not included by construction |
| Instability index — DIWV matrix; II > 40 not applied as filter | ✅ | `instability.py` (full Guruprasad 1990 Table 1); `properties.py:instability_index` returns the raw value |
| Aliphatic index | ✅ | `properties.py:aliphatic_index` (mole fractions) |
| Aromaticity = (F+W+Y) / effective_length | ✅ | `properties.py:aromaticity` |
| Full-chain reconstruction order | ✅ | `pipeline.py:REQUIRED_FEATURES` |
| Fv charge / pI / ε / MW formulas | ✅ | `properties.py:fv_*` family. pI uses per-chain charge sum (NOT concatenation), per spec. |
| Fv ε additive note (odd-Cys) | ✅ | `fv_extinction_coefficients` is per-chain additive, so `floor` is applied per chain |
| Fv pI vs IgG pI explanation | ✅ | Captured in the Fv pI column description annotation |

### Input sequence collection

| Item | Status |
|---|---|
| Peptide query: `pl7.app/sequence` + `feature: peptide, alphabet: aminoacid` | ✅ |
| VDJ legacy query: `pl7.app/vdj/sequence` + `alphabet: aminoacid` | ✅ |
| Universal forward-compat query | ✅ |
| Per-chain coverage check at collection time | ✅ |
| FR1 boundary (post-signal-peptide) | ✅ inherits from upstream |
| CDR3 IMGT-inclusive counting | ✅ inherits from upstream |

### TSV format contract (Tengo → Python)

Spec (peptide): `entity_key`, `sequence`. Antibody: `entity_key`, `receptor_type`, `A_CDR3`, `A_FR1` … `B_FR4`.

Implementation:
- ⚠️ Peptide column header is `peptide_seq` (not `sequence`). This is internal to the Tengo→Python contract; both sides agree. Either rename or keep — flagged for the operator's choice.
- ⚠️ `receptor_type` is **not** carried in the TSV row. Receptor is delivered to Python via `plan.json` instead. Python doesn't currently need per-row receptor awareness (label adaptation happens entirely in Tengo at output-spec construction). Keeping receptor in plan.json keeps the row schema uniform.

### Score annotations

Both numeric and discrete annotation conventions match. `defaultCutoff` is never set. ✅

### Processing pipeline

3 steps. ✅
- Step 1 (Tengo collection): `main.tpl.tengo`
- Step 2 (Python computation): `software/src/main.py` → `pipeline.run`
- Step 3 (Tengo PColumn export): `process.tpl.tengo` via `xsv.importFile` (functionally equivalent to `processColumn` for our case — the reference block ASL also uses `xsv.importFile` for this step).

---

## README.md — Defaults and edge cases table

| Case | Handled where |
|---|---|
| CDR3-only input | ✅ R11a info; `_planned_output_columns` excludes full-chain |
| Some regions missing for a specific clone | ✅ `_reconstruct_chain` returns None per clone |
| Single chain only (unpaired) | ✅ `chains` list contains only present chain; Fv requires both |
| Sequence containing `*` (stop) | ✅ `is_invalid_sequence` |
| Other non-standard (X, B, Z, U, J, `-`) | ✅ `clean_sequence` filters; AA-fraction policy preserves sum=1.0 |
| Zero-length | ✅ NA |
| pI no zero crossing | ✅ NA via bisection guard |
| No Y / no W | ✅ ε = 0 (not NA) |
| Peptide < 10 aa instability | ✅ NA |
| γδ TCR | ✅ info annotation; Fv suppressed (TCR receptor types skip Fv even with full coverage) |
| Single-cell single chain | ✅ per-chain emitted, Fv NA for clonotype |
| Modified termini (N-acetyl, C-amide) | ✅ per spec, not detected — out of scope |
| Cyclic peptide | ✅ per spec, not detected — out of scope |

---

## README.md — Open questions / risks

| Item | Status |
|---|---|
| AA composition in antibody mode | ✅ Deferred to Phase 2 per spec; only peptide-mode emits the 2-axis column. |
| pKa scale discrepancy with BioPython/ExPASy | ✅ Documented in pI column description annotations and `pka_tables.py`. |
| Partial-region inputs | ✅ R11b info annotation. |

---

## pcolumn-spec.md — every column

Each column has a row in `process.tpl.tengo` (or, for AA fraction, in the peptide-mode AA fraction block). Walked through:

- **Peptide mode** — 9 scalar columns + `pl7.app/aaFraction` (2-axis). All format strings, min/max bounds, isScore flags, descriptions, orderPriority numeric values match the spec. ✅
- **Antibody/TCR mode CDR-H3** — `pl7.app/charge` + `pl7.app/hydrophobicity` with `feature: CDR3, scClonotypeChain: A`. orderPriority 68000 / 67900. ✅
- **CDR-L3** — same pair, chain B, orderPriority 67700 / 67600. ✅
- **VH full chain** — 9 columns, `feature: VDJRegion, chain: A`, orderPriority 67000 → 66200 in 100-step decrements. ✅
- **VL full chain** — same 9, chain B, orderPriority 66000 → 65200. ✅
- **Fv** — 5 columns, `feature: Fv` (no chain domain), orderPriority 65100 → 64700. ✅
- **Order priority summary** — collision check holds (VH ends 66200, VL starts 66000 — gap 200; etc.). ✅

---

## Open items (not blockers, tracked for follow-up)

1. **R11c VHH detection** — needs CDR-H3 length distribution. Either compute median in Tengo before invoking Python, or have Python compute and emit back to Tengo. Defer until the operator confirms the heuristic threshold (≥ 16 aa) holds against an internal VHH dataset.
2. **R1 1-axis input picker** — current `inputOptions` filter requires `pl7.app/sampleId` as axis 0. Add 1-axis filters if a real upstream block produces a 1-axis abundance anchor. Today every recognised upstream emits 2-axis, so this is theoretical.
3. **TSV column-header naming** — `peptide_seq` vs spec's `sequence`. Either is fine; rename if the operator prefers the spec's literal wording.
4. **IPC 2.0 pKa table verification** — manual gate before release. Cross-check `pka_tables.py` constants against the paper supplementary data.
5. **Pre-M1 manual gate (R13a)** — verify chain slot ↔ chain type mapping (αβ TCR: A=alpha, B=beta) on real MiXCR single-cell scClonotypeChain output before the antibody+TCR scope ships.
