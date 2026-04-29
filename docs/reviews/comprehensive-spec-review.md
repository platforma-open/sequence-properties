# Comprehensive Spec Review — Sequence Properties Block

**Date:** 2026-04-29
**Reviewer:** Claude (block-dev skill, Opus 4.7 1M)
**Scope:** Full audit of `blocks/sequence-properties` against `docs/text/work/projects/sequence-properties` (`README.md` + `pcolumn-spec.md`).
**Method:** Line-by-line check of every Requirement, every PColumn spec, every milestone acceptance, every Python formula. UI / pframe layout / block layout / upgrade audit follows.

Files reviewed:

- `model/src/{index.ts, dataModel.ts, types.ts}`
- `workflow/src/{main.tpl.tengo, process.tpl.tengo}`
- `ui/src/{app.ts, main.ts, pages/MainPage.vue}`
- `software/src/{main.py, pipeline.py, properties.py, instability.py, pka_tables.py, aa_tables.py, io_layer.py}`
- `software/tests/unit/*.py`, `software/tests/integration/test_corpus_e2e.py`
- `software/tests/data/corpus/*` (manifest + inputs)
- `block/package.json`, `package.json`, `pnpm-workspace.yaml`, `.changeset/icy-walls-flow.md`

Assertion key:

- **PASS** — implementation matches spec.
- **PARTIAL** — implementation matches spec intent; spec wording/strategy not literally followed.
- **DEVIATION** — implementation diverges from spec; may or may not be intentional.
- **MISSING** — required item is absent.

---

## Section 1. README — Overview, Concept, Properties

### Overview Section (L15–32)

| Spec Claim | Status | Notes |
|---|---|---|
| Block accepts any peptide or antibody/TCR sequence dataset | PASS | `inputAnchorSpecs` covers `pl7.app/variantKey`, `pl7.app/vdj/cloneId`, `pl7.app/vdj/clonotypeKey`, `pl7.app/vdj/scClonotypeKey`. |
| Outputs standardised PColumns flowing into Lead Selection | PASS | Score columns carry `pl7.app/isScore: "true"` and stamp `pl7.app/blockId` on the export pframe (`process.tpl.tengo:251–259`). |
| Connects to `lead-selection` | PASS by construction | Score annotations + result-pool export. |
| Useful standalone | PASS | UI exposes a `PlAgDataTableV2` browse view of the property table. |

### Concept Section (L34–69)

Modality auto-detection — spec says: `pl7.app/variantKey` with peptide/extractionRunId domain → peptide; with vdj/clonotypingRunId → antibody_tcr_universal; `pl7.app/vdj/cloneId` → legacy bulk; `pl7.app/vdj/scClonotypeKey` → legacy single-cell.

| Spec mode | Workflow detect path | Status |
|---|---|---|
| `pl7.app/variantKey` + peptide/extractionRunId | `main.tpl.tengo:detectMode` line returning `"peptide"` | PASS |
| `pl7.app/variantKey` + vdj/clonotypingRunId | returns `"antibody_tcr_universal"` | PASS |
| `pl7.app/vdj/cloneId` | returns `"antibody_tcr_legacy_bulk"` | PASS |
| `pl7.app/vdj/clonotypeKey` (NOT in spec) | also returns `"antibody_tcr_legacy_bulk"` | DEVIATION (extension) |
| `pl7.app/vdj/scClonotypeKey` | returns `"antibody_tcr_legacy_sc"` | PASS |

> The `clonotypeKey` alias is documented as "alias seen in real MiXCR output" in `test/src/wf.test.ts`. Confirms an undocumented MiXCR-side fact and looks like a legitimate operational extension.

**DEVIATION (axis selection):** `main.tpl.tengo:115` does `keyAxisIdx := len(axes) - 1` — picks the **last** axis. Spec R1a says "scanning `inputSpec.axesSpec` for the **first** axis whose name + domain match". For 2-axis inputs `[sampleId, variantKey]` the result is identical (key axis is last). For higher-axis inputs (e.g. `[sampleId, replicateId, variantKey]`) "first matching axis" and "last axis" diverge. **Functionally correct for current upstream blocks; deviates from spec wording.**

Output granularity ("All output columns are 1-axis (per-sequence, not per-sample)") — PASS. Verified per-column construction in `process.tpl.tengo`.

pH reference 7.0 — `pipeline.py:PH = 7.0`. PASS.

### Properties Section (L71–112) — Shared Table

Properties listed in the shared table (L77–88):

| Property | Method spec | Implementation | Status |
|---|---|---|---|
| Net charge (pH 7) | Henderson-Hasselbalch + IPC 2.0 pKa | `properties.py:charge_at_ph` | PASS |
| Hydrophobicity (GRAVY) | Kyte-Doolittle | `properties.py:gravy` over `KD_SCALE` | PASS |
| Molecular weight | Σ avg residue mass + H₂O | `properties.py:molecular_weight`, `H2O_AVG_MASS = 18.0153` | PASS |
| Isoelectric point | Bisection at net charge = 0 | `properties.py:isoelectric_point` + `_bisect_zero` | PASS |
| Extinction coeff oxidised | Y·1490 + W·5500 + floor(C/2)·125 | `properties.py:extinction_coefficients` | PASS |
| Extinction coeff reduced | Y·1490 + W·5500 | same fn, second tuple element | PASS |
| Instability index | Guruprasad 1990 dipeptide composition | `properties.py:instability_index` + DIWV in `instability.py` | PASS |
| Aliphatic index | Ala + 2.9·Val + 3.9·(Ile+Leu) fractions × 100 | `properties.py:aliphatic_index` | PASS |
| Aromaticity | Fraction (F+W+Y) | `properties.py:aromaticity` | PASS |
| AA composition (peptide only, 2-axis) | per residue fraction | `properties.py:aa_fractions` + `process.tpl.tengo` aaSpecs | PASS |

### Antibody/TCR-mode-only Table (L93–110)

| Spec column | Implementation | Status |
|---|---|---|
| CDR-H3/α3 net charge | `_compute_cdr3_row → charge_at_ph(IPC2_PEPTIDE, include_cys=True)` | PASS |
| CDR-H3/α3 hydrophobicity | `_compute_cdr3_row → gravy` | PASS |
| CDR-L3/β3 net charge / hydrophobicity | same per-chain loop in `_compute_row_for` | PASS |
| Fv net charge | `properties.py:fv_charge` (per-chain sum, not concatenation) | PASS |
| Fv pI | `fv_isoelectric_point` (bisects per-chain charge sum) | PASS |
| Fv ε oxidised | `fv_extinction_coefficients` (additive) | PASS |
| Fv ε reduced | same | PASS |
| Fv MW | `fv_molecular_weight` (additive, two H₂O terms) | PASS |

scFv exclusion (L111) — PASS by virtue of nothing detecting scFv; falls through the VHH heuristic per spec direction.

TCR vs antibody naming (L113–116) — PASS. `process.tpl.tengo:labelFragments` correctly maps `IG/A → CDR-H3, VH`, `TCRAB/A → CDR-α3, Vα`, `TCRGD/A → CDR-γ3, Vγ`, etc.

---

## Section 2. README — Requirements (R1–R17)

### R1 — Block accepts any input column with a recognised sequence-key axis

**Status: PASS.**

`inputAnchorSpecs` in `model/src/index.ts:11–32` lists 4 axis patterns; each entry has 2 axes `[sampleId, key]`. `bb.addAnchor("main", args.inputAnchor)` and downstream column queries use `axes: [{anchor: "main", idx: 1}]` so the block targets the second axis (the entity key). No minimum axis count violation.

### R1a — Modality determined by scanning axesSpec for first matching axis

**Status: DEVIATION (selection direction).** **Status: PASS (matching pattern).**

Code (`main.tpl.tengo:113–118`):

```go
keyAxisIdx := len(axes) - 1
keyAxisSpec := axes[keyAxisIdx]
mode := detectMode(keyAxisSpec)
```

The implementation picks the last axis instead of scanning forward. Equivalent for 2-axis inputs, divergent for ≥3-axis. The error message wording on no match (`main.tpl.tengo:121`) matches the spec verbatim: `"no recognized sequence key axis found; connect a peptide extraction or MiXCR clonotyping dataset"`.

The 4 supported axis patterns match the spec; one extra (`pl7.app/vdj/clonotypeKey`) is an extension covering current MiXCR aliases.

### R1b — Detected input axis emitted verbatim on all output columns

**Status: PASS.**

`process.tpl.tengo:53` builds `axes := [{ column: "entity_key", spec: keyAxisSpec }]` from the resolved input axis spec. The same `keyAxisSpec` is reused for the AA-fraction PColumn's first axis (`process.tpl.tengo:303`). No renaming or re-domaining anywhere.

### R2 — Detected modality + coverage tier recorded in run report

**Status: PASS.**

`main.tpl.tengo:235–241` emits `infoBlob` with `mode`, `receptor`, `coverageTier ∈ {peptide, full_chain, cdr3_only, partial}`, and `messages: [...]`. Surfaced to the UI through `outputs.info` and rendered as `PlAlert` rows in `MainPage.vue:31–36`.

### R3 — Single inputAnchor reference

**Status: PASS.**

`BlockArgs` in `model/src/types.ts:11–14`:

```ts
export type BlockArgs = {
  inputAnchor: PlRef;
  defaultBlockLabel: string;
};
```

Single anchor per spec.

### R4 — Collect aa sequence columns with `pl7.app/alphabet: "aminoacid"`

**Status: PASS.**

`main.tpl.tengo:79–104`:
- peptide path: name `pl7.app/sequence`, domain `{pl7.app/feature: "peptide", pl7.app/alphabet: "aminoacid"}`.
- VDJ legacy path: name `pl7.app/vdj/sequence`, domain `{pl7.app/alphabet: "aminoacid"}` — does **not** filter on `isAssemblingFeature` (R4 explicit).
- Universal forward-compat path: name `pl7.app/sequence`, domain `{pl7.app/alphabet: "aminoacid"}`.

Per-chain coverage assessed independently in `main.tpl.tengo:178–198` by scanning `chainsFound[chain]` for the 7 required features. Looks correct.

### R5 — Full-chain reconstruction by FR1→CDR1→FR2→CDR2→FR3→CDR3→FR4 concatenation

**Status: PASS.**

`pipeline.py:_reconstruct_chain` walks `REQUIRED_FEATURES = ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")` and concatenates. Returns `None` if any region is empty for that clone — properly NA-flags downstream computation.

### R6 — 9 scalar peptide PColumns + AA composition; no length

**Status: PASS.**

`process.tpl.tengo:65–162` emits exactly the 9 columns (`charge_peptide`, `gravy_peptide`, `mw_peptide`, `pi_peptide`, `eox_peptide`, `ered_peptide`, `instability_peptide`, `aliphatic_peptide`, `aromaticity_peptide`) plus the AA fraction PColumn. No `pl7.app/sequenceLength` is emitted. Verified by inspection of all emitted column names.

### R7 — AA composition: 2-axis, 20 standard AAs, non-standard excluded from numerator and denominator

**Status: PASS.**

- 2-axis: `process.tpl.tengo:303–315` declares `[entity_key, pl7.app/aminoAcid]`.
- 20 standard codes: `aa_tables.py:STANDARD_AAS = "ACDEFGHIKLMNPQRSTVWY"`.
- Non-standard exclusion: `properties.py:aa_counts` filters by `STANDARD_AA_SET`; `aa_fractions` uses the count sum as denominator. Verified by `tests/unit/test_properties.py::TestAaFractions::test_nonstandard_dropped` and `test_fractions_sum_to_one`.
- Pipeline emits one row per (entity, aa) pair so each variantKey gets exactly 20 rows (`pipeline.py:run_peptide`, lines 119–134); enforced by `tests/integration/test_corpus_e2e.py::test_peptide_aa_fraction_shape`.

### R8 — Extinction coefficient: oxidised + reduced (two columns)

**Status: PASS.**

Both emitted in peptide and full-chain paths, with the disulfide bonus = `floor(C/2) × 125` only on oxidised (`properties.py:extinction_coefficients`).

### R9 — Instability index < 10 aa: NA

**Status: PASS.**

`properties.py:instability_index` line: `if cleaned is None or len(cleaned) < 10: return None`. **Effective length** is used (post-non-standard-strip) — verified by `test_floor_uses_effective_length`. Correct spec interpretation.

### R10 — CDR3 charge + hydrophobicity for both chains when present

**Status: PASS.**

`process.tpl.tengo:161–187` iterates `chains` (sorted) and emits `charge_<chain>_CDR3` + `gravy_<chain>_CDR3`. CDR3 length is **not** re-emitted — spec compliance.

### R11 — Full-chain emitted per chain only when all 7 regions available

**Status: PASS.**

Workflow gates: `present == len(REQUIRED_FEATURES)` → `fullChain[chain] = true` (`main.tpl.tengo:184–186`). `process.tpl.tengo:189–264` iterates only `fullChains`. Per-clone reconstruction in `_compute_row_for` returns NA for that clone when `_reconstruct_chain` returns `None`.

### R11a — CDR3-only info annotation, exact wording

**Status: PASS (exact).**

`main.tpl.tengo:200–202`:
> "CDR3-only input detected — full-chain properties not computed. To enable them, use a MiXCR preset that exports all VDJ regions."

Verbatim with spec.

### R11b — Partial-region info annotation, "[heavy/light]" chain label

**Status: PARTIAL.**

Spec wording: `"…of 7 required regions found for [heavy/light] chain — full-chain properties not computed…"`. Code (`main.tpl.tengo:188–190`) substitutes the receptor-aware chain label via `chainLabel(chain)` (heavy/light, alpha/beta, gamma/delta) and otherwise mirrors the spec string. **Better than spec literal** — spec only listed the antibody case; code generalises to TCR.

**MINOR DEVIATION:** When a chain has 1–6 regions but **CDR3 is absent**, no info message is emitted (`main.tpl.tengo:185–192` else-if chain triggers only when `feats["CDR3"]` is truthy). Spec wording says "1–6 of the 7 required regions are present"; CDR3-absent partial coverage is silently dropped. Edge case unlikely with current MiXCR presets but worth noting.

### R11c — VHH/single-domain antibody detection (≥16 aa median CDR-H3, IG, chain A only)

**Status: MISSING.**

No code path computes CDR-H3 length distributions or emits the VHH info annotation. The Tengo workflow does not pass CDR-H3 length data to a downstream summarisation step; the Python pipeline doesn't compute medians. Confirmed by full grep of `main.tpl.tengo`, `process.tpl.tengo`, `pipeline.py` for "VHH", "16", "median", "single-domain" — **zero matches**.

This is a user-facing requirement per spec ("This is a user-facing requirement, not optional" wording in R11a applies to R11a; R11c is described as "after CDR-H3 length computation (Step 2)…emit a block-level info annotation"). **Not implemented.**

### R12 — Fv: charge, pI, ε ox/red, MW; antibody only; both chains full

**Status: PASS.**

`main.tpl.tengo:213`:

```go
hasFv := mode != "peptide" && receptor == "IG" && fullChain["A"] && fullChain["B"]
```

Excludes TCR (correct), gates on both chains full (correct). Fv columns emitted in `process.tpl.tengo:266–306`. Five columns: charge, pi, eox, ered, mw. Additivity in Python: `fv_charge`, `fv_extinction_coefficients`, `fv_molecular_weight` are explicit additive functions; `fv_isoelectric_point` bisects on per-chain sum (NOT concatenated string) — matches spec.

### R13 — Chain identity from `pl7.app/vdj/scClonotypeChain`

**Status: PASS.**

`main.tpl.tengo:154`: `chain := d["pl7.app/vdj/scClonotypeChain"]`. With one extension: when chain is undefined / empty (bulk MiXCR data without chain annotation), assumes chain "A" (`main.tpl.tengo:155–157`). Reasonable extension but not in spec.

### R13a — Labels adapt to receptor at runtime; column name + chain domain unchanged

**Status: PASS.**

`process.tpl.tengo:13–28` `labelFragments` returns the right (cdr3, fullChain) pair for each (receptor, chain) tuple:

| Receptor | A | B |
|---|---|---|
| IG | (CDR-H3, VH) | (CDR-L3, VL) |
| TCRAB | (CDR-α3, Vα) | (CDR-β3, Vβ) |
| TCRGD | (CDR-γ3, Vγ) | (CDR-δ3, Vδ) |

Column `name` stays `pl7.app/charge` etc. Domain `pl7.app/vdj/scClonotypeChain` stays "A" / "B".

γδ TCR info annotation (`main.tpl.tengo:206–208`):
> "γδ TCR input detected — displaying with γδ-specific labels; Fv columns are not computed for TCR inputs."

Matches spec exactly.

**Note:** Spec calls out a CDR3 length description that's modality-aware. Since this block does **not** emit CDR3 length (MiXCR upstream does), this requirement does not apply here.

### R13b — Receptor type from `pl7.app/vdj/receptor` domain; default IG + warning

**Status: PASS, with extension SD-003 noted in deviations.**

Two paths:
1. Read from input anchor's axis domain (`main.tpl.tengo:127–134`) — extension noted in code as SD-003. MiXCR puts receptor on the clonotypeKey axis domain.
2. Fallback: per-column scan inside the VDJ loop (`main.tpl.tengo:166–172`).

Default to "IG" if unseen. Warning info message (`main.tpl.tengo:209–211`):
> "Receptor type not found on input columns (pl7.app/vdj/receptor); defaulting to antibody (IG) labels. Use a MiXCR preset that emits the receptor annotation if this is a TCR dataset."

Reasonable extension of "log a warning" from spec.

### R14 — `isScore: "true"` on the spec-listed columns; `rankingOrder: "increasing"` only on hydrophobicity

**Status: PASS.**

Verified column-by-column from `process.tpl.tengo`:

| Column | Expected isScore | Actual | Expected rankingOrder | Actual |
|---|---|---|---|---|
| Peptide charge | yes | yes | — | — |
| Peptide hydrophobicity | yes | yes | increasing | increasing |
| Peptide MW / pI / eox / ered / instability / aliphatic / aromaticity | no | no | — | — |
| CDR3 charge (A,B) | yes | yes | — | — |
| CDR3 hydrophobicity (A,B) | yes | yes | increasing | increasing |
| VH/VL charge | yes | yes | — | — |
| VH/VL pI | yes | yes | — | — |
| VH/VL hydrophobicity | no | no | — | — |
| VH/VL mw / eox / ered / instability / aliphatic / aromaticity | no | no | — | — |
| Fv charge | yes | yes | — | — |
| Fv pI | yes | yes | — | — |
| Fv eox / ered / mw | no | no | — | — |
| AA fraction | no | no | — | — |

Min/max bounds: pI columns have `min=0 max=14`, MW / EC / aliphatic / aromaticity have `min=0`, aromaticity also `max=1`. Matches `pcolumn-spec.md`.

### R15 — No `defaultCutoff` on any column

**Status: PASS.**

Grep of `process.tpl.tengo` for `defaultCutoff` returns zero matches.

### R16 — AA composition columns NOT marked isScore

**Status: PASS.**

`process.tpl.tengo:303–323`'s aaCols entry has no `pl7.app/isScore` key.

### R17 — IPC 2.0 peptide pKa for peptide + CDR3; protein pKa for VH/VL

**Status: PASS.**

`pipeline.py:_compute_peptide_row` uses `IPC2_PEPTIDE`. `_compute_cdr3_row` uses `IPC2_PEPTIDE`. `_compute_full_chain_row` uses `IPC2_PROTEIN`. `_compute_fv_row` uses `IPC2_PROTEIN`. `pka_tables.py` defines both sets with `n_terminus`, `c_terminus`, and per-residue pKa. Cys handling is correct: `include_cys=True` for peptide and CDR3, `include_cys=False` for full-chain and Fv.

**Caveat — pKa value verification:** The pKa tables say "NOTE FOR IMPLEMENTORS: these constants were transcribed from the published IPC 2.0 supplementary data. Verify each value against the paper before shipping a release." This release gate is documented but I cannot independently verify the transcription against the paper from inside this review.

---

## Section 3. README — Technical Specification

### Modality detection block (L213–242)

Pseudocode in spec uses **"break"** (first-match-then-break). Implementation does last-axis selection (described above as DEVIATION on R1a). For the supported 2-axis upstream blocks the result is identical.

### Key formulas (L242–340)

#### Charge formula (L243–262)

Acid: `−1/(1 + 10^(pKa − pH))`. Base: `+1/(1 + 10^(pH − pKa))`. **Code (`properties.py:_residue_charge`) matches both signs and forms exactly. PASS.**

Independent verification (matching the spec example for glycine at pH 7 with peptide pKa):
- N-term: `1/(1+10^(7−9.094)) = 0.99205`
- C-term: `−1/(1+10^(2.869−7)) = −0.99993`
- Net: ~−0.008. The block's `charge_at_ph("G", 7.0, IPC2_PEPTIDE)` returns `−0.007915` — matches.

Charge formula preconditions documented in spec (linear, free termini, etc.) — implementation does not detect modifications, consistent with spec acknowledging this as out of scope.

#### Cysteine handling (L249–254)

Three regimes, all correctly mapped:
- Peptide: include_cys=True (free thiol) — `_compute_peptide_row`. PASS.
- Full chain: include_cys=False — `_compute_full_chain_row`. PASS.
- CDR3: include_cys=True — `_compute_cdr3_row`. PASS.

#### Isoelectric point (L266–267)

Spec says: bisection [0, 14], tol 0.001, NA guard via same-sign check before bisecting. **PASS — `_bisect_zero` checks `math.copysign(1.0, f_lo) == math.copysign(1.0, f_hi)` and returns None.** Tolerance is the default `0.001`.

**Verification corner case (spec L266):** "polybasic synthetic" example "RRRRRRRRRR" emitting NA — **technically incorrect example.** With IPC 2.0 peptide pKa_R=11.84 and free C-terminus, RRRRRRRRRR has pI ≈ 12.79 (charge crosses zero in [0,14]). The implementation correctly returns ~12.79 for this sequence, not None. The NA path **is** verified by `test_pi_polybasic_no_zero_crossing_returns_none` using `R*50` with the protein pKa set, which does saturate. The feature works; the spec's M1 acceptance example is a misnomer (see Milestone audit below).

#### Molecular weight (L269–276)

`MW = Σ avg_residue_mass + 18.02`. Implementation: `sum(AVG_RESIDUE_MASS[c] for c in cleaned) + H2O_AVG_MASS`. `H2O_AVG_MASS = 18.0153` (more precise than spec's 18.02 — fine, more accurate). PASS.

Spot-check: GGGGGG = 6 × 57.0519 + 18.0153 = 360.3267 — block returns 360.3267. Match.

#### GRAVY + non-standard handling (L278–289)

Both numerator and denominator exclude non-standard residues — matches BioPython behaviour. `properties.py:gravy` operates on `cleaned` (post `_prepare`). PASS.

Spec verbatim: stop codon `*` invalidates whole sequence; non-standard residues `B, Z, X, U, J, -` excluded from sum and effective_length. Verified by `test_clean_sequence` and `test_invalid_yields_none`.

#### Non-standard residues in charge / pI (L291–293)

"Apply the same exclusion policy" — code does. `_residue_charge` returns 0.0 when `pka_set.get(aa) is None` (which is true for non-standard codes), so they don't contribute. AND `clean_sequence` strips them before iteration. PASS.

#### Extinction coefficient (L295–305)

```
ε_ox  = Y·1490 + W·5500 + floor(C/2)·125
ε_red = Y·1490 + W·5500
```

Code: `aa_tables.py:EC_TYR=1490, EC_TRP=5500, EC_DISULFIDE=125`. `extinction_coefficients` line: `eps_ox = eps_red + (c // 2) * EC_DISULFIDE`. PASS.

Sequences with no Y or W → 0 (not NA): `test_no_aromatic_yields_zero_not_na`. PASS.

Cys scope note: code counts only Cys residues in the sequence passed in; for reconstructed VH this excludes the boundary CH1 Cys (out of `FR1+...+FR4`). Spec calls for verification on a known VH sequence — the golden test (VH=…) pins ε_ox at 22460 (`test_golden_values.py:VH eox`). This is a regression pin, not an external reference cross-check (which spec defers to M3).

#### Instability index (L307–319)

`II = (10/L) × Σ DIWV(aa_i, aa_{i+1})` for i in 1..L−1. Code (`properties.py:instability_index`):

```python
for i in range(n - 1):
    v = diwv(cleaned[i], cleaned[i + 1])
    if v is None: continue
    total += v
return (10.0 / n) * total
```

`n = len(cleaned)`, range is L−1 dipeptides. Formula matches. PASS.

Closed-form verification: poly-A 10 → every dipeptide AA = 1.0 → sum = 9 → II = (10/10) × 9 = 9.0. Test `test_polyalanine_known_value` asserts 9.0 — passes.

#### Aliphatic index (L321–326)

`AI = 100 × (X_A + 2.9·X_V + 3.9·(X_I + X_L))`. Code (`aliphatic_index`):

```python
return 100.0 * (x_a + 2.9 * x_v + 3.9 * (x_i + x_l))
```

PASS.

Verification: `test_alvi_equal_fractions` — ALVI = 100 × (0.25 + 2.9·0.25 + 3.9·0.5) = 292.5 — passes.

#### Aromaticity (L328–331)

Denominator = `effective_length`. Code: `arom / n` where `n = len(cleaned)`. PASS.

#### Full-chain reconstruction (L333–337)

`VH = FR1+CDR1+FR2+CDR2+FR3+CDR3+FR4`. Code: `_reconstruct_chain` walks `REQUIRED_FEATURES` in that order. PASS.

#### Fv properties (L339–365)

- `fv_charge = charge(VH) + charge(VL)` — `properties.py:fv_charge`. PASS.
- `fv_pI = bisect[charge(VH, pH) + charge(VL, pH) = 0]` — `fv_isoelectric_point` bisects on per-chain sum. **NOT concatenated.** PASS.
- `fv_ε_ox = ε_ox(VH) + ε_ox(VL)` — additive. PASS.
- `fv_MW = MW(VH) + MW(VL)` — each chain contributes one H₂O. PASS.

The spec's Fv-pI rationale (per-chain sum counts 2 N-termini and 2 C-termini correctly) is exactly the implemented semantics. Verified by direct calculation: `fv_charge_perchainsum @ pH=7` ≠ `concatenated chain charge @ pH=7` for a synthetic "DEKR/DEKR" Fv (~0.003 difference at pH 7), confirming the distinction is not vacuous.

### Input sequence collection (L367–399)

Peptide path query — PASS (covered above in R4).
VDJ legacy path query — PASS.
FR1 boundary note — verified by spec and golden test (VH MW = 6050.7302); not externally cross-validated.
CDR3 length convention — IMGT inclusive of anchor residues — block does not emit length, so out of scope here.

### PColumn output summary (L401–445)

Tables in this section duplicate the shape in `pcolumn-spec.md`. Cross-checked column-by-column against `process.tpl.tengo`. Match.

### Score annotations for Lead Selection (L447–462)

`isScore`, `rankingOrder` semantics, no `defaultCutoff` rule — all verified in R14, R15.

### Processing pipeline (L464–495)

**Step 1 (Tengo):** anchor resolution + modality detection + per-entity sequence TSV. Implemented (`main.tpl.tengo:wf.body`). PASS.

**TSV format contract** — peptide columns `entity_key, peptide_seq` (note: spec says `entity_key, sequence`; implementation says `peptide_seq`). DEVIATION (column name).

> Spec L457: "Peptide mode columns: `entity_key`, `sequence`"
> Implementation `main.tpl.tengo:131`: header `peptide_seq`
> Python `pipeline.py:run_peptide`: `reads["peptide_seq"]`

The Tengo and Python sides are consistent. Spec's wording differs. Functionally fine but minor spec drift.

Antibody columns `entity_key, A_FR1, A_CDR1, A_FR2, A_CDR2, A_FR3, A_CDR3, A_FR4, B_FR1, …` — implementation matches. Note: spec also lists `receptor_type` as a column; implementation does **not** emit `receptor_type` as a TSV column — instead passes receptor in `plan.json`. **DEVIATION (TSV schema).** Functionally equivalent (Python pipeline reads receptor from plan, not TSV), but spec was explicit.

> Spec L460: "`receptor_type`: literal string `"IG"`, `"TCRAB"`, or `"TCRGD"` — taken from the `pl7.app/vdj/receptor` domain annotation"
> Implementation: `plan.json` carries receptor; TSV does not.

Missing region empty-string convention — PASS. Polars reads empty as null, treated as absent.

**Step 2 (Python):**

- "Read sequence TSV" — PASS (`io_layer.read_input_tsv`).
- "Normalise all sequences to uppercase" — PASS (`properties.py:_prepare → clean_sequence` upper-cases).
- "Compute all 9 scalar properties + AA composition per row" (peptide) — PASS.
- "Antibody mode: compute CDR-H3 properties for all clones; attempt full-chain reconstruction per clone; compute full-chain and Fv properties where reconstruction succeeds" — PASS.
- **"Vectorize using NumPy: represent pKa lookups as arrays, compute charge over the full sequence batch in parallel, run bisection on a vector of sequences simultaneously."**

> **DEVIATION (performance approach).** `pipeline.py` uses pure-Python list comprehensions and per-row dict operations. `polars` is used for IO/schema management, not for vectorised arithmetic. Property functions (`properties.py`) operate on single strings — not array-batched.
>
> Spec performance estimate: "<60 seconds with vectorised implementation" for typical single-cell. Current scalar implementation is likely well above that on 10⁶-clone datasets. May or may not be in scope of an existing perf-review (that file is in `docs/reviews/` but excluded from this review per user instructions).

**Step 3 (Tengo):**

Spec: "`pframes.processColumn` per output column with score annotations". Implementation: uses `xsv.importFile` + `pframes.pFrameBuilder`. **PARTIAL/CLARIFICATION.** `processColumn` is for per-group computation on an existing PColumn (not what this block does — Python writes a TSV, Tengo imports it). `xsv.importFile` is the correct primitive. Spec wording is imprecise; implementation is correct.

### Defaults and edge cases (L479–498)

| Case | Spec | Implementation | Status |
|---|---|---|---|
| CDR3-only input | CDR-H3 emitted; full-chain absent | `fullChain` empty → no full-chain columns; CDR3 emitted | PASS |
| Some regions missing for a clone | Full-chain NA per-clone | `_reconstruct_chain` returns None → NA | PASS |
| Single chain only | Per-chain emitted; Fv absent | `hasFv` requires both A and B full | PASS |
| `*` stop codon | All NA | `is_invalid_sequence` covers | PASS |
| Non-standard `X B Z U J -` | Excluded from numerator + denominator | `clean_sequence` strips; `aa_counts` filters | PASS |
| Zero-length | All NA | `_prepare` returns None | PASS |
| pI no zero crossing | NA | `_bisect_zero` same-sign check | PASS |
| No Tyr/Trp | ε = 0 (not NA) | `extinction_coefficients` returns (0, 0) | PASS |
| Peptide < 10 aa instability | NA | `instability_index` length floor | PASS |
| γδ TCR | γδ labels, no Fv | TCRGD path in `labelFragments`, hasFv stays false | PASS |
| Single-cell single chain | per-chain only, Fv NA | hasFv false → no Fv columns | PASS |
| Modified termini | Not detected (acknowledged) | Not detected | PASS (per spec) |
| Cyclic peptide | Not detected (acknowledged) | Not detected | PASS (per spec) |

---

## Section 4. PColumn Spec — Per-Column Audit

### Order priorities (`pcolumn-spec.md` L549–559)

| Range | Spec | Implementation |
|---|---|---|
| Peptide 70000–69200 | charge 70000, hydrophobicity 69900, MW 69800, pI 69700, eox 69600, ered 69500, instability 69400, aliphatic 69300, aromaticity 69200, AA fraction 69000 | All match (`process.tpl.tengo:65–162`, AA fraction at 69000 in `process.tpl.tengo:317`) |
| CDR-H3 68000/67900 (A) | Charge 68000, hydrophobicity 67900 | Match (`cdr3OrderA = 68000; gravyOrder = chargeOrder − 100`) |
| CDR-L3 67700/67600 (B) | Charge 67700, hydrophobicity 67600 | Match (`cdr3OrderB = 67700`) |
| VH 67000–66200 | charge 67000, pI 66900, gravy 66800, mw 66700, eox 66600, ered 66500, instability 66400, aliphatic 66300, aromaticity 66200 | Match (`fcOrderBaseA = 67000; base − {0,100,…,800}`) |
| VL 66000–65200 | All −1000 from VH | Match (`fcOrderBaseB = 66000`) |
| Fv 65100–64700 | charge 65100, pI 65000, eox 64900, ered 64800, mw 64700 | Match (`process.tpl.tengo:266–305`) |

**No collisions or missing slots.** Phase-2 reservations 64200–64600 are not used (correct per spec).

### Annotation parity per column (peptide)

| Column | Spec annotations | Code annotations | Diff |
|---|---|---|---|
| `pl7.app/charge` | label, format, isScore, **description**, visibility, orderPriority | label, format, isScore, **shorter description**, visibility, orderPriority | description abbreviated |
| `pl7.app/hydrophobicity` | label, format, isScore, rankingOrder, **description**, visibility, orderPriority | label, format, isScore, rankingOrder, **shorter description**, visibility, orderPriority | description abbreviated |
| `pl7.app/molecularWeight` | label, format, **min**, visibility, orderPriority | label, format, min, visibility, orderPriority | identical |
| `pl7.app/isoelectricPoint` | label, format, **min, max**, visibility, orderPriority | label, format, min, max, visibility, orderPriority | identical |
| `pl7.app/extinctionCoefficientOx` | label, format, min, **description**, visibility, orderPriority | label, format, min, **shorter description**, visibility, orderPriority | description abbreviated |
| `pl7.app/extinctionCoefficientRed` | label, format, min, **description**, visibility, orderPriority | label, format, min, visibility, orderPriority | **MISSING description** |
| `pl7.app/instabilityIndex` | label, format, **description**, visibility, orderPriority | label, format, **shorter description**, visibility, orderPriority | description abbreviated |
| `pl7.app/aliphaticIndex` | label, format, min, **description**, visibility, orderPriority | label, format, min, visibility, orderPriority | **MISSING description** |
| `pl7.app/aromaticity` | label, format, min, max, **description (long)**, visibility, orderPriority | label, format, min, max, **shorter description**, visibility, orderPriority | description abbreviated |
| `pl7.app/aaFraction` | label, format, min, max, visibility, orderPriority (no description in spec) | label, format, min, max, visibility, orderPriority | identical |

**MISSING:** Peptide ε reduced description and peptide aliphatic description. Both are user-facing tooltip content.

### Annotation parity per column (CDR-H3 / CDR-L3)

| Column | Spec annotation diff | Status |
|---|---|---|
| `pl7.app/charge` (CDR-H3, chain A) | spec has long description; CDR-L3 (chain B) has a *different* description per spec | Code uses one description for both A and B — **DEVIATION (description per chain not differentiated)** |
| `pl7.app/hydrophobicity` (CDR-H3) | spec has long description; CDR-L3 (chain B) has *different* description per spec | Same — **DEVIATION** |

The spec table at L268–272 lists CDR-L3-specific descriptions ("Strongly positive CDR-L3 charge contributes to paratope polyreactivity. …" / "Same aggregation and polyreactivity signal as CDR-H3 hydrophobicity; lower independent predictive weight…"). Code emits an identical description for both chains in the loop body (`process.tpl.tengo:166–187`). Minor.

### Annotation parity per column (VH/VL VDJRegion)

All 9 columns × 2 chains. Spec says "All annotations carry over unchanged from the VH versions" for VL — **PASS** (loop emits identical annotations, only label differs via `frag.fullChain`).

VH annotation diffs vs spec:

| Column | Status |
|---|---|
| charge | description abbreviated |
| isoelectricPoint | spec has no description; code has no description ✓ |
| hydrophobicity | description abbreviated |
| molecularWeight | description matches semantically |
| extinctionCoefficientOx | spec has no description; code has no description ✓ |
| extinctionCoefficientRed | spec has no description; code has no description ✓ |
| instabilityIndex | description abbreviated |
| aliphaticIndex | description abbreviated |
| aromaticity | description abbreviated |

### Annotation parity per column (Fv)

| Column | Spec | Code | Status |
|---|---|---|---|
| charge | label, format, isScore, visibility, orderPriority (no description) | matches | PASS |
| isoelectricPoint | label, format, isScore, min, max, **description**, visibility, orderPriority | matches with shorter description | PASS |
| extinctionCoefficientOx | label, format, min, **description**, visibility, orderPriority | matches | PASS |
| extinctionCoefficientRed | label, format, min, **description**, visibility, orderPriority | label, format, min, visibility, orderPriority | **MISSING description** |
| molecularWeight | label, format, min, **description**, visibility, orderPriority | matches | PASS |

### `pl7.app/aminoAcid` axis values

Spec L196: "20 standard single-letter codes" listed. Implementation: `aa_tables.py:STANDARD_AAS = "ACDEFGHIKLMNPQRSTVWY"`. PASS. Pipeline emits 20 rows per entity (`run_peptide` iterates `STANDARD_AAS`). Verified by `test_peptide_aa_fraction_shape`.

### TCR-mode label substitutions

Verified above in R13a. All 4 chain combinations (TCRAB-A, TCRAB-B, TCRGD-A, TCRGD-B) covered in `labelFragments`.

---

## Section 5. Milestones — Acceptance Criteria

### M1 (Peptide mode)

| Criterion | Status |
|---|---|
| 9 scalar property columns + AA composition emitted | PASS — corpus test `pep_all20` covers schema |
| All isScore columns visible as Lead Selection ranking criteria | NOT VERIFIED end-to-end (`test/src/wf.test.ts` workflow tests are all `it.todo`) |
| Charge / GRAVY / MW / EC verified vs BioPython ProtParam (≥3 sequences, 2 dp) | PARTIAL — golden values pinned at 1e-6 internally but no external BioPython cross-check committed |
| pI bisection NA guard verified for polybasic | PASS (uses `R*50` not `R*10`; spec example was inexact) |
| Pre-M1: IPC 2.0 pKa tables embedded (peptide + protein), verified vs paper | PARTIAL — embedded; "verify before release" gate is documented but external verification not auditable from this review |

Implementation strategy spec (L520–522): "Use BioPython's `IsoelectricPoint` with a custom IPC 2.0 pKa dictionary (override `Bio.SeqUtils.IsoelectricPoint.pKa`) rather than reimplementing Henderson-Hasselbalch from scratch." → **DEVIATION.** Code reimplements from scratch (`properties.py`). Edge cases (ambiguity codes, zero-length, stop codons) are explicitly handled in custom code; correctness preserved. Trade-off: looser dependency surface, but M3 validation must check formula correctness AND pKa correctness rather than just pKa correctness.

### M2 (Antibody/TCR mode)

| Criterion | Status |
|---|---|
| CDR-H3 charge + hydrophobicity emitted for CDR3-only MiXCR dataset | PASS — corpus test `ab_cdr3_only` |
| Full-chain VH and VL emitted for whole-region MiXCR dataset | PASS — corpus test `ab_full_paired` |
| Fv charge and pI emitted for paired-chain whole-region | PASS — corpus test + golden values |
| Clone with any missing region: NA for full-chain, no block failure | PASS — `_reconstruct_chain` + per-clone null fill |
| CDR3-only input: no full-chain columns (absent, not NA-filled) | PASS — workflow gates emission, not just data |
| All isScore columns visible in Lead Selection ranking panel | NOT VERIFIED end-to-end |

### M3 (Validation)

| Criterion | Status |
|---|---|
| Charge/GRAVY/MW for ≥5 VH sequences vs BioPython, 2 dp | NOT COMMITTED — no external cross-validation script |
| Extinction coefficient: exact integer match vs BioPython | NOT COMMITTED |
| Instability index for ≥3 peptides vs ProteinAnalysis.instability_index() | NOT COMMITTED |
| pI for ≥5 VH sequences within 0.1 pH units of IPC 2.0 reference | NOT COMMITTED |
| VL pI on ≥2 sequences | PARTIAL — golden test pins one value |
| Fv charge / Fv pI manually verified on ≥2 paired sequences | PARTIAL — golden values pinned (`fv_charge ≈ 1.995494`, `fv_pi ≈ 9.330627`) but no external manual computation log |
| CDR-L3 charge by manual H-H on ≥3 sequences | NOT COMMITTED |
| CDR-H3 charge by manual computation on ≥10 sequences | NOT COMMITTED |
| Aliphatic for ≥3 VH sequences with closed-form check | PARTIAL — closed-form on ALVI in unit tests; one VH golden value |
| CDR-H3 length counting convention vs MiXCR IMGT | NOT VERIFIED (block does not emit length) |
| AstraZeneca/BMS dataset end-to-end | NOT COMMITTED — out-of-tree client validation |

> M3 is largely undischarged in committed artifacts. The corpus tests provide drift detection but are not external-tool cross-validation. This is a meaningful gap to flag for release readiness, though it may be tracked as future work.

---

## Section 6. Other Spec Items

### Score annotations rankingOrder semantics (L437)

Spec: "increasing = lower values ranked first". Code: `pl7.app/score/rankingOrder: "increasing"` is set on hydrophobicity columns only. PASS. Cross-reference annotation in spec: confirmed from `antibody-tcr-lead-selection/model/src/util.ts:247`.

### Open questions (L498)

"AA composition in antibody mode … default: defer to Phase 2." Implementation: `run_antibody_tcr` returns an empty AA-fraction frame (`pipeline.py:265`). **PASS** (deferred per spec default).

### Risks section (L566–572)

- pKa scale discrepancy with partner tools — partly addressed (`pl7.app/description` annotations on peptide charge mention IPC 2.0).
- Partial-region inputs — addressed via R11b info annotation (PASS).
- AA composition antibody mode — deferred per spec default (PASS).

---

## Section 7. Python Calculations — Formula Audit

Independent verification (re-derived formulas from spec; computed in this review session against reference values):

| Formula | Spec | Code | Independent reference | Status |
|---|---|---|---|---|
| Net charge `G` peptide pH 7 | H-H w/ IPC2_PEPTIDE | `−0.007915` | `0.99205 + (−0.99993) = −0.00788` (3 sig fig: −0.008) | PASS |
| MW(GGGGGG) | 6×57.0519 + 18.0153 | `360.3267` | 360.3267 | PASS |
| GRAVY(AAAA) | 4×1.8/4 | `1.800` | 1.8 | PASS |
| EC(YYWCCCC) ε_red | 2×1490 + 1×5500 = 8480 | `8480.0` | 8480 | PASS |
| EC(YYWCCCC) ε_ox | 8480 + (4//2)×125 = 8730 | `8730.0` | 8730 | PASS |
| AI(AAAVL) | 100×(0.6 + 2.9·0.2 + 3.9·(0+0.2)) = 196.0 | `196.0` | 196.0 | PASS |
| Aromaticity(WFY) | 3/3 = 1.0 | `1.000` | 1.0 | PASS |
| GRAVY(BXZJ) — all non-standard | None | `None` | NA expected | PASS |
| pI(GGGG) — termini-only | between pKa_n and pKa_c (~5.98) | covered by `test_pi_polyglycine_between_termini` | 4.5 < pI < 7.5 | PASS |
| pI(R*50) protein pKa — saturated polybasic | None (no zero crossing) | `None` | confirmed by guard logic | PASS |
| Fv(per-chain sum) ≠ concatenated charge for `DEKR/DEKR` | should differ (different terminus counts) | `−0.0003` vs `+0.0024` | difference confirmed | PASS |

Polyalanine instability closed-form (II = 9.0): PASS.
Glycine peptide pI between termini (~5.98): PASS.
Polyacidic 30 D pI low: PASS.
Cys include vs exclude differs by single Cys side-chain charge contribution: PASS (spot-checked in `test_cys_include_versus_exclude`).

`_bisect_zero` correctness:
- Same-sign guard at endpoints — PASS.
- Tolerance 0.001 — PASS.
- Returns midpoint when interval narrows below tolerance — PASS.
- Iterative narrowing: maintains sign-preserving substitution — PASS.

`_residue_charge` sign convention:
- Acid: `−1/(1 + 10^(pKa − pH))` — matches.
- Base: `+1/(1 + 10^(pH − pKa))` — matches.
- Cys: tracked separately via `include_cys` — matches.
- Non-ionizable (`pka_set.get(aa) is None`): returns 0.0 — matches.

`molecular_weight`:
- `Σ AVG_RESIDUE_MASS[c] for c in cleaned + H2O_AVG_MASS`
- AVG_RESIDUE_MASS values: NIST/UniMod condensed form — sample check: G=57.0519, A=71.0788, W=186.2132 — match published tables.
- H2O_AVG_MASS = 18.0153 (more precise than spec's 18.02 — improves accuracy).

`extinction_coefficients`:
- Operates on **raw** input (`aa_counts(seq)` is called on `seq` directly, not `cleaned`). `aa_counts` filters internally to standard residues. Equivalent to operating on cleaned sequence for counting purposes. PASS.
- floor division integer truncation correct.

`instability_index`:
- Length floor: `len(cleaned) < 10 → None`. PASS.
- Sum over `range(n - 1)`: correct dipeptide window.
- Skips dipeptides where DIWV is None (defensive, said to be impossible after `_prepare`).

`aliphatic_index`:
- Mole fractions over **cleaned** length (uses `aa_counts(seq)` then divides by `n = len(cleaned)`).
- Subtle: `aa_counts(seq)` filters non-standard internally so counts are consistent with `n`.

`aromaticity`:
- (F + W + Y) / cleaned_length — PASS.

`aa_fractions`:
- Returns dict over all 20 standard AAs always.
- Sum to 1.0 verified by `test_fractions_sum_to_one`.
- Returns None for invalid sequence — caller (`run_peptide`) writes 20 NA rows for the entity (uniform shape across entities).

`fv_*`:
- All four functions correctly delegate to per-chain `properties` and add. Fv pI bisection uses `lambda ph: charge_at_ph(VH) + charge_at_ph(VL)` — correct.

**Test results:** 96 unit tests + 24 integration corpus tests pass cleanly (`uv run pytest`). Corpus tests cover the spec's edge cases (`pep_with_stop`, `pep_polybasic_150`, `pep_no_aromatic`, `ab_cdr3_only`, `ab_h_only`, etc.) with manifest-driven expectations.

**Conclusion: Python computations are mathematically correct against the spec.** The reimplementation deviation (vs BioPython) is a code-level choice that does not affect output correctness; it shifts the M3 validation burden onto formula correctness as well as pKa correctness.

---

## Section 8. Summary Table — All Findings

| Severity | Item | Section |
|---|---|---|
| MISSING | R11c VHH/single-domain detection (≥16 aa median CDR-H3, IG, chain A only) | R11c |
| MISSING (M3) | External BioPython / IPC 2.0 cross-validation suite | M3 |
| MISSING | Peptide ε reduced description annotation | PColumn audit |
| MISSING | Peptide aliphatic index description annotation | PColumn audit |
| MISSING | Fv ε reduced description annotation | PColumn audit |
| DEVIATION | R1a: last-axis selection vs spec's "first-matching axis" scan | R1a |
| DEVIATION | Extra axis pattern `pl7.app/vdj/clonotypeKey` (alias for legacy bulk) | R1a |
| DEVIATION | Implementation reimplements H-H instead of leveraging BioPython | M1 strategy |
| DEVIATION | Python pipeline not vectorised (NumPy) — affects performance, not correctness | Step 2 |
| DEVIATION | TSV peptide column named `peptide_seq` not `sequence` | Step 1 contract |
| DEVIATION | TSV antibody schema does not include `receptor_type` column (carried via plan.json) | Step 1 contract |
| DEVIATION | CDR-L3 description not differentiated from CDR-H3 description | PColumn audit |
| DEVIATION | Spec example "RRRRRRRRRR emits NA" is mathematically incorrect; feature works | M1 acceptance |
| DEVIATION (subtle) | Per-chain partial coverage with no CDR3 silently emits no info message | R11b edge case |
| DEVIATION (subtle) | Mutation of shared spec object stamps blockId on both `propertiesPf` and `exports.properties` (comment says "export-only") | process.tpl.tengo:251–259 |
| EXTENSION | Bulk MiXCR fallback: chain undefined → assume "A" | R13 |
| EXTENSION | SD-003 receptor read from axis domain in addition to per-column | R13b |
| PASS (everything else) | All formulas, all column specs, all order priorities, all min/max bounds, all isScore + rankingOrder, all edge case handling | — |

This concludes the line-by-line spec compliance review. Idiomatic-ness audit, deviations cross-check, and upgrade audit follow in subsequent sections after `spec-deviations.md` is read (per review protocol).

---

## Section 9. Cross-Check Against `spec-deviations.md`

Read after the spec compliance review concluded. The block's `docs/spec-deviations.md` file documents three deviations:

### SD-001 — Skip Secondary Alleles in Single-Cell Paired Data

**Documented:** Yes, with full root-cause + decision rationale.
**Deviation type:** Real-world MiXCR data shape not anticipated by spec (single-cell secondary alleles). Code skips columns where `pl7.app/vdj/scClonotypeChain/index != "primary"`.
**My review status before reading:** I had not flagged this. The implementation is correct given the chosen filter; the deviation is from an undocumented MiXCR schema detail, not from the spec text. **Now confirmed: properly documented and justified.**

### SD-002 — Treat `FR4InFrame` as `FR4`

**Documented:** Yes, with root-cause + alternative options + decision.
**Deviation type:** MiXCR's actual feature label is `FR4InFrame`, not `FR4`. Code normalises before matching against `REQUIRED_FEATURES`.
**My review status before reading:** Not flagged. The spec used the literal `"FR4"`; reality requires the normalisation. **Properly documented.**

### SD-003 — Read Receptor From `clonotypeKey` Axis Domain

**Documented:** Yes, with root-cause + alternative options + decision.
**Deviation type:** MiXCR places `pl7.app/vdj/receptor` on the input anchor's secondary axis domain rather than per-column. Code reads the axis first, falls back to per-column lookup.
**My review status before reading:** I noted this as "extension SD-003" because the inline workflow comment referenced it. **Properly documented.**

### Cross-Check Verdict — What's Documented vs What's Not

`spec-deviations.md` covers the **3 MiXCR data-shape** deviations. The other deviations / missing items I flagged above are **not** documented in `spec-deviations.md`:

| Finding | In spec-deviations.md? |
|---|---|
| **R11c VHH detection** (MISSING) | NOT documented |
| **M3 external validation** (NOT COMMITTED) | NOT documented |
| Peptide ε reduced description (MISSING) | NOT documented |
| Peptide aliphatic index description (MISSING) | NOT documented |
| Fv ε reduced description (MISSING) | NOT documented |
| **R1a last-axis vs first-matching axis** (DEVIATION) | NOT documented |
| **Extra axis pattern `pl7.app/vdj/clonotypeKey`** (DEVIATION/extension) | NOT documented |
| **BioPython implementation strategy not followed** (DEVIATION) | NOT documented |
| **Python pipeline not vectorised** (DEVIATION) | NOT documented |
| TSV column name `peptide_seq` vs `sequence` | NOT documented |
| TSV antibody schema lacks `receptor_type` column | NOT documented |
| CDR-L3 description not differentiated from CDR-H3 | NOT documented |
| Spec example "RRRRRRRRRR emits NA" inaccurate (spec issue, not code) | NOT documented |
| Per-chain partial coverage with no CDR3 silently no message (R11b edge) | NOT documented |
| `blockId` mutation affects both `propertiesPf` and `exports.properties` (subtle bug) | NOT documented |
| Bulk MiXCR fallback chain → "A" (extension) | NOT documented |

**Recommendation:** The undocumented items break into three classes:

1. **Spec-level corrections** (e.g., RRRRRRRRRR example, axis pattern naming, TSV schema details, spec performance approach assumptions) — should be either fixed in spec or recorded in `spec-deviations.md`.

2. **Real omissions** (R11c VHH detection, missing description annotations, M3 external validation harness) — should either be implemented or moved to a "deferred to follow-up" entry in `spec-deviations.md` with a justification.

3. **Subtle bugs** (`blockId` mutation, R11b CDR3-absent silent path) — fix in code, no deviation entry needed.

The rule of thumb from the harness is: deviations recorded in code comments AND in `spec-deviations.md`. The 3 MiXCR-shape items meet that bar; the rest do not.

---

## Section 10. Idiomatic-ness Audit

### 10.1 UI

#### What's Idiomatic (PASS)

- **`PlBlockPage` with title + append slot** for the Logs button (`MainPage.vue:23–29`). Standard pattern.
- **`PlDropdownRef` with `v-model` on `app.model.data.inputAnchor`** and options from `app.model.outputs.inputOptions` (`MainPage.vue:31–37`). Canonical input selection pattern.
- **`PlAlert` for info messages** rendered from `app.model.outputs.info?.messages` array (`MainPage.vue:39–46`). Reactive, idiomatic.
- **`PlAgDataTableV2` with `usePlDataTableSettingsV2`** binding to `app.model.outputs.propertiesTable` and state at `app.model.data.tableState` (`MainPage.vue:18–22, 48`). Idiomatic V2 table.
- **`PlSlideModal` + `PlLogView`** for log inspection (`MainPage.vue:50–55`). Standard.
- **Local Vue ref for modal open state** (`const logOpen = ref(false)`). Correct — modal open/close belongs in local component state, not the model.
- **`defineAppV3(platforma, ...)`** entry point in `app.ts`. Correct V3 app shape.
- **No hand-rolled HTML table, no direct DOM, no bypass of standard SDK components.**

#### Hairpin (block-label exception — currently TOLERATED but DEAD)

`app.ts:7–13`:

```ts
watchEffect(() => {
  const anchor = app.model.data.inputAnchor;
  const opts = app.model.outputs.inputOptions ?? [];
  const match = anchor
    ? opts.find((o) => o.ref?.blockId === anchor.blockId && o.ref?.name === anchor.name)
    : undefined;
  app.model.data.defaultBlockLabel = match?.label ?? "";
});
```

This is the hairpin pattern (output → data write). Per `harness/hairpin.md` it is a tolerated exception **for the block-label use case only.** However:

**The `defaultBlockLabel` is dead code:**
- `model/src/index.ts` `.title()` callback hardcodes `"Sequence Properties"` (line `.title(() => "Sequence Properties")`). It does **not** read `ctx.data.defaultBlockLabel`.
- `BlockArgs.defaultBlockLabel` is projected to args (line `defaultBlockLabel: data.defaultBlockLabel ?? "Sequence Properties"`) but the workflow (`main.tpl.tengo`) does **not** consume `args.defaultBlockLabel` anywhere.

**Net effect:**
1. The watcher fires on dataset changes, mutating `data.defaultBlockLabel`.
2. The args projection picks up the mutation, **marking the block stale** (per `BlockModelV3` semantics: any args change fires the stale gate).
3. The workflow runs (or prompts user to Run). Workflow consumes `args.inputAnchor` and ignores `args.defaultBlockLabel`.
4. The title in the project listing remains "Sequence Properties" regardless.

**FINDING (Idiomatic): Dead-code hairpin.** Either:
- Remove the watcher + `defaultBlockLabel` field entirely, OR
- Wire `.title()` to read `ctx.data.defaultBlockLabel` and remove `defaultBlockLabel` from `args` projection (it's UI-only state, not workflow input).

The current state is the worst of both worlds: incurs the stale-gate cost without delivering the dynamic title benefit.

#### Other UI Concerns

- `usePlDataTableSettingsV2` is correct V2 settings hook (matches V3 model's `createPlDataTableV3` selector mode = "enrichment").
- Outputs null-checks — `app.model.outputs.info?.messages ?? []` on the alerts iteration — correct defensive pattern.
- No localisation, but block strings are consistent with sibling blocks. Acceptable.

### 10.2 PFrame Layout

#### What's Idiomatic (PASS)

- **Universal value column names** (`pl7.app/charge`, `pl7.app/hydrophobicity`, etc.) — matches the cross-block convention introduced in peptide-discovery.
- **Domain disambiguation by `pl7.app/feature`** (`"peptide"`, `"CDR3"`, `"VDJRegion"`, `"Fv"`) — exactly the pattern the harness recommends.
- **Domain disambiguation by `pl7.app/vdj/scClonotypeChain`** (`"A"`, `"B"`) — chain identity preserved.
- **Trace via `pSpec.makeTrace` + `trace.inject`** — never hand-built. Done at `process.tpl.tengo:286–290`. Correct.
- **`pframes.pFrameBuilder()` then `.build()`** — canonical construction.
- **Sorted iteration** of result columns via `maps.getKeys(scalarOut)` before `fb.add(...)` — preserves dedup-stable ordering. Correct.
- **`pframes.exportFrame(...)` for the model output** — correct V3 contract.
- **`exports: { properties }` in body return** — populated. Required `exports` key present.
- **Axes pass-through** — output `axes[0]` uses input `keyAxisSpec` verbatim. Cross-block joins preserved.
- **Per-block-instance disambiguation via `pl7.app/blockId` domain key** on isScore export columns — supports multiple instances of the block in the same project.

#### Subtle Bug — `blockId` Mutation Leak

`process.tpl.tengo:251–259`:

```go
exportColumns := []
for col in columns {
  if col.spec.annotations && col.spec.annotations["pl7.app/isScore"] == "true" {
    if !col.spec.domain { col.spec.domain = {} }
    col.spec.domain["pl7.app/blockId"] = blockId
    exportColumns += [col]
  }
}
```

The same `col` reference is in **both** `columns` (used in `outputSpecs`, fed to `xsv.importFile` for `scalarOut`) **and** `exportColumns` (used in `exportSpecs`, fed to `xsv.importFile` for `exportOut`). Mutating `col.spec.domain` in place stamps `blockId` on the *internal* output spec as well as the export spec.

The inline comment says "Stamp blockId on the export-only columns" — but the implementation stamps it on both. In practice this is mostly benign (the internal `propertiesPf` is presence-checked but the displayed table uses result-pool enrichment from `exports.properties`), so the user sees correctly-labelled data either way. But:

- If a future debugging session asserts on the `propertiesPf` spec shape, the unexpected `pl7.app/blockId` domain key will surprise.
- If the model ever switches from result-pool enrichment to direct `propertiesPf` consumption, downstream queries that don't include `blockId` in their spec pattern still match (subset matching) — no breakage. But the comment is misleading.

**Recommendation:** Either deep-copy the spec before mutation, or update the comment to reflect that both outputs receive the stamp.

#### PColumn Spec Construction

- **Inherit when possible:** input axis verbatim, AA-fraction second axis introduces a new dimension `pl7.app/aminoAcid` (correct namespacing). Domain on the new axis is empty (no need for one — axis values are the 20 standard codes). Annotations include `pl7.app/label: "Amino Acid"`. Idiomatic.
- **Domain construction for new columns:** `{pl7.app/feature: <…>}` with optional `pl7.app/vdj/scClonotypeChain` for VDJ. Domain is used for disambiguation, not arbitrary tagging. Idiomatic.
- **Annotation set:** label, format, isScore, score/rankingOrder (where applicable), description, table/visibility, table/orderPriority, min/max — all standard annotations. No off-spec annotations.

### 10.3 Standard Block Layout

#### Workspace Structure (PASS)

```
sequence-properties/
├── block/
│   ├── package.json         # block manifest with components map
│   └── block-pack/          # build artifact
├── model/                    # TS BlockModelV3
│   └── src/{index.ts, dataModel.ts, types.ts}
├── ui/                       # Vue 3 + ui-vue SDK
│   └── src/{app.ts, main.ts, pages/MainPage.vue}
├── workflow/                 # Tengo
│   └── src/{main.tpl.tengo, process.tpl.tengo, wf.test.ts}
├── software/                 # Python
│   ├── pyproject.toml
│   ├── src/{main.py, pipeline.py, properties.py, instability.py, pka_tables.py, aa_tables.py, io_layer.py}
│   └── tests/{unit, integration, data/corpus}
├── test/                     # TS integration test (currently all .todo)
│   └── src/wf.test.ts
├── logos/
├── package.json              # turbo orchestrator + tooling
├── pnpm-workspace.yaml
├── pnpm-lock.yaml
├── .changeset/icy-walls-flow.md
└── .github/workflows/{build, mark-stable}.yaml
```

This matches the canonical block layout (`harness/block-anatomy.md`). Components map in `block/package.json` correctly references the workflow + model + ui packages.

#### Tooling

- pnpm workspace + turbo — standard.
- `@platforma-sdk/block-tools` for pack / publish — standard.
- `@platforma-sdk/eslint-config` not present in workspace deps; instead uses `oxlint`/`oxfmt` (per `model/.oxlintrc.json` etc.). This is non-standard for blocks per the harness recommendation but consistent with internal tooling — not a breaking issue.
- Changeset present (`icy-walls-flow.md`) — covers software, workflow, root, model, test, ui packages with `minor` bumps. Well-formed.

#### `index.html` for UI (REQUIRED)

Verified present at `ui/index.html` (349 bytes). PASS.

#### Workflow Test Coverage

`test/src/wf.test.ts` has a comprehensive R-by-R structure but **all tests are `.todo`**:

- 6 modality detection tests — all `.todo`
- 4 peptide-mode tests — all `.todo`
- 11 antibody/TCR tests — all `.todo`
- 5 edge-case tests — all `.todo`
- 3 downstream-consumption tests — all `.todo`

**FINDING:** No end-to-end integration coverage of the workflow + model + UI wiring. Unit tests (Python) and corpus tests (Python pipeline) do cover the property computation. But **no test verifies** that the workflow produces the expected PColumn shapes in a real platforma backend, that score columns are discoverable by Lead Selection, that trace annotations land correctly, or that the info messages flow to UI.

This is a meaningful gap before marking the block stable — the .todo skeleton lays out exactly what should be checked but the checks are not wired up.

### 10.4 Block Upgrade Story

#### Current State

`model/src/dataModel.ts`:

```ts
export const blockDataModel = new DataModelBuilder()
  .from<BlockData>("Ver_2026_04_28")
  .init(() => ({
    tableState: createPlDataTableStateV2(),
  }));
```

**Properties of this setup:**

- **Single named version** `"Ver_2026_04_28"`. First and only release version.
- **No `.upgradeLegacy(...)`** — correct, this is a fresh V3 block, not a V1 migration.
- **No `.migrate(...)` chain yet** — also correct, no prior version to migrate from.
- **`init()` provides `tableState` only.** `inputAnchor` and `defaultBlockLabel` start `undefined`.

**Version-naming convention.** Date-based (`Ver_2026_04_28`) is unusual — sibling V3 blocks more commonly use sequential names like `"v1"`, `"v2"`. Date-based names work but reduce ergonomics: `.migrate("Ver_2026_04_28", "Ver_2026_05_15", fn)` is harder to read and order at a glance than `.migrate("v1", "v2", fn)`. **Not a bug, just a convention call.**

#### Upgrade Scenarios — Will customers upgrade cleanly?

**Scenario A: Patch / minor with no schema change** (e.g., add a description annotation, fix a formula):
- Old `BlockData` shape unchanged.
- Customer's persisted `data` (still `Ver_2026_04_28`) is read with no migration.
- Workflow re-runs against same args: produces a new CID (different code → different output bytes if anything changed).
- Old result is replaced in the result graph.
- **Clean upgrade.** No re-run required if formula didn't change; re-run required if it did.

**Scenario B: Schema change adding a new optional field** (e.g., adding a `filterMode?: string` to `BlockData`):
- Need to bump version, e.g., `from<BlockData>("Ver_2026_06_01")` and chain `.migrate("Ver_2026_04_28", "Ver_2026_06_01", (old) => ({ ...old, filterMode: undefined }))`.
- If forgotten, customers' old projects break (`data.filterMode` is `undefined`, but if downstream code assumes a field shape, that breaks).
- **Requires care from author.**

**Scenario C: PColumn schema change** (e.g., new column added):
- New CIDs everywhere — re-run required for the new column.
- Old downstream consumers that match by spec patterns continue working (subset matching).
- **Clean upgrade for downstream**, but customer must press Run to populate the new column.

**Scenario D: PColumn schema change** (e.g., domain key renamed):
- Major risk of breaking downstream blocks (Lead Selection won't find the score columns under the new domain).
- Requires careful coordination — should be a `major` bump with release notes.
- **Not graceful** — but this is a fundamental constraint of any cross-block schema change.

#### CID Conflict Story (Block Re-Use, Multi-Instance)

The block stamps `pl7.app/blockId` on isScore export columns. Each block instance has a unique blockId. So two `sequence-properties` blocks in the same project producing the same property columns can coexist in the result pool without spec collision.

- Lead Selection picker queries by `pl7.app/isScore: "true"` and matching feature/chain domain — **gets both instances**, distinguished by blockId.
- The user can pick which instance to use as the source of scores.

**Verdict:** Good multi-instance design.

#### Dedup Story

The block goes to extra effort for content-addressed dedup (per the changeset `icy-walls-flow.md`):

1. **Canonical JSON resources** — `canonicalJsonResource` wraps `canonical.encode(...)` for `plan.json` (Python step input), `params` (process template input), `infoBlob`. This sorts keys at every nesting level so resource bytes are stable.
2. **Sorted Tengo map iteration** — `maps.getKeys(...)` everywhere a map is walked into a pframe builder or output structure.
3. **Sorted Python TSV output** — `write_output_tsv(..., sort_keys=["entity_key"])` and `["entity_key", "aminoAcid"]`. Resource bytes hash identically across runs.
4. **Quantization for ULP stability** — `_quantize_for_cid` rounds `charge_*` and `pi_*` columns to 3 decimals before TSV write. The pipeline comment correctly identifies these as the only outputs whose CID can vary across libm/numpy SIMD changes (transcendental `10**x`).

**This is exemplary dedup engineering.** Two runs of the same input on the same code produce byte-identical output → identical CID → result resource is shared, not recomputed.

The 3-decimal quantization is well-justified: `_bisect_zero` tolerance is 0.001, display format is `.2f` (so `.001` is below display precision). Correctness preserved; non-determinism in sub-display ULPs is squashed.

#### Pre-PR Sanity (per harness)

- `pnpm-lock.yaml` is committed alongside `pnpm-workspace.yaml`. Should be checked at PR time.
- `@platforma-sdk/block-tools` version: latest catalog entry per pnpm-workspace; can't verify against npm latest from inside this review.
- Changeset present and well-formed.

### 10.5 Summary — Idiomatic-ness

| Layer | Verdict |
|---|---|
| UI components and patterns | PASS, with one dead-code hairpin (defaultBlockLabel) |
| PFrame construction (specs, axes, domain, trace, exports) | PASS, with one mutation leak (blockId stamping) |
| Block layout (workspace, packages, components map) | PASS |
| Workflow test coverage | INCOMPLETE (all `.todo`) |
| Upgrade path (DataModel versioning) | READY for future migrations; first version named oddly |
| CID dedup engineering | EXCELLENT |
| Multi-instance support (blockId stamping) | PASS |
| `ui/index.html` | PASS (verified present) |

---

## Section 11. Top Recommendations (Actionable)

Ordered by impact:

1. **Implement R11c (VHH detection info annotation)** OR add to `spec-deviations.md` with a deferred-to-Phase-2 justification. Currently entirely missing.
2. **Wire up the workflow integration tests** (`test/src/wf.test.ts` is currently 100% `.todo`). The skeleton enumerates exactly the right checks (R1–R16); just needs implementation against `blockTest`.
3. **Fix the `defaultBlockLabel` dead code path:**
   - Either: remove the field, the watcher, and the args projection (simplest);
   - Or: read `ctx.data.defaultBlockLabel` from `.title(...)` and remove the field from args projection (UI-only state).
4. **Fix the `blockId` mutation leak** in `process.tpl.tengo:251–259` — deep-copy the spec or document that both outputs receive the stamp.
5. **Add missing description annotations** on:
   - peptide ε reduced
   - peptide aliphatic index
   - Fv ε reduced
6. **Differentiate CDR-L3 description** from CDR-H3 description per spec table at `pcolumn-spec.md` L268–272.
7. **Document the spec / implementation drift** in `spec-deviations.md`:
   - R1a "first-matching axis" vs implementation last-axis
   - Extra axis pattern `pl7.app/vdj/clonotypeKey`
   - BioPython implementation strategy not followed
   - Pipeline not vectorised (with link to perf plan if it exists)
   - TSV column name `peptide_seq` vs spec's `sequence`
   - TSV antibody schema lacks `receptor_type` column
8. **Drive M3 external validation** to a committed state — at minimum, a runnable script comparing block outputs to BioPython on a fixed dataset; pin tolerances per spec criteria. Currently undischarged.
9. **Edge case for R11b**: if a chain has 1–6 regions but **no CDR3**, no info message fires. Consider adding a parallel branch.
10. **Spec correction (out of code's hands):** `README.md:528` cites `RRRRRRRRRR` as triggering pI NA. Mathematically, this sequence has pI ≈ 12.79 with the IPC2 peptide pKa set; the example needs updating to a properly polybasic sequence (e.g., R*150) or the protein pKa set needs to be specified.

End of review.

---

## Amendment — SD-001 / SD-002 Upstream Verification

**Date added:** 2026-04-29
**Scope:** Verify that the deviation documents accurately describe upstream block behavior, and check whether any other upstream block in the workspace would break the current deviation handling.

### Method

Direct grep against the source files of all relevant upstream blocks in the workspace. Files examined:

- `blocks/mixcr-clonotyping/workflow/src/{process.tpl.tengo, calculate-export-specs.lib.tengo}`
- `blocks/import-vdj-data/workflow/src/{infer-columns-cellranger.lib.tengo, infer-columns-mixcr.lib.tengo, process-bulk.tpl.tengo, process-single-cell.tpl.tengo, process-utils.lib.tengo}`
- `blocks/antibody-sequence-liabilities/workflow/src/{main.tpl.tengo, process.tpl.tengo}`
- `blocks/cellecta-drivermap-air-mixcr-clonotyping/workflow/src/get-export-params.lib.tengo`

### SD-001 — Skip Secondary Alleles — Verification

#### Claim 1 (deviation doc): MiXCR clonotyping emits `pl7.app/vdj/scClonotypeChain/index` with values `"primary"` / `"secondary"` on single-cell paired data.

**CONFIRMED.** Source: `blocks/mixcr-clonotyping/workflow/src/process.tpl.tengo:728–784`.

```go
for isPrimary in [true, false] {
    pPrefixU := isPrimary ? "Primary" : "Secondary"
    pPrefixL := text.to_lower(pPrefixU)
    …
    domain: {
        "pl7.app/vdj/receptor": receptor,
        "pl7.app/vdj/scClonotypeChain": chainLetterU,
        "pl7.app/vdj/scClonotypeChain/index": pPrefixL,
        …
    }
}
```

The loop runs over `[true, false]` so every chain emits both a `Primary` and `Secondary` variant in single-cell mode. The lower-cased value (`"primary"` / `"secondary"`) goes into the domain key — exact match for the SD-001 filter (`if idx != undefined && idx != "primary" { continue }`).

**Additional finding:** MiXCR's secondary chain emits **only `aaSeqCDR3`** (`process.tpl.tengo:687–693`):

```go
//removing columns from secondary chain except aaSeqCDR3
columnsSpecPerClonotypeSecondary := []
for col in columnsSpecPerClonotypeNoAggregates {
    if col.column == "aaSeqCDR3" {
        columnsSpecPerClonotypeSecondary += [ col ]
    }
}
```

So secondary alleles never carry FR1–FR4 columns. Without the SD-001 filter, the block sees two `A_CDR3` columns and panics on the duplicate header — exactly the symptom recorded. The filter is the minimum viable fix.

#### Claim 2 (deviation doc): `pl7.app/vdj/scClonotypeChain/index` is single-cell-specific.

**CONFIRMED.** The `for isPrimary` loop in MiXCR's `process.tpl.tengo` is inside the single-cell branch (axis = `axisByScClonotypeKeyGen(receptor)`, line 766). Bulk MiXCR output never sets the `/index` domain key. Bulk MiXCR data is unaffected by the filter — the `idx != undefined` guard is correct.

#### Claim 3 (deviation doc): "Pattern block also lacks this filter, so the same bug likely affects `antibody-sequence-liabilities` on paired single-cell data."

**INCORRECT — `antibody-sequence-liabilities` already filters for primary.** Source: `blocks/antibody-sequence-liabilities/workflow/src/main.tpl.tengo:112, 127, 136`:

```go
if additionalSeq.spec.domain["pl7.app/vdj/scClonotypeChain/index"] == "primary" {
    …
}
…
} else if seq_entry.spec.domain["pl7.app/vdj/scClonotypeChain/index"] == "primary" {
    …
}
…
if ann_entry.spec.domain["pl7.app/vdj/scClonotypeChain/index"] == "primary" {
    …
}
```

The pattern block uses an inverse predicate — *include* only when `index == "primary"`. SD-001 in sequence-properties uses *exclude* when `index != "primary"`. They are equivalent in effect. The pattern block does **not** carry the bug; the deviation doc's parenthetical claim is wrong.

**Recommendation:** Update the deviation doc to remove the antibody-sequence-liabilities reference, OR clarify that the pattern block uses an inclusive filter that the sequence-properties block adopted as an exclusive filter.

#### Claim 4 (additional sources): Are there other upstream blocks emitting `scClonotypeChain/index`?

**YES — `import-vdj-data` cellranger single-cell path.** Source: `blocks/import-vdj-data/workflow/src/process-single-cell.tpl.tengo:317–319`:

```go
…
"pl7.app/vdj/scClonotypeChain": text.to_upper(chainLetter),
"pl7.app/vdj/scClonotypeChain/index": pPrefixL
```

Same domain key, same primary/secondary values. SD-001's filter handles this correctly without modification.

**Cellecta-DriverMap-AIR-MiXCR-Clonotyping** uses MiXCR via composition (`get-export-params.lib.tengo:25`); the same emission path. Also handled.

#### SD-001 Verdict

- The technical claim (MiXCR emits primary/secondary on single-cell) is correct.
- The fix is correct.
- The pattern-block bug-propagation note in the deviation doc is wrong — that block already filters.
- No additional upstream blocks would break the SD-001 filter.

### SD-002 — Treat `FR4InFrame` as `FR4` — Verification

#### Claim 1 (deviation doc): MiXCR exports FR4 with `pl7.app/vdj/feature: "FR4InFrame"` not `"FR4"` for amino acid alphabet.

**CONFIRMED.** Source: `blocks/mixcr-clonotyping/workflow/src/calculate-export-specs.lib.tengo:18–21, 505`:

```go
inFrameFeatures := {
    "FR4": "FR4InFrame",
    "VDJRegion": "VDJRegionInFrame"
}
…
featureInFrameU := isAminoAcid ? inFrameFeatures[featureU] : featureU
if is_undefined(featureInFrameU) {
    featureInFrameU = featureU
}
…
domain: {
    "pl7.app/vdj/feature": featureInFrameU,
    "pl7.app/alphabet": alphabet
}
```

The substitution is **AA-only**. Confirmed by the `isAminoAcid ? inFrameFeatures[featureU] : featureU` ternary. Implications:

- aa FR4 column → `pl7.app/vdj/feature: "FR4InFrame"` — needs SD-002 normalization.
- nucleotide FR4 column → `pl7.app/vdj/feature: "FR4"` — but the block only queries aa alphabet, so this is irrelevant.
- aa VDJRegion column → `pl7.app/vdj/feature: "VDJRegionInFrame"` — see additional finding below.

#### Claim 2 (additional finding the deviation doc does NOT mention): MiXCR also emits `VDJRegionInFrame`.

**CONFIRMED** by the same `inFrameFeatures` map. Sequence-properties does not handle this name — the loop check `if !contains(REQUIRED_FEATURES, feat) { continue }` silently drops it because `REQUIRED_FEATURES` lists FR1, CDR1, FR2, CDR2, FR3, CDR3, FR4 only.

**Implication:** The block does **not** miss any required feature, because VDJRegion is not in the 7-region reconstruction list (the block reconstructs from individual FR/CDR regions, not from the pre-assembled VDJRegion column). However, this is a missed optimisation per SD-002 option B in the deviation doc — using the pre-assembled `VDJRegionInFrame` directly would skip the 7-region concat. The deviation doc already records this as deferred.

#### Claim 3 (deviation doc): Spec did not anticipate the `InFrame` suffix on FR4.

**CONFIRMED.** `README.md` lists `"FR4"` in R4 verbatim, with no mention of `InFrame` variants. The spec drift is real.

#### Claim 4 (additional sources): Do any other upstream blocks emit FR4 with a non-standard name?

**No.** Source: `blocks/import-vdj-data/workflow/src/infer-columns-cellranger.lib.tengo:345`:

```go
domain: {
    "pl7.app/alphabet": "aminoacid",
    "pl7.app/vdj/feature": "FR4"
},
```

Cellranger emits aa FR4 as plain `"FR4"`. Other import paths (`infer-columns-immunoSeq`, `infer-columns-qiagen`, `infer-columns-custom`, `infer-columns-mixcr`) — searched and none emit `InFrame` variants.

`infer-columns-mixcr.lib.tengo:35` does have:

```go
"aaSeqFR4": "fr4-aa",
```

But this is the column **name** mapping, not the `pl7.app/vdj/feature` domain value. The domain still ends up as `"FR4"` for cellranger / immunoSeq / qiagen / custom imports.

So the InFrame substitution is **MiXCR-specific** (and Cellecta-DriverMap-AIR-MiXCR-Clonotyping by composition). All other upstream blocks emit `"FR4"` directly, which already matches `REQUIRED_FEATURES`. SD-002's normalization works for both worlds:

- MiXCR / Cellecta path: `"FR4InFrame"` → normalized → `"FR4"` → matches.
- import-vdj-data path: `"FR4"` → no normalization → matches directly.

#### Claim 5 (additional concern): Will future MiXCR additions to `inFrameFeatures` break the block?

The map currently has:

```go
inFrameFeatures := {
    "FR4": "FR4InFrame",
    "VDJRegion": "VDJRegionInFrame"
}
```

If MiXCR adds entries (e.g., `"FR1": "FR1InFrame"` to filter productive FR1 sequences specifically), sequence-properties would silently drop those columns (`!contains(REQUIRED_FEATURES, feat)` is true for `"FR1InFrame"`), reverting to CDR3-only output for affected chains. The user would see the partial-region info message but not understand the cause.

**Mitigation options (not currently implemented):**
- Maintain a parallel `IN_FRAME_NORMALIZATIONS` map mirroring MiXCR's `inFrameFeatures`.
- Strip a trailing `"InFrame"` substring before checking against `REQUIRED_FEATURES`.

The risk is low (the `inFrameFeatures` map has been stable in MiXCR), but worth flagging.

#### SD-002 Verdict

- The technical claim (MiXCR uses `FR4InFrame`) is correct.
- The fix is correct and works for both MiXCR and import-vdj-data.
- VDJRegionInFrame is silently skipped — not a bug under the current 7-region reconstruction approach.
- Future-proofing: a generalised `*InFrame` strip would be more robust than the current single-feature special-case.

### Other Upstream Considerations Discovered During Verification

#### A. Bulk `import-vdj-data` Uses `pl7.app/vdj/chain` Instead of `pl7.app/vdj/scClonotypeChain`

Source: `blocks/import-vdj-data/workflow/src/process-bulk.tpl.tengo:128, 150, 252` and `process-utils.lib.tengo:21, 42`:

```go
domain: {
    "pl7.app/vdj/chain": chain,
    "pl7.app/vdj/clonotypingRunId": blockId,
    …
}
```

The bulk path uses `pl7.app/vdj/chain` (different domain key) and emits **separate `clonotypeKey` axes per chain**. Each axis is single-chain.

The sequence-properties block reads chain via `d["pl7.app/vdj/scClonotypeChain"]` and falls back to `"A"` when undefined (`main.tpl.tengo:155–157`). For bulk import-vdj-data:

- Input anchor is one of the per-chain `clonotypeKey` axes (single-chain).
- All sequence columns under that axis carry only one chain's data.
- All columns get labeled `"A"` (correct — single-chain dataset).
- No Fv emitted (correct — `hasFv` requires both `fullChain["A"]` and `fullChain["B"]`).
- Labels would default to antibody (heavy/VH) regardless of actual chain (heavy/light/alpha/beta/etc.) because the receptor heuristic could read receptor from axis domain (SD-003) but the chain identity is lost.

**Concrete user-facing concern:** If a user imports light-chain bulk data via `import-vdj-data`, sequence-properties would emit `"VH Net Charge (pH 7)"` labels (since chain "A" maps to "Heavy" / "VH"), but the data is actually light-chain. **Labels misled.** No filter or detection differentiates heavy from light when the chain key is missing.

**Severity:** Moderate. Affects a real-world flow (importing single-chain bulk data). Properties are computed correctly; only the labels are wrong.

**Recommended action:** Either:
1. Detect `pl7.app/vdj/chain` on the input anchor's axis domain (similar to the SD-003 receptor read) and use it to set the chain letter, OR
2. Read `pl7.app/vdj/chain` per-column as a fallback when `scClonotypeChain` is absent, mapping (heavy/IGH→A, light/IGK/IGL→B, alpha/TRA→A, beta/TRB→B, etc.).

This is a **NEW finding**, not in `spec-deviations.md`, and not in the original review's section 8. Worth opening as a follow-up.

#### B. Spec Wording vs Actual MiXCR Axis Name

The spec lists `pl7.app/vdj/cloneId` as the legacy MiXCR bulk axis. The original review noted this as "extension" because the block also accepts `pl7.app/vdj/clonotypeKey`. **Verification confirms `pl7.app/vdj/clonotypeKey` is the actual MiXCR axis name** — there is no `pl7.app/vdj/cloneId` emission anywhere in `mixcr-clonotyping/workflow/src/`. The spec's R1a name is fictional / aspirational.

**The block's "extension" axis pattern is, in fact, the only working pattern.** The `pl7.app/vdj/cloneId` entry in `inputAnchorSpecs` is dead code — no upstream block emits that axis. Could be removed without functional impact, OR the spec should be updated to reflect reality.

### Amendment Summary

| Finding | Status |
|---|---|
| SD-001 technical claim (MiXCR primary/secondary) | CONFIRMED |
| SD-001 antibody-sequence-liabilities bug-propagation note | INCORRECT — pattern block already filters |
| SD-001 also affects `import-vdj-data` cellranger SC path | NEW SOURCE (handled by current filter) |
| SD-002 technical claim (FR4InFrame) | CONFIRMED |
| SD-002 covers VDJRegionInFrame implication | DEFERRED (per option B in deviation doc) |
| SD-002 generalised *InFrame strip would be more robust | OPEN — not currently implemented |
| Bulk `import-vdj-data` uses `pl7.app/vdj/chain` not `scClonotypeChain` — labels would be wrong | NEW FINDING — not in deviations doc |
| `pl7.app/vdj/cloneId` axis pattern in spec / model is fictional | NOTED — `clonotypeKey` is the real name |

### Recommendations Added to Section 11

12. **Update SD-001 deviation doc:** remove the antibody-sequence-liabilities bug-propagation claim — that block already filters on primary.
13. **Document import-vdj-data cellranger SC support:** add a note in SD-001 that the same filter handles `import-vdj-data` cellranger output.
14. **Investigate bulk `import-vdj-data` chain-label drift:** light-chain or non-A-chain bulk imports get heavy-chain labels because `scClonotypeChain` is absent. Add chain-letter resolution from `pl7.app/vdj/chain` axis-domain.
15. **Consider generalising SD-002:** strip a trailing `"InFrame"` from feature names rather than special-casing FR4. Future-proofs against MiXCR adding new InFrame variants.
16. **Reconcile spec / reality on bulk MiXCR axis name:** spec says `pl7.app/vdj/cloneId`, MiXCR emits `pl7.app/vdj/clonotypeKey`. Either update the spec or document the gap in `spec-deviations.md`.
