# Punchlist — Sequence Properties Block

**Date:** 2026-04-29
**Branch:** `paulnewling/MILAB-6140_init-block`
**Last commits:** `b9407ae` (R11c + BioPython migration) → `5902213` (test coverage)

**Sources audited for this list:**
- `docs/text/work/projects/sequence-properties/README.md` (spec)
- `docs/text/work/projects/sequence-properties/pcolumn-spec.md` (column contracts)
- `docs/reviews/comprehensive-spec-review.md` (original audit, 17 findings)
- `docs/reviews/post-changes-review.md` (this session's work)
- `docs/reviews/biopython-tradeoffs.md` (migration rationale)
- Current code on this branch

## 1. Resolved this session

| Severity | Item | Verification |
|---|---|---|
| MISSING | R11c VHH/single-domain antibody detection | `pipeline._median_cdr3_length_by_chain` + `process.tpl.tengo` VHH rule + `TestR11cStats` |
| MISSING | Peptide ε reduced description | `process.tpl.tengo:ered_peptide` |
| MISSING | Peptide aliphatic index description | `process.tpl.tengo:aliphatic_peptide` |
| MISSING | Fv ε reduced description | `process.tpl.tengo:ered_Fv` |
| DEVIATION | Implementation reimplemented H-H instead of using BioPython | `properties.py` rewritten + `biopython-tradeoffs.md` |
| DEVIATION | CDR-L3 description not differentiated from CDR-H3 | `process.tpl.tengo:cdr3ChargeDesc` / `cdr3GravyDesc` per-chain maps |

127 / 127 Python tests pass. Coverage 99%. Block builds clean (`pnpm run build:dev` 7/7).

## 2. CRITICAL — release-blocking

None. All MISSING items from the comprehensive review are resolved.

## 3. HIGH — spec compliance gaps to document or fix

### 3.1 `spec-deviations.md` is incomplete

`docs/spec-deviations.md` records SD-001 / SD-002 / SD-003 only. The
deviations below are real-but-functional and should be documented as
SD-004…SD-006 (or fixed):

| Proposed | Item | Source |
|---|---|---|
| SD-004 | TSV peptide column named `peptide_seq`, spec says `sequence` | review §3 |
| SD-005 | TSV antibody schema omits `receptor_type` column (carried via `plan.json`) | review §3 |
| SD-006 | Key axis selected by last-axis index, spec says first-matching axis scan (R1a) | review §2 |
| SD-007 | Extra axis pattern `pl7.app/vdj/clonotypeKey` accepted alongside spec's named anchors | review §2 |

**Action:** Either add SD-004…SD-007 entries to `spec-deviations.md`
matching the SD-001…003 format, or fix the implementation. Document
preferred since these are functionally correct.

### 3.2 R11b edge case — partial coverage with no CDR3 silently emits no info

`main.tpl.tengo:215-247`: when a chain has 1-6 of 7 required regions but
NO CDR3, it falls into none of the message paths (the `else if feats["CDR3"]`
guards both partial and CDR3-only branches). Result: silent fallthrough,
no info annotation.

**Action:** Add a guard for "partial regions seen but no CDR3" — emit a
"CDR3 absent for $chain — no per-chain properties computed" message.

### 3.3 `blockId` mutation stamps both `propertiesPf` and `exports.properties`

`process.tpl.tengo:307-315` mutates `col.spec.domain` in place inside the
loop that builds `exportColumns`. Because `outputSpecs` and `exportSpecs`
share references to the same column dicts, the `pl7.app/blockId` annotation
that should be export-only ends up on `propertiesPf` columns too.

**Action:** Deep-clone before mutating, or build two separate column lists.
Subtle correctness issue — flagged in review §10.

## 4. MEDIUM — deviations and cleanup

### 4.1 Pipeline not vectorised (NumPy)

`pipeline.py:_compute_row_for` runs row-by-row through pure-Python
BioPython calls. Spec L520-525 mentions vectorisation as a goal but does
not block on it. Performance only — not correctness.

**Action:** Defer until measured perf becomes a problem. Vectorising
BioPython calls is non-trivial since `ProteinAnalysis` operates on single
sequences.

### 4.2 Spec example "RRRRRRRRRR emits NA" mathematically incorrect

Spec L535 example claims a 10-Arg sequence emits NA. With IPC 2.0 peptide
pKa, charge crosses zero around pH 13.5 — `R*10` does have a defined pI.
Implementation correctly computes it; the spec example is wrong.

**Action:** Spec issue, not implementation. Surface to spec author for
correction; recommend updating the example to `R*50` or `R*150` (which
genuinely lacks a zero crossing in [0, 14]).

### 4.3 Dead lookup tables in `aa_tables.py` and `pka_tables.py`

After the BioPython migration, the following constants are unused:

- `aa_tables.py`: `KD_SCALE`, `AVG_RESIDUE_MASS`, `H2O_AVG_MASS`,
  `AROMATIC_AAS`, `EC_TYR`, `EC_TRP`, `EC_DISULFIDE`
- `pka_tables.py`: `ACIDIC_AAS`, `BASIC_AAS`

Only `STANDARD_AAS` / `STANDARD_AA_SET` from `aa_tables.py` and `PKaSet` /
`IPC2_PEPTIDE` / `IPC2_PROTEIN` from `pka_tables.py` are still used.

**Action:** Delete unused constants. Low-risk cleanup — 1 commit.

### 4.4 `defaultBlockLabel` dead code path

`model/src/types.ts` defines `defaultBlockLabel` on `BlockData` and
`BlockArgs`, but `index.ts` never reads it from `ctx.data.defaultBlockLabel`
in `.title()`. Field is plumbed but inert.

**Action:** Either:
- Remove the field, the watcher (if any), and the args projection (simplest); or
- Read `ctx.data.defaultBlockLabel` from `.title(...)` and remove the field
  from args projection (UI-only state).

Review §11 recommendation #3.

## 5. LOW — extensions / minor items

| Item | Status | Notes |
|---|---|---|
| Bulk MiXCR fallback chain undefined → "A" | EXTENSION (kept) | Compatibility with non-scClonotype data; document in spec or deviations |
| SD-003 receptor read from axis domain | EXTENSION (already documented as SD-003) | OK |
| MW values shifted ~0.07-0.13 Da per chain on BioPython migration | Documented | Below `.1f` display precision |
| `aa_tables.py` and `pka_tables.py` constants dead | See 4.3 | Cleanup only |

## 6. Test infrastructure outstanding

### 6.1 Workflow integration tests are all `.todo`

`test/src/wf.test.ts` enumerates the right checks (R1-R16) but every
`it.todo(...)` — none execute. Per review §11 #2 this is the highest-
leverage testing gap because the pframe / xsv / model integration is the
most surface area not covered by Python or corpus tests.

**Action:** Implement the `.todo` cases against `blockTest` (the SDK's
template runner). Reference: existing patterns in
`blocks/clonotype-clustering/test/`.

**Estimated effort:** 1-2 days. Each R item is a small input fixture +
output assertion.

### 6.2 R11c integration test not yet covered

The new `info.messages` output now contains the VHH message under spec
conditions. No integration test exercises it.

**Action:** Add a `wf.test.ts` case feeding a heavy-only-IG dataset with
median CDR-H3 ≥ 16 and asserting `info.messages` includes the VHH
substring. Block on 6.1.

## 7. M3 validation status

Spec M3 (`README.md` L545-575) calls for external cross-validation:

| Criterion | Status | Notes |
|---|---|---|
| Charge / GRAVY / MW / EC for ≥5 VH vs BioPython, 2 dp | INHERENT | Implementation now IS BioPython — no need to validate against itself |
| Instability for ≥3 peptides vs `ProteinAnalysis.instability_index()` | INHERENT | Same |
| pI for ≥5 VH within 0.1 pH units of IPC 2.0 reference | NOT VERIFIED | IPC 2.0 reference is the paper's supplementary data; need to compare against published values |
| VL pI on ≥2 sequences | PARTIAL | One golden value pinned in `test_golden_values.py` |
| Fv charge / Fv pI manually verified ≥2 paired sequences | PARTIAL | Pinned values; no external manual computation log |
| CDR-L3 charge by manual H-H ≥3 sequences | NOT VERIFIED | |
| CDR-H3 charge by manual computation ≥10 sequences | NOT VERIFIED | |
| Aliphatic for ≥3 VH with closed-form check | PARTIAL | Closed-form on `ALVI`; one VH pinned |
| CDR-H3 length counting convention vs MiXCR IMGT | NOT VERIFIED | Block doesn't emit CDR3 length as a column (only median in `stats.json`) |
| AstraZeneca / BMS dataset end-to-end | NOT IN SCOPE | External client validation |

The BioPython migration collapses many M3 criteria — the implementation
*is* BioPython, so most of M3 reduces to: validate IPC 2.0 pKa values
against the paper, and validate the wrapper layer (Cys exclusion, [0, 14]
range, instability floor, Fv pairing).

**Action:**
1. Cross-check IPC 2.0 pKa values in `pka_tables.py` against the Kozlowski
   2021 supplementary data (one-time, ~1 hour).
2. Add a small validation script under
   `software/tests/validation/test_ipc2_paper_values.py` that pins ≥5 VH
   pI values from a published IPC 2.0 reference set if available.
3. Defer AstraZeneca / BMS external dataset validation to a follow-up — no
   internal artifact required.

## 8. Recommended priority order for next session

1. **3.1 SD-004…SD-007** — quick documentation, closes the spec-deviation
   audit trail. ~30 min.
2. **3.3 `blockId` mutation** — subtle correctness fix. ~1 hour including
   test.
3. **3.2 R11b no-CDR3 edge case** — small Tengo guard. ~30 min.
4. **4.3 dead-code cleanup** — `aa_tables.py` / `pka_tables.py`. ~15 min.
5. **6.1 wire up `wf.test.ts`** — biggest leverage; closes the integration
   test gap and de-risks future refactors. 1-2 days.
6. **7 IPC 2.0 paper cross-check** — release readiness. ~1 hour.
7. **4.4 `defaultBlockLabel` cleanup** — low priority; either path is
   clean. ~30 min.

The first six items together would put the block in a release-ready
state vs the comprehensive review's recommendations.

## 9. Items intentionally not on this list

- Spec text issues (4.2 RRRRRRRRRR example) — flagged for spec author, not
  implementation work.
- "Pipeline not vectorised" (4.1) — performance, not correctness; no
  measured pain.
- TCR mode label coverage — already PASS in the comprehensive review.
- V3 model migration — already complete on this branch (prior commit).
- Dedup byte-stability — already complete on this branch (prior commits).
