# Punchlist — Sequence Properties Block

**Date:** 2026-04-29
**Branch:** `paulnewling/MILAB-6140_init-block`
**Sessions:**
- Session 1: `b9407ae` (R11c + BioPython migration) → `5902213` (test coverage) → `78b4990` (initial punchlist) → `2e7599c` (principled-coding refactor)
- Session 2: this update — punchlist items 1-4, 6, 7 (skipping #5 wf.test.ts wire-up per user direction)

**Sources audited:**
- `docs/text/work/projects/sequence-properties/README.md` (spec)
- `docs/text/work/projects/sequence-properties/pcolumn-spec.md` (column contracts)
- `docs/reviews/comprehensive-spec-review.md` (original audit, 17 findings)
- `docs/reviews/post-changes-review.md` (session 1 work)
- `docs/reviews/biopython-tradeoffs.md` (migration rationale)
- Current code on this branch

## 1. Resolved across both sessions

| Severity | Item | Verification |
|---|---|---|
| MISSING | R11c VHH/single-domain antibody detection | `pipeline._median_cdr3_length_by_chain` + `process.tpl.tengo` VHH rule + `TestR11cStats` |
| MISSING | Peptide ε reduced description | `process.tpl.tengo:ered_peptide` |
| MISSING | Peptide aliphatic index description | `process.tpl.tengo:aliphatic_peptide` |
| MISSING | Fv ε reduced description | `process.tpl.tengo:ered_Fv` |
| DEVIATION | Implementation reimplemented H-H instead of using BioPython | `properties.py` rewritten + `biopython-tradeoffs.md` |
| DEVIATION | CDR-L3 description not differentiated from CDR-H3 | `process.tpl.tengo:cdr3ChargeDesc` / `cdr3GravyDesc` per-chain maps |
| **CRITICAL** | **IPC 2.0 pKa values were IPC 1.0 with labels swapped** | `pka_tables.py` rewritten with verified IPC 2.0 values from `http://ipc2-isoelectric-point.org/theory.html`. Golden tests + corpus expectations updated to new pI / charge values. |
| HIGH | `spec-deviations.md` missing entries SD-004…SD-007 | Added: peptide_seq column name, omitted receptor_type column, last-axis selection, clonotypeKey alias |
| HIGH | `blockId` mutation leaking to `propertiesPf` | `process.tpl.tengo` now builds export columns via fresh dict construction; `outputSpecs.columns` no longer carries `pl7.app/blockId` |
| HIGH | R11b silent fallthrough when partial regions seen but no CDR3 | New `else` branch in `main.tpl.tengo` chain-collection loop emits a "CDR3 absent" message |
| MEDIUM | Dead constants in `aa_tables.py` / `pka_tables.py` | Purged in session 1; `aa_tables.py` shrank from 78 → 16 lines |
| MEDIUM | `defaultBlockLabel` plumbed but inert | `model/src/index.ts`: `.title()` now reads `ctx.data.defaultBlockLabel`; field removed from `BlockArgs` (UI-only state) |

127 / 127 Python tests pass after IPC 2.0 fix. Coverage 99%. Block builds clean (`pnpm run build:dev` 7/7).

## 2. CRITICAL note — IPC 2.0 pKa value error

**This is the load-bearing finding from session 2.** The spec README L514-525 explicitly warns:

> "Verify each value against the paper before committing. Wrong pKa values
> will quietly produce systematically biased pI / charge values."

Verification was not done before the original commit. The values stored as
`IPC2_PEPTIDE` and `IPC2_PROTEIN` in `pka_tables.py` were actually the IPC 1.0
peptide / protein scales (Kozlowski 2016), with the labels swapped between
the two contexts. Every pI and charge result emitted to date has been
miscalibrated against the spec's intended pKa reference.

**Verification source:** the official IPC 2.0 site
(`http://ipc2-isoelectric-point.org/theory.html`) tabulates IPC2_protein and
IPC2_peptide pKa scales linked directly to the paper's DOI. Our values now
match the official table verbatim.

**Impact of the fix on numerics:**
- VH pI: 7.018 → 6.007 (Δ ≈ 1 unit)
- VL pI: 9.799 → 9.166 (Δ ≈ 0.6 unit)
- Fv pI: 9.331 → 7.634 (Δ ≈ 1.7 units)
- VH charge at pH 7: +0.001 → −0.837 (sign flip)

The polybasic-NA case (`R*150 → NA`) no longer holds under IPC 2.0
because the protein C-terminus pKa (6.065) generates a real zero crossing
even for pure-Arg sequences; expectation updated to a numeric range.

**Cache impact:** all CIDs for `charge_*` and `pi_*` columns will change on
next run. Documented in the changeset.

## 3. Outstanding — open items

### 3.1 Workflow integration tests are 100% `.todo`

`test/src/wf.test.ts` enumerates the right checks (R1–R16) but every case
is `it.todo(...)`. Per the comprehensive review §11 #2 this is the highest-
leverage testing gap because the pframe / xsv / model integration is the
largest surface area not exercised by Python or corpus tests.

**Action:** Implement the `.todo` cases against `blockTest`. Reference
patterns: `blocks/clonotype-clustering/test/`.

**Estimated effort:** 1-2 days. Each R item is a small input fixture +
output assertion.

**Status:** Deferred per user direction in session 2. Tracked here for
follow-up.

### 3.2 R11c integration test not yet covered

The `info.messages` output now contains the VHH message under spec
conditions. No integration test exercises it. Block on 3.1.

### 3.3 M3 external validation status

The BioPython migration collapses many M3 criteria — the implementation
*is* BioPython for the closed-form properties (charge / GRAVY / MW / EC /
instability / aromaticity / AA fractions). For pI and charge, the IPC 2.0
fix in this session brings the pKa values into alignment with the spec
reference.

| Criterion | Status |
|---|---|
| Charge / GRAVY / MW / EC for ≥5 VH vs BioPython, 2 dp | INHERENT (impl IS BioPython) |
| Instability for ≥3 peptides vs `ProteinAnalysis.instability_index()` | INHERENT |
| pI for ≥5 VH within 0.1 pH units of IPC 2.0 reference | NEEDS WEBSERVER CROSS-CHECK — submit ≥5 VH sequences to `ipc2-isoelectric-point.org` and compare against our IPC2_protein output |
| VL pI on ≥2 sequences | PARTIAL — one golden value pinned in `test_golden_values.py` |
| Fv charge / Fv pI manually verified ≥2 paired sequences | PARTIAL |
| CDR-L3 charge by manual H-H ≥3 sequences | NOT VERIFIED |
| CDR-H3 charge by manual computation ≥10 sequences | NOT VERIFIED |
| Aliphatic for ≥3 VH with closed-form check | PARTIAL |
| CDR-H3 length counting convention vs MiXCR IMGT | NOT VERIFIED — block emits median in `stats.json` not as a column |
| AstraZeneca / BMS dataset end-to-end | OUT OF SCOPE (external client validation) |

**Action:** Submit ≥5 VH sequences from `test_golden_values.py` to the IPC
2.0 webserver and pin the published pI values as a separate test (or
update existing golden values to match). ~1 hour.

### 3.4 Spec text issue — `RRRRRRRRRR emits NA`

Spec L535 example claims `R*10` emits NA. Under any reasonable pKa scale
including IPC 2.0 peptide, `R*10` has a defined pI. Spec issue, not
implementation. Surface to spec author.

## 4. Items intentionally not on this list

- "Pipeline not vectorised" — performance, not correctness; no measured
  pain. Defer until measured.
- TCR mode label coverage — already PASS.
- V3 model migration — complete on this branch.
- Dedup byte-stability — complete on this branch.
- TSV peptide column name `peptide_seq` (vs spec's `sequence`) — documented
  as SD-004; no further action.
- `receptor_type` carried via plan.json — documented as SD-005; no further
  action.
- Last-axis selection — documented as SD-006; no further action.
- `pl7.app/vdj/clonotypeKey` alias — documented as SD-007; no further
  action.

## 5. Recommended priority order for next session

1. **3.1 wire up `wf.test.ts`** — biggest leverage; closes the integration
   test gap, de-risks the R11b / R11c / blockId fixes from this session.
2. **3.3 IPC 2.0 webserver cross-check** — release readiness; pin reference
   pI values from the authoritative source.
3. **3.4** — spec author conversation, not implementation work.

After 3.1 and 3.3 the block is in a release-ready state vs the comprehensive
review's recommendations.
