# Punchlist — Sequence Properties Block

**Date:** 2026-04-30
**Branch:** `paulnewling/MILAB-6140_init-block`
**Inputs reviewed:**
- `docs/reviews/punchlist-2026-04-29.md` (prior session's open items)
- `docs/reviews/comprehensive-spec-review.md` (17-finding audit)
- `docs/reviews/post-changes-review.md` (session 1 closeout)
- `docs/reviews/perf-review.md` (Tier A/B/C performance items)
- `docs/reviews/biopython-tradeoffs.md` (migration rationale)
- `docs/plans/columnar-and-quantization-plan.md` (Polars rewrite plan)
- `docs/plans/dedup-plan.md` (canonical encoding, applied)
- `docs/investigations/dedup-investigation.md` (background for dedup-plan)
- `docs/spec-deviations.md` (SD-001…SD-007, all still active)
- Current code on this branch

## 1. Status Snapshot

| Area | Status |
|---|---|
| Spec deviations (SD-001…SD-007) | All 7 still in code — none reconciled, all still need to live in `spec-deviations.md` |
| R11c VHH detection | Implemented in code; NOT covered by integration test |
| IPC 2.0 pKa values | Corrected in `pka_tables.py` (session 2); cross-checked against the IPC 2.0 webserver — 5 VH + 2 VL match within 0.05 pH (well under spec's 0.1) |
| M3 cross-validation | Discharged for charge / pI / Fv / aliphatic via `test_m3_validation.py` (38 tests, all pass). CDR-H3 length conv. and AZ/BMS dataset still out-of-tree |
| BioPython migration | Done for closed-form properties (charge / GRAVY / MW / EC / instability / aromaticity / AA fractions) |
| Dedup byte-stability | Applied (`canonical.encode` + `smart.createValueResource` + `maps.getKeys`) |
| `wf.test.ts` integration tests | 41/41 still `.todo` — zero actual tests |
| Polars columnar pipeline | DEFERRED (perf, not correctness; no measured pain) |
| Python tests | 165/165 pass (was 127; +38 M3 cross-checks) |
| Build | Clean (`pnpm run build:dev` 7/7) |

## 2. Outstanding — Open Items

### 2.1 Wire up `test/src/wf.test.ts` (highest leverage)

All 41 cases are `it.todo(...)` — zero actual workflow integration coverage. This is the largest untested surface (pframe / xsv / model integration).

- **Reference pattern:** `blocks/clonotype-clustering/test/`.
- **Effort:** 1-2 days. Each R item is a small input fixture + output assertion.
- **Unblocks:** §2.2 (R11c integration coverage).

### 2.2 R11c VHH integration test

`info.messages` now carries the VHH message (verified via `TestR11cStats`), but no end-to-end test exercises the workflow→model→info path under VHH conditions. Blocked on §2.1.

### 2.3 IPC 2.0 webserver cross-check — DONE (2026-05-01)

5 VH and 2 VL sequences submitted to the IPC 2.0 webserver (`https://ipc2.mimuw.edu.pl/`); published `IPC2_protein` consensus values pinned in `software/tests/unit/test_m3_validation.py` (`TestVHPIWebserverCrosscheck`, `TestVLPIWebserverCrosscheck`). All 7 sequences match within 0.05 pH (spec tolerance was 0.1).

Cross-check is run with `include_cys=True` to match the webserver convention; production VH/VL pI uses `include_cys=False` per spec L535-540 (disulfide-bonded). The spec-convention values are separately pinned in `TestVHVLPISpecConvention` so refactors that change Cys handling fail loudly.

### 2.4 M3 external validation — DONE (2026-05-01)

`software/tests/unit/test_m3_validation.py` discharges the M3 acceptance criteria via two cross-check shapes:

1. **Pinned IPC 2.0 webserver values** for VH/VL pI (§2.3).
2. **Independent textbook Henderson-Hasselbalch reference** (`_ref_charge_hh`, pure-Python, no BioPython) for charge cross-checks. Production code IS BioPython under the hood, so this catches a BioPython formula bug.

Status now:

| Criterion | Status |
|---|---|
| Charge / GRAVY / MW / EC for ≥5 VH vs BioPython, 2 dp | INHERENT (impl IS BioPython) |
| Instability for ≥3 peptides vs `ProteinAnalysis.instability_index()` | INHERENT |
| pI for ≥5 VH within 0.1 pH units of IPC 2.0 reference | DONE — 5 VH pinned in `TestVHPIWebserverCrosscheck`, all within 0.05 pH |
| VL pI on ≥2 sequences | DONE — 2 VL pinned in `TestVLPIWebserverCrosscheck` and `TestVHVLPISpecConvention` |
| Fv charge / Fv pI manually verified ≥2 paired sequences | DONE — synthetic + trastuzumab pairs in `TestFvManualVerification` against `_manual_fv_charge` / `_manual_fv_pi` |
| CDR-L3 charge by manual H-H ≥3 sequences | DONE — 5 sequences in `TestCDRL3ChargeManualHH` |
| CDR-H3 charge by manual computation ≥10 sequences | DONE — 11 sequences in `TestCDRH3ChargeManualHH` |
| Aliphatic for ≥3 VH with closed-form check | DONE — 4 VH in `TestAliphaticIndexClosedForm` against `_closed_form_aliphatic` |
| CDR-H3 length counting convention vs MiXCR IMGT | NOT VERIFIED — block emits median in `stats.json`, not as a column |
| AstraZeneca / BMS dataset end-to-end | OUT OF SCOPE (external client validation) |

### 2.5 Spec correction — `RRRRRRRRRR` example

Spec README L535 claims `R*10` emits NA. Under IPC 2.0 peptide pKa, `R*10` has a defined pI (~12.79). Spec issue, not implementation. Surface to spec author for L535 update (recommend `R*150` or specify protein pKa scale).

### 2.6 Performance — Tier A (perf-review §1, 2)

Pipeline currently per-row; no measured pain reported, but two cheap wins identified:

1. **`Counter`-based `aa_counts`** — replace per-char Python loop with `collections.Counter(seq.upper())` mapped to standard codes. ~30 min, **~2-3× on `aa_counts`**, benefits `extinction_coefficients` / `aliphatic_index` / `aromaticity` / `aa_fractions`. Trivial.
2. **`SequenceContext` per row** — compute upper-cased + cleaned + counts ONCE per input sequence; pass into property fns. ~3 hours, **~3-5× on peptide-mode pipeline; ~2× on antibody**. Touches every property fn signature; golden values + corpus tests pin behaviour.

Combined: ~1 hour for #1 alone; ~3.5 hours for both → **~3× end-to-end** (perf-review §1.).

### 2.7 Performance — Tier B / C (Polars columnar plan)

`docs/plans/columnar-and-quantization-plan.md` is a 7-phase, ~17-hour plan to rewrite the pipeline as count-matrix + polars expressions + numpy bisection. Rationale (Tier B/C in `perf-review.md`):
- ~5-10× when datasets reach ~10⁵ entities.
- ~50-100× when single-cell datasets land.

**Trigger to start:** measured perf pain at scale, or prior to processing a multi-day SC dataset. Until then, keep deferred.

### 2.8 Spec deviations still live

All 7 SD entries remain in `docs/spec-deviations.md`. Cross-reference:

| SD | Subject | Why kept |
|---|---|---|
| SD-001 | Skip secondary alleles in single-cell paired data | Smallest fix; no schema change. Revisit if customer requests secondary-allele properties. |
| SD-002 | Treat `FR4InFrame` as `FR4` | One-line normalisation. Long-term: switch to `VDJRegionInFrame` (Option B in SD-002). |
| SD-003 | Read receptor from clonotypeKey axis domain | Two-source fallback; correct on real MiXCR data. Spec wording incomplete. |
| SD-004 | TSV peptide column named `peptide_seq` (vs spec `sequence`) | Internal Tengo↔Python contract; not exposed downstream. No further action. |
| SD-005 | TSV antibody schema omits `receptor_type` (carried via `plan.json`) | Correct shape for per-run constants. No further action. |
| SD-006 | Last-axis selection vs first-matching-axis scan | Correct for every observed input shape. Revisit if upstream emits non-last entity-key axis. |
| SD-007 | Accept `clonotypeKey` alongside `cloneId` | Forward+backward MiXCR compatibility. No further action. |

## 3. Recommended Priority Order

1. **§2.1 wire up `wf.test.ts`** — biggest leverage; closes the integration test gap; unblocks §2.2.
2. **§2.5 spec author conversation** — L535 example correction (out of code's hands).
3. **§2.4 row "CDR-H3 length counting convention vs MiXCR IMGT"** — only undischarged item in the M3 table; block emits median in `stats.json` rather than as a column, so verification needs to run against MiXCR IMGT output and compare counting.

After §2.1 and §2.2 the block is in a release-ready state vs the comprehensive review's recommendations.

§2.6 (Tier-A perf) and §2.7 (Polars rewrite) stay deferred — perf, not correctness; no measured pain.

## 4. Items Intentionally Not on This List

- **Pipeline not vectorised at the polars level** — perf, not correctness; deferred per §2.7 trigger condition.
- **Tier-A perf wins** (§2.6) — same; ~3× speedup on the table when measured pain emerges.
- **TCR mode label coverage** — already PASS.
- **V3 model migration** — complete on this branch.
- **Dedup byte-stability** — complete on this branch (canonical encoding, sorted-key map iteration, `infoBlob`).
- **R11c VHH detection (impl)** — complete on this branch; only the integration test remains (§2.2).
- **IPC 2.0 webserver cross-check** — DONE 2026-05-01, see §2.3.
- **M3 manual cross-checks** (CDR-H3, CDR-L3, Fv, aliphatic) — DONE 2026-05-01, see §2.4.
- **AA composition antibody mode** — deferred per spec L498 default.
- **`blockId` mutation** — fixed on this branch.
- **R11b silent fallthrough** — fixed on this branch.
- **CDR-L3 description differentiation** — fixed on this branch.
- **Spec deviations SD-001…SD-007** — all documented, all intentional; see §2.8.
