# Post-Changes Review — Sequence Properties Block

**Date:** 2026-04-29
**Scope:** Verify adherence to spec + comprehensive review after this
session's changes (R11c, missing description annotations, CDR-L3
differentiation, BioPython migration).

## Summary

| Severity (pre-session) | Item | Status now |
|---|---|---|
| MISSING | R11c VHH/single-domain detection | RESOLVED |
| MISSING | Peptide ε reduced description annotation | RESOLVED |
| MISSING | Peptide aliphatic index description annotation | RESOLVED |
| MISSING | Fv ε reduced description annotation | RESOLVED |
| DEVIATION | Implementation reimplements H-H instead of BioPython | RESOLVED |
| DEVIATION | CDR-L3 description not differentiated from CDR-H3 | RESOLVED |

Tests: **120 / 120 passing** after migration. `pnpm run build:dev` succeeds.

## Per-item adherence

### R11c — VHH detection (RESOLVED)

**Spec (R11c):** "if receptor type is `"IG"` and only chain `"A"` is present
(no chain `"B"` CDR3) and the median CDR-H3 length across the dataset is
≥ 16 aa, emit a block-level info annotation: '...'."

**Implementation:**
- `pipeline.py:_median_cdr3_length_by_chain` computes per-chain median CDR-H3
  length from the input frame.
- `pipeline.run_antibody_tcr` returns `stats: {medianCdr3Length: {chain: median}}`.
- `main.py` writes `stats.json` (sorted-key, separators=(",",":") for
  byte-stability).
- `main.tpl.tengo` adds `--stats stats.json` + `saveFileContent("stats.json")`,
  computes `chainsWithCdr3`, passes both to `process.tpl.tengo` along with
  Tengo-side `infoMessages`.
- `process.tpl.tengo` reads `args.stats.getDataAsJson()`, applies the VHH
  rule (receptor IG + only chain A has CDR3 + median ≥ 16), appends the
  exact spec message to `infoMessages`, builds the info JSON resource.
- `main.tpl.tengo` retrieves `info` from `processResult.output("info", ...)`.

**Message text** (verified verbatim against spec):
> "Possible VHH/single-domain antibody input detected (heavy chain only;
> CDR-H3 length distribution consistent with VHH). IgG-calibrated CDR-H3
> length thresholds (>15 aa elevated risk, >20 aa high risk) do not apply
> to VHH — disregard these thresholds for nanobody libraries."

### Missing description annotations (RESOLVED)

`process.tpl.tengo`:
- `ered_peptide` description added — verbatim from `pcolumn-spec.md`.
- `aliphatic_peptide` description added — verbatim from `pcolumn-spec.md`.
- `ered_Fv` description added — verbatim from `pcolumn-spec.md`.

### CDR-L3 description differentiation (RESOLVED)

`process.tpl.tengo` introduces per-chain description maps `cdr3ChargeDesc`
and `cdr3GravyDesc`, indexed by chain key (`A` / `B`). Chain A retains the
prior CDR-H3 wording; chain B uses the spec's CDR-L3 wording.

### BioPython migration (RESOLVED)

`properties.py` rewritten as thin wrappers over `Bio.SeqUtils.ProtParam` and
`Bio.SeqUtils.IsoelectricPoint`:
- IPC 2.0 pKa overrides applied per-instance (no module-global monkey-patching).
- Cys-include flag controls whether Cys is in `neg_pKs`.
- pI bisection range broadened to [0, 14] locally — BioPython's `pi()` is
  hard-bracketed [4.05, 12].
- 10-residue floor for instability_index preserved per spec.
- Aliphatic index kept custom (BioPython has no equivalent).
- Fv pI uses BioPython's `charge_at_pH` per chain, with local bisection
  on the chain-pair charge sum.

Detail: `docs/reviews/biopython-tradeoffs.md`.

`instability.py` deleted (BioPython has the same DIWV table).

## Items NOT addressed in this session

These were flagged by the comprehensive review but are out of scope for the
current task ("Python work + missing items") or unaffected by the changes.
They remain valid follow-ups:

| Item | Source | Why deferred |
|---|---|---|
| M3 external validation suite | Review §5, §8 | Separate test infrastructure; spec direction now satisfied (BioPython is reference). The `test_golden_values.py` characterization tests now pin BioPython output, providing a partial substitute. |
| R1a last-axis vs first-matching axis scan | Review §2 | Tengo workflow logic; not flagged as MISSING. Should be added to `spec-deviations.md`. |
| Pipeline not vectorised (NumPy) | Review §3 | Performance — not correctness; deferred per comprehensive review's classification. |
| TSV column name `peptide_seq` vs spec's `sequence` | Review §3 | Cosmetic; Tengo and Python sides consistent. Should be added to `spec-deviations.md` or fix both sides. |
| TSV antibody schema lacks `receptor_type` column | Review §3 | Functionally equivalent (carried via plan.json). Should be added to `spec-deviations.md`. |
| `defaultBlockLabel` dead code path | Review §11 | Model layer cleanup; outside this session's scope. |
| Spec example "RRRRRRRRRR emits NA" mathematically incorrect | Review §5 | Spec issue, not implementation. |
| Per-chain partial coverage with no CDR3 silently emits no message | Review §2 | Subtle edge case; recommend separate fix or doc. |
| `blockId` mutation affects both `propertiesPf` and exports | Review §10 | Subtle correctness issue; unrelated to this session. |

## Recommended follow-up — `spec-deviations.md` updates

The block's `docs/spec-deviations.md` currently records SD-001 / SD-002 /
SD-003. After this session, the BioPython deviation is now resolved, but
the cosmetic / behavioural deviations above should be added. Suggested
new entries:
- SD-004: TSV peptide column named `peptide_seq` (not `sequence`).
- SD-005: TSV antibody schema does not emit `receptor_type` column
  (carried via plan.json instead).
- SD-006: R1a uses last-axis selection rather than first-matching axis scan
  for the key axis. (If kept; otherwise fix the implementation.)

## Files changed in this session

| File | Reason |
|---|---|
| `software/src/properties.py` | Rewrote to use BioPython. |
| `software/src/pipeline.py` | Added `_median_cdr3_length_by_chain`, threaded `stats` dict through `run` / `run_antibody_tcr` / `run_peptide`. |
| `software/src/main.py` | Added `--stats` arg, writes stats JSON. |
| `software/src/instability.py` | Deleted (BioPython has DIWV). |
| `software/pyproject.toml` | Added `biopython>=1.84` to dev group. |
| `software/src/requirements.txt` | Added `biopython==1.84` to runtime deps. |
| `software/tests/unit/test_golden_values.py` | Updated VH MW + Fv MW golden values to BioPython output; updated docstring. |
| `software/tests/unit/test_properties.py` | Updated `test_three_glycines` MW golden value. |
| `software/tests/integration/test_cli.py` | Added `--stats` to all CLI invocations. |
| `workflow/src/main.tpl.tengo` | Added VHH chain detection prep, `--stats` wiring, moved info-blob assembly into process template. |
| `workflow/src/process.tpl.tengo` | Reads stats, applies VHH rule, builds info resource as third output; added missing description annotations; added per-chain CDR3 description maps. |
| `docs/reviews/biopython-tradeoffs.md` | New — documents the BioPython migration's wrapper layer. |
| `docs/reviews/post-changes-review.md` | This file. |

## Why BioPython is preferable here

Beyond spec compliance: BioPython is the canonical reference for ProtParam
and IsoelectricPoint computations. Using it directly:
- Eliminates a class of bugs (DIWV table transcription, Pace constants
  drift, terminus charge sign errors) that custom code would need to
  re-prove correct.
- Gives the M3 validation criteria a clear meaning ("pKa-only differences"
  rather than "formula + pKa differences").
- Reduces our maintenance surface — when BioPython is patched, we get
  the patch.

The only correctness-sensitive places where BioPython could not be used
directly (pI range, Cys exclusion, instability floor, Fv pairing,
empty/stop handling, EC tuple order, AA fraction unit) are documented in
`biopython-tradeoffs.md` as boundary glue, not reimplementation.
