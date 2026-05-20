/**
 * Block integration tests for sequence-properties.
 *
 * Tests run on a real platforma backend (`PL_ADDRESS` configured via
 * `.test_auth.json`). Each `blockTest` spins up a fresh project; teardown
 * is automatic.
 *
 * Two upstream paths (see helpers.ts):
 *   - Synthetic: xsv-import publishes a hand-crafted anchor PColumn from a
 *     committed TSV fixture. Fast, deterministic, covers any axis pattern.
 *   - MiXCR canary: samples-and-data + mixcr-clonotyping-2 against fastq
 *     fixtures. Slow; guards the synthetic axis specs against MiXCR drift.
 *
 * Status: scaffolding + idle smoke land in this slice. Synthetic /
 * canary execution paths are wired (helpers, fixtures, deps) but the
 * upstream `xsv-import` block did not finish in 120s during this
 * session's first dry-run — root cause not yet identified (file-handle
 * import, local backend file-driver registration, or our spec). Tests
 * targeting those upstream runs are kept as `it.todo` with the existing
 * spec-rule comments. Fill them in once the upstream timing or wiring is
 * resolved — the helpers + scenarios in `./fixtures/scenarios.ts` are
 * ready.
 *
 * Running: `cd blocks/sequence-properties && pnpm --filter
 * @platforma-open/milaboratories.sequence-properties.test test`
 */

import type { platforma } from '@platforma-open/milaboratories.sequence-properties.model';
import type { InferBlockState } from '@platforma-sdk/model';
import { wrapOutputs } from '@platforma-sdk/model';
import { awaitStableState, blockTest } from '@platforma-sdk/test';
import { blockSpec as seqPropsBlockSpec } from 'this-block';
import { existsSync } from 'node:fs';
import { resolve } from 'node:path';
import { describe, it } from 'vitest';
import { setupTwoSeqPropsCoInstances } from './helpers';

const FIXTURE_R1 = resolve(__dirname, '../assets/canary-sc-ig.R1.fastq.gz');
const FIXTURE_R2 = resolve(__dirname, '../assets/canary-sc-ig.R2.fastq.gz');
const HAS_CANARY_FIXTURES = existsSync(FIXTURE_R1) && existsSync(FIXTURE_R2);

// ---------------------------------------------------------------------------
// Cheap idle-state smoke (model + UI baseline)
// ---------------------------------------------------------------------------

// Idle block has no inputAnchor selected and no upstream. inputOptions must
// be an empty array; workflow-driven outputs absent. The model's args lambda
// throws on missing inputAnchor (caught by the runtime) — observable as
// argsValid === false / no run. This is the universal smoke test —
// confirms test scaffolding + backend connection + dev-block load all work.
blockTest('empty inputs', { timeout: 30000 }, async ({ rawPrj, expect }) => {
  const blockId = await rawPrj.addBlock('Sequence Properties', seqPropsBlockSpec);
  const stableState = (await awaitStableState(
    rawPrj.getBlockState(blockId),
    20000,
  )) as InferBlockState<typeof platforma>;

  expect(stableState.outputs).toMatchObject({
    inputOptions: { ok: true, value: [] },
  });
  const outputs = wrapOutputs(stableState.outputs);
  expect(outputs.propertiesTable).toBeUndefined();
  expect(outputs.info).toBeUndefined();
  expect(outputs.isRunning).toBeFalsy();
});

// ---------------------------------------------------------------------------
// Modality detection (R1, R1a, R1b)
// ---------------------------------------------------------------------------

describe('modality detection', () => {
  it.todo('detects peptide mode from variantKey + extractionRunId domain'); // R1, peptideScenario
  it.todo('detects antibody/TCR single-cell from vdj/scClonotypeKey axis'); // R1, scIgFullScenario
  it.todo('detects antibody/TCR legacy bulk from vdj/cloneId axis'); // R1
  it.todo('detects antibody/TCR legacy bulk from vdj/clonotypeKey axis'); // R1
  it.todo('panics with R1a message when input has no recognised sequence key'); // R1a, bogusAxisScenario
  it.todo('emits output PColumns with input entity axis verbatim'); // R1b
});

// ---------------------------------------------------------------------------
// Peptide mode (R6, R7, R8, R9, R16) — column-shape assertions need
// `ml.driverKit.pFrameDriver.listColumns(handle)` to inspect the
// fullPframeHandle. Pattern from titeseq + clonotype-browser tests.
// ---------------------------------------------------------------------------

describe('peptide mode output', () => {
  it.todo('emits the 9 peptide scalar PColumns with feature=peptide domain'); // R6
  it.todo('emits 2-axis pl7.app/aaFraction PColumn with aminoAcid axis'); // R7
  it.todo('AA fraction PColumn has no isScore annotation'); // R16
  it.todo('extinction coefficients oxidised and reduced both emitted'); // R8
  it.todo('instability NA for peptides under 10 aa (per-row, not column-wide)'); // R9
});

// ---------------------------------------------------------------------------
// Antibody / TCR mode (R10, R11, R12, R13, R14)
// ---------------------------------------------------------------------------

describe('antibody/TCR coverage handling', () => {
  it.todo('emits CDR3 charge + hydrophobicity per chain when CDR3 present'); // R10, scIgFullScenario
  it.todo('emits full-chain + Fv columns on full IG paired coverage'); // R11+R12, scIgFullScenario
  it.todo('CDR3-only input emits R11a info message and no full-chain columns'); // R11a, scIgCdr3OnlyScenario
  it.todo('partial-region info message uses receptor-aware chain label'); // R11b
  it.todo('missing region for one clone yields NA for that row only');
  it.todo('TCR receptor never emits Fv columns regardless of coverage'); // R12, scTcrAbFullScenario
  it.todo('γδ TCR input emits the γδ-labels info message');
  it.todo('chain labels adapt to receptor while name + chain domain stay constant'); // R13
  it.todo('bulk MiXCR IGHeavy axis derives receptor=IG from pl7.app/vdj/chain, no R13b warning'); // SD-008, bulkIgHeavyScenario
  it.todo('bulk MiXCR TCRAlpha axis derives receptor=TCRAB from pl7.app/vdj/chain, no R13b warning'); // SD-008, bulkTcrAlphaScenario
  it.todo('R13b warning fires only when neither receptor nor recognised chain key is present'); // R13b post-SD-008
  it.todo('isScore annotation on the spec-mandated columns only'); // R14
  it.todo('rankingOrder=increasing only on hydrophobicity columns'); // R14
  it.todo('no defaultCutoff annotation anywhere in the output'); // R15
});

// ---------------------------------------------------------------------------
// Edge cases per spec "Defaults and edge cases" table — synthetic-corpus
// route required (real fastq can't easily encode stop codons / polybasic
// peptides). Use xsv-import scenarios with crafted TSVs once that path is
// verified.
// ---------------------------------------------------------------------------

describe('edge cases', () => {
  it.todo('stop codon in any region NAs the whole reconstructed chain');
  it.todo('non-standard residues do not affect fractions / pI / GRAVY denominator');
  it.todo('polybasic with no zero crossing emits NA pI without error');
  it.todo('no aromatic residues emits ε = 0 (not NA)');
  it.todo('single-cell single chain emits available chain, Fv NA');
});

// ---------------------------------------------------------------------------
// Result-pool integration & enrichment (Lead Selection downstream)
// ---------------------------------------------------------------------------

describe('downstream consumption', () => {
  it.todo('score columns are discoverable by a downstream block via the result pool');
  it.todo('output PColumns carry pl7.app/trace stamped with this block\'s label');
  it.todo('export PFrame stamps blockId on score column domains');
});

// ---------------------------------------------------------------------------
// Model / UI behaviors
// ---------------------------------------------------------------------------

describe('model + UI', () => {
  it.todo('inputOptions surfaces a published vdj/scClonotypeKey anchor'); // synthetic
  it.todo('run button disabled when inputAnchor not set'); // covered by 'empty inputs' assertions
  it.todo('inputOptions filters to the 4 supported abundance anchor patterns');
  it.todo('info output exposes mode, receptor, coverageTier, messages array');
  it.todo('changing inputAnchor invalidates and re-runs cleanly');
  it.todo('title is the configured static label');
});

// ---------------------------------------------------------------------------
// Dedup wiring (see docs/dedup-plan.md)
// ---------------------------------------------------------------------------

describe('dedup', () => {
  it.todo('second project on identical upstream lands on Done via dedup');
  it.todo('changed upstream input breaks dedup and triggers fresh run');

  // Regression for PR #9. Two co-instances on identical input must:
  //   (1) reach Done with no CIDConflictError on any output, and
  //   (2) emit per-instance pl7.app/trace.id so downstream pickers
  //       disambiguate.
  //
  // Runs whenever the SC IG fastq fixtures exist at
  // test/assets/canary-sc-ig.R{1,2}.fastq.gz (staged from
  // mixcr-clonotyping's SRR11233625 slices). The synthetic xsv-import
  // route remains blocked on the empty-tarball upstream documented in
  // helpers.ts — the MiXCR canary is the only working path today.
  blockTest.skipIf(!HAS_CANARY_FIXTURES)(
    'two co-instances on identical upstream run without CID conflicts',
    { timeout: 600_000 },
    async ({ expect, rawPrj, ml, helpers }) => {
      const ctx = { expect, rawPrj, ml, helpers };
      const { clonotypingBlockId, seqPropsBlockIdA, seqPropsBlockIdB } =
        await setupTwoSeqPropsCoInstances(ctx, { r1Path: FIXTURE_R1, r2Path: FIXTURE_R2 });

      // Run MiXCR first so its outputs publish into the result pool and
      // become discoverable as seqProps inputOptions.
      await rawPrj.runBlock(clonotypingBlockId);

      const seqPropsAStateWithOpts = (await awaitStableState(
        rawPrj.getBlockState(seqPropsBlockIdA),
        500_000,
      )) as InferBlockState<typeof platforma>;
      const opts = (seqPropsAStateWithOpts.outputs.inputOptions as { value: { ref: unknown; label: string }[] }).value;
      expect(opts?.length).toBeGreaterThan(0);
      const anchorRef = opts[0].ref;

      await rawPrj.setBlockArgs(seqPropsBlockIdA, { inputAnchor: anchorRef });
      await rawPrj.setBlockArgs(seqPropsBlockIdB, { inputAnchor: anchorRef });

      await rawPrj.runBlock(seqPropsBlockIdA);
      await rawPrj.runBlock(seqPropsBlockIdB);

      const stateA = (await awaitStableState(
        rawPrj.getBlockState(seqPropsBlockIdA),
        500_000,
      )) as InferBlockState<typeof platforma>;
      const stateB = (await awaitStableState(
        rawPrj.getBlockState(seqPropsBlockIdB),
        500_000,
      )) as InferBlockState<typeof platforma>;

      // Neither instance should hit a CIDConflictError.
      for (const [name, state] of [['A', stateA], ['B', stateB]] as const) {
        for (const [key, val] of Object.entries(state.outputs)) {
          const v = val as { ok?: boolean; errors?: { message?: string }[] };
          if (v.ok === false) {
            for (const err of v.errors ?? []) {
              expect(
                err.message ?? '',
                `instance ${name} output ${key} carries an error`,
              ).not.toMatch(/CIDConflict|CID conflict/);
            }
          }
        }
      }

      // Both must reach a state where propertiesPfCols resolves — proves the
      // pure-template + xsv.importFile pipeline ran end-to-end for both.
      const aCols = wrapOutputs(stateA.outputs).propertiesPfCols;
      const bCols = wrapOutputs(stateB.outputs).propertiesPfCols;
      expect(aCols).toBeTruthy();
      expect(bCols).toBeTruthy();
      expect(Array.isArray(aCols)).toBe(true);
      expect((aCols as unknown[]).length).toBeGreaterThan(0);
      expect((aCols as unknown[]).length).toBe((bCols as unknown[]).length);

      // Per-instance trace.id is what lets downstream pickers disambiguate
      // co-instances. The label comes from resolveTraceLabel(data) in
      // model/src/label.ts; this test leaves customBlockLabel unset so both
      // labels collapse to the same defaultBlockLabel — assertion is on id,
      // not label.
      const aTrace = ((aCols as { spec: { annotations?: Record<string, string> } }[])[0]
        .spec.annotations ?? {})['pl7.app/trace'];
      const bTrace = ((bCols as { spec: { annotations?: Record<string, string> } }[])[0]
        .spec.annotations ?? {})['pl7.app/trace'];
      expect(aTrace).toBeTruthy();
      expect(bTrace).toBeTruthy();
      const aSelf = JSON.parse(aTrace as string).at(-1) as { type: string; id: string };
      const bSelf = JSON.parse(bTrace as string).at(-1) as { type: string; id: string };
      expect(aSelf.type).toBe('milaboratories.sequence-properties');
      expect(bSelf.type).toBe('milaboratories.sequence-properties');
      expect(aSelf.id).not.toBe(bSelf.id);
    },
  );
});

// ---------------------------------------------------------------------------
// Cross-block composition
// ---------------------------------------------------------------------------

describe('cross-block composition', () => {
  it.todo('canary: real MiXCR fastq → sequence-properties detects sc IG and emits scores');
  it.todo('end-to-end: peptide-extraction → sequence-properties → lead-selection');
  it.todo('two sequence-properties instances in one project disambiguate via blockId');
  it.todo('VDJ chain: full-coverage MiXCR → seq-properties → lead-selection');
});
