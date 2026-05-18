/**
 * Block-test helpers for sequence-properties.
 *
 * Two upstream paths planned:
 *
 *   1. Synthetic — xsv-import publishes a hand-crafted anchor PColumn from a
 *      TSV fixture. Fast, deterministic, supports any axis pattern + edge
 *      case via crafted TSVs (see ./fixtures/scenarios.ts). Currently
 *      blocked: @platforma-open/milaboratories.xsv-import@1.0.5 ships an
 *      empty tarball (LICENSE + package.json only), so its blockSpec can't
 *      be imported. Resurrect once a fixed version is published or another
 *      synthetic-publisher block is identified.
 *
 *   2. MiXCR canary — `setupMixcrAnchor` runs the real samples-and-data +
 *      mixcr-clonotyping-2 pipeline against fastq fixtures. Slow, but is
 *      the only path currently runnable end-to-end. Used by the canary
 *      test (kept as `it.todo` until the upstream awaitBlockDone timing is
 *      confirmed against this workspace's local platforma).
 *
 * Each blockTest spins up a fresh platforma container, so co-locating
 * multiple assertions per test is the standard cost optimization. Helpers
 * here are deliberately small composable steps so individual tests can
 * inline the parts they need to assert on.
 */

import { blockSpec as samplesAndDataBlockSpec } from '@platforma-open/milaboratories.samples-and-data';
import { blockSpec as mixcrClonotypingBlockSpec } from '@platforma-open/milaboratories.mixcr-clonotyping-2';
import { blockSpec as seqPropsBlockSpec } from 'this-block';
import { uniquePlId } from '@platforma-sdk/model';
import type { ML, RawHelpers } from '@platforma-sdk/test';
import type { expect as vitestExpect } from 'vitest';

export type TestCtx = {
  rawPrj: ML.Project;
  ml: ML.MiddleLayer;
  helpers: RawHelpers;
  expect: typeof vitestExpect;
};

/**
 * Add the sequence-properties block under test.
 */
export async function addSequenceProperties(ctx: TestCtx, label = 'Sequence Properties'): Promise<string> {
  return await ctx.rawPrj.addBlock(label, seqPropsBlockSpec);
}

/**
 * Add samples-and-data + mixcr-clonotyping-2 wired to fastq fixtures, plus
 * two sequence-properties blocks pending an `inputAnchor` from the caller.
 *
 * Locks PR #9's canonical-specs-in / stamp-out invariant. If a future change
 * leaks blockId back into a pure template's inputs (or into the column specs
 * passed to `xsv.importFile`), two co-instances will trigger CIDConflictError
 * — this helper makes that regression surface as a test failure.
 */
export async function setupTwoSeqPropsCoInstances(
  ctx: TestCtx,
  opts: {
    preset?: string;
    chains?: string[];
    r1Path: string;
    r2Path: string;
  },
): Promise<{
  sndBlockId: string;
  clonotypingBlockId: string;
  seqPropsBlockIdA: string;
  seqPropsBlockIdB: string;
}> {
  const preset = opts.preset ?? '10x-sc-xcr-vdj-rhapsody';
  const chains = opts.chains ?? ['IGHeavy', 'IGLight'];

  const { rawPrj, helpers } = ctx;
  const sndBlockId = await rawPrj.addBlock('Samples & Data', samplesAndDataBlockSpec);
  const clonotypingBlockId = await rawPrj.addBlock('MiXCR Clonotyping', mixcrClonotypingBlockSpec);
  const seqPropsBlockIdA = await addSequenceProperties(ctx, 'Sequence Properties A');
  const seqPropsBlockIdB = await addSequenceProperties(ctx, 'Sequence Properties B');

  const sample1Id = uniquePlId();
  const dataset1Id = uniquePlId();
  const r1Handle = await helpers.getLocalFileHandle(opts.r1Path);
  const r2Handle = await helpers.getLocalFileHandle(opts.r2Path);

  await rawPrj.setBlockArgs(sndBlockId, {
    metadata: [],
    sampleIds: [sample1Id],
    sampleLabelColumnLabel: 'Sample Name',
    sampleLabels: { [sample1Id]: 'Sample 1' },
    datasets: [
      {
        id: dataset1Id,
        label: 'Dataset 1',
        content: {
          type: 'Fastq',
          readIndices: ['R1', 'R2'],
          gzipped: true,
          data: { [sample1Id]: { R1: r1Handle, R2: r2Handle } },
        },
      },
    ],
  });

  await rawPrj.setBlockArgs(clonotypingBlockId, {
    preset: { type: 'name', name: preset },
    chains,
  });

  return { sndBlockId, clonotypingBlockId, seqPropsBlockIdA, seqPropsBlockIdB };
}

/**
 * Add samples-and-data + mixcr-clonotyping-2 wired to fastq fixtures.
 * Single-cell IG preset by default — for the MiXCR canary test.
 *
 * Fastq fixtures are NOT bundled — drop appropriate small paired-end
 * fastq.gz files into `test/assets/` (e.g. a single-cell IG slice from
 * SRA) and pass their paths via `opts.r1Path` / `opts.r2Path`. Sibling
 * blocks like mixcr-clonotyping vendor SRR-prefixed fixtures in their
 * own test/assets — those can be reused here when the canary lands.
 *
 * Returns block ids; caller drives runs + assertions.
 */
export async function setupMixcrAnchor(
  ctx: TestCtx,
  opts: {
    preset?: string;
    chains?: string[];
    r1Path: string;
    r2Path: string;
  },
): Promise<{ sndBlockId: string; clonotypingBlockId: string; seqPropsBlockId: string }> {
  const preset = opts.preset ?? '10x-sc-xcr-vdj-rhapsody';
  const chains = opts.chains ?? ['IGHeavy', 'IGLight'];

  const { rawPrj, helpers } = ctx;
  const sndBlockId = await rawPrj.addBlock('Samples & Data', samplesAndDataBlockSpec);
  const clonotypingBlockId = await rawPrj.addBlock('MiXCR Clonotyping', mixcrClonotypingBlockSpec);
  const seqPropsBlockId = await addSequenceProperties(ctx);

  const sample1Id = uniquePlId();
  const dataset1Id = uniquePlId();
  const r1Handle = await helpers.getLocalFileHandle(opts.r1Path);
  const r2Handle = await helpers.getLocalFileHandle(opts.r2Path);

  await rawPrj.setBlockArgs(sndBlockId, {
    metadata: [],
    sampleIds: [sample1Id],
    sampleLabelColumnLabel: 'Sample Name',
    sampleLabels: { [sample1Id]: 'Sample 1' },
    datasets: [
      {
        id: dataset1Id,
        label: 'Dataset 1',
        content: {
          type: 'Fastq',
          readIndices: ['R1', 'R2'],
          gzipped: true,
          data: { [sample1Id]: { R1: r1Handle, R2: r2Handle } },
        },
      },
    ],
  });

  await rawPrj.setBlockArgs(clonotypingBlockId, {
    preset: { type: 'name', name: preset },
    chains,
  });

  return { sndBlockId, clonotypingBlockId, seqPropsBlockId };
}
