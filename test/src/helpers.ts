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
 *   2. MiXCR canary — runs the real samples-and-data + mixcr-clonotyping-2
 *      pipeline against fastq fixtures. Slow, but is the only path currently
 *      runnable end-to-end. Used by the canary test and by
 *      `setupTwoSeqPropsCoInstances` to seed the two-instance dedup test.
 *
 * V3 plumbing: both upstreams (`samples-and-data@^1.20` and
 * `mixcr-clonotyping-2@^2.22`) are `PlatformaV3` blocks and reject
 * `setBlockArgs` with `ModelAPIVersionMismatchError`. The helpers below use
 * `mutateBlockStorage({ operation: 'update-block-data', value: <BlockData> })`
 * which is the V3 update path.
 *
 * Both upstreams are slim facades, so a block is added by its exported
 * `<Block>BlockPointer` and not by the retired `blockSpec`. The block-data
 * literals below carry no `satisfies <Block>BlockData`: a facade bundles its
 * own copy of the brand symbols behind `PlId` and `CanonicalizedJson`, so
 * those types never unify with our catalog SDK's. `mutateBlockStorage` takes
 * the value as `unknown` and the brands are erased at runtime.
 *
 * Each blockTest spins up a fresh platforma container, so co-locating
 * multiple assertions per test is the standard cost optimization.
 */

import { MixcrClonotyping2BlockPointer } from "@platforma-open/milaboratories.mixcr-clonotyping-2";
import { SamplesAndDataBlockPointer } from "@platforma-open/milaboratories.samples-and-data";
import { createPlDataTableStateV2, uniquePlId } from "@platforma-sdk/model";
import { SequencePropertiesBlockPointer } from "this-block";
import type { ML, RawHelpers } from "@platforma-sdk/test";
import { awaitStableState } from "@platforma-sdk/test";
import type { expect as vitestExpect } from "vitest";

export type TestCtx = {
  rawPrj: ML.Project;
  ml: ML.MiddleLayer;
  helpers: RawHelpers;
  expect: typeof vitestExpect;
};

/**
 * Add the sequence-properties block under test.
 */
export async function addSequenceProperties(
  ctx: TestCtx,
  label = "Sequence Properties",
): Promise<string> {
  return await ctx.rawPrj.addBlock(label, SequencePropertiesBlockPointer);
}

/**
 * Configure samples-and-data with a one-sample fastq dataset. Uses V3
 * mutateBlockStorage with the full BlockData payload (V1 setBlockArgs is
 * rejected by samples-and-data@^1.20 — `PlatformaV3` model).
 */
async function configureSamplesAndData(
  ctx: TestCtx,
  sndBlockId: string,
  opts: { r1Path: string; r2Path: string },
): Promise<void> {
  const { rawPrj, helpers } = ctx;
  const sample1Id = uniquePlId();
  const dataset1Id = uniquePlId();
  const r1Handle = await helpers.getLocalFileHandle(opts.r1Path);
  const r2Handle = await helpers.getLocalFileHandle(opts.r2Path);

  await rawPrj.mutateBlockStorage(sndBlockId, {
    operation: "update-block-data",
    value: {
      metadata: [],
      sampleIds: [sample1Id],
      sampleLabelColumnLabel: "Sample Name",
      sampleLabels: { [sample1Id]: "Sample 1" },
      datasets: [
        {
          id: dataset1Id,
          label: "Dataset 1",
          content: {
            type: "Fastq",
            readIndices: ["R1", "R2"],
            gzipped: true,
            data: { [sample1Id]: { R1: r1Handle, R2: r2Handle } },
          },
        },
      ],
      h5adFilesToPreprocess: [],
      seuratFilesToPreprocess: [],
      suggestedImport: false,
    },
  });
}

/**
 * Configure mixcr-clonotyping-2 with the preset + chains, wired to
 * samples-and-data's published output. V3 — requires the upstream input
 * ref, which is read from clonotyping's inputOptions after samples-and-data
 * has run to Done.
 */
async function configureMixcrClonotyping(
  ctx: TestCtx,
  clonotypingBlockId: string,
  preset: string,
  chains: string[],
): Promise<void> {
  const { rawPrj } = ctx;

  // After samples-and-data Done, clonotyping's inputOptions populates from
  // the result pool. Wait for it to stabilize, then take the first option.
  type ClonotypingInputOption = { ref: { __isRef: true; blockId: string; name: string } };
  const clonotypingState = await awaitStableState(rawPrj.getBlockState(clonotypingBlockId), 25000);
  const inputOptions = (
    clonotypingState.outputs as Record<string, { value?: ClonotypingInputOption[] }>
  ).inputOptions?.value;
  if (!inputOptions || inputOptions.length === 0) {
    throw new Error("mixcr-clonotyping-2 inputOptions did not populate after samples-and-data");
  }

  // mixcr-clonotyping-2@^2.22 is built against a current SDK, so its
  // tableState is the same PlDataTableStateV2 shape our catalog SDK emits.
  // The hand-rolled version-5 literal the old pin needed is gone.
  const tableState = createPlDataTableStateV2();

  await rawPrj.mutateBlockStorage(clonotypingBlockId, {
    operation: "update-block-data",
    value: {
      defaultBlockLabel: "",
      customBlockLabel: "",
      input: inputOptions[0].ref,
      preset: { type: "name", name: preset },
      chains,
      tableState,
      runMode: "full",
    },
  });
}

/**
 * Older helper retained for the MiXCR canary test scaffolding. Same V3
 * configuration flow as setupTwoSeqPropsCoInstances but with a single
 * sequence-properties block.
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
  const preset = opts.preset ?? "10x-sc-xcr-vdj-rhapsody";
  const chains = opts.chains ?? ["IGHeavy", "IGLight"];

  const { rawPrj, helpers } = ctx;
  const sndBlockId = await rawPrj.addBlock("Samples & Data", SamplesAndDataBlockPointer);
  const clonotypingBlockId = await rawPrj.addBlock(
    "MiXCR Clonotyping",
    MixcrClonotyping2BlockPointer,
  );
  const seqPropsBlockId = await addSequenceProperties(ctx);

  await configureSamplesAndData(ctx, sndBlockId, opts);
  await rawPrj.runBlock(sndBlockId);
  await helpers.awaitBlockDone(sndBlockId, 30000);

  await configureMixcrClonotyping(ctx, clonotypingBlockId, preset, chains);

  return { sndBlockId, clonotypingBlockId, seqPropsBlockId };
}
