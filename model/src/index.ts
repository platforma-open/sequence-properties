import type { InferOutputsType, PColumnIdAndSpec, PFrameHandle } from "@platforma-sdk/model";
import {
  Annotation,
  BlockModelV3,
  createPFrameForGraphs,
  createPlDataTableV2,
} from "@platforma-sdk/model";
import { blockDataModel } from "./dataModel";
import type { BlockArgs, WorkflowInfo } from "./types";

export type * from "@milaboratories/helpers";
export { blockDataModel } from "./dataModel";
export type { BlockArgs, BlockData, WorkflowInfo, WorkflowMode, WorkflowReceptor } from "./types";

const inputAnchorSpecs = [
  // Peptide mode — universal naming
  {
    axes: [{ name: "pl7.app/sampleId" }, { name: "pl7.app/variantKey" }],
    annotations: { "pl7.app/isAnchor": "true" },
  },
  // Antibody/TCR — legacy MiXCR bulk (cloneId per spec; clonotypeKey per current MiXCR output)
  {
    axes: [{ name: "pl7.app/sampleId" }, { name: "pl7.app/vdj/cloneId" }],
    annotations: { "pl7.app/isAnchor": "true" },
  },
  {
    axes: [{ name: "pl7.app/sampleId" }, { name: "pl7.app/vdj/clonotypeKey" }],
    annotations: { "pl7.app/isAnchor": "true" },
  },
  // Antibody/TCR — legacy MiXCR single-cell
  {
    axes: [{ name: "pl7.app/sampleId" }, { name: "pl7.app/vdj/scClonotypeKey" }],
    annotations: { "pl7.app/isAnchor": "true" },
  },
];

export const platforma = BlockModelV3.create(blockDataModel)
  .args<BlockArgs>((data) => {
    if (data.inputAnchor === undefined) {
      throw new Error("Select an input dataset");
    }
    return {
      inputAnchor: data.inputAnchor,
    };
  })
  .output("inputOptions", (ctx) => ctx.resultPool.getOptions(inputAnchorSpecs))
  .output("inputSpec", (ctx) =>
    ctx.data.inputAnchor ? ctx.resultPool.getPColumnSpecByRef(ctx.data.inputAnchor) : undefined,
  )
  .output("info", (ctx) => ctx.outputs?.resolve("info")?.getDataAsJson<WorkflowInfo>())
  .output("isRunning", (ctx) => ctx.outputs?.getIsReadyOrError() === false)
  .output("processingLog", (ctx) => ctx.outputs?.resolve("processingLog")?.getLogHandle())
  .outputWithStatus("propertiesTable", (ctx) => {
    if (ctx.data.inputAnchor === undefined) return undefined;
    const ownCols = ctx.outputs?.resolve("propertiesPf")?.getPColumns();
    if (ownCols === undefined) return undefined;
    // Temporary downgrade to V2 with a fixed, propertiesPf-only column list.
    // V3 default-visibility rules and the column picker are not yet stable for
    // the multi-source layout this block needs (upstream sequence column +
    // own scalar properties). Once V3 supports it natively, rewire across
    // blocks. aaFraction is 2-axis (variantKey × aminoAcid) — already filtered
    // from the graph pFrame; filtered here too so it doesn't widen the table.
    const tableCols = ownCols.filter((c) => c.spec.axesSpec.length === 1);
    return createPlDataTableV2(ctx, tableCols, ctx.data.tableState);
  })
  .outputWithStatus("propertiesPfHandle", (ctx): PFrameHandle | undefined => {
    const allPCols = ctx.outputs?.resolve("propertiesPf")?.getPColumns();
    if (allPCols === undefined) return undefined;
    // Drop the AA fraction column from the pframe entirely. Two-axis
    // (variantKey × aminoAcid), at 50k peptides ~1M cells — enough to trip
    // graph-maker's cell-count guard on its own. The picker already excludes
    // it via `isNumericScalar` (axesSpec.length === 1), so the data was
    // pure overhead.
    const pCols = allPCols.filter((c) => c.spec.name !== "pl7.app/aaFraction");
    // Use `ctx.createPFrame` instead of `createPFrameForGraphs`. The latter
    // walks the result pool and pulls in this block's `exports.properties`
    // — a `trace.inject`-stamped re-emission of every column already in
    // `propertiesPf`, published for Lead Selection — so axis dropdowns
    // show e.g. "Net Charge (pH7) / IG" twice. Same workaround chosen by
    // cdr3-spectratype, batch-correction, cell-type-annotation, and
    // dimensionality-reduction.
    //
    // Pull single-axis metadata anchored to the input dataset's two axes
    // (idx 0 = sample, idx 1 = entity key) so sample groups / patient IDs /
    // peptide abundance and similar cols remain available for grouping and
    // filtering. Drop self-trace to keep our own exports out.
    const inputAnchor = ctx.data.inputAnchor;
    const upstreamMeta =
      inputAnchor !== undefined
        ? (
            ctx.resultPool.getAnchoredPColumns({ main: inputAnchor }, [
              { axes: [{ anchor: "main", idx: 0 }] },
              { axes: [{ anchor: "main", idx: 1 }] },
            ]) ?? []
          ).filter(
            (c) =>
              !c.spec.annotations?.[Annotation.Trace]?.includes(
                "milaboratories.sequence-properties",
              ),
          )
        : [];
    return createPFrameForGraphs(ctx, [...pCols, ...upstreamMeta]);
  })
  .output("propertiesPfCols", (ctx): PColumnIdAndSpec[] | undefined => {
    const pCols = ctx.outputs?.resolve("propertiesPf")?.getPColumns();
    if (pCols === undefined) return undefined;
    return pCols.map((c) => ({ columnId: c.id, spec: c.spec }) satisfies PColumnIdAndSpec);
  })
  .title(() => "Sequence Properties")
  .subtitle((ctx) => ctx.data.defaultBlockLabel ?? "")
  .sections(() => [
    { type: "link", href: "/", label: "Main" },
    { type: "link", href: "/scatter", label: "Property Relationships" },
    { type: "link", href: "/histogram", label: "Property Distribution" },
  ])
  .done();

export type BlockOutputs = InferOutputsType<typeof platforma>;
