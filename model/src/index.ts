import type { InferOutputsType } from "@platforma-sdk/model";
import { BlockModelV3, createPlDataTableV3 } from "@platforma-sdk/model";
import { blockDataModel } from "./dataModel";
import type { BlockArgs, WorkflowInfo } from "./types";

export type { BlockArgs, BlockData, WorkflowInfo, WorkflowMode, WorkflowReceptor } from "./types";
export { blockDataModel } from "./dataModel";

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
  .output("inputOptions", (ctx) =>
    ctx.resultPool.getOptions(inputAnchorSpecs, { refsWithEnrichments: true }),
  )
  .output("inputSpec", (ctx) =>
    ctx.data.inputAnchor ? ctx.resultPool.getPColumnSpecByRef(ctx.data.inputAnchor) : undefined,
  )
  .output("info", (ctx) => ctx.outputs?.resolve("info")?.getDataAsJson<WorkflowInfo>())
  .output("isRunning", (ctx) => ctx.outputs?.getIsReadyOrError() === false)
  .output("processingLog", (ctx) => ctx.outputs?.resolve("processingLog")?.getLogHandle())
  .outputWithStatus("propertiesTable", (ctx) => {
    if (ctx.data.inputAnchor === undefined) return undefined;
    const propertiesPf = ctx.outputs?.resolve("propertiesPf");
    if (propertiesPf === undefined) return undefined;
    return createPlDataTableV3(ctx, {
      tableState: ctx.data.tableState,
      columns: {
        anchors: { main: ctx.data.inputAnchor },
        selector: { mode: "enrichment", maxHops: 0 },
      },
    });
  })
  .title((ctx) => ctx.data.defaultBlockLabel || "Sequence Properties")
  .sections(() => [{ type: "link" as const, href: "/" as const, label: "Main" }])
  .done();

export type BlockOutputs = InferOutputsType<typeof platforma>;
