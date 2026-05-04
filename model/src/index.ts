import type { ColumnSource, InferOutputsType } from "@platforma-sdk/model";
import {
  Annotation,
  ArrayColumnProvider,
  BlockModelV3,
  createPlDataTableV3,
} from "@platforma-sdk/model";
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
    const ownCols = ctx.outputs?.resolve("propertiesPf")?.getPColumns();
    if (ownCols === undefined) return undefined;

    // Build sources explicitly: upstream cols from the result pool minus
    // anything traced back to this block, plus this block's own cols from
    // `propertiesPf`. The workflow also publishes `exports.properties` —
    // a blockId-stamped score-only variant for downstream consumers like
    // Lead Selection — into the result pool. Filtering by trace excludes
    // it here so score cols don't duplicate the propertiesPf variant.
    const upstreamCols = ctx.resultPool.selectColumns(
      (spec) =>
        !spec.annotations?.[Annotation.Trace]?.includes("milaboratories.sequence-properties"),
    );
    const sources: ColumnSource[] = [
      new ArrayColumnProvider(upstreamCols),
      new ArrayColumnProvider(ownCols),
    ];

    return createPlDataTableV3(ctx, {
      tableState: ctx.data.tableState,
      columns: {
        sources,
        anchors: { main: ctx.data.inputAnchor },
        selector: { mode: "enrichment" },
      },
      // UX policy (not in spec): default-show only this block's columns;
      // demote upstream cols to optional to keep the table uncluttered.
      // This block's cols fall through unmatched and keep the visibility
      // stamped at workflow-time (`pl7.app/table/visibility`).
      displayOptions: {
        visibility: [
          {
            match: (spec) =>
              !spec.annotations?.[Annotation.Trace]?.includes(
                "milaboratories.sequence-properties",
              ) && spec.annotations?.["pl7.app/isLinkerColumn"] !== "true",
            visibility: "optional",
          },
        ],
      },
    });
  })
  .title(() => "Sequence Properties")
  .subtitle((ctx) => ctx.data.defaultBlockLabel ?? "")
  .sections(() => [{ type: "link" as const, href: "/" as const, label: "Main" }])
  .done();

export type BlockOutputs = InferOutputsType<typeof platforma>;
