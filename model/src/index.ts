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
    // `coverageTier` is set in workflow/main.tpl.tengo and surfaced via the
    // `info` JSON resource. Allowed values are defined in types.ts::WorkflowInfo.
    // Gate on `info` so the table renders consistently with the chosen aa
    // column rather than briefly without it while `info` is still resolving.
    const info = ctx.outputs?.resolve("info")?.getDataAsJson<WorkflowInfo>();
    if (info === undefined) return undefined;
    const tier = info.coverageTier;

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
      // Default-visible: this block's columns + a single source amino-acid
      // sequence column matching the analysed coverage tier. Reviewer asked
      // for one sequence next to the properties — full-chain VDJRegion when
      // available (it contains the CDR3); CDR3 alone when that is all the
      // input has; peptide for peptide mode. Chain A (heavy / alpha / gamma)
      // only — chain B stays available via the column picker. Other upstream
      // cols → optional. This block's cols fall through unmatched and keep
      // their workflow-time `pl7.app/table/visibility` annotation.
      displayOptions: {
        visibility: [
          {
            match: (spec) => {
              if (spec.domain?.["pl7.app/vdj/scClonotypeChain/index"] === "secondary") {
                return false;
              }
              if (spec.domain?.["pl7.app/alphabet"] !== "aminoacid") return false;

              const isVdj = spec.name === "pl7.app/vdj/sequence";
              const isUniversal = spec.name === "pl7.app/sequence";
              if (!isVdj && !isUniversal) return false;

              const feature = isVdj
                ? spec.domain?.["pl7.app/vdj/feature"]
                : spec.domain?.["pl7.app/feature"];

              if (tier === "peptide") {
                return isUniversal && feature === "peptide";
              }

              const chain = spec.domain?.["pl7.app/vdj/scClonotypeChain"];
              if (chain !== undefined && chain !== "A") return false;

              if (tier === "full_chain") {
                return feature === "VDJRegion" || feature === "VDJRegionInFrame";
              }
              if (tier === "cdr3_only" || tier === "partial") {
                return feature === "CDR3";
              }
              return false;
            },
            visibility: "default",
          },
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
