import { kind } from "@platforma-open/milaboratories.sequence-properties.kind";
import type {
  DatasetOption,
  InferOutputsType,
  PColumnIdAndSpec,
  PFrameHandle,
} from "@platforma-sdk/model";
import {
  Annotation,
  BlockModelV3,
  buildDatasetOptions,
  createGlobalPObjectId,
  createPlDataTableV2,
  getRelatedColumns,
  isHiddenFromGraphColumn,
  isHiddenFromUIColumn,
  isPColumnSpec,
} from "@platforma-sdk/model";
import { blockDataModel } from "./dataModel";
import { resolveSubtitle, resolveTraceLabel } from "./label";
import type { BlockArgs, WorkflowInfo } from "./types";

export type * from "@milaboratories/helpers";
export { blockDataModel } from "./dataModel";
export type { BlockArgs, BlockData, WorkflowInfo, WorkflowMode, WorkflowReceptor } from "./types";

// Key axes of the datasets this block accepts: peptides and variants, legacy bulk clones, bulk
// and single-cell clonotypes.
const DATASET_KEY_AXES: ReadonlySet<string> = new Set([
  "pl7.app/variantKey",
  "pl7.app/vdj/cloneId",
  "pl7.app/vdj/clonotypeKey",
  "pl7.app/vdj/scClonotypeKey",
]);

export const platforma = BlockModelV3.create({ dataModel: blockDataModel, kind })
  .args<BlockArgs>((data) => {
    if (data.inputAnchor === undefined) {
      throw new Error("Select an input dataset");
    }
    return {
      inputAnchor: data.inputAnchor,
      traceLabel: resolveTraceLabel(data),
      // The subset's canonical column id, which the workflow stamps as the outputs'
      // `pl7.app/inputSubset`. Absent without a filter, so unfiltered args are unchanged.
      ...(data.filterRef !== undefined && {
        inputFilter: createGlobalPObjectId(data.filterRef.blockId, data.filterRef.name),
      }),
    };
  })
  // Inverse of the kind's init-params contract: the same fields `init`
  // consumes. The label fields and the table / plot states stay out — they are
  // view state, not configuration a template carries.
  .templateParams((data) => ({ inputAnchor: data.inputAnchor, filterRef: data.filterRef }))
  // Datasets, each with its subset columns (e.g. Repertoire Labeling labels) as filters.
  .output(
    "inputOptions",
    (ctx): DatasetOption[] =>
      buildDatasetOptions(ctx, {
        primary: (spec) =>
          isPColumnSpec(spec) &&
          spec.annotations?.["pl7.app/isAnchor"] === "true" &&
          spec.axesSpec.length === 2 &&
          spec.axesSpec[0].name === "pl7.app/sampleId" &&
          DATASET_KEY_AXES.has(spec.axesSpec[1].name),
        // Only subsets keyed by the entity axis alone: properties are per entity.
        filter: (spec) =>
          isPColumnSpec(spec) &&
          spec.axesSpec.length === 1 &&
          DATASET_KEY_AXES.has(spec.axesSpec[0].name),
        requireEnrichments: false,
      }) ?? [],
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
    // Same pool walk as `createPFrameForGraphs` (upstream and linker columns such as cluster size),
    // minus this block's `exports.properties`, which repeat every column in `propertiesPf`.
    const related = getRelatedColumns(ctx, {
      columns: pCols,
      predicate: (spec) =>
        !isHiddenFromUIColumn(spec) &&
        !isHiddenFromGraphColumn(spec) &&
        !spec.annotations?.[Annotation.Trace]?.includes("milaboratories.sequence-properties"),
    });
    const ownIds = new Set(pCols.map((c) => c.id));
    return ctx.createPFrame([...pCols, ...related.filter((c) => !ownIds.has(c.id))]);
  })
  .output("propertiesPfCols", (ctx): PColumnIdAndSpec[] | undefined => {
    const pCols = ctx.outputs?.resolve("propertiesPf")?.getPColumns();
    if (pCols === undefined) return undefined;
    return pCols.map((c) => ({ columnId: c.id, spec: c.spec }) satisfies PColumnIdAndSpec);
  })
  .title(() => "Sequence Properties")
  .subtitle((ctx) => resolveSubtitle(ctx.data))
  .sections(() => [
    { type: "link", href: "/", label: "Main" },
    { type: "link", href: "/scatter", label: "Property Relationships" },
    { type: "link", href: "/histogram", label: "Property Distribution" },
  ])
  .done();

export type BlockOutputs = InferOutputsType<typeof platforma>;
