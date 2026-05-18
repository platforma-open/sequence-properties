import type { GraphMakerState } from "@milaboratories/graph-maker";
import { createPlDataTableStateV2, DataModelBuilder } from "@platforma-sdk/model";
import type { BlockData, BlockDataV1, BlockDataV2 } from "./types";

const DEFAULT_SCATTER_STATE: GraphMakerState = {
  title: "Property Relationships",
  template: "dots",
  currentTab: null,
};

const DEFAULT_HISTOGRAM_STATE: GraphMakerState = {
  title: "Property Distribution",
  template: "bins",
  currentTab: null,
  layersSettings: {
    bins: { fillColor: "#99e099" },
  },
};

export const blockDataModel = new DataModelBuilder()
  .from<BlockDataV1>("Ver_2026_04_28")
  // Already-deployed V2 step. Body must stay byte-equivalent to what shipped
  // — projects already tagged Ver_2026_05_05 skip this step and rely on the
  // next migration below to backfill the label fields.
  .migrate<BlockDataV2>("Ver_2026_05_05", (v1) => ({
    ...v1,
    graphStateScatter: { ...DEFAULT_SCATTER_STATE },
    graphStateHistogram: { ...DEFAULT_HISTOGRAM_STATE },
  }))
  // Backfills label fields onto projects already tagged at the deployed V2
  // version. Both reads preserve any pre-existing value an interim
  // deployment may have written; new V2 projects coalesce undefined to "".
  // Keeps the workflow's `traceLabel := args.customBlockLabel ||
  // args.defaultBlockLabel` wiring receiving strings, never undefined.
  .migrate<BlockData>("Ver_2026_05_18", (v2) => ({
    ...v2,
    defaultBlockLabel: v2.defaultBlockLabel ?? "",
    customBlockLabel: v2.customBlockLabel ?? "",
  }))
  .init(() => ({
    tableState: createPlDataTableStateV2(),
    defaultBlockLabel: "",
    customBlockLabel: "",
    graphStateScatter: { ...DEFAULT_SCATTER_STATE },
    graphStateHistogram: { ...DEFAULT_HISTOGRAM_STATE },
  }));
