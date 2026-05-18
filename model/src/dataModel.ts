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
  // Already-deployed step. Future field additions must go into a new step
  // below — editing this body has no effect on projects already tagged
  // Ver_2026_05_05 (DataModelBuilder skips matching-version migrations).
  .migrate<BlockDataV2>("Ver_2026_05_05", (v1) => ({
    ...v1,
    graphStateScatter: { ...DEFAULT_SCATTER_STATE },
    graphStateHistogram: { ...DEFAULT_HISTOGRAM_STATE },
  }))
  // Backfills label fields onto V2-tagged projects. `?? ""` preserves any
  // interim-deployed value; missing fields default to "". The workflow's
  // `traceLabel := args.customBlockLabel || args.defaultBlockLabel`
  // requires both args to be strings, never undefined.
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
