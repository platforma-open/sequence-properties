import type { GraphMakerState } from "@milaboratories/graph-maker";
import { createPlDataTableStateV2, DataModelBuilder } from "@platforma-sdk/model";
import type { BlockData, BlockDataV1 } from "./types";

const DEFAULT_SCATTER_STATE: GraphMakerState = {
  title: "Properties Scatter Plot",
  template: "dots",
  currentTab: null,
};

const DEFAULT_HISTOGRAM_STATE: GraphMakerState = {
  title: "Property Histogram",
  template: "bins",
  currentTab: null,
  layersSettings: {
    bins: { fillColor: "#99e099" },
  },
};

export const blockDataModel = new DataModelBuilder()
  .from<BlockDataV1>("Ver_2026_04_28")
  .migrate<BlockData>("Ver_2026_05_05", (v1) => ({
    ...v1,
    graphStateScatter: { ...DEFAULT_SCATTER_STATE },
    graphStateHistogram: { ...DEFAULT_HISTOGRAM_STATE },
  }))
  .init(() => ({
    tableState: createPlDataTableStateV2(),
    graphStateScatter: { ...DEFAULT_SCATTER_STATE },
    graphStateHistogram: { ...DEFAULT_HISTOGRAM_STATE },
  }));
