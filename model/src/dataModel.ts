import type { GraphMakerState } from "@milaboratories/graph-maker";
import { kind } from "@platforma-open/milaboratories.sequence-properties.kind";
import { createPlDataTableStateV2, DataModelBuilder } from "@platforma-sdk/model";
import type { BlockData, BlockDataV1, BlockDataV2, BlockDataV2_1 } from "./types";

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

export const migrateV1toV2 = (v1: BlockDataV1): BlockDataV2 => ({
  ...v1,
  graphStateScatter: { ...DEFAULT_SCATTER_STATE },
  graphStateHistogram: { ...DEFAULT_HISTOGRAM_STATE },
});

export const migrateV2toV2_1 = (v2: BlockDataV2): BlockDataV2_1 => ({
  ...v2,
  defaultBlockLabel: v2.defaultBlockLabel ?? "",
  customBlockLabel: v2.customBlockLabel ?? "",
});

// The old label was exactly the chosen input's (or subset's) option label, so
// it carries over as the subtitle text unchanged.
export const migrateV2_1toV3 = ({ defaultBlockLabel, ...rest }: BlockDataV2_1): BlockData => ({
  ...rest,
  subtitleLabel: defaultBlockLabel,
});

export const blockDataModel = new DataModelBuilder({ kind })
  .from<BlockDataV1>("Ver_2026_04_28")
  // Already-deployed step. Future field additions must go into a new step
  // below — editing this body has no effect on projects already tagged
  // Ver_2026_05_05 (DataModelBuilder skips matching-version migrations).
  .migrate<BlockDataV2>("Ver_2026_05_05", migrateV1toV2)
  // Backfills label fields onto V2-tagged projects. `?? ""` preserves any
  // interim-deployed value; missing fields default to "".
  .migrate<BlockDataV2_1>("Ver_2026_05_18", migrateV2toV2_1)
  .migrate<BlockData>("Ver_2026_11_06", migrateV2_1toV3)
  // `params` carries the kind's init-params contract, and is undefined when the
  // block is created outside a template. `inputAnchor` and `filterRef` are the
  // only fields a creator supplies; everything else is view state that always defaults.
  .init(({ params }) => ({
    inputAnchor: params?.inputAnchor,
    filterRef: params?.filterRef,
    tableState: createPlDataTableStateV2(),
    customBlockLabel: "",
    subtitleLabel: "",
    graphStateScatter: { ...DEFAULT_SCATTER_STATE },
    graphStateHistogram: { ...DEFAULT_HISTOGRAM_STATE },
  }));
