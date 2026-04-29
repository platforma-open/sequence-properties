import { createPlDataTableStateV2, DataModelBuilder } from "@platforma-sdk/model";
import type { BlockData } from "./types";

export const blockDataModel = new DataModelBuilder().from<BlockData>("Ver_2026_04_28").init(() => ({
  tableState: createPlDataTableStateV2(),
}));
