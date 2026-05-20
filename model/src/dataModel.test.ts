// Tests drive the migration callbacks directly because @platforma-sdk/model's
// DataModelBuilder exposes no introspection / external-apply API. The builder
// wiring at dataModel.ts (.migrate("Ver_2026_05_05", migrateV1toV2) and
// .migrate("Ver_2026_05_18", migrateV2toV2_1)) is covered by integration use
// — if someone removes a .migrate(...) line but keeps the named callback
// exported, these unit tests will still pass. Live block load is the only
// signal that catches that class of regression today.

import { describe, expect, it } from "vitest";
import type { BlockData, BlockDataV1, BlockDataV2 } from "./types";
import { migrateV1toV2, migrateV2toV2_1 } from "./dataModel";

const tableState = {
  pTableParams: { defaultFilters: null, filters: null, hiddenColIds: null, sorting: [], sourceId: null },
  stateCache: [],
  version: 6,
} as BlockDataV2["tableState"];

const v2Graph: Pick<BlockDataV2, "graphStateScatter" | "graphStateHistogram"> = {
  graphStateScatter: { currentTab: null, template: "dots", title: "Property Relationships" },
  graphStateHistogram: { currentTab: null, layersSettings: { bins: { fillColor: "#99e099" } }, template: "bins", title: "Property Distribution" },
};

describe("blockDataModel Ver_2026_05_18 backfill", () => {
  it("backfills both label fields to '' on a bare V2 payload", () => {
    const v2: BlockDataV2 = { tableState, ...v2Graph };
    const upgraded: BlockData = migrateV2toV2_1(v2);
    expect(upgraded.customBlockLabel).toBe("");
    expect(upgraded.defaultBlockLabel).toBe("");
  });

  it("preserves an interim-deployed customBlockLabel", () => {
    const v2: BlockDataV2 = { tableState, ...v2Graph, customBlockLabel: "X" };
    const upgraded: BlockData = migrateV2toV2_1(v2);
    expect(upgraded.customBlockLabel).toBe("X");
    expect(upgraded.defaultBlockLabel).toBe("");
  });

  it("preserves an interim-deployed defaultBlockLabel", () => {
    const v2: BlockDataV2 = { tableState, ...v2Graph, defaultBlockLabel: "Y" };
    const upgraded: BlockData = migrateV2toV2_1(v2);
    expect(upgraded.defaultBlockLabel).toBe("Y");
    expect(upgraded.customBlockLabel).toBe("");
  });

  it("runs the full V1 → V2 → V2.1 chain on legacy data", () => {
    const v1: BlockDataV1 = { tableState, defaultBlockLabel: "Old" };
    const upgraded: BlockData = migrateV2toV2_1(migrateV1toV2(v1));
    expect(upgraded.defaultBlockLabel).toBe("Old");
    expect(upgraded.customBlockLabel).toBe("");
    expect(upgraded.graphStateScatter).toBeDefined();
    expect(upgraded.graphStateHistogram).toBeDefined();
  });
});
