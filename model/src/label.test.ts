import { describe, expect, it } from "vitest";
import type { BlockData } from "./types";
import { resolveSubtitle, resolveTraceLabel } from "./label";

const base: Omit<BlockData, "customBlockLabel" | "defaultBlockLabel"> = {
  tableState: {
    pTableParams: {
      defaultFilters: null,
      filters: null,
      hiddenColIds: null,
      sorting: [],
      sourceId: null,
    },
    stateCache: [],
    version: 6,
  } as BlockData["tableState"],
  graphStateScatter: {
    currentTab: null,
    template: "dots",
    title: "Property Relationships",
  } as BlockData["graphStateScatter"],
  graphStateHistogram: {
    currentTab: null,
    layersSettings: { bins: { fillColor: "#99e099" } },
    template: "bins",
    title: "Property Distribution",
  } as BlockData["graphStateHistogram"],
};

const make = (custom: string, def: string): BlockData => ({
  ...(base as BlockData),
  customBlockLabel: custom,
  defaultBlockLabel: def,
});

describe("resolveSubtitle", () => {
  it("uses customBlockLabel when set", () => {
    expect(resolveSubtitle(make("My label", "Dataset"))).toBe("My label");
  });

  it("falls back to defaultBlockLabel when customBlockLabel is empty", () => {
    expect(resolveSubtitle(make("", "Dataset"))).toBe("Dataset");
  });

  it("returns empty string when both label fields are empty", () => {
    expect(resolveSubtitle(make("", ""))).toBe("");
  });
});

describe("resolveTraceLabel", () => {
  it("uses customBlockLabel when set", () => {
    expect(resolveTraceLabel(make("My label", "Dataset"))).toBe("My label");
  });

  it("falls back to defaultBlockLabel when customBlockLabel is empty", () => {
    expect(resolveTraceLabel(make("", "Dataset"))).toBe("Dataset");
  });

  it("falls back to the block-type default when both label fields are empty", () => {
    expect(resolveTraceLabel(make("", ""))).toBe("Sequence Properties");
  });
});
