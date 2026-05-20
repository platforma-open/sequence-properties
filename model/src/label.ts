import type { BlockData } from "./types";

const STATIC_FALLBACK = "Sequence Properties";

// Subtitle resolution for the PlBlockPage. Empty string is a valid result
// (the page renders the title alone).
export function resolveSubtitle(data: BlockData): string {
  return data.customBlockLabel || data.defaultBlockLabel;
}

// Trace-label resolution for the workflow's pl7.app/trace.label. Adds a
// static block-type last-resort over resolveSubtitle so automated pipelines
// that run before the UI watchEffect populates defaultBlockLabel never see
// an empty label downstream.
export function resolveTraceLabel(data: BlockData): string {
  return data.customBlockLabel || data.defaultBlockLabel || STATIC_FALLBACK;
}
