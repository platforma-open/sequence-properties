import type { BlockData } from "./types";

const BLOCK_NAME = "Sequence Properties";

// Subtitle resolution for the PlBlockPage. Empty string is a valid result
// (the page renders the title alone).
export function resolveSubtitle(data: BlockData): string {
  return data.customBlockLabel || data.subtitleLabel;
}

// Trace-label resolution for the workflow's pl7.app/trace.label. A pure
// function of block data, so a block that was never opened gets the same
// label. The block has no settings that tell two instances apart; the input
// dataset (or subset) is already part of the inherited upstream trace.
export function resolveTraceLabel(data: BlockData): string {
  return data.customBlockLabel || BLOCK_NAME;
}
