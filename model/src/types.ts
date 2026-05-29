// canary: coverage-gap red-check test (throwaway, do not merge)
import type { GraphMakerState } from "@milaboratories/graph-maker";
import type { PlDataTableStateV2, PlRef } from "@platforma-sdk/model";

// Pre-Visualizations shape. Retained as the v1 type so the data-model
// migration chain stays well-typed.
export type BlockDataV1 = {
  inputAnchor?: PlRef;
  tableState: PlDataTableStateV2;
  // Historically optional UI-only state. Required in the current BlockData
  // shape; consumed by the label helpers in label.ts (resolveSubtitle for
  // the PlBlockPage subtitle, resolveTraceLabel for the workflow trace).
  defaultBlockLabel?: string;
};

// V2 shape — what the deployed Ver_2026_05_05 migration produces. Input to
// the new Ver_2026_05_18 step that backfills the label fields. Both label
// fields are optional here so the V2→V2.1 migration can read-or-default any
// value an interim deployment may have written.
export type BlockDataV2 = Omit<BlockDataV1, "defaultBlockLabel"> & {
  defaultBlockLabel?: string;
  customBlockLabel?: string;
  graphStateScatter: GraphMakerState;
  graphStateHistogram: GraphMakerState;
};

export type BlockData = Omit<BlockDataV2, "defaultBlockLabel"> & {
  defaultBlockLabel: string;
  customBlockLabel: string;
};

export type BlockArgs = {
  inputAnchor: PlRef;
  traceLabel: string;
};

export type WorkflowMode =
  | "peptide"
  | "antibody_tcr_universal"
  | "antibody_tcr_legacy_bulk"
  | "antibody_tcr_legacy_sc";

export type WorkflowReceptor = "IG" | "TCRAB" | "TCRGD";

export type WorkflowInfo = {
  mode?: WorkflowMode;
  receptor?: WorkflowReceptor;
  coverageTier?: "full_chain" | "cdr3_only" | "partial" | "peptide";
  messages: string[];
};
