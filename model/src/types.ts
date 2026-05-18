import type { GraphMakerState } from "@milaboratories/graph-maker";
import type { PlDataTableStateV2, PlRef } from "@platforma-sdk/model";

// Pre-Visualizations shape. Retained as the v1 type so the data-model
// migration chain stays well-typed.
export type BlockDataV1 = {
  inputAnchor?: PlRef;
  tableState: PlDataTableStateV2;
  // Historically optional UI-only state. Required in the current BlockData
  // shape (V2) and projected into BlockArgs so the workflow trace label can
  // fall back to it when the user has not typed a custom subtitle.
  defaultBlockLabel?: string;
};

export type BlockData = Omit<BlockDataV1, "defaultBlockLabel"> & {
  defaultBlockLabel: string;
  customBlockLabel: string;
  graphStateScatter: GraphMakerState;
  graphStateHistogram: GraphMakerState;
};

export type BlockArgs = {
  inputAnchor: PlRef;
  defaultBlockLabel: string;
  customBlockLabel: string;
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
