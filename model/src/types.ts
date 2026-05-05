import type { GraphMakerState } from "@milaboratories/graph-maker";
import type { PlDataTableStateV2, PlRef } from "@platforma-sdk/model";

// Pre-Visualizations shape. Retained as the v1 type so the data-model
// migration chain stays well-typed.
export type BlockDataV1 = {
  inputAnchor?: PlRef;
  tableState: PlDataTableStateV2;
  // UI-only state. Tracks the selected input dataset's label so the block
  // subtitle can reflect it — populated by the UI watcher in app.ts. Not
  // projected into BlockArgs because the workflow does not consume it.
  defaultBlockLabel?: string;
};

export type BlockData = BlockDataV1 & {
  graphStateScatter: GraphMakerState;
  graphStateHistogram: GraphMakerState;
};

export type BlockArgs = {
  inputAnchor: PlRef;
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
