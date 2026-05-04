import type { PlDataTableStateV2, PlRef } from "@platforma-sdk/model";

export type BlockData = {
  inputAnchor?: PlRef;
  tableState: PlDataTableStateV2;
  // UI-only state. Tracks the selected input dataset's label so the block
  // subtitle can reflect it — populated by the UI watcher in app.ts. Not
  // projected into BlockArgs because the workflow does not consume it.
  defaultBlockLabel?: string;
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
