import type { PlDataTableStateV2, PlRef } from "@platforma-sdk/model";

export type BlockData = {
  inputAnchor?: PlRef;
  tableState: PlDataTableStateV2;
  defaultBlockLabel?: string;
};

export type BlockArgs = {
  inputAnchor: PlRef;
  defaultBlockLabel: string;
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
