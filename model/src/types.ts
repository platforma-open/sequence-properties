import type { GraphMakerState } from "@milaboratories/graph-maker";
import type { PlDataTableStateV2, PlRef } from "@platforma-sdk/model";

// Pre-Visualizations shape. Retained as the v1 type so the data-model
// migration chain stays well-typed.
export type BlockDataV1 = {
  inputAnchor?: PlRef;
  tableState: PlDataTableStateV2;
  // Historically optional UI-only state, written by the UI from the chosen
  // input's option label. Carried into `subtitleLabel` by Ver_2026_11_06.
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

// V2.1 shape — what the Ver_2026_05_18 migration produces. Input to the
// Ver_2026_11_06 step that moves the input label into `subtitleLabel`.
export type BlockDataV2_1 = Omit<BlockDataV2, "defaultBlockLabel"> & {
  defaultBlockLabel: string;
  customBlockLabel: string;
  /** Subset column picked with the input in PlDatasetSelector (e.g. a Repertoire Labeling label). */
  filterRef?: PlRef;
};

export type BlockData = Omit<BlockDataV2_1, "defaultBlockLabel"> & {
  /** Label of the chosen input (or subset), shown in the subtitle; written when the user picks the input. */
  subtitleLabel: string;
};

export type BlockArgs = {
  inputAnchor: PlRef;
  traceLabel: string;
  inputFilter?: string;
};

export type WorkflowMode =
  | "peptide"
  // synthetic-repertoire-profiler amplicon variants: whole-sequence properties
  // on the variant's amino-acid sequence (same math as peptide mode, labeled
  // with the "amplicon-sequence" feature instead of "peptide").
  | "amplicon"
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
