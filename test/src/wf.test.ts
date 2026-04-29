/**
 * Block integration tests for sequence-properties.
 *
 * Tests live here (per harness `testing.md`): real platforma backend, real
 * project, the block under test added with `blockTest` from
 * `@platforma-sdk/test`. Each test:
 *
 *   1. Creates / opens a project.
 *   2. Adds an upstream block that publishes an abundance anchor PColumn
 *      shaped like one of our four supported axis patterns
 *      (peptide variantKey, vdj cloneId, vdj clonotypeKey, vdj scClonotypeKey).
 *   3. Adds the sequence-properties block, sets `inputAnchor` to the upstream
 *      ref, runs the workflow, awaits stable outputs.
 *   4. Asserts on the *observable* contract — emitted PColumn names, axis
 *      passthrough, info messages, score annotations — never on internal
 *      Tengo or Python machinery.
 *
 * No mocks: this is the integration layer. Unit-level Python behavior is
 * already covered by `software/tests/`. These tests guard the workflow +
 * model + UI wiring, plus the spec-mandated cross-cutting behaviors that
 * unit tests can't reach.
 *
 * Running: `cd blocks/sequence-properties && pnpm --filter
 * @platforma-open/MiLaboratories.sequence-properties.test test`
 *
 * A local platforma backend must be running — see the platforma-dev skill's
 * `local-server.md`.
 */

import { describe, it } from "vitest";

// ---------------------------------------------------------------------------
// Modality detection (R1, R1a, R1b)
// ---------------------------------------------------------------------------

describe("modality detection", () => {
  // Peptide variantKey + extractionRunId → peptide mode. Outputs expose
  // peptide property columns (charge_peptide, gravy_peptide, ...) and the
  // 2-axis AA fraction PColumn. info.mode == "peptide".
  it.todo("detects peptide mode from variantKey + extractionRunId domain");

  // Legacy MiXCR bulk: vdj/cloneId axis → antibody_tcr_legacy_bulk.
  it.todo("detects antibody/TCR legacy bulk from vdj/cloneId axis");

  // Legacy MiXCR alias seen in real MiXCR output: vdj/clonotypeKey axis
  // also resolves to antibody_tcr_legacy_bulk (current platform emits this
  // alias even though spec says cloneId).
  it.todo("detects antibody/TCR legacy bulk from vdj/clonotypeKey axis");

  // Single-cell: vdj/scClonotypeKey → antibody_tcr_legacy_sc.
  it.todo("detects antibody/TCR single-cell from vdj/scClonotypeKey axis");

  // No recognised axis → workflow panics with the spec-mandated message:
  // "no recognized sequence key axis found; connect a peptide extraction or
  // MiXCR clonotyping dataset."
  it.todo("panics with R1a message when input has no recognised sequence key");

  // R1b: every output column's axis spec must be byte-identical to the
  // input's entity axis (name + domain + type + annotations). Renaming or
  // re-domaining the axis breaks downstream Lead Selection joins.
  it.todo("emits output PColumns with input entity axis verbatim");
});

// ---------------------------------------------------------------------------
// Peptide mode (R6, R7, R8, R9, R16)
// ---------------------------------------------------------------------------

describe("peptide mode output", () => {
  // R6: nine scalar PColumn names with the universal `pl7.app/*` value
  // names (pl7.app/charge, pl7.app/hydrophobicity, etc.) and the
  // distinguishing domain `{pl7.app/feature: "peptide"}`. No
  // pl7.app/sequenceLength column emitted (that comes from
  // peptide-extraction upstream).
  it.todo("emits the 9 peptide scalar PColumns with feature=peptide domain");

  // R7: the AA fraction PColumn is 2-axis — primary axis matches the
  // peptide variantKey input, second axis name == "pl7.app/aminoAcid".
  // Each variantKey has exactly 20 rows (one per std AA).
  it.todo("emits 2-axis pl7.app/aaFraction PColumn with aminoAcid axis");

  // R16: AA fraction must NOT carry pl7.app/isScore — display only.
  it.todo("AA fraction PColumn has no isScore annotation");

  // R8: ε oxidised carries floor(C/2)·125 disulfide bonus; ε reduced
  // omits Cys term. Both columns have valueType Double, min=0.
  it.todo("extinction coefficients oxidised and reduced both emitted");

  // R9: peptides shorter than 10 aa receive NA for instability_peptide,
  // not 0 and not a column-level warning that would taint mixed-length
  // datasets.
  it.todo("instability NA for peptides under 10 aa (per-row, not column-wide)");
});

// ---------------------------------------------------------------------------
// Antibody / TCR mode (R10, R11, R12, R13, R14)
// ---------------------------------------------------------------------------

describe("antibody/TCR coverage handling", () => {
  // R10: CDR3 charge + hydrophobicity emitted for each chain present in
  // the input; CDR3 sequence length NOT re-emitted (MiXCR upstream).
  it.todo("emits CDR3 charge + hydrophobicity per chain when CDR3 present");

  // R11 + R12: when both chains carry full 7 regions and receptor is IG,
  // full-chain VH + VL columns appear AND the 5 Fv columns appear.
  it.todo("emits full-chain + Fv columns on full IG paired coverage");

  // R11a: CDR3-only input (no FR/CDR1/CDR2 columns in the result pool)
  // — full-chain columns absent from output schema entirely; info message
  // exact wording: "CDR3-only input detected — full-chain properties not
  // computed. To enable them, use a MiXCR preset that exports all VDJ
  // regions."
  it.todo("CDR3-only input emits R11a info message and no full-chain columns");

  // R11b: 1-6 of 7 regions present for a chain — partial-region message
  // names the chain by anatomical label ("heavy"/"light" for IG,
  // "alpha"/"beta" for TCRAB, "gamma"/"delta" for TCRGD), not by chain id.
  it.todo("partial-region info message uses receptor-aware chain label");

  // Per-clone partial coverage: full-chain columns are emitted at schema
  // level (chain has all 7 region columns dataset-wide) but NA for the
  // specific clone whose region was empty. Other clones unaffected.
  it.todo("missing region for one clone yields NA for that row only");

  // R12: TCR inputs (TCRAB / TCRGD) MUST NOT emit Fv columns even when
  // both chains have full coverage. The hasFv flag in the workflow plan
  // gates this.
  it.todo("TCR receptor never emits Fv columns regardless of coverage");

  // γδ TCR: emits an info message — "γδ TCR input detected — displaying
  // with γδ-specific labels; Fv columns are not computed for TCR inputs."
  it.todo("γδ TCR input emits the γδ-labels info message");

  // R13: chain "A" → heavy (IG) / alpha (TCRAB) / gamma (TCRGD) labels in
  // the column `pl7.app/label` annotation. Column name and chain domain
  // value ("A"/"B") stay constant — only the label varies.
  it.todo("chain labels adapt to receptor while name + chain domain stay constant");

  // R13b: when no `pl7.app/vdj/receptor` annotation is present on collected
  // sequence columns, block defaults to IG and emits a warning info
  // message.
  it.todo("missing receptor annotation defaults to IG with warning info message");

  // R14: scalar score columns carry isScore=true exactly on the spec-
  // listed set (CDR3 charge+gravy, peptide charge+gravy, full-chain charge
  // +pi, Fv charge+pi). All non-score numeric columns omit isScore.
  it.todo("isScore annotation on the spec-mandated columns only");

  // R14: hydrophobicity columns (peptide, CDR3) carry rankingOrder
  // "increasing"; charge / pi columns leave rankingOrder unset for user
  // choice in Lead Selection.
  it.todo("rankingOrder=increasing only on hydrophobicity columns");

  // R15: no column anywhere carries pl7.app/score/defaultCutoff.
  it.todo("no defaultCutoff annotation anywhere in the output");
});

// ---------------------------------------------------------------------------
// Edge cases per spec "Defaults and edge cases" table
// ---------------------------------------------------------------------------

describe("edge cases", () => {
  // Stop codon `*` anywhere in a region invalidates the whole sequence —
  // every property NA for that entity (including reconstructed full chain
  // when stop is in any region).
  it.todo("stop codon in any region NAs the whole reconstructed chain");

  // Non-standard residues (X, B, Z, U, J, gap '-') drop from numerator
  // AND denominator — fractions still sum to 1.0 across present residues.
  it.todo("non-standard residues do not affect fractions / pI / GRAVY denominator");

  // Polybasic synthetic with no zero crossing in [0, 14] — pI emits NA,
  // not an exception, not 14.0 clamped.
  it.todo("polybasic with no zero crossing emits NA pI without error");

  // Sequence with no Tyr or Trp → ε = 0 (defined zero, not NA).
  it.todo("no aromatic residues emits ε = 0 (not NA)");

  // Single-cell single chain (chain dropout) — per-chain columns emitted
  // for the present chain; Fv column NA for that clonotype.
  it.todo("single-cell single chain emits available chain, Fv NA");
});

// ---------------------------------------------------------------------------
// Result-pool integration & enrichment (Lead Selection downstream)
// ---------------------------------------------------------------------------

describe("downstream consumption", () => {
  // Lead Selection picker queries the result pool for `isScore: "true"`
  // numeric columns. Verify our emitted score columns are visible to
  // a downstream block via getOptions / getAnchoredPColumns.
  it.todo("score columns are discoverable by a downstream block via the result pool");

  // Trace annotation: every output PColumn's `pl7.app/trace` includes our
  // block's id + label, prepended to the upstream's trace chain.
  it.todo("output PColumns carry pl7.app/trace stamped with this block's label");

  // The export PFrame (the `exports.properties` key returned by the
  // workflow body) carries the same score columns with `pl7.app/blockId`
  // stamped into the domain — disambiguates outputs from multiple
  // sequence-properties block instances in the same project.
  it.todo("export PFrame stamps blockId on score column domains");
});

// ---------------------------------------------------------------------------
// Model / UI behaviors
// ---------------------------------------------------------------------------

describe("model + UI", () => {
  // Args validation: with no inputAnchor selected, the run button is
  // disabled — model `.args(data => ...)` lambda throws and the runtime
  // catches it.
  it.todo("run button disabled when inputAnchor not set");

  // inputOptions output advertises the 4 supported abundance-anchor
  // patterns. Adding an upstream that emits one of them makes that ref
  // appear in the dropdown; an unrelated abundance anchor (e.g. with a
  // sampleId+geneSymbol axis pair) does not.
  it.todo("inputOptions filters to the 4 supported abundance anchor patterns");

  // info output is consumed by the UI's PlAlert list — confirm it's a
  // structured object with mode / receptor / coverageTier / messages keys
  // and that messages is a non-null array.
  it.todo("info output exposes mode, receptor, coverageTier, messages array");

  // Re-running with a different inputAnchor triggers a fresh workflow —
  // outputs from the prior run do not leak into the new run's PFrame.
  it.todo("changing inputAnchor invalidates and re-runs cleanly");

  // Title is the static "Sequence Properties" string (per current
  // implementation). If the operator later adds dynamic title from
  // detected modality (e.g. "Sequence Properties — peptide"), update this
  // test.
  it.todo("title is the configured static label");
});

// ---------------------------------------------------------------------------
// Dedup wiring (see docs/dedup-plan.md)
// ---------------------------------------------------------------------------

describe("dedup", () => {
  // Two projects pointing at the same upstream peptide-extraction (or same
  // MiXCR clonotyping) — the second project's sequence-properties block must
  // land on Done immediately without spawning a fresh Python step. Verified
  // by inspecting that the resource handle returned for `exports.properties`
  // points at the same resource-pool entry across both projects.
  //
  // Python-side byte-stability of properties.tsv / aa_fraction.tsv is already
  // covered by `software/tests/integration/test_cli.py` (sha256 across two
  // runs). This test guards the Tengo-side wiring: sorted map iteration in
  // main.tpl.tengo + process.tpl.tengo and canonical JSON resource encoding
  // for plan.json / params / infoBlob.
  it.todo("second project on identical upstream lands on Done via dedup");

  // Negative: change a single byte in the upstream input — the second
  // project's block must run a fresh Python step. Proves the CID is
  // genuinely distinguishing inputs and not always cache-hitting.
  it.todo("changed upstream input breaks dedup and triggers fresh run");
});

// ---------------------------------------------------------------------------
// Cross-block composition
// ---------------------------------------------------------------------------

describe("cross-block composition", () => {
  // Realistic chain: peptide-extraction → sequence-properties →
  // lead-selection. Each block in turn sees the upstream's outputs via
  // result-pool queries; no manual ref wiring beyond the inputAnchor on
  // sequence-properties.
  it.todo("end-to-end: peptide-extraction → sequence-properties → lead-selection");

  // Same block instance twice in one project (e.g. same dataset analysed
  // with different parameter sets in the future): blockId domain stamping
  // keeps the two output sets distinguishable in the result pool.
  it.todo("two sequence-properties instances in one project disambiguate via blockId");

  // VDJ chain: mixcr-clonotyping (full preset, all 7 regions) →
  // sequence-properties → antibody-tcr-lead-selection. Verify the Fv
  // score columns appear in lead-selection's ranking criteria UI.
  it.todo("VDJ chain: full-coverage MiXCR → seq-properties → lead-selection");
});
