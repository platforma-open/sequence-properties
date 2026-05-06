import type { PColumnIdAndSpec, PColumnSpec } from "@platforma-sdk/model";
import type {
  WorkflowInfo,
  WorkflowMode,
} from "@platforma-open/milaboratories.sequence-properties.model";

const NUMERIC_VALUE_TYPES = new Set(["Int", "Long", "Float", "Double"]);

// R18a: numeric scalar pcolumns shown in axis pickers. The 2-axis check
// excludes the AA fraction column (R7); the explicit name check is
// belt-and-braces against any future scalar carrying the same name.
export function isNumericScalar(spec: PColumnSpec): boolean {
  return (
    NUMERIC_VALUE_TYPES.has(spec.valueType) &&
    spec.axesSpec.length === 1 &&
    spec.name !== "pl7.app/aaFraction"
  );
}

function isPeptideMode(mode: WorkflowMode | undefined): boolean {
  return mode === "peptide";
}

function matchesDomain(spec: PColumnSpec, requiredDomain: Record<string, string>): boolean {
  const domain = spec.domain;
  if (domain === undefined) return false;
  for (const key of Object.keys(requiredDomain)) {
    if (domain[key] !== requiredDomain[key]) return false;
  }
  return true;
}

function pickByName(
  cols: PColumnIdAndSpec[],
  name: string,
  domain: Record<string, string>,
): PColumnSpec | undefined {
  return cols.find((c) => c.spec.name === name && matchesDomain(c.spec, domain))?.spec;
}

// R19: default scatter X = charge, Y = hydrophobicity (per modality).
export function defaultScatterAxes(
  cols: PColumnIdAndSpec[],
  info: WorkflowInfo | undefined,
): { x?: PColumnSpec; y?: PColumnSpec } {
  if (info?.mode === undefined) return {};
  if (isPeptideMode(info.mode)) {
    const domain = { "pl7.app/feature": "peptide" };
    return {
      x: pickByName(cols, "pl7.app/charge", domain),
      y: pickByName(cols, "pl7.app/hydrophobicity", domain),
    };
  }
  // Antibody/TCR — IG/TCRAB/TCRGD, bulk or single-cell. Domain is identical
  // across receptors; the user-facing label (CDR-H3 / CDR-α3 / CDR-γ3) is
  // already encoded in spec annotations by the workflow (R13a).
  const domain = {
    "pl7.app/feature": "CDR3",
    "pl7.app/vdj/scClonotypeChain": "A",
  };
  return {
    x: pickByName(cols, "pl7.app/charge", domain),
    y: pickByName(cols, "pl7.app/hydrophobicity", domain),
  };
}

// R20: default histogram metric = hydrophobicity (per modality).
export function defaultHistogramMetric(
  cols: PColumnIdAndSpec[],
  info: WorkflowInfo | undefined,
): PColumnSpec | undefined {
  if (info?.mode === undefined) return undefined;
  if (isPeptideMode(info.mode)) {
    return pickByName(cols, "pl7.app/hydrophobicity", { "pl7.app/feature": "peptide" });
  }
  return pickByName(cols, "pl7.app/hydrophobicity", {
    "pl7.app/feature": "CDR3",
    "pl7.app/vdj/scClonotypeChain": "A",
  });
}

// R19a / R20a: emission-order fallback when the modality default is absent.
export function numericScalarsInOrder(cols: PColumnIdAndSpec[]): PColumnSpec[] {
  return cols.filter((c) => isNumericScalar(c.spec)).map((c) => c.spec);
}
