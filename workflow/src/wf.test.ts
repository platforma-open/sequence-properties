// Workflow integration tests will be added once Python property computation lands.
// Placeholder kept so vitest discovers the suite without failures.
import { describe, it } from "vitest";

describe("sequence-properties workflow", () => {
  it.skip("computes peptide properties end-to-end", () => {
    // TODO: implement once compute_properties.py is real.
  });
});
