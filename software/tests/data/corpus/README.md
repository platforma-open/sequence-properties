# Sequence-properties E2E corpus

Static input TSVs + a `manifest.json` describing the expected pipeline behavior
per case. Tests in `tests/integration/test_corpus_e2e.py` load these and run
`pipeline.run` against the committed manifest.

## Files

- `peptide_input.tsv` — peptide-mode entity TSV.
- `peptide_plan.json` — peptide-mode plan.
- `antibody_input.tsv` — antibody/TCR-mode entity TSV (covers IG full-coverage,
  CDR3-only, partial-region, single-chain, and TCRAB cases).
- `antibody_plan.json` — antibody/TCR-mode plan (legacy MiXCR bulk, IG, both
  chains full → Fv emitted).
- `tcr_plan.json` — same input rows, plan switched to TCRAB (no Fv).
- `manifest.json` — per-case expected behavior. See "Manifest structure".

## Manifest structure

```jsonc
{
  "peptide": {
    "cases": {
      "<entity_key>": {
        "notes": "why this case exists",
        "expected": {
          "<column>": "na" | <number> | { "min": x, "max": y }
        }
      }
    }
  },
  "antibody": { ... },
  "tcr": { ... }
}
```

`"na"` asserts the cell is NA (None / null). A bare number asserts an exact
match (use only when the value is closed-form — e.g. instability index of a
pure-A 10-mer is 9.0). A `{min, max}` range guards against silent drift in
formulas without locking to a specific numeric.

## Updating

The corpus is hand-crafted, not generated. Add new cases by:
1. Append a row to the relevant TSV.
2. Append an entry to `manifest.json`.
3. Run `uv run pytest tests/integration/test_corpus_e2e.py` and confirm pass.

Commit input + manifest changes together with any code change that required
a new case.
