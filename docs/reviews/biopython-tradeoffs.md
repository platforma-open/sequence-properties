# BioPython Migration — Tradeoffs and Wrapper Layer

**Date:** 2026-04-29
**Context:** Migration from custom Henderson-Hasselbalch / DIWV / Pace
implementations to `Bio.SeqUtils.ProtParam.ProteinAnalysis` and
`Bio.SeqUtils.IsoelectricPoint.IsoelectricPoint`, per spec M1 strategy
direction (`docs/text/work/projects/sequence-properties/README.md` L520-522).

## Verdict

BioPython is the correct strategy. Every property the spec requires either
maps directly onto BioPython's API or is preserved by a thin wrapper that
adapts BioPython behaviour to spec contracts. **No case was found where
reimplementing the property from scratch is preferable.**

The wrapper layer below is glue, not duplication of BioPython math.

## Where BioPython is used directly

| Property | BioPython entry point |
|---|---|
| GRAVY | `ProteinAnalysis.gravy()` |
| Molecular weight | `ProteinAnalysis.molecular_weight()` |
| Aromaticity | `ProteinAnalysis.aromaticity()` |
| Instability index | `ProteinAnalysis.instability_index()` |
| Extinction coefficient | `ProteinAnalysis.molar_extinction_coefficient()` |
| AA percent / fraction | `ProteinAnalysis.amino_acids_percent` |
| Charge at pH | `IsoelectricPoint.charge_at_pH(pH)` |
| pI (default range) | `IsoelectricPoint.pi()` (overridden — see below) |

Underlying tables (residue masses, Kyte-Doolittle hydropathy, DIWV matrix,
Pace ε constants) all come from BioPython after this migration. Our local
`aa_tables.py` constants for these are now unused; only `STANDARD_AAS` /
`STANDARD_AA_SET` are still referenced (sequence cleanup).

## Wrapper layer — where in-house code remains

These wrappers are necessary because BioPython's defaults do not match the
spec. Each is a small, well-bounded adaptation, not a reimplementation.

### 1. pI bisection range — must be [0, 14], not BioPython's [4.05, 12]

**BioPython behaviour:** `IsoelectricPoint.pi()` brackets the bisection
between pH 4.05 and 12.0 internally. Sequences with no charge zero-crossing
in that range (polybasic synthetic sequences with median pI > 12, polyacidic
with median pI < 4) get clamped to a boundary value with a large residual
charge.

**Spec requirement:** "NA when the sequence has no zero crossing in [0, 14]."

**Fix:** Reuse BioPython's `charge_at_pH` (with our IPC 2.0 pKa overrides),
but drive bisection locally over [0, 14] (`_bisect_charge_zero` in
`properties.py`). Same-sign endpoints → return None.

**Verified by:** `pep_polybasic_150` and `pep_polyacidic_30` corpus cases
both correctly return NA.

### 2. IPC 2.0 pKa override

**BioPython default:** Bjellqvist pKa values (its module-level
`positive_pKs` / `negative_pKs` / `pKnterminal` / `pKcterminal`).

**Spec requirement:** IPC 2.0 pKa values (Kozlowski 2021), peptide set for
peptide / CDR3 sequences, protein set for full-chain.

**Fix:** Construct an `IsoelectricPoint(seq)` instance, then overwrite
`ip.pos_pKs` and `ip.neg_pKs` with our IPC 2.0 dicts before calling
`charge_at_pH` or `pi`. Per-instance override avoids the thread-safety
hazard of monkey-patching module globals.

### 3. Cys exclusion for full-chain mode

**BioPython behaviour:** Cys is always in `neg_pKs` (free thiol assumption).

**Spec requirement:** Full-chain mode treats Cys as disulfide-bonded — Cys
must NOT contribute to the ionizable sum. CDR3 and peptide modes keep Cys
ionizable.

**Fix:** `_ipc2_isoelectric_point(seq, pka_set, include_cys)` conditionally
omits the Cys entry from `neg_pKs`. BioPython's `charge_at_pH` only iterates
keys in `neg_pKs`, so an absent entry produces no contribution — clean.

### 4. Instability index 10-residue floor

**BioPython behaviour:** No length check — computes a value for any
sequence ≥ 2 residues.

**Spec requirement:** NA when effective length < 10.

**Fix:** `instability_index` checks `len(cleaned) < 10` before calling
BioPython. Spec floor preserved.

### 5. Empty / stop-codon / ambiguity-code sequences

**BioPython behaviour:** Raises `ValueError` on `'X'`, `'B'`, `'Z'`, `'*'`,
etc.; raises `IndexError` on empty sequences; counts ambiguity codes
incorrectly in some methods.

**Spec requirement:** Empty + stop codon → NA for all properties; ambiguity
codes filtered residue-by-residue, NOT invalidating the whole sequence.

**Fix:** `_prepare(seq)` cleans (uppercase + filter to standard AAs) and
returns None for empty/stop. Every BioPython call is gated by
`_prepare`. BioPython operates only on sanitised input.

### 6. Aliphatic index

**BioPython behaviour:** Not provided. BioPython has `flexibility`,
`secondary_structure_fraction`, and `protein_scale`, but no Ikai aliphatic
index helper.

**Fix:** Retain the custom Ikai 1980 closed-form computation:
`AI = 100 × (X_A + 2.9·X_V + 3.9·(X_I + X_L))`. No BioPython equivalent
exists.

### 7. Fv pI — paired-chain charge sum bisection

**BioPython behaviour:** `IsoelectricPoint` operates on a single sequence.
There is no paired-chain charge sum primitive.

**Spec requirement:** Fv pI is the pH where `charge(VH, pH) + charge(VL, pH)
= 0`, NOT the pI of the concatenated string (concatenation loses one
terminus pair and adds a fake peptide bond).

**Fix:** Construct two `IsoelectricPoint` instances (VH, VL with IPC 2.0
overrides + Cys exclusion) and bisect their charge sum locally. Reuses
`_bisect_charge_zero`.

### 8. Extinction coefficient tuple order

**BioPython behaviour:** `molar_extinction_coefficient()` returns
`(reduced, oxidized)`.

**Spec / TSV contract:** Output columns `eox_*` (oxidized) come before
`ered_*` (reduced); functions return `(eox, ered)`.

**Fix:** `extinction_coefficients` swaps the BioPython tuple to
`(oxidized, reduced)`.

### 9. AA fraction unit (percent → fraction)

**BioPython behaviour:** `amino_acids_percent` returns values in PERCENT
(sum = 100.0).

**Spec / PColumn contract:** AA fraction PColumn values are FRACTIONS
(sum = 1.0).

**Fix:** `aa_fractions` divides each percent by 100.

### 10. Mass-value drift vs prior NIST/UniMod table

**Observation:** BioPython's residue mass table differs from the NIST /
UniMod table previously used by ~0.001 Da per residue, accumulating to
~0.07 Da on a 54-residue chain and ~0.13 Da on a paired Fv (~104 residues).

**Spec direction:** "MW for ≥5 VH sequences agree to 2 decimal places"
against BioPython — so BioPython is the reference. Our display format is
`.1f`, so the drift is below the displayed precision.

**Fix:** Pinned golden values in `test_golden_values.py` and
`test_three_glycines` updated to BioPython output.

## Cases evaluated and rejected

| Question | Answer |
|---|---|
| Could we use `IsoelectricPoint.pi()` and skip the local bisect? | No — its [4.05, 12] bracket fails spec's [0, 14] NA rule. |
| Could we monkey-patch `IPMod.positive_pKs` globally for IPC 2.0? | No — process-wide and not thread-safe. Per-instance override is safer. |
| Could we drop the 10-residue floor for instability_index? | No — spec mandates NA below 10. |
| Could we drop sequence cleanup and let BioPython raise? | No — the spec convention is "ambiguity codes filtered, not raised". |
| Could we use BioPython for aliphatic index? | No — BioPython does not provide it. |
| Could we represent Fv pI via concatenation? | No — explicitly forbidden by spec; loses one terminus pair. |

## Net code shape after migration

- `properties.py`: thin wrappers over BioPython + 1 local bisection helper +
  custom aliphatic_index. ~250 lines (down from ~360 lines pre-migration).
- `instability.py`: deleted (BioPython has the same DIWV table).
- `aa_tables.py`: `STANDARD_AAS` / `STANDARD_AA_SET` still used; mass /
  hydropathy / EC tables now unreferenced (kept for now, not load-bearing).
- `pka_tables.py`: `PKaSet` and IPC2 sets still load-bearing for the pKa
  override layer.
