# Performance review — `software/`

After the refactor pass (104 tests, 97% coverage). Profiled on a 2025 Apple
Silicon dev box — single-threaded, no warmup tricks.

## Current numbers

| Workload | Per-row | At 10⁴ rows | At 10⁵ rows | At 10⁶ rows |
|---|---:|---:|---:|---:|
| Peptide (~15-20 aa) full panel | 49 µs | 0.5 s | 5 s | 49 s |
| Antibody chain (~120 aa) full panel | 243 µs | 2.4 s | 24 s | 4 min |
| pI alone (per ~120 aa chain) | 190 µs | 1.9 s | 19 s | 3 min |

**Where the time goes:** pI bisection dominates antibody chain cost (78%
of the per-chain budget). The driver is ~40 calls to `charge_at_ph` per pI;
each call re-prepares (validates + cleans) the already-cleaned sequence
and re-iterates the Python residue loop.

## Optimisations — effort vs payoff

Sorted by `payoff / effort`. "Payoff" is the speedup factor on the relevant
hot path; "Effective at" is the dataset size where it starts to matter.

### Tier A — cheap and impactful

| # | Change | Effort | Payoff | Effective at | Risk |
|---|---|---|---:|---|---|
| 1 | **Cache cleaned sequence in `isoelectric_point`** — pass cleaned into a `_charge_at_ph_cleaned` private variant; bisection lambda calls that. Eliminates 40 redundant `_prepare` calls per pI. | 30 min | **~1.7x** on pI (≈ 1.5x on antibody pipeline) | every dataset | trivial — golden values stay stable |
| 2 | **`Counter`-based `aa_counts`** — replace per-char Python loop with `collections.Counter(seq.upper())`, then map standard codes. C-implemented inner loop. | 30 min | **~2-3x** on `aa_counts`; `extinction_coefficients` / `aliphatic_index` / `aromaticity` / `aa_fractions` all benefit | every dataset | trivial |
| 3 | **`SequenceContext` per row** — compute upper-cased + cleaned + counts ONCE per input sequence; pass into property fns instead of letting each fn re-derive. | 3 h | **~3-5x** on peptide-mode pipeline; ~2x on antibody | 10³+ | medium — touches every property fn signature, but golden values + corpus tests pin behaviour |

If you do nothing else, do #1+#2 — under an hour, and the antibody pipeline
gets ~3× faster end-to-end with zero test churn.

### Tier B — medium effort, scale-dependent

| # | Change | Effort | Payoff | Effective at |
|---|---|---|---:|---|
| 4 | **Numpy DIWV lookup** — convert the 20×20 DIWV dict-of-dicts to a numpy `int8`-keyed array; map sequence chars to ints once, do `arr[ints[:-1], ints[1:]].sum()`. | 1 h | ~5-10× on `instability_index` alone (~10% of antibody chain budget today, so ~5% end-to-end) | 10⁴+ |
| 5 | **Numpy residue-charge vector** — precompute pKa / acid-mask / base-mask arrays indexed by `ord(c)`; `_residue_charge` becomes a vectorised lookup. Especially impactful inside the pI inner loop. | 2 h | ~2-3× on `charge_at_ph`, ~2× on pI | 10⁴+ |
| 6 | **Build columns directly, skip per-row dicts** — pipeline currently builds `[{col1: v, col2: v, …}, …]`. Build `{col1: [v, v, …], col2: …}` instead; less dict allocation, faster polars schema inference. | 1 h | ~1.3-1.5× on the pipeline tail | 10⁵+ |

### Tier C — high effort, only worth it for scale

| # | Change | Effort | Payoff | Effective at |
|---|---|---|---:|---|
| 7 | **Vectorise charge over sequence batches** — represent N peptides as an `(N, 20)` count matrix; charge per pH = matrix · `(20,)` pKa-derived vector — single numpy op for the whole batch. | 1 day | ~10-100× at scale | 10⁵+ |
| 8 | **Vectorise pI bisection over the batch** — same bisection algorithm, but evaluating charge on the whole vector each iteration. Couples to #7. | 1 day | ~50-200× at scale | 10⁵+ |
| 9 | **Polars expressions for AA fraction** — push the long-form `(entity, aa, fraction)` build into native polars via `str.count_match` per AA; avoids the Python row-construction loop entirely. | 1 day | ~5-10× on AA fraction generation alone | 10⁶+ |

The Tier-C set is what the spec was hinting at with "Typical single-cell
dataset should complete in <60 seconds with vectorized implementation."
Without it, 10⁶ antibody clones is firmly minutes-long; with it, sub-minute.

## Counter-incentives — things NOT to do

- **Don't shrink the bisection tolerance** (currently 0.001 → ~40 iterations).
  Reporting precision is `.2f`, so 0.0005 is overkill — but golden-value
  tests pin the current value to 1e-6, so loosening the tolerance breaks
  the test suite without a payoff anywhere users see. Bundle this with
  Tier A/B if at all.
- **Don't cache properties at the workflow / Tengo level.** The dataset
  shape varies between runs; cache miss-rate would be ~100%. Caching
  belongs inside Python, on per-sequence work, not across runs.

## Recommended sequencing

1. **Now (1 hour, ~3× end-to-end):** apply Tier-A items #1 and #2.
2. **When dataset hits ~10⁵ entities (~5 hours, ~5-10× end-to-end):** add #3
   and #4.
3. **When SC datasets land (multi-day, ~50-100×):** Tier C — #7 + #8 — as a
   single rewrite. Will require a new internal API (batch property
   functions); golden-value tests still pin per-sequence correctness.

The 104-test suite + 16 golden values + 18 corpus cases is enough safety
net to do any of these without behaviour drift.
