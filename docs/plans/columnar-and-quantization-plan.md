# Polars columnar implementation plan + quantization re-evaluation

Two related items, kept in one doc because the polars conversion changes
which floating-point paths matter and therefore which columns need
quantization for CID stability.

---

## Part A — Polars columnar implementation plan

### Goals

1. Add a **polars-native fast path** for peptide mode (everything except pI
   and instability does well as columnar work).
2. **Vectorise pI bisection** across the batch using numpy on a pre-built
   residue-count matrix.
3. **Vectorise instability** using a numpy gather on a 20×20 DIWV array.
4. Keep antibody/TCR mode's outer loop in Python (per-clone region
   reconstruction + NA logic doesn't suit columnar) — but lift the
   inner property panels onto the same batch primitives so VH / VL /
   Fv panels get the same speedup as peptide mode.
5. Preserve every value pinned by `tests/unit/test_golden_values.py`,
   `test_properties.py`, `test_pipeline.py`, and the manifest-driven
   corpus. Sequence-level correctness has 104 tests guarding it.

### Two-path architecture

```
                ┌─────────────────────────────────────┐
                │           pipeline.run()            │
                └────────────────┬────────────────────┘
              peptide mode       │      antibody/TCR mode
            ┌──────────┴──────────┐         │
            │                     │         │
   N < THRESHOLD             N ≥ THRESHOLD  │
   single-row Python         columnar batch │
   (existing path)                          │
            │                     │         │
            └──────────┬──────────┘         ▼
                       ▼              row-by-row Python loop
                                      over reads.iter_rows()
                                              │
                                              ▼
                              for each chain reconstruction:
                              call batch primitive on a small
                              list of strings (still wins via
                              numpy gather + dot products)
```

THRESHOLD ≈ 1 000 rows. Below it, polars expression-graph overhead
dominates the actual work; above it, columnar wins by 5-100×.

### Module layout

| New / changed | Role |
|---|---|
| `src/properties.py` | unchanged — single-sequence functions stay as the fallback path and as the spec for behaviour |
| `src/properties_batch.py` | **new** — batch primitives: residue-count matrix, vectorised charge / pI / instability; closed-form scalar properties as polars expressions |
| `src/aa_tables.py` | add precomputed numpy lookup tables: `KD_VECTOR[20]`, `MASS_VECTOR[20]`, `ALIPHATIC_WEIGHTS[20]`, `AROMATIC_MASK[20]`, `EC_VECTOR[20]`, `CHAR_TO_IDX[256]`, `DIWV_MATRIX[20, 20]`. All built once at import time. |
| `src/pipeline.py` | dispatch at top of `run_peptide` / `run_antibody_tcr` based on row count; call into `properties_batch` for the fast path |
| `tests/unit/test_properties_batch.py` | **new** — verify each batch primitive matches the per-row property fns to 1e-12 |
| `tests/unit/test_pipeline.py` | extend with size-parametrised cases (N=1, N=10, N=THRESHOLD-1, N=THRESHOLD, N=10×THRESHOLD) — both paths must produce identical output |
| existing tests | unchanged; pass on both paths |

### Phase 1 — count-matrix primitive (~3 hours)

New `properties_batch.residue_count_matrix(seqs: pl.Series)`:

1. **Validity column** via polars expressions:
   ```python
   is_invalid = seqs.is_null() | (seqs.str.len_chars() == 0) | seqs.str.contains(r"\*", literal=False)
   ```
2. **Cleanup** via polars regex:
   ```python
   cleaned = seqs.str.to_uppercase().str.replace_all(r"[^ACDEFGHIKLMNPQRSTVWY]", "")
   effective_length = cleaned.str.len_chars()
   ```
3. **Counts** — 20 polars columns via `str.count_matches`:
   ```python
   counts_df = pl.DataFrame({aa: cleaned.str.count_matches(aa).cast(pl.Int32) for aa in STANDARD_AAS})
   ```
4. **Convert to numpy** (`(N, 20)` int32 array) for downstream gather / dot operations.

Returns `(is_invalid: pl.Series, effective_length: pl.Series, counts_df: pl.DataFrame, counts_np: np.ndarray)`.

Behavioural test: for a 100-row sample including invalid / non-standard /
empty rows, the per-row `aa_counts(seq)` and `effective_length(seq)`
must agree with the matrix output. Add to `test_properties_batch.py`.

### Phase 2 — closed-form properties as polars expressions (~3 hours)

Each property becomes a single polars expression on the counts_df:

```python
def gravy_expr(counts_df: pl.DataFrame, eff_len: pl.Series) -> pl.Series:
    weighted = sum(pl.col(aa) * KD_SCALE[aa] for aa in STANDARD_AAS)
    return (weighted / eff_len).fill_nan(None)
```

Similar:
- `mw_expr`: `sum(count * mass) + 18.0153`
- `extinction_ox_expr`: `n_Y * 1490 + n_W * 5500 + (n_C // 2) * 125`
- `extinction_red_expr`: `n_Y * 1490 + n_W * 5500`
- `aliphatic_expr`: closed-form on n_A, n_V, n_I, n_L
- `aromaticity_expr`: `(n_F + n_W + n_Y) / eff_len`

NA propagation: every expression wrapped in
`pl.when(is_invalid | (eff_len == 0)).then(None).otherwise(expr)`.

Behavioural test: column-by-column equality vs the per-row functions on
the 100-row sample, abs tolerance 1e-12 (these are exact arithmetic — no
transcendentals).

### Phase 3 — vectorised charge + pI bisection (numpy) (~4 hours)

Charge depends on `10**x` and is the only place transcendentals appear.
Implementing it in polars expressions is awkward; numpy is the right
tool.

**Charge primitive** (`(N,)` charge array given counts and pH):

```python
def batch_charge(counts: np.ndarray, ph: float, pka_set: PKaSet, include_cys: bool) -> np.ndarray:
    # Precompute coefficients for the 7 ionizable residues:
    #   acid_coef[a] = -1 / (1 + 10**(pKa[a] - pH))
    #   base_coef[b] = +1 / (1 + 10**(pH - pKa[b]))
    # Apply Cys mask if include_cys is False.
    coefs = _ionic_coefs(pka_set, ph, include_cys)  # (20,) float64
    n_term = 1.0 / (1.0 + 10.0 ** (ph - pka_set.n_terminus))
    c_term = -1.0 / (1.0 + 10.0 ** (pka_set.c_terminus - ph))
    return counts @ coefs + n_term + c_term  # shape (N,)
```

One matrix-vector multiply replaces N × O(L) Python loops.

**Vectorised pI** (`(N,)` pI array via batch bisection):

```python
def batch_pi(counts: np.ndarray, pka_set: PKaSet, include_cys: bool, tol: float = 1e-3) -> np.ndarray:
    lo = np.zeros(N); hi = np.full(N, 14.0)
    f_lo = batch_charge(counts, 0.0, pka_set, include_cys)
    f_hi = batch_charge(counts, 14.0, pka_set, include_cys)
    # Same-sign rows → NA (no zero crossing). Mark and skip.
    no_crossing = np.sign(f_lo) == np.sign(f_hi)
    # Bisection: vectorised midpoint update
    for _ in range(_iterations_for(tol)):
        mid = 0.5 * (lo + hi)
        f_mid = batch_charge(counts, mid, pka_set, include_cys)  # but this is per-pH, not per-row pH!
        ...
```

**Subtlety**: `batch_charge` evaluates one pH for all rows. Bisection wants
*different pH per row*. So either:

- Outer loop iterates pH-rows-in-lockstep (every row gets the same number
  of bisection steps; this is the standard approach), OR
- Use scipy.optimize.brentq per row (per-row but C-implemented; surprisingly
  competitive for N up to ~10⁴).

Standard approach is the lockstep bisection, accepting that some rows
"converge" earlier but we still iterate to global tolerance. Same number
of inner ops as today (~40 charge calls × N rows) but each op is a
single matrix-vector multiply instead of N Python loops.

Behavioural test: pI matches single-row computation to abs tolerance
1e-3 (matches bisection tolerance). May need 1e-4 if numpy `10**x`
produces ULP-different intermediates from Python's `**`.

### Phase 4 — vectorised instability (~2 hours)

DIWV as a `(20, 20)` numpy array. CHAR_TO_IDX as a 256-int LUT
(`CHAR_TO_IDX[ord('A')] = 0`, etc.; non-standard ords map to -1).

Per-batch:

```python
def batch_instability(cleaned_seqs: list[str], eff_lens: np.ndarray) -> np.ndarray:
    out = np.full(len(cleaned_seqs), np.nan)
    for i, (s, L) in enumerate(zip(cleaned_seqs, eff_lens)):
        if L < 10: continue
        ints = np.frombuffer(s.encode("ascii"), dtype=np.uint8)
        idx = CHAR_TO_IDX[ints]              # (L,)
        out[i] = (10.0 / L) * DIWV_MATRIX[idx[:-1], idx[1:]].sum()
    return out
```

Per-row Python (the outer loop) but per-row work is one numpy gather +
sum. ~5-10× faster than the current Python dipeptide loop. Pure-numpy
batching across rows isn't practical because dipeptide indices have
variable length per row.

(If we *really* want to push: pad sequences to a uniform length, mask
the pad positions in the gather. Adds complexity for diminishing
returns. Defer.)

### Phase 5 — AA fraction long-form via polars melt (~1 hour)

After Phase 1 we have `counts_df` (N rows × 20 columns). For peptide
mode:

```python
fractions_df = counts_df.with_columns([pl.col(aa) / eff_len for aa in STANDARD_AAS])
long_df = fractions_df.with_columns(entity_key=keys).melt(
    id_vars="entity_key",
    value_vars=list(STANDARD_AAS),
    variable_name="aminoAcid",
    value_name="value",
)
```

Native polars unpivot. ~10-20× faster than the Python row-by-row build
on 10⁴+ entities.

Invalid-row handling: replace fraction columns with None where
`is_invalid` before melting, so each invalid entity emits 20 NA rows
(maintaining the uniform 20-rows-per-entity shape required by the
2-axis PColumn).

### Phase 6 — pipeline glue (~2 hours)

```python
BATCH_THRESHOLD = 1000

def run_peptide(reads: pl.DataFrame) -> dict[str, pl.DataFrame]:
    if reads.height < BATCH_THRESHOLD:
        return _run_peptide_python(reads)   # current implementation
    return _run_peptide_batch(reads)

def _run_peptide_batch(reads):
    # Phase 1: count matrix + cleaned seqs
    is_invalid, eff_len, counts_df, counts_np = residue_count_matrix(reads["peptide_seq"])
    # Phase 2: closed-form props as polars cols
    properties = pl.DataFrame({
        "entity_key": reads["entity_key"],
        "charge_peptide": batch_charge(counts_np, 7.0, IPC2_PEPTIDE, True),  # numpy
        "gravy_peptide": gravy_expr(...).to_list(),
        ...
    })
    # Phase 3: vectorised pI on the count matrix
    properties = properties.with_columns(pi_peptide=batch_pi(...))
    # Phase 4: instability per-row numpy gather
    properties = properties.with_columns(instability_peptide=batch_instability(...))
    # NA propagation: where is_invalid, replace every property with None
    properties = _apply_invalid_mask(properties, is_invalid)
    # Phase 5: AA fraction long-format via melt
    aa_fraction = _aa_fraction_long(reads["entity_key"], counts_df, eff_len, is_invalid)
    return {"properties": properties, "aa_fraction": aa_fraction}
```

For antibody mode, the outer Python loop over clones stays. Inside,
each chain reconstruction produces a single string; we batch the
reconstructed-chain strings into one `_run_chain_batch(seqs, pka_set)`
call per chain group:

```python
def run_antibody_tcr(reads, plan):
    # Outer Python: per-clone reconstruction + NA logic.
    clones = list(reads.iter_rows(named=True))
    cdr3s_A = [c.get("A_CDR3") or "" for c in clones]
    cdr3s_B = [c.get("B_CDR3") or "" for c in clones]
    chain_A = [_reconstruct_chain(c, "A") for c in clones]  # list of str|None
    chain_B = [_reconstruct_chain(c, "B") for c in clones]
    # Inner batch: each list of N strings goes through the batch primitives
    cdr3_A_props = batch_compute(cdr3s_A, IPC2_PEPTIDE, include_cys=True)
    full_A_props = batch_compute([c or "" for c in chain_A], IPC2_PROTEIN, include_cys=False)
    ... # similarly for B and Fv (Fv via per-chain results)
    # Assemble final DataFrame in one go (no per-row dict construction).
```

### Phase 7 — tests (~2 hours)

1. **Path equivalence**: parametrise key tests with `[N=1, 10, 100, 999, 1000, 10000]` so the dispatch flips paths inside the test. Both paths must produce the same output (column-by-column abs tol 1e-3 for charge/pi where libm vs numpy differ; 1e-12 for closed-form properties).
2. **Batch primitive tests**: each batch fn vs single-row baseline on 100-row randoms.
3. **Performance regression**: a `@pytest.mark.slow` test that asserts the batch path runs in under (current_python_path_time / 5) at N=10 000.

### Effort summary

| Phase | Effort |
|---|---|
| 1: count-matrix primitive | 3 h |
| 2: closed-form polars expressions | 3 h |
| 3: vectorised charge + pI | 4 h |
| 4: vectorised instability | 2 h |
| 5: AA fraction melt | 1 h |
| 6: pipeline glue + dispatch | 2 h |
| 7: tests | 2 h |
| **Total** | **17 h ≈ 2.5 days** |

### Risks

- **ULP differences between Python `**` and numpy `**`** — already on the
  table; goldens at 1e-6 may need to relax to 1e-4 OR we add a "batch
  path" set of goldens with their own pinned values.
- **Polars version pin tightening** — minor expression API changes
  between polars versions are common. Pin `polars-lts-cpu==1.33.1`
  exactly in `requirements.txt` (already done).
- **Numpy reduction order** — `np.sum` on small arrays uses different
  strategies than on large arrays (chunked vs simple). Same-version
  same-machine deterministic; cross-version possibly not.

### Out of scope for this plan

- GPU acceleration (numpy CPU is plenty up to 10⁷ rows).
- Any change to `process.tpl.tengo` — the column contract is unchanged;
  this is an internal compute substitution.
- Any change to single-sequence `properties.py` API; it remains the
  fallback path and the documentation of intent.

---

## Part B — Quantization re-evaluation under the CID-mismatch lens

### What CID mismatch means here (block-dev)

Platforma uses content-addressable storage. A workflow output's CID is
derived from its bytes. If two runs of the same block on the same input
produce byte-identical outputs, downstream cache hits fire. Mismatches
trigger redundant downstream recomputation.

The operator's frame: **same computer, same input, similar runs — CIDs
should match.** That's the bar. Cross-machine CID stability is a
different (harder) problem and not what we're solving here.

### What's actually deterministic on a single machine today

Pure Python with no randomness, fixed iteration order, IEEE-754 FP:

| Operation | Same-machine same-version determinism |
|---|---|
| `int / int` (e.g. AA fraction = `n / total`) | bit-exact, every run |
| `sum(constants)` (e.g. MW) | bit-exact, fixed iteration order |
| `sum(constants) * scalar` (instability, aliphatic) | bit-exact |
| Boolean / integer counts (ε from Y / W / C counts) | bit-exact |
| `aromaticity = n / L` | bit-exact |
| `gravy = sum(KD[c]) / L` | bit-exact (KD constants are floats but sum order is fixed) |
| `10**x` via `pow` libm call (charge, pI inner loop) | **deterministic on same machine + libm**; varies across libm versions |
| pI bisection over the libm-charge | **deterministic** while libm doesn't change |

So on a single machine, with a fixed Python interpreter and a fixed
libm, the current code already produces CID-stable output.

### Where same-machine CID could still drift

1. **Library upgrade in the dev environment**. A `uv sync` after a
   minor numpy / polars / Python patch release could change `10**x`
   for charge / pI by a ULP. Two runs separated by a `uv sync` are
   different bytes, different CIDs.
2. **Code path change**. Today's pure Python charge vs tomorrow's
   numpy charge produce values that differ by ~1 ULP for some inputs.
   The output_value at .2f format is identical, but the un-rounded
   stored bytes differ.
3. **Runtime environment**. Same dev computer, but the block runs
   in a Docker image whose libm is different from the host's. CIDs
   diverge.

(1) and (3) are real same-machine concerns. (2) becomes real the moment
the columnar plan above lands.

### Where quantization helps for CID

Only on columns whose values depend on a transcendental (`10**x`) or
on a sum-order-sensitive reduction. Everything else is bit-exact under
IEEE-754 closed-form arithmetic and quantizing it is wasted churn.

| Column | FP path | CID risk | Quantize? |
|---|---|---|---|
| `charge_peptide`, `charge_<chain>_CDR3`, `charge_<chain>_VDJRegion`, `charge_Fv` | `10**x` via libm; ULP variance across libm versions and Python→numpy code-path swaps | **real** | **YES → 2 decimals** (matches `.2f` display, well within bisection tolerance) |
| `pi_peptide`, `pi_<chain>_VDJRegion`, `pi_Fv` | bisection on libm-charge → propagates the same risk | **real** | **YES → 2 decimals** (matches `.2f`; bisection tolerance is 0.001 anyway) |
| `gravy_*` | `sum(KD_constants) / int` — closed-form, no transcendentals | none on same machine | **NO** |
| `mw_*` | `sum(mass_constants) + 18.0153` — closed-form | none | **NO** |
| `eox_*`, `ered_*` | integer arithmetic cast to float | none | **NO** (already exact) |
| `instability_*` | `(10/L) * sum(table_lookups)` — closed-form | none | **NO** |
| `aliphatic_*` | mole fractions × scalar coefficients — closed-form | none | **NO** |
| `aromaticity_*` | `int / int` | none | **NO** |
| `aa_fraction.value` | `int / int` — and the sum-to-1 invariant breaks if quantized per-cell | none + would harm | **NO** (operator direction; also FP-exact) |

So under the CID-mismatch lens: **quantize only `charge_*` and `pi_*` columns to 2 decimals.**

### Updated recommendation

Drop in a `_quantize_for_cid(properties: pl.DataFrame) -> pl.DataFrame`
helper at the tail of `pipeline.run` (both peptide and antibody
branches). Single polars `with_columns` operation:

```python
CID_QUANTIZED_PREFIXES = ("charge_", "pi_")

def _quantize_for_cid(df: pl.DataFrame) -> pl.DataFrame:
    cols_to_round = [
        c for c in df.columns
        if any(c.startswith(p) for p in CID_QUANTIZED_PREFIXES)
    ]
    if not cols_to_round:
        return df
    return df.with_columns([pl.col(c).round(2) for c in cols_to_round])
```

That's the entire change. Two columns of behaviour — both already at
`.2f` display precision, both with ULP-variance risk under the polars
plan, both inside bisection tolerance.

### Test impact

- Golden values currently pinned at 1e-6 abs tolerance. After
  quantization:
  - `pi VH = 7.018372` becomes `7.02` in the *output* — but the test asserts on `isoelectric_point()` directly, which is the un-quantized internal function. **Goldens unchanged.**
  - The corpus test reads the *output DataFrame* via `pipeline.run`, so its assertions on `pi_*` see quantized values. The corpus manifest currently has `pi` as `"not_na"` or `{"min": 1.0, "max": 4.5}` ranges — quantization to 2 decimals doesn't break either form.
- One new test in `test_pipeline.py`: assert that `output["properties"]["charge_peptide"]` and `pi_peptide` are rounded to 2 decimals (column-level invariant).

### Effort

30 minutes. Three lines of helper code, one new test, one test re-run
to confirm everything green.

### What we explicitly do NOT do

- Quantize MW / GRAVY / aliphatic / instability / aromaticity / ε / aa_fraction. They are FP-deterministic on closed-form arithmetic and quantization would be noise.
- Add a config flag for the quantization. It's fixed behaviour. Anyone wanting full precision computes from the input sequences directly — that's what `properties.py` provides as the public API.
- Round inside the single-row functions. Quantization is a *boundary* concern (inter-process / CID stability). Internal computations stay full-precision so golden tests stay sharp.
