# Sequence-Properties Parallel Compute Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers-extended-cc:subagent-driven-development (recommended) or superpowers-extended-cc:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let the `sequence-properties` block compute a huge dataset on the first run by parallelizing the Python property kernel across CPU cores, behind a frozen characterization net that proves the output is byte-identical to today's.

**Architecture:** Three phases. **Phase A (characterization, do first):** freeze a full-output golden master + structural-invariant scale test from the *current* sequential code, so any behavior change in later phases fails loudly. **Phase B (parallelism, option #1):** map the existing pure per-row compute functions across a `ProcessPoolExecutor`, assemble results in input order, keep the existing sort + CID-quantization boundary — so worker count never changes a byte; bump the workflow's `.cpu()`/`.mem()`. **Phase C (streaming, option #7, OPTIONAL):** stream sorted input in batches to bound peak memory; gated on a measurement decision after Phase B.

**Tech Stack:** Python 3.12, `polars` (TSV I/O + frame assembly), `biopython` (property math), `concurrent.futures.ProcessPoolExecutor` (parallelism), `pytest` + `uv` + `ruff` (test/lint, run from `software/`), Tengo workflow (`exec.builder()` resource requests).

---

## Before You Start (read once, do not skip)

1. **Worktree + branch.** Do not work on `main` (workspace rule). Create a worktree for the block repo under `worktrees/sequence-properties/<branch>/` (use the `worktree` skill). Ask the operator whether this is tied to a Notion ticket (`MILAB-NNNN`); if yes, name the branch `MILAB-NNNN_seqprops-parallel-compute` and prefix every commit subject `MILAB-NNNN: …`. If no ticket, use `feat/seqprops-parallel-compute`.
2. **Move this plan into the worktree** and commit it as the first commit on the branch, so it travels with the work.
3. **Baseline green run** — from `blocks/sequence-properties/software/`:
   ```bash
   uv sync --locked
   uv run pytest
   uv run ruff check
   uv run ruff format --check
   ```
   All must pass before you change anything. If they don't, stop and tell the operator — you cannot characterize behavior on a red baseline.
4. **All Python commands run from `blocks/sequence-properties/software/`** (the dir holding `pyproject.toml`). All `pnpm` build commands run from the block root `blocks/sequence-properties/`.

## The Determinism Contract (the invariant every task must preserve)

The block's output must be **byte-identical run-to-run and machine-to-machine** so the platform's content-addressed cache dedups it. The current code earns this with two mechanisms — do not break either:

- **Stable row order:** output rows are sorted by key before writing (`io_layer.py:48-50`, `main.py:47-48`).
- **FP reproducibility:** `charge_*`, `chargeShift_*`, `pi_*` touch the `10^x` transcendental and are rounded to 3 decimals at the output boundary (`pipeline.py:55-63`); every other property is bit-exact integer/constant arithmetic given a **fixed summation order**.

The governing rule for this plan: **parallelism across rows is safe** (each row is an independent, pure computation; results are reassembled in input order; the final sort restores byte order). **Parallelism inside a single float reduction is not** (it reorders a sum). Phase B parallelizes only across rows — it never changes how any single property is summed. The golden master in Phase A is what proves this held.

## File Structure

| Path | Responsibility | Touched in |
|---|---|---|
| `software/tests/regression/__init__.py` | New regression-suite package marker | Task 1 |
| `software/tests/regression/golden_cases.py` | Single source of truth for golden input rows + plans (DRY: shared by regen + test) | Task 1 |
| `software/tests/regression/regen_golden.py` | Deliberate regeneration of committed golden outputs from current code | Task 1 |
| `software/tests/regression/test_golden_master.py` | Asserts full-output bytes (sha256) match committed goldens | Task 1 |
| `software/tests/data/golden/<case>/{properties.tsv,aa_fraction.tsv,stats.json}` | Committed frozen outputs (generated artifacts) | Task 1 |
| `software/tests/regression/test_scale_invariants.py` | Structural invariants at scale: count/keys preserved, byte-stable | Task 2 |
| `software/src/pipeline.py` | Add worker pool + `resolve_workers`; thread `workers` through `run` | Task 3 |
| `software/tests/unit/test_parallel_invariance.py` | Worker-count invariance: `workers=1` ≡ `workers=N`, byte-identical | Task 3 |
| `software/src/main.py` | Resolve worker count from env (argv unchanged), pass to `run`, log it | Task 4 |
| `workflow/src/main.tpl.tengo` | Request more cores/mem on the Python exec step | Task 5 |
| `software/.changeset/*.md`, `workflow/.changeset` (root `.changeset/`) | Version bump for software + workflow | Task 5 |
| `software/src/io_layer.py`, `software/src/pipeline.py` | (OPTIONAL) streaming sorted-batch I/O | Task 6 |

---

## Phase A — Characterization (do first)

### Task 1: Freeze a full-output golden master

**Goal:** A committed, byte-exact snapshot of the full CLI output (`properties.tsv`, `aa_fraction.tsv`, `stats.json`) for every mode + coverage tier + edge case, generated from the *current* (pre-refactor) code, asserted via sha256. This is the safety net for Phases B and C.

**Files:**
- Create: `software/tests/regression/__init__.py`
- Create: `software/tests/regression/golden_cases.py`
- Create: `software/tests/regression/regen_golden.py`
- Create: `software/tests/regression/test_golden_master.py`
- Create (generated): `software/tests/data/golden/<case>/{properties.tsv,aa_fraction.tsv,stats.json}`

**Acceptance Criteria:**
- [ ] Cases cover: peptide (valid, <10aa→II NA, no-aromatic ε=0, paired-Cys ε, non-standard residue excluded, all-non-standard→all NA, stop-codon→all NA, empty→all NA), antibody full+Fv (incl. a clone with one region missing → full-chain NA for that clone, and a clone with empty CDR3), antibody CDR3-only, antibody partial (`fullChains:["A"]`), TCR (no Fv columns).
- [ ] `test_golden_master.py` runs `main()` per case into a tmp dir and asserts sha256 of all three output files equals the committed golden.
- [ ] A non-vacuous guard: at least one case asserts a known column set is present (proves the snapshot isn't empty).
- [ ] Goldens are generated by `regen_golden.py` from current code and committed.

**Verify:** `uv run pytest tests/regression/test_golden_master.py -v` → all cases PASS.

**Steps:**

- [ ] **Step 1: Create the regression package marker**

Create `software/tests/regression/__init__.py` (empty file).

- [ ] **Step 2: Define the golden cases (single source of truth)**

Create `software/tests/regression/golden_cases.py`:

```python
"""Golden-master case definitions — the single source of truth for both the
regeneration script and the byte-compare test. Each case is (input rows,
column order, plan dict). Inputs live in code; only the *outputs* are
committed under tests/data/golden/<name>/.

Edge cases are drawn from the spec's documented behaviour tables (NA
propagation, non-standard residues, stop codon, per-clone missing region,
TCR has no Fv).
"""

from __future__ import annotations

_AB_COLS = (
    ["entity_key"]
    + [f"A_{f}" for f in ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")]
    + [f"B_{f}" for f in ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")]
)

_AB_FULL_ROW = {
    "entity_key": "c1",
    "A_FR1": "EVQLVES", "A_CDR1": "GFTFSSY", "A_FR2": "AMSWVRQ", "A_CDR2": "ISGSGGS",
    "A_FR3": "TYYAESVKGRFTI", "A_CDR3": "CARDYW", "A_FR4": "WGQGTLV",
    "B_FR1": "DIQMTQS", "B_CDR1": "QSISSY", "B_FR2": "LNWYQQK", "B_CDR2": "AASSLQS",
    "B_FR3": "GVPSRFSGSG", "B_CDR3": "CQQYNS", "B_FR4": "FGQGTKV",
}
# Clone with one heavy region missing → full-chain A is NA for this clone only.
_AB_MISSING_REGION = {**_AB_FULL_ROW, "entity_key": "c2", "A_FR3": "", "A_CDR3": "CARGFW", "B_CDR3": "CQHFSS"}
# Clone with empty heavy CDR3 → CDR3-A NA for this clone only.
_AB_EMPTY_CDR3 = {**_AB_FULL_ROW, "entity_key": "c3", "A_CDR3": ""}

# CASES: name -> (rows, columns, plan)
CASES: dict[str, tuple[list[dict], list[str], dict]] = {
    "peptide": (
        [
            {"entity_key": "p_valid", "sequence": "ACDEFGHIKL"},
            {"entity_key": "p_basic", "sequence": "KKKKHHHHHH"},
            {"entity_key": "p_acidic", "sequence": "DDDDEEEEEE"},
            {"entity_key": "p_short", "sequence": "RPPGFSPF"},          # 8 aa -> instability NA
            {"entity_key": "p_no_aromatic", "sequence": "AAAAAAAAAA"},  # eox/ered = 0
            {"entity_key": "p_paired_cys", "sequence": "CYIQNCPLG"},    # ε floor(C/2)*125
            {"entity_key": "p_nonstd", "sequence": "ACDXEFGHIK"},       # X excluded
            {"entity_key": "p_all_nonstd", "sequence": "XXXXX"},        # all NA
            {"entity_key": "p_stop", "sequence": "ACDE*GHIK"},          # stop codon -> all NA
            {"entity_key": "p_empty", "sequence": ""},                  # empty -> all NA
        ],
        ["entity_key", "sequence"],
        {"mode": "peptide"},
    ),
    "antibody_full": (
        [_AB_FULL_ROW, _AB_MISSING_REGION, _AB_EMPTY_CDR3],
        _AB_COLS,
        {"mode": "antibody_tcr_legacy_bulk", "receptor": "IG",
         "chains": ["A", "B"], "fullChains": ["A", "B"], "hasFv": True},
    ),
    "antibody_cdr3_only": (
        [
            {"entity_key": "c1", "A_CDR3": "CARDYW", "B_CDR3": "CQQYNS"},
            {"entity_key": "c2", "A_CDR3": "CARGFW", "B_CDR3": "CQHFSS"},
        ],
        ["entity_key", "A_CDR3", "B_CDR3"],
        {"mode": "antibody_tcr_legacy_sc", "receptor": "IG",
         "chains": ["A", "B"], "fullChains": [], "hasFv": False},
    ),
    "antibody_partial": (
        [_AB_FULL_ROW],
        _AB_COLS,
        {"mode": "antibody_tcr_legacy_bulk", "receptor": "IG",
         "chains": ["A", "B"], "fullChains": ["A"], "hasFv": False},
    ),
    "tcr": (
        [_AB_FULL_ROW],
        _AB_COLS,
        {"mode": "antibody_tcr_legacy_bulk", "receptor": "TCRAB",
         "chains": ["A", "B"], "fullChains": ["A", "B"], "hasFv": False},
    ),
}
```

- [ ] **Step 3: Write the regeneration helper**

Create `software/tests/regression/regen_golden.py`:

```python
"""Regenerate committed golden outputs from the CURRENT code.

Run ONLY deliberately:
    uv run python tests/regression/regen_golden.py

DO NOT run this during the parallelism/streaming refactor — if the golden
needs regenerating to make test_golden_master pass, the refactor changed
observable behaviour and that is a bug, not a snapshot update. Regenerate
only on an intentional behaviour change (BioPython bump, pKa update, spec
change) and call it out in the PR.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from main import main  # noqa: E402

from golden_cases import CASES  # noqa: E402

GOLDEN_DIR = Path(__file__).resolve().parents[1] / "data" / "golden"


def _write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join(row.get(c, "") for c in columns))
    path.write_text("\n".join(lines) + "\n")


def regen() -> None:
    for name, (rows, columns, plan) in CASES.items():
        case_dir = GOLDEN_DIR / name
        case_dir.mkdir(parents=True, exist_ok=True)
        in_tsv = case_dir / "input.tsv"
        plan_json = case_dir / "plan.json"
        _write_tsv(in_tsv, rows, columns)
        plan_json.write_text(json.dumps(plan))
        rc = main([
            "--input", str(in_tsv),
            "--plan", str(plan_json),
            "--output", str(case_dir / "properties.tsv"),
            "--aa-fraction", str(case_dir / "aa_fraction.tsv"),
            "--stats", str(case_dir / "stats.json"),
        ])
        assert rc == 0, f"regen failed for case {name}"
        print(f"regenerated golden/{name}")


if __name__ == "__main__":
    regen()
```

- [ ] **Step 4: Write the byte-compare test (red — goldens don't exist yet)**

Create `software/tests/regression/test_golden_master.py`:

```python
"""Full-output byte snapshot — the refactor safety net.

Each case runs main() into a tmp dir and asserts every output file is
byte-identical (sha256) to the committed golden produced by regen_golden.py
from pre-refactor code. If any of Phase B / Phase C changes a byte, this
fails — which is the entire point.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from main import main
from golden_cases import CASES

GOLDEN_DIR = Path(__file__).resolve().parents[1] / "data" / "golden"
_OUTPUTS = ("properties.tsv", "aa_fraction.tsv", "stats.json")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join(row.get(c, "") for c in columns))
    path.write_text("\n".join(lines) + "\n")


@pytest.mark.parametrize("name", sorted(CASES), ids=sorted(CASES))
def test_output_matches_golden(tmp_path: Path, name: str):
    rows, columns, plan = CASES[name]
    in_tsv = tmp_path / "input.tsv"
    plan_json = tmp_path / "plan.json"
    _write_tsv(in_tsv, rows, columns)
    plan_json.write_text(json.dumps(plan))

    out = {n: tmp_path / n for n in _OUTPUTS}
    rc = main([
        "--input", str(in_tsv),
        "--plan", str(plan_json),
        "--output", str(out["properties.tsv"]),
        "--aa-fraction", str(out["aa_fraction.tsv"]),
        "--stats", str(out["stats.json"]),
    ])
    assert rc == 0

    golden = GOLDEN_DIR / name
    assert golden.is_dir(), f"missing golden dir for {name}; run regen_golden.py"
    for n in _OUTPUTS:
        assert _sha256(out[n]) == _sha256(golden / n), (
            f"{name}/{n} diverged from golden — refactor changed behaviour"
        )


# Non-vacuous guard: the antibody_full golden actually carries the Fv columns,
# and the tcr golden does NOT (R12). Proves the snapshot is meaningful, not empty.
def test_golden_column_presence():
    import polars as pl

    ab = pl.read_csv(GOLDEN_DIR / "antibody_full" / "properties.tsv", separator="\t")
    assert {"charge_Fv", "pi_Fv", "charge_A_VDJRegion"} <= set(ab.columns)
    tcr = pl.read_csv(GOLDEN_DIR / "tcr" / "properties.tsv", separator="\t")
    assert not any(c.endswith("_Fv") for c in tcr.columns)
```

- [ ] **Step 5: Run the test — confirm it fails for the right reason**

Run: `uv run pytest tests/regression/test_golden_master.py -v`
Expected: FAIL — `missing golden dir for <name>; run regen_golden.py`.

- [ ] **Step 6: Generate goldens from current code, then go green**

```bash
uv run python tests/regression/regen_golden.py
uv run pytest tests/regression/test_golden_master.py -v
```
Expected: regen prints one line per case; pytest PASSES all cases + `test_golden_column_presence`.

- [ ] **Step 7: Lint, then commit (goldens included)**

```bash
uv run ruff format
uv run ruff check
git add tests/regression/ tests/data/golden/
git commit -m "test: freeze full-output golden master for sequence-properties refactor"
```

---

### Task 2: Structural invariants at scale

**Goal:** Pin the behavioral invariants that the parallel/streaming refactor must preserve — output has exactly one row per input entity (no loss/duplication), the same set of keys, the AA-fraction frame has 20 rows per peptide, and two in-process runs are byte-identical — at a size (5000 rows) that will span multiple worker chunks in Phase B.

**Files:**
- Create: `software/tests/regression/test_scale_invariants.py`

**Acceptance Criteria:**
- [ ] Deterministic synthetic input (seeded `random.Random(0)` — no unseeded RNG, no clock).
- [ ] Peptide: output row count == input count; output `entity_key` set == input set; AA-fraction height == 20 × count.
- [ ] Antibody full: output row count == input count; key set preserved; `charge_Fv` column present.
- [ ] Two runs on the same input produce byte-identical `properties.tsv`.
- [ ] Marked `@pytest.mark.slow` so it deselects during fast iteration.

**Verify:** `uv run pytest tests/regression/test_scale_invariants.py -v` → PASS (and `-m "not slow"` deselects it).

**Steps:**

- [ ] **Step 1: Write the invariant test (passes on current sequential code)**

Create `software/tests/regression/test_scale_invariants.py`:

```python
"""Structural invariants at scale — the per-row-map contract the parallel and
streaming refactors must preserve: one output row per input entity, same key
set, AA-fraction 20 rows/peptide, byte-stable across runs.

Synthetic data is built from a SEEDED RNG (random.Random(0)) so the suite is
deterministic. Marked slow: 5000 rows is enough to exercise >1 worker chunk
once Phase B lands, but too slow for every iteration.
"""

from __future__ import annotations

import hashlib
import random
from pathlib import Path

import polars as pl
import pytest

from io_layer import write_output_tsv
from pipeline import run

_AAS = "ACDEFGHIKLMNPQRSTVWY"
_REGIONS = ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")


def _rand_seq(rng: random.Random, lo: int, hi: int) -> str:
    return "".join(rng.choice(_AAS) for _ in range(rng.randint(lo, hi)))


def _peptide_frame(n: int) -> pl.DataFrame:
    rng = random.Random(0)
    return pl.DataFrame(
        {"entity_key": [f"p{i}" for i in range(n)],
         "sequence": [_rand_seq(rng, 5, 25) for _ in range(n)]},
        schema={"entity_key": pl.Utf8, "sequence": pl.Utf8},
    )


def _antibody_frame(n: int) -> pl.DataFrame:
    rng = random.Random(1)
    cols: dict[str, list[str]] = {"entity_key": [f"c{i}" for i in range(n)]}
    for ch in ("A", "B"):
        for feat in _REGIONS:
            cols[f"{ch}_{feat}"] = [_rand_seq(rng, 6, 14) for _ in range(n)]
    schema = {k: pl.Utf8 for k in cols}
    return pl.DataFrame(cols, schema=schema)


_AB_PLAN = {"mode": "antibody_tcr_legacy_bulk", "receptor": "IG",
            "chains": ["A", "B"], "fullChains": ["A", "B"], "hasFv": True}

N = 5000


@pytest.mark.slow
def test_peptide_row_and_key_preservation():
    reads = _peptide_frame(N)
    out = run(reads, {"mode": "peptide"})
    props = out["properties"]
    assert props.height == N
    assert set(props["entity_key"].to_list()) == set(reads["entity_key"].to_list())
    assert out["aa_fraction"].height == 20 * N


@pytest.mark.slow
def test_antibody_row_and_key_preservation():
    reads = _antibody_frame(N)
    out = run(reads, _AB_PLAN)
    props = out["properties"]
    assert props.height == N
    assert set(props["entity_key"].to_list()) == set(reads["entity_key"].to_list())
    assert "charge_Fv" in props.columns


@pytest.mark.slow
def test_peptide_byte_stable_two_runs(tmp_path: Path):
    reads = _peptide_frame(N)
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    write_output_tsv(run(reads, {"mode": "peptide"})["properties"], a, sort_keys=["entity_key"])
    write_output_tsv(run(reads, {"mode": "peptide"})["properties"], b, sort_keys=["entity_key"])
    assert hashlib.sha256(a.read_bytes()).hexdigest() == hashlib.sha256(b.read_bytes()).hexdigest()
```

- [ ] **Step 2: Run it (and confirm slow-deselect works)**

```bash
uv run pytest tests/regression/test_scale_invariants.py -v
uv run pytest tests/regression/test_scale_invariants.py -m "not slow"   # collects 0
```
Expected: first PASSES (3 tests); second deselects all 3.

- [ ] **Step 3: Lint and commit**

```bash
uv run ruff format
uv run ruff check
git add tests/regression/test_scale_invariants.py
git commit -m "test: pin row/key-preservation + byte-stability invariants at scale"
```

---

## Phase B — Option #1: Parallelize across cores

### Task 3: Map the per-row compute across a process pool

**Goal:** Run the existing pure per-row compute functions across a `ProcessPoolExecutor`, assembling results in input order, so the output is byte-identical regardless of worker count. Keep a sequential fallback for small inputs and `workers=1`.

**Files:**
- Modify: `software/src/pipeline.py` (add imports, `resolve_workers`, `_pmap`, worker functions; thread `workers` through `run` / `run_peptide` / `run_antibody_tcr` — `pipeline.py:71-91`, `159-221`, `414-439`)
- Test: `software/tests/unit/test_parallel_invariance.py` (create)

**Acceptance Criteria:**
- [ ] `run(reads, plan, workers=None)` accepts a worker count; `run_peptide` / `run_antibody_tcr` accept and use it.
- [ ] `workers=1` and `workers=4` produce byte-identical `properties.tsv`, `aa_fraction.tsv`, and identical `stats` for peptide, antibody-full, and the 5000-row scale frame.
- [ ] Worker functions are module-level and picklable (spawn-safe): they take only strings / plain dicts.
- [ ] Results are assembled in input order (no reliance on completion order, no `set`/`dict`-order leakage into output).
- [ ] The Phase A golden master and scale invariants still pass unchanged.

**Verify:** `uv run pytest tests/unit/test_parallel_invariance.py tests/regression/ -v` → PASS.

**Steps:**

- [ ] **Step 1: Write the invariance test first (red)**

Create `software/tests/unit/test_parallel_invariance.py`:

```python
"""Worker-count invariance — the core safety property for parallelism.

The output of run() must not depend on how many workers compute it. We
compare workers=1 against workers=4 (and a larger frame that actually spills
into the pool) and assert byte-identical serialized output + identical stats.
"""

from __future__ import annotations

import random
from pathlib import Path

import polars as pl
import pytest

from io_layer import write_output_tsv
from pipeline import run

_AAS = "ACDEFGHIKLMNPQRSTVWY"
_REGIONS = ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")
_AB_PLAN = {"mode": "antibody_tcr_legacy_bulk", "receptor": "IG",
            "chains": ["A", "B"], "fullChains": ["A", "B"], "hasFv": True}


def _peptide_frame(n: int, seed: int) -> pl.DataFrame:
    rng = random.Random(seed)
    seqs = ["".join(rng.choice(_AAS) for _ in range(rng.randint(5, 25))) for _ in range(n)]
    return pl.DataFrame({"entity_key": [f"p{i}" for i in range(n)], "sequence": seqs},
                        schema={"entity_key": pl.Utf8, "sequence": pl.Utf8})


def _antibody_frame(n: int, seed: int) -> pl.DataFrame:
    rng = random.Random(seed)
    cols = {"entity_key": [f"c{i}" for i in range(n)]}
    for ch in ("A", "B"):
        for feat in _REGIONS:
            cols[f"{ch}_{feat}"] = ["".join(rng.choice(_AAS) for _ in range(rng.randint(6, 14))) for _ in range(n)]
    return pl.DataFrame(cols, schema={k: pl.Utf8 for k in cols})


def _serialize(out: dict, tmp: Path, tag: str) -> tuple[bytes, bytes, dict]:
    p = tmp / f"props_{tag}.tsv"
    a = tmp / f"aa_{tag}.tsv"
    write_output_tsv(out["properties"], p, sort_keys=["entity_key"])
    write_output_tsv(out["aa_fraction"], a, sort_keys=["entity_key", "aminoAcid"])
    return p.read_bytes(), a.read_bytes(), out["stats"]


@pytest.mark.parametrize("n", [3, 3000], ids=["sequential-path", "pool-path"])
def test_peptide_invariant_to_workers(tmp_path: Path, n: int):
    reads = _peptide_frame(n, seed=0)
    one = _serialize(run(reads, {"mode": "peptide"}, workers=1), tmp_path, "1")
    four = _serialize(run(reads, {"mode": "peptide"}, workers=4), tmp_path, "4")
    assert one == four


@pytest.mark.parametrize("n", [3, 3000], ids=["sequential-path", "pool-path"])
def test_antibody_invariant_to_workers(tmp_path: Path, n: int):
    reads = _antibody_frame(n, seed=1)
    one = _serialize(run(reads, _AB_PLAN, workers=1), tmp_path, "1")
    four = _serialize(run(reads, _AB_PLAN, workers=4), tmp_path, "4")
    assert one == four
```

- [ ] **Step 2: Run it — confirm it fails (run() has no `workers` param yet)**

Run: `uv run pytest tests/unit/test_parallel_invariance.py -v`
Expected: FAIL — `TypeError: run() got an unexpected keyword argument 'workers'`.

- [ ] **Step 3: Add the pool infrastructure to `pipeline.py`**

At the top of `software/src/pipeline.py`, add to the imports block (after `import logging`):

```python
import functools
import os
from concurrent.futures import ProcessPoolExecutor
```

Then add, just below the `PH = 7.0` line (`pipeline.py:34`):

```python
# Parallelism. Below this row count the pool's process startup + pickle cost
# outweighs the benefit, so we stay in-process. The threshold also keeps the
# whole existing unit-test suite on the fast sequential path.
_PARALLEL_MIN_ROWS = 2000


def resolve_workers(workers: int | None) -> int:
    """How many worker processes to use. Explicit arg wins (used by tests and
    by main.py once it reads the platform's CPU allocation). Falls back to the
    PL_COMPUTE_WORKERS env var, then os.cpu_count(). The RESULT never depends on
    this number — only the wall-clock does — so an over- or under-estimate is a
    speed concern, never a correctness one.
    """
    if workers is not None:
        return max(1, int(workers))
    env = os.environ.get("PL_COMPUTE_WORKERS")
    if env and env.isdigit() and int(env) > 0:
        return int(env)
    return max(1, os.cpu_count() or 1)


def _pmap(fn, items: list, workers: int, chunksize: int = 256) -> list:
    """Map fn over items, preserving input order. Sequential below the
    threshold or when workers<=1; otherwise a process pool. ProcessPoolExecutor
    .map() preserves input order, so results align with items by index — the
    property the byte-stable output depends on.
    """
    if workers <= 1 or len(items) < _PARALLEL_MIN_ROWS:
        return [fn(x) for x in items]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        return list(ex.map(fn, items, chunksize=chunksize))
```

- [ ] **Step 4: Add module-level, picklable worker functions**

In `software/src/pipeline.py`, add a peptide worker just above `run_peptide` (`pipeline.py:159`):

```python
def _peptide_worker(seq: str) -> tuple[dict[str, float | None], list[float | None] | None]:
    """Picklable per-peptide unit: (properties row, 20 AA fractions in
    STANDARD_AAS order | None). One SequenceContext per sequence, shared across
    all 11 reads — same sharing the old inline loop relied on.
    """
    ctx = SequenceContext.from_seq(seq)
    if ctx is None:
        return (dict(_NA_PEPTIDE_ROW), None)
    props = _compute_peptide_row_from_ctx(ctx)
    fr = ctx.aa_fractions()
    return (props, [fr[aa] for aa in STANDARD_AAS])
```

And an antibody worker just above `run_antibody_tcr` (`pipeline.py:414`):

```python
def _antibody_worker(record: dict, plan: dict[str, Any]) -> dict[str, Any]:
    """Picklable per-clone unit. `plan` is bound per-call via functools.partial;
    it is a plain dict and pickles cleanly. Delegates to the existing
    _compute_row_for so the computation is unchanged.
    """
    return _compute_row_for(record, plan)
```

- [ ] **Step 5: Rewrite `run` to thread `workers` through**

Replace `run` (`pipeline.py:71-91`) with:

```python
def run(reads: pl.DataFrame, plan: dict[str, Any], workers: int | None = None) -> dict[str, Any]:
    """Dispatch by mode. `workers` controls parallelism only — output is
    identical for any value (see test_parallel_invariance). Returns a dict with
    `properties`, `aa_fraction`, and `stats` entries (unchanged contract).
    """
    n_workers = resolve_workers(workers)
    mode = plan["mode"]
    if mode == "peptide":
        log.info("Running peptide mode (%d entities, %d workers)", reads.height, n_workers)
        return run_peptide(reads, n_workers)
    log.info(
        "Running antibody/TCR mode (receptor=%s, %d clones, %d workers)",
        plan.get("receptor", "IG"), reads.height, n_workers,
    )
    return run_antibody_tcr(reads, plan, n_workers)
```

- [ ] **Step 6: Rewrite `run_peptide`'s loop to use the pool**

Replace the body of `run_peptide` (`pipeline.py:159-221`) — keep the signature change and the DataFrame/stats construction identical, only the per-row work moves into `_pmap`:

```python
def run_peptide(reads: pl.DataFrame, workers: int = 1) -> dict[str, Any]:
    """Compute peptide-mode outputs. Per-sequence work runs through `_pmap`
    (sequential or pooled); results are reassembled in input order so the
    serialized output is byte-identical regardless of worker count.
    """
    keys = reads["entity_key"].to_list()
    seqs = reads["sequence"].to_list()
    n = len(seqs)

    log.info("Computing peptide properties + AA fractions (%d sequences)", n)
    results = _pmap(_peptide_worker, seqs, workers)

    prop_cols: dict[str, list[Any]] = {"entity_key": list(keys), **{c: [] for c in PEPTIDE_PROPERTY_COLUMNS}}
    aa_entity: list[str] = []
    aa_amino: list[str] = []
    aa_value: list[float | None] = []
    for k, (props, fractions) in zip(keys, results):
        for c in PEPTIDE_PROPERTY_COLUMNS:
            prop_cols[c].append(props[c])
        if fractions is None:
            for aa in STANDARD_AAS:
                aa_entity.append(k)
                aa_amino.append(aa)
                aa_value.append(None)
        else:
            for aa, val in zip(STANDARD_AAS, fractions):
                aa_entity.append(k)
                aa_amino.append(aa)
                aa_value.append(val)
    properties = pl.DataFrame(
        prop_cols,
        schema={"entity_key": pl.Utf8, **{c: pl.Float64 for c in PEPTIDE_PROPERTY_COLUMNS}},
    )
    aa_fraction = pl.DataFrame(
        {"entity_key": aa_entity, "aminoAcid": aa_amino, "value": aa_value},
        schema={"entity_key": pl.Utf8, "aminoAcid": pl.Utf8, "value": pl.Float64},
    )

    has_below_floor = any(0 < effective_length(s) < INSTABILITY_MIN_LENGTH for s in seqs if s)
    stats = {
        "medianCdr3Length": {},
        "hasPeptideBelowInstabilityFloor": has_below_floor,
    }

    return {
        "properties": _quantize_for_cid(properties),
        "aa_fraction": aa_fraction,
        "stats": stats,
    }
```

- [ ] **Step 7: Rewrite `run_antibody_tcr`'s loop to use the pool**

Replace the body of `run_antibody_tcr` (`pipeline.py:414-439`):

```python
def run_antibody_tcr(reads: pl.DataFrame, plan: dict[str, Any], workers: int = 1) -> dict[str, Any]:
    chains = plan.get("chains", [])
    full_chains = plan.get("fullChains", [])
    n = reads.height
    if chains:
        log.info("Computing CDR3 properties for chains %s (%d clones)", list(chains), n)
    if full_chains:
        log.info("Reconstructing full chains %s and computing full-chain properties", list(full_chains))
    if plan.get("hasFv"):
        log.info("Computing Fv properties (paired VH+VL)")

    out_cols = _planned_output_columns(plan)
    records = list(reads.iter_rows(named=True))
    worker = functools.partial(_antibody_worker, plan=plan)
    rows = _pmap(worker, records, workers)

    columns: dict[str, list[Any]] = {"entity_key": [], **{c: [] for c in out_cols}}
    for row in rows:
        columns["entity_key"].append(row["entity_key"])
        for c in out_cols:
            columns[c].append(row[c])
    schema = {"entity_key": pl.Utf8, **{c: pl.Float64 for c in out_cols}}
    properties = pl.DataFrame(columns, schema=schema)
    aa_fraction = pl.DataFrame(schema={"entity_key": pl.Utf8, "aminoAcid": pl.Utf8, "value": pl.Float64})
    stats = {"medianCdr3Length": _median_cdr3_length_by_chain(reads, chains)}
    return {
        "properties": _quantize_for_cid(properties),
        "aa_fraction": aa_fraction,
        "stats": stats,
    }
```

- [ ] **Step 8: Run the invariance test + the full safety net**

```bash
uv run pytest tests/unit/test_parallel_invariance.py -v
uv run pytest tests/regression/ -v
uv run pytest
```
Expected: invariance PASSES (both sequential-path and pool-path ids); golden master + scale invariants PASS unchanged; full suite green. If the golden master fails, **do not regenerate it** — the parallel path changed a byte; debug ordering/assembly until it matches.

- [ ] **Step 9: Lint and commit**

```bash
uv run ruff format
uv run ruff check
git add src/pipeline.py tests/unit/test_parallel_invariance.py
git commit -m "feat: parallelize per-row property compute across a process pool"
```

---

### Task 4: Resolve worker count from the environment in `main.py`

**Goal:** Make the CLI use all allocated cores at runtime without adding anything to argv — so the workflow exec node's command line stays constant and its content-addressed identity does not change.

**Files:**
- Modify: `software/src/main.py` (`main.py:18`, `main.py:44-46`)
- Test: `software/tests/integration/test_cli.py` (add one test)

**Acceptance Criteria:**
- [ ] `main()` resolves workers via `pipeline.resolve_workers(None)` (env → `os.cpu_count()`), passes to `run`, and logs the count to stderr.
- [ ] **No new CLI argument** — argv is unchanged from today (preserves the exec CID).
- [ ] Output bytes are identical with `PL_COMPUTE_WORKERS=1` vs unset (env changes speed, not output).

**Verify:** `uv run pytest tests/integration/test_cli.py -v` → PASS.

**Steps:**

- [ ] **Step 1: Add the env-invariance test (red)**

Append to `software/tests/integration/test_cli.py`:

```python
# Worker count comes from the environment, not argv, and must not change output
# bytes — only wall-clock. Guards the CID-stability property at the CLI layer.
def test_env_worker_count_does_not_change_output(tmp_path: Path, monkeypatch):
    rows = [
        {"entity_key": "p1", "sequence": "ACDEFGHIKL"},
        {"entity_key": "p2", "sequence": "MNPQRSTVWY"},
    ]
    monkeypatch.setenv("PL_COMPUTE_WORKERS", "1")
    out_a, _ = _run_peptide(tmp_path, "_w1", rows)
    monkeypatch.setenv("PL_COMPUTE_WORKERS", "4")
    out_b, _ = _run_peptide(tmp_path, "_w4", rows)
    assert _sha256(out_a) == _sha256(out_b)
```

- [ ] **Step 2: Run it — confirm it passes already (worker count not yet wired, but defaults are equivalent on tiny input)**

Run: `uv run pytest tests/integration/test_cli.py::test_env_worker_count_does_not_change_output -v`
Expected: PASS (tiny input stays sequential either way). This test locks the contract; the wiring below makes it meaningful at scale.

- [ ] **Step 3: Wire `resolve_workers` into `main.py`**

In `software/src/main.py`, change the import (`main.py:19`):

```python
from pipeline import resolve_workers, run
```

Then in `main()`, replace the `outputs = run(reads, plan)` line (`main.py:46`) with:

```python
    workers = resolve_workers(None)
    logging.getLogger(__name__).info("Using %d compute worker(s)", workers)
    outputs = run(reads, plan, workers=workers)
```

- [ ] **Step 4: Run the CLI suite + full suite**

```bash
uv run pytest tests/integration/test_cli.py -v
uv run pytest
```
Expected: all PASS, including the determinism subprocess test (`test_determinism.py`) which now also exercises env-resolved workers in a fresh process.

- [ ] **Step 5: Lint and commit**

```bash
uv run ruff format
uv run ruff check
git add src/main.py tests/integration/test_cli.py
git commit -m "feat: resolve compute worker count from environment in CLI"
```

---

### Task 5: Request more cores + memory on the workflow exec step

**Goal:** Tell the platform to schedule the Python step with multiple cores and more memory, so the parallel kernel actually has cores to use and a large input doesn't OOM. Keep the exec node's identity stable.

**Files:**
- Modify: `workflow/src/main.tpl.tengo` (`main.tpl.tengo:344-360`)
- Create: `.changeset/<name>.md` (block root)
- Verify against SDK: `core/platforma/sdk/workflow-tengo/src/exec/index.lib.tengo`

**Acceptance Criteria:**
- [ ] The `pyRun` exec builder requests multiple cores (`.cpu(8)`) and more memory (`.mem("8GiB")`).
- [ ] **No `--workers` (or any machine-dependent value) added to argv** — Python reads cores from the environment at runtime.
- [ ] Confirmed (by reading `exec/index.lib.tengo`) whether `cpu()`/`mem()` enter the exec node's content-addressed id. Document the finding in the PR description: if they do, the bump is a one-time cache invalidation and the requested values are fixed constants in code, so the CID stays stable run-to-run.
- [ ] `pnpm run build:dev` succeeds.
- [ ] A changeset covers `.workflow` and `.software`.

**Verify:** `pnpm run build:dev` (from block root) → build succeeds; `grep -n "cpu(" workflow/src/main.tpl.tengo` shows the new value on the `pyRun` step.

**Steps:**

- [ ] **Step 1: Verify whether cpu()/mem() are part of the exec CID**

Read `core/platforma/sdk/workflow-tengo/src/exec/index.lib.tengo` and confirm whether `.cpu()` / `.mem()` feed the resource's content-addressed identity or are scheduling-only hints. Record the answer in your working notes — it determines whether bumping them invalidates existing cached results (a one-time, acceptable cost) and confirms the CID stays stable across machines (the requested values are fixed code constants, not machine-derived).

- [ ] **Step 2: Bump cpu + mem on the Python step**

In `workflow/src/main.tpl.tengo`, in the `pyRun := exec.builder()...` chain (`main.tpl.tengo:345-360`), change:

```tengo
		software(soft).
		mem("4GiB").
		cpu(1).
```
to:
```tengo
		software(soft).
		mem("8GiB").
		cpu(8).
```

Leave the `arg("--input")...arg("--stats")` chain exactly as-is — **do not add a `--workers` arg** (argv must stay constant for CID stability; the Python step reads its worker count from the environment via `resolve_workers`).

- [ ] **Step 3: Add a changeset**

Create `.changeset/seqprops-parallel-compute.md` (from block root):

```markdown
---
'@platforma-open/milaboratories.sequence-properties.software': patch
'@platforma-open/milaboratories.sequence-properties.workflow': patch
---

Parallelize the property computation across CPU cores and raise the compute step's CPU/memory request, so large datasets compute on the first run without the previous single-core bottleneck. Output is byte-identical to the prior version (worker count never affects results).
```

Do not hand-edit any `package.json` `version` field — changesets + CI handle the bump.

- [ ] **Step 4: Build the block (dev)**

From the block root `blocks/sequence-properties/`:
```bash
pnpm install
pnpm run build:dev
```
Expected: turbo builds workflow + model + ui + software; build succeeds. If `pnpm install` modified `pnpm-lock.yaml`, stage it.

- [ ] **Step 5: Commit**

```bash
git add workflow/src/main.tpl.tengo .changeset/seqprops-parallel-compute.md
git add -A pnpm-lock.yaml 2>/dev/null || true
git commit -m "feat: request 8 cores + 8GiB for the sequence-properties compute step"
```

- [ ] **Step 6: (Operator / optional) live verification**

The block-level workflow test (`workflow/src/wf.test.ts`) needs a live backend (`run-platforma` skill). If a backend is available, run the block test suite per the `block-dev` skill (`pnpm test --filter=...`). Otherwise the Python golden master + invariance tests are the primary correctness gate; flag to the operator that an end-to-end run on a representative large dataset should be done before marking the PR ready.

---

## Phase C — Option #7: Streaming I/O (OPTIONAL — decide at end)

> **OPTIONAL.** Do not implement until Phase B is merged and measured. Decision gate at Step 0 below. Phase B's `mem("8GiB")` bump may already make this unnecessary — streaming is only worth its complexity if a target dataset still risks OOM at the raised memory ceiling.

### Task 6 (OPTIONAL): Stream sorted input in batches to bound peak memory

**Goal:** Process the input in sorted batches and write output incrementally, so peak memory is O(batch) instead of O(dataset). Output stays byte-identical (the golden master is the guard).

**Files:**
- Modify: `software/src/io_layer.py`, `software/src/pipeline.py`
- Test: `software/tests/regression/test_streaming_equivalence.py` (create)

**Acceptance Criteria:**
- [ ] Decision recorded (see Step 0): implement or skip, with the measurement that justified it.
- [ ] If implemented: a streaming path that (a) sorts input by `entity_key` once up front so output needs no global re-sort, (b) computes + writes in batches of `B` rows (header on first batch only), (c) computes `stats` (median CDR3) without materializing the whole dataset.
- [ ] If implemented: `test_streaming_equivalence.py` asserts the streaming output is byte-identical to the non-streaming output on the 5000-row frame; the Phase A golden master still passes.

**Steps:**

- [ ] **Step 0: Decision gate (do this first, record the outcome in the PR)**

Measure peak RSS of the Phase B build on the largest realistic input (or a synthetic frame at the agreed upper bound — confirm the target row count with the operator: 10⁶? 10⁷?). If peak memory stays comfortably under the `8GiB` request, **skip this task** and record "streaming not needed: peak RSS = X GiB at N rows < 8 GiB". Only proceed if it approaches or exceeds the ceiling.

- [ ] **Step 1: Write the equivalence test (red)**

Create `software/tests/regression/test_streaming_equivalence.py`:

```python
"""Streaming must not change a byte. Compares the streaming CLI path against
the in-memory path on a 5000-row frame, and relies on the golden master for
the edge-case modes.
"""

from __future__ import annotations

import hashlib
import json
import random
from pathlib import Path

import pytest

from main import main

_AAS = "ACDEFGHIKLMNPQRSTVWY"


def _write_peptides(path: Path, n: int) -> None:
    rng = random.Random(0)
    lines = ["entity_key\tsequence"]
    for i in range(n):
        seq = "".join(rng.choice(_AAS) for _ in range(rng.randint(5, 25)))
        lines.append(f"p{i}\t{seq}")
    path.write_text("\n".join(lines) + "\n")


def _run(tmp: Path, tag: str, in_tsv: Path, env_streaming: str, monkeypatch) -> str:
    monkeypatch.setenv("PL_STREAM_BATCH", env_streaming)  # "0" disables, ">0" sets batch size
    plan = tmp / f"plan_{tag}.json"
    plan.write_text(json.dumps({"mode": "peptide"}))
    out = tmp / f"out_{tag}.tsv"
    rc = main(["--input", str(in_tsv), "--plan", str(plan), "--output", str(out),
               "--aa-fraction", str(tmp / f"aa_{tag}.tsv"), "--stats", str(tmp / f"stats_{tag}.json")])
    assert rc == 0
    return hashlib.sha256(out.read_bytes()).hexdigest()


@pytest.mark.slow
def test_streaming_matches_in_memory(tmp_path: Path, monkeypatch):
    in_tsv = tmp_path / "input.tsv"
    _write_peptides(in_tsv, 5000)
    in_memory = _run(tmp_path, "mem", in_tsv, "0", monkeypatch)
    streamed = _run(tmp_path, "stream", in_tsv, "1000", monkeypatch)
    assert in_memory == streamed
```

- [ ] **Step 2: Implement the streaming path**

In `software/src/pipeline.py`, add a streaming entry that the CLI selects when `PL_STREAM_BATCH` > 0:
- Sort `reads` by `entity_key` once (`reads.sort("entity_key")`).
- Iterate in slices of `B` rows; for each slice call the existing `run_peptide` / `run_antibody_tcr` (which already parallelize via Phase B), apply `_quantize_for_cid`, and append to the output TSV — header only on the first slice (`write_csv(..., include_header=(first_slice))`, append mode).
- Compute `stats` separately: `hasPeptideBelowInstabilityFloor` is a boolean OR across slices; `medianCdr3Length` is computed once over just the CDR3 columns (lengths are small ints — a single `effective_length` pass, no full compute).
- Because the input is pre-sorted, the appended output is already in `entity_key` order — no global re-sort, no full-dataset materialization.

In `software/src/main.py`, read `int(os.environ.get("PL_STREAM_BATCH", "0"))` and dispatch to the streaming entry when > 0; otherwise the in-memory `run`. (Keep argv unchanged — streaming is env-gated, like worker count.)

- [ ] **Step 3: Verify byte-equivalence + golden master**

```bash
uv run pytest tests/regression/test_streaming_equivalence.py -v
uv run pytest tests/regression/test_golden_master.py -v
uv run pytest
```
Expected: streaming output == in-memory output; golden master unchanged; full suite green.

- [ ] **Step 4: Lint, changeset, commit**

```bash
uv run ruff format
uv run ruff check
```
Add a `patch` changeset entry for `.software`, then:
```bash
git add src/pipeline.py src/main.py tests/regression/test_streaming_equivalence.py .changeset/
git commit -m "feat: optional streaming compute path to bound peak memory"
```

---

## Self-Review

- **Spec coverage / scope:** This plan does not change any property formula, column, axis, domain, or annotation in the `sequence-properties` spec — it is a pure performance change behind a behavior-freeze. The determinism contract from `docs/text/work/projects/sequence-properties` (byte-stable output, CID quantization) is preserved by construction and verified by Task 1.
- **Characterization-first:** Phase A (Tasks 1–2) is committed before any compute change; the golden master is generated from pre-refactor code and must not be regenerated during Phases B/C.
- **Type/name consistency:** `resolve_workers`, `_pmap`, `_peptide_worker`, `_antibody_worker` are defined in Task 3 and consumed by Tasks 3–4; `run(reads, plan, workers=None)`, `run_peptide(reads, workers)`, `run_antibody_tcr(reads, plan, workers)` signatures are consistent across tasks; `PL_COMPUTE_WORKERS` (Task 3/4) and `PL_STREAM_BATCH` (Task 6) are the only new env knobs, both output-invariant.
- **CID stability:** every speed knob (worker count, streaming batch) is read from the environment at runtime, never added to argv — so the exec node's identity is unchanged. Task 5 Step 1 verifies the cpu/mem question against the SDK rather than assuming.
- **Open items to confirm with the operator at execution time:** the Notion ticket ID (branch naming), the realistic upper-bound row count (sizes the Task 6 decision gate), and the exact platform env var for allocated CPUs (Task 4 — `os.cpu_count()` is the safe default if none exists, since the result is worker-count-invariant).
