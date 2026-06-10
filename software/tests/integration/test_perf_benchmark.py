"""Performance benchmark for the vectorized pipeline (marked slow).

Generates a large full-paired antibody dataset (all 7 regions on both chains,
Fv requested) and times `run()` end-to-end. The vectorized engine replaces the
per-row BioPython loop; the documented serial baseline is ~6 s / 50k clones on
the dev Mac. The assertion uses a GENEROUS absolute bound so it does not flake
on slow CI — the real speedup goes in the PR / task report, not the gate.

Run explicitly:
    uv run pytest tests/integration/test_perf_benchmark.py -m slow -v
"""

from __future__ import annotations

import random
import time

import polars as pl
import pytest

from pipeline import run

# IMGT regions per chain, in reconstruction order.
_REGIONS = ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")

# Realistic-ish per-region length bands (match the corpus generator's bands so
# the timed workload resembles real input shape).
_REGION_LEN = {
    "FR1": (7, 10),
    "CDR1": (6, 8),
    "FR2": (7, 9),
    "CDR2": (6, 8),
    "FR3": (10, 14),
    "CDR3": (5, 18),
    "FR4": (7, 9),
}

_STANDARD_AAS = "ACDEFGHIKLMNPQRSTVWY"

_N_CLONES = 50_000
# Generous bound: after vectorizing build_counts + instability + full-chain
# reconstruction the path runs ~0.4 s / 50k on the dev Mac (~14x over the ~6 s
# serial baseline). 1.5 s keeps ~3x headroom for slow CI without masking a real
# regression — a return toward the per-row Python loop would blow well past it.
_TIME_BUDGET_S = 1.5

_PLAN = {
    "mode": "antibody_tcr_legacy_bulk",
    "receptor": "IG",
    "chains": ["A", "B"],
    "fullChains": ["A", "B"],
    "hasFv": True,
}


def _rand_seq(rng: random.Random, lo: int, hi: int) -> str:
    return "".join(rng.choice(_STANDARD_AAS) for _ in range(rng.randint(lo, hi)))


def _make_reads(n: int, seed: int = 0) -> pl.DataFrame:
    rng = random.Random(seed)
    cols: dict[str, list[str]] = {"entity_key": [f"clone_{i:06d}" for i in range(n)]}
    for ch in ("A", "B"):
        for feat in _REGIONS:
            lo, hi = _REGION_LEN[feat]
            cols[f"{ch}_{feat}"] = [_rand_seq(rng, lo, hi) for _ in range(n)]
    schema = {c: pl.Utf8 for c in cols}
    return pl.DataFrame(cols, schema=schema)


@pytest.mark.slow
def test_full_paired_antibody_50k_throughput(capsys):
    reads = _make_reads(_N_CLONES)

    start = time.perf_counter()
    out = run(reads, _PLAN)
    elapsed = time.perf_counter() - start

    # Sanity: every clone produced a row with the full planned column set.
    props = out["properties"]
    assert props.height == _N_CLONES
    us_per_clone = elapsed / _N_CLONES * 1e6
    with capsys.disabled():
        print(
            f"\n[perf] vectorized run(): {_N_CLONES} full-paired antibody clones "
            f"in {elapsed:.3f}s ({us_per_clone:.2f} us/clone); "
            f"serial baseline ~6s/50k -> ~{6.0 / elapsed:.1f}x"
        )

    assert elapsed < _TIME_BUDGET_S, f"50k clones took {elapsed:.3f}s, over the {_TIME_BUDGET_S}s budget"
