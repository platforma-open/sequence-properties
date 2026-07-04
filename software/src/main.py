"""CLI entry point — reads input TSV + plan, dispatches pipeline, writes outputs.

Invoked from the Tengo workflow as:

    python main.py --input input.tsv --plan plan.json
                   --output properties.tsv --aa-fraction aa_fraction.tsv
                   --stats stats.json
"""

from __future__ import annotations

import os

# Thread configuration MUST be set before numpy / polars are imported — both
# read their thread-pool size once, at import time.
#
# BLAS / OpenMP intra-op threads -> 1 (forced). The vectorized engine's only
# order-sensitive step is the `counts @ weights` matvec (vectorized.py): if a
# single dot-product reduction is split across threads, the float summation
# order — and therefore the emitted TSV bytes and the resource CID — becomes
# thread-count dependent. Pinning every BLAS backend to one intra-op thread
# keeps that reduction bit-identical. This is the load-bearing half of the
# determinism contract, so it is forced, not setdefault.
for _thread_var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ[_thread_var] = "1"

# Polars parallelism is safe to scale across cores: every polars op the pipeline
# uses (CSV read, concat_str chain reconstruction, the aa_fraction explode, the
# unique-key output sort, CSV write) is deterministic regardless of thread count
# — none is a cross-row float reduction, and the sorts are total orders on unique
# keys. This value also sizes the numpy pI-bisection worker pool (see
# vectorized._n_workers). The workflow sets POLARS_MAX_THREADS to the cores the
# backend actually grants (the {system.cpu} command expression = the resolved
# cpuFormula value); setdefault to 1 keeps local / test runs single-threaded
# unless they opt in.
os.environ.setdefault("POLARS_MAX_THREADS", "1")

import argparse
import json
import logging
import sys
from pathlib import Path

from io_layer import read_input_tsv, read_plan, write_output_tsv
from pipeline import run


def _configure_logging() -> None:
    # Pipeline milestones go to stderr so the Tengo workflow's stderr stream
    # captures them. force=True lets repeated test invocations re-bind handlers.
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
        force=True,
    )


def main(argv: list[str] | None = None) -> int:
    _configure_logging()
    logging.getLogger(__name__).info(
        "thread config: POLARS_MAX_THREADS=%s, BLAS intra-op threads pinned to 1",
        os.environ["POLARS_MAX_THREADS"],
    )
    parser = argparse.ArgumentParser(prog="compute-properties")
    parser.add_argument("--input", required=True, help="path to input entity TSV")
    parser.add_argument("--plan", required=True, help="path to plan JSON")
    parser.add_argument("--output", required=True, help="path to write properties TSV")
    parser.add_argument("--aa-fraction", required=True, help="path to write AA fraction TSV")
    parser.add_argument("--stats", required=True, help="path to write dataset stats JSON")
    args = parser.parse_args(argv)

    reads = read_input_tsv(args.input)
    plan = read_plan(args.plan)
    outputs = run(reads, plan)
    write_output_tsv(outputs["properties"], args.output, sort_keys=["entity_key"])
    write_output_tsv(outputs["aa_fraction"], args.aa_fraction, sort_keys=["entity_key", "aminoAcid"])
    Path(args.stats).write_text(json.dumps(outputs["stats"], sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
