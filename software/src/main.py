"""CLI entry point — reads input TSV + plan, dispatches pipeline, writes outputs.

Invoked from the Tengo workflow as:

    python main.py --input input.tsv --plan plan.json
                   --output properties.tsv --aa-fraction aa_fraction.tsv
                   --stats stats.json
"""

from __future__ import annotations

import os

# Pin the polars thread pool to 1 BEFORE polars (or anything importing it,
# including `pipeline`) is loaded — POLARS_MAX_THREADS is read once at import.
# The block reserves cpu(1), and single-threaded execution keeps output
# byte-stable (no thread-count-dependent reduction order). setdefault so an
# explicit override from the environment still wins.
os.environ.setdefault("POLARS_MAX_THREADS", "1")

import argparse
import json
import logging
import sys
from pathlib import Path

from io_layer import TsvChunkWriter, read_input_batches, read_input_tsv, read_plan
from pipeline import (
    InputNotSorted,
    Reads,
    checked_key_order,
    has_integer_keys,
    key_order,
    run_streaming,
)

log = logging.getLogger(__name__)


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
    parser = argparse.ArgumentParser(prog="compute-properties")
    parser.add_argument("--input", required=True, help="path to input entity TSV")
    parser.add_argument("--plan", required=True, help="path to plan JSON")
    parser.add_argument("--output", required=True, help="path to write properties TSV")
    parser.add_argument("--aa-fraction", required=True, help="path to write AA fraction TSV")
    parser.add_argument("--stats", required=True, help="path to write dataset stats JSON")
    args = parser.parse_args(argv)

    plan = read_plan(args.plan)

    def compute(reads: Reads) -> dict:
        with TsvChunkWriter(args.output) as properties, TsvChunkWriter(args.aa_fraction) as aa_fraction:

            def emit(props, aa):
                properties.write(props)
                aa_fraction.write(aa)

            return run_streaming(reads, plan, emit)

    # Output rows are in ascending entity_key order, so the output bytes do not
    # depend on the input row order. An input already in that order is computed
    # batch by batch; any other input is read whole and sorted first.
    integer_keys = has_integer_keys(plan)
    try:
        stats = compute(checked_key_order(read_input_batches(args.input), integer_keys))
    except InputNotSorted:
        log.warning("Input rows are not in entity_key order: reading the whole table to sort it")
        stats = compute(read_input_tsv(args.input).sort(key_order(integer_keys), maintain_order=True))
    Path(args.stats).write_text(json.dumps(stats, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
