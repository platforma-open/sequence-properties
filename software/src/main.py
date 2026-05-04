"""CLI entry point — reads input TSV + plan, dispatches pipeline, writes outputs.

Invoked from the Tengo workflow as:

    python main.py --input input.tsv --plan plan.json
                   --output properties.tsv --aa-fraction aa_fraction.tsv
                   --stats stats.json
"""

from __future__ import annotations

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
