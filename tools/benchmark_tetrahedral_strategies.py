#!/usr/bin/env python3
"""Compare MorseFrames strategies on injective tetrahedral volumes."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from benchmark_simplicial_strategies import (
    DEFAULT_ALGORITHMS,
    benchmark_volume,
    write_rows,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", type=int, nargs="+", default=(4, 8, 12))
    parser.add_argument("--seeds", type=int, nargs="+", default=(0, 1, 2))
    parser.add_argument("--algorithms", nargs="+", default=DEFAULT_ALGORITHMS)
    parser.add_argument("--parallel-workers", type=int, default=8)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--format", choices=("csv", "json"), default="csv")
    parser.add_argument("--output", type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    rows = [
        row
        for grid_size in args.sizes
        for seed in args.seeds
        for row in benchmark_volume(
            seed,
            grid_size,
            algorithms=args.algorithms,
            parallel_workers=args.parallel_workers,
            repeats=args.repeats,
            warmups=args.warmups,
        )
    ]
    if args.output is None:
        write_rows(rows, sys.stdout, args.format)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", newline="") as output:
            write_rows(rows, output, args.format)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
