#!/usr/bin/env python3
"""Benchmark worker scaling on injective tetrahedral volumes."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TextIO


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "tools"))

import benchmark_simplicial_strategies as strategy_benchmark  # noqa: E402
import morseframes as mp  # noqa: E402


STRATEGIES = {
    "process-lower-stars": (
        mp.PROCESS_LOWER_STARS_SEQUENCE,
        mp.PROCESS_LOWER_STARS_PARALLEL_SEQUENCE,
    ),
    "reduction-kernel": (
        mp.FLOODING_REDUCTION_KERNEL_SEQUENCE,
        mp.FLOODING_REDUCTION_KERNEL_PARALLEL_SEQUENCE,
    ),
}


@dataclass(frozen=True)
class TetrahedralWorkerScalingRow:
    family: str
    name: str
    seed: int
    grid_size: int
    strategy: str
    algorithm: str
    max_workers: int
    repeats: int
    warmups: int
    cpp_backend: bool
    num_simplices: int
    num_levels: int
    num_critical_simplices: int
    critical_simplices_by_dimension: tuple[int, ...]
    sequence_matches_sequential: bool
    barcode_matches_standard: bool
    sequence_seconds: float
    total_seconds: float
    sequence_speedup_vs_one_worker: float
    total_speedup_vs_one_worker: float
    sequence_parallel_efficiency: float
    total_parallel_efficiency: float
    total_time_ratio_vs_f_max: float


def benchmark_scaling(
    seed: int,
    grid_size: int,
    *,
    workers: tuple[int, ...] = (1, 2, 4, 8),
    strategies: tuple[str, ...] = tuple(STRATEGIES),
    repeats: int = 3,
    warmups: int = 1,
) -> list[TetrahedralWorkerScalingRow]:
    selected_workers = tuple(dict.fromkeys(workers))
    selected_strategies = tuple(dict.fromkeys(strategies))
    if not selected_workers or any(worker < 1 for worker in selected_workers):
        raise ValueError("workers must contain positive integers")
    if 1 not in selected_workers:
        raise ValueError("workers must include the one-worker baseline")
    unknown = set(selected_strategies) - STRATEGIES.keys()
    if unknown:
        raise ValueError(f"Unknown strategies: {', '.join(sorted(unknown))}")

    complex_ = strategy_benchmark.make_injective_volume(seed, grid_size)
    sequential = {
        strategy: mp.compute_morse_sequence(complex_, algorithm=STRATEGIES[strategy][0])
        for strategy in selected_strategies
    }
    measured: dict[
        str, list[tuple[int, strategy_benchmark.SimplicialStrategyBenchmarkRow]]
    ] = {strategy: [] for strategy in selected_strategies}

    for strategy in selected_strategies:
        _, parallel_algorithm = STRATEGIES[strategy]
        for worker_count in selected_workers:
            parallel = mp.compute_morse_sequence(
                complex_, algorithm=parallel_algorithm, max_workers=worker_count
            )
            if parallel.steps != sequential[strategy].steps:
                raise AssertionError(
                    f"{worker_count}-worker {strategy} differs from sequential"
                )
            rows = strategy_benchmark.benchmark_volume(
                seed,
                grid_size,
                algorithms=(parallel_algorithm, mp.F_MAX_SEQUENCE),
                parallel_workers=worker_count,
                repeats=repeats,
                warmups=warmups,
            )
            measured[strategy].append((worker_count, rows[0]))

    output: list[TetrahedralWorkerScalingRow] = []
    for strategy in selected_strategies:
        baseline = next(
            row for worker_count, row in measured[strategy] if worker_count == 1
        )
        for worker_count, row in measured[strategy]:
            sequence_speedup = baseline.sequence_seconds / row.sequence_seconds
            total_speedup = baseline.total_seconds / row.total_seconds
            output.append(
                TetrahedralWorkerScalingRow(
                    family="injective-volume",
                    name=row.name,
                    seed=seed,
                    grid_size=grid_size,
                    strategy=strategy,
                    algorithm=row.algorithm,
                    max_workers=worker_count,
                    repeats=repeats,
                    warmups=warmups,
                    cpp_backend=row.cpp_backend,
                    num_simplices=row.num_simplices,
                    num_levels=row.num_levels,
                    num_critical_simplices=row.num_critical_simplices,
                    critical_simplices_by_dimension=(
                        row.critical_simplices_by_dimension
                    ),
                    sequence_matches_sequential=True,
                    barcode_matches_standard=row.barcode_matches_standard,
                    sequence_seconds=row.sequence_seconds,
                    total_seconds=row.total_seconds,
                    sequence_speedup_vs_one_worker=sequence_speedup,
                    total_speedup_vs_one_worker=total_speedup,
                    sequence_parallel_efficiency=sequence_speedup / worker_count,
                    total_parallel_efficiency=total_speedup / worker_count,
                    total_time_ratio_vs_f_max=1.0 / row.total_speedup_vs_f_max,
                )
            )
    return output


def write_rows(
    rows: list[TetrahedralWorkerScalingRow], output: TextIO, output_format: str
) -> None:
    if not rows:
        raise ValueError("No benchmark rows were generated")
    if output_format == "csv":
        writer = csv.DictWriter(output, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)
        return
    if output_format == "json":
        json.dump([asdict(row) for row in rows], output, indent=2)
        output.write("\n")
        return
    raise ValueError(f"Unknown output format: {output_format}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", type=int, nargs="+", default=(4, 8, 12))
    parser.add_argument("--seeds", type=int, nargs="+", default=(0, 1, 2))
    parser.add_argument("--workers", type=int, nargs="+", default=(1, 2, 4, 8))
    parser.add_argument(
        "--strategies", nargs="+", choices=tuple(STRATEGIES), default=tuple(STRATEGIES)
    )
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
        for row in benchmark_scaling(
            seed,
            grid_size,
            workers=tuple(args.workers),
            strategies=tuple(args.strategies),
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
