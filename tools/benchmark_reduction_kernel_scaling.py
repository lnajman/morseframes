#!/usr/bin/env python3
"""Benchmark deterministic reduction-kernel scaling on injective terrains."""

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

import benchmark_simplicial_strategies as terrain_benchmark  # noqa: E402
import morseframes as mp  # noqa: E402


@dataclass(frozen=True)
class ReductionKernelScalingRow:
    family: str
    name: str
    seed: int
    grid_size: int
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
    parallel_level_batches: int
    max_parallel_levels: int
    parallel_facet_batches: int
    max_parallel_facets: int
    facet_discovery_parallel_tasks: int
    essential_parallel_tasks: int
    aggregation_parallel_tasks: int


def benchmark_scaling(
    seed: int,
    grid_size: int,
    *,
    workers: tuple[int, ...] = (1, 2, 4, 8),
    repeats: int = 3,
    warmups: int = 1,
) -> list[ReductionKernelScalingRow]:
    selected_workers = tuple(dict.fromkeys(workers))
    if not selected_workers or any(worker < 1 for worker in selected_workers):
        raise ValueError("workers must contain positive integers")
    if 1 not in selected_workers:
        raise ValueError("workers must include the one-worker baseline")

    complex_ = terrain_benchmark.make_injective_terrain(seed, grid_size)
    sequential = mp.compute_morse_sequence(
        complex_, algorithm=mp.FLOODING_REDUCTION_KERNEL_SEQUENCE
    )
    measured: list[tuple[int, terrain_benchmark.SimplicialStrategyBenchmarkRow]] = []
    for worker_count in selected_workers:
        parallel = mp.compute_morse_sequence(
            complex_,
            algorithm=mp.FLOODING_REDUCTION_KERNEL_PARALLEL_SEQUENCE,
            max_workers=worker_count,
        )
        if parallel.steps != sequential.steps:
            raise AssertionError(
                f"{worker_count}-worker reduction kernel differs from sequential"
            )
        rows = terrain_benchmark.benchmark_terrain(
            seed,
            grid_size,
            algorithms=(
                mp.FLOODING_REDUCTION_KERNEL_PARALLEL_SEQUENCE,
                mp.F_MAX_SEQUENCE,
            ),
            parallel_workers=worker_count,
            repeats=repeats,
            warmups=warmups,
        )
        measured.append((worker_count, rows[0]))

    baseline = next(row for worker_count, row in measured if worker_count == 1)
    output: list[ReductionKernelScalingRow] = []
    for worker_count, row in measured:
        sequence_speedup = baseline.sequence_seconds / row.sequence_seconds
        total_speedup = baseline.total_seconds / row.total_seconds
        output.append(
            ReductionKernelScalingRow(
                family="injective-terrain",
                name=row.name,
                seed=seed,
                grid_size=grid_size,
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
                parallel_level_batches=row.reduction_kernel_parallel_level_batches,
                max_parallel_levels=row.reduction_kernel_max_parallel_levels,
                parallel_facet_batches=row.reduction_kernel_parallel_batches,
                max_parallel_facets=row.reduction_kernel_max_parallel_facets,
                facet_discovery_parallel_tasks=(
                    row.reduction_kernel_facet_discovery_parallel_tasks
                ),
                essential_parallel_tasks=(
                    row.reduction_kernel_essential_parallel_tasks
                ),
                aggregation_parallel_tasks=(
                    row.reduction_kernel_aggregation_parallel_tasks
                ),
            )
        )
    return output


def write_rows(
    rows: list[ReductionKernelScalingRow], output: TextIO, output_format: str
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
    parser.add_argument("--sizes", type=int, nargs="+", default=(16, 32, 64))
    parser.add_argument("--seeds", type=int, nargs="+", default=(0, 1, 2))
    parser.add_argument("--workers", type=int, nargs="+", default=(1, 2, 4, 8))
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
