#!/usr/bin/env python3
"""Compare gradient construction on injective 2D and 3D lower-star complexes."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TextIO


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "tools"))

import benchmark_simplicial_strategies as generators  # noqa: E402
import morseframes as mp  # noqa: E402


SEQUENTIAL_ALGORITHMS = (
    mp.F_MAX_SEQUENCE,
    mp.PROCESS_LOWER_STARS_SEQUENCE,
    mp.FLOODING_REDUCTION_KERNEL_SEQUENCE,
)
PARALLEL_ALGORITHMS = (
    mp.PROCESS_LOWER_STARS_PARALLEL_SEQUENCE,
    mp.FLOODING_REDUCTION_KERNEL_PARALLEL_SEQUENCE,
)
SEQUENTIAL_COUNTERPART = {
    mp.PROCESS_LOWER_STARS_PARALLEL_SEQUENCE: mp.PROCESS_LOWER_STARS_SEQUENCE,
    mp.FLOODING_REDUCTION_KERNEL_PARALLEL_SEQUENCE: (
        mp.FLOODING_REDUCTION_KERNEL_SEQUENCE
    ),
}


@dataclass(frozen=True)
class GradientStrategyRow:
    family: str
    dimension: int
    name: str
    seed: int
    grid_size: int
    algorithm: str
    max_workers: int
    repeats: int
    warmups: int
    cpp_backend: bool
    num_simplices: int
    num_vertices: int
    num_levels: int
    num_critical_simplices: int
    critical_simplices_by_dimension: tuple[int, ...]
    num_regular_pairs: int
    matches_sequential: bool
    critical_count_delta_vs_f_max: int
    critical_count_ratio_vs_f_max: float
    construction_seconds: float
    construction_time_ratio_vs_f_max: float
    construction_speedup_vs_f_max: float
    simplices_per_second: float
    worker_speedup_vs_one: float
    parallel_efficiency: float


def _best_profile(
    complex_: mp.FilteredComplex,
    algorithm: str,
    max_workers: int | None,
    *,
    repeats: int,
    warmups: int,
) -> mp.MorseSequenceProfile:
    for _ in range(warmups):
        mp.profile_morse_sequence(
            complex_, algorithm=algorithm, max_workers=max_workers
        )
    profiles = [
        mp.profile_morse_sequence(
            complex_, algorithm=algorithm, max_workers=max_workers
        )
        for _ in range(repeats)
    ]
    return min(profiles, key=lambda profile: profile.construction_seconds)


def benchmark_complex(
    complex_: mp.FilteredComplex,
    *,
    family: str,
    dimension: int,
    name: str,
    seed: int,
    grid_size: int,
    num_vertices: int,
    workers: tuple[int, ...] = (1, 2, 4, 8),
    repeats: int = 5,
    warmups: int = 1,
) -> list[GradientStrategyRow]:
    selected_workers = tuple(dict.fromkeys(workers))
    if not selected_workers or any(worker < 1 for worker in selected_workers):
        raise ValueError("workers must contain positive integers")
    if 1 not in selected_workers:
        raise ValueError("workers must include the one-worker baseline")
    if repeats < 1:
        raise ValueError("repeats must be positive")
    if warmups < 0:
        raise ValueError("warmups must be non-negative")

    sequences = {
        algorithm: mp.compute_morse_sequence(complex_, algorithm=algorithm)
        for algorithm in SEQUENTIAL_ALGORITHMS
    }
    profiles: dict[tuple[str, int], mp.MorseSequenceProfile] = {}
    matches: dict[tuple[str, int], bool] = {}
    for algorithm in SEQUENTIAL_ALGORITHMS:
        profiles[(algorithm, 1)] = _best_profile(
            complex_, algorithm, None, repeats=repeats, warmups=warmups
        )
        matches[(algorithm, 1)] = True
    for algorithm in PARALLEL_ALGORITHMS:
        sequential = sequences[SEQUENTIAL_COUNTERPART[algorithm]]
        for worker_count in selected_workers:
            parallel = mp.compute_morse_sequence(
                complex_, algorithm=algorithm, max_workers=worker_count
            )
            matches[(algorithm, worker_count)] = parallel.steps == sequential.steps
            if not matches[(algorithm, worker_count)]:
                raise AssertionError(
                    f"{worker_count}-worker {algorithm} differs from sequential"
                )
            profiles[(algorithm, worker_count)] = _best_profile(
                complex_,
                algorithm,
                worker_count,
                repeats=repeats,
                warmups=warmups,
            )

    f_max = profiles[(mp.F_MAX_SEQUENCE, 1)]
    rows: list[GradientStrategyRow] = []
    for (algorithm, worker_count), profile in profiles.items():
        one_worker = profiles.get((algorithm, 1), profile)
        worker_speedup = (
            one_worker.construction_seconds / profile.construction_seconds
        )
        rows.append(
            GradientStrategyRow(
                family=family,
                dimension=dimension,
                name=name,
                seed=seed,
                grid_size=grid_size,
                algorithm=algorithm,
                max_workers=worker_count,
                repeats=repeats,
                warmups=warmups,
                cpp_backend=mp.cpp_backend_active(complex_),
                num_simplices=profile.num_simplices,
                num_vertices=num_vertices,
                num_levels=profile.num_levels,
                num_critical_simplices=profile.num_critical_simplices,
                critical_simplices_by_dimension=(
                    profile.critical_simplices_by_dimension
                ),
                num_regular_pairs=profile.num_regular_pairs,
                matches_sequential=matches[(algorithm, worker_count)],
                critical_count_delta_vs_f_max=(
                    profile.num_critical_simplices
                    - f_max.num_critical_simplices
                ),
                critical_count_ratio_vs_f_max=(
                    profile.num_critical_simplices
                    / f_max.num_critical_simplices
                    if f_max.num_critical_simplices
                    else math.inf
                ),
                construction_seconds=profile.construction_seconds,
                construction_time_ratio_vs_f_max=(
                    profile.construction_seconds / f_max.construction_seconds
                ),
                construction_speedup_vs_f_max=(
                    f_max.construction_seconds / profile.construction_seconds
                ),
                simplices_per_second=profile.simplices_per_second,
                worker_speedup_vs_one=worker_speedup,
                parallel_efficiency=worker_speedup / worker_count,
            )
        )
    return rows


def benchmark_terrain(
    seed: int,
    grid_size: int,
    *,
    workers: tuple[int, ...] = (1, 2, 4, 8),
    repeats: int = 5,
    warmups: int = 1,
) -> list[GradientStrategyRow]:
    return benchmark_complex(
        generators.make_injective_terrain(seed, grid_size),
        family="injective-terrain",
        dimension=2,
        name=f"injective-terrain-n{grid_size}-seed{seed}",
        seed=seed,
        grid_size=grid_size,
        num_vertices=grid_size**2,
        workers=workers,
        repeats=repeats,
        warmups=warmups,
    )


def benchmark_volume(
    seed: int,
    grid_size: int,
    *,
    workers: tuple[int, ...] = (1, 2, 4, 8),
    repeats: int = 5,
    warmups: int = 1,
) -> list[GradientStrategyRow]:
    return benchmark_complex(
        generators.make_injective_volume(seed, grid_size),
        family="injective-volume",
        dimension=3,
        name=f"injective-volume-n{grid_size}-seed{seed}",
        seed=seed,
        grid_size=grid_size,
        num_vertices=grid_size**3,
        workers=workers,
        repeats=repeats,
        warmups=warmups,
    )


def write_rows(
    rows: list[GradientStrategyRow], output: TextIO, output_format: str
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
    parser.add_argument(
        "--families",
        nargs="+",
        choices=("terrain", "volume"),
        default=("terrain", "volume"),
    )
    parser.add_argument("--terrain-sizes", type=int, nargs="+", default=(16, 32, 64))
    parser.add_argument("--volume-sizes", type=int, nargs="+", default=(4, 8, 12))
    parser.add_argument("--seeds", type=int, nargs="+", default=(0, 1, 2))
    parser.add_argument("--workers", type=int, nargs="+", default=(1, 2, 4, 8))
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--format", choices=("csv", "json"), default="csv")
    parser.add_argument("--output", type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    options = {
        "workers": tuple(args.workers),
        "repeats": args.repeats,
        "warmups": args.warmups,
    }
    rows: list[GradientStrategyRow] = []
    if "terrain" in args.families:
        rows.extend(
            row
            for grid_size in args.terrain_sizes
            for seed in args.seeds
            for row in benchmark_terrain(seed, grid_size, **options)
        )
    if "volume" in args.families:
        rows.extend(
            row
            for grid_size in args.volume_sizes
            for seed in args.seeds
            for row in benchmark_volume(seed, grid_size, **options)
        )
    if args.output is None:
        write_rows(rows, sys.stdout, args.format)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", newline="") as output:
            write_rows(rows, output, args.format)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
