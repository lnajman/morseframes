#!/usr/bin/env python3
"""Compare MorseFrames strategies on connected injective lower-star terrains."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import random
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from time import perf_counter
from typing import Iterable, TextIO


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

import morseframes as mp  # noqa: E402


DEFAULT_ALGORITHMS = (
    mp.PROCESS_LOWER_STARS_SEQUENCE,
    mp.PROCESS_LOWER_STARS_PARALLEL_SEQUENCE,
    mp.FLOODING_REDUCTION_KERNEL_SEQUENCE,
    mp.FLOODING_REDUCTION_KERNEL_PARALLEL_SEQUENCE,
    mp.F_MAX_SEQUENCE,
    mp.F_MIN_SEQUENCE,
    mp.COREDUCTION_SEQUENCE,
    mp.SATURATED_SEQUENCE,
)


@dataclass(frozen=True)
class SimplicialStrategyBenchmarkRow:
    family: str
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
    max_dimension: int
    num_critical_simplices: int
    critical_simplices_by_dimension: tuple[int, ...]
    num_regular_pairs: int
    critical_ratio: float
    critical_count_delta_vs_f_max: int
    critical_count_ratio_vs_f_max: float
    sequence_matches_process_lower_stars: bool
    barcode_matches_standard: bool
    sequence_seconds: float
    persistence_seconds: float
    total_seconds: float
    sequence_seconds_per_eliminated_simplex: float
    sequence_speedup_vs_f_max: float
    total_speedup_vs_f_max: float
    reduction_kernel_metrics_available: bool
    reduction_kernel_levels: int
    reduction_kernel_rounds: int
    reduction_kernel_facet_kernels: int
    reduction_kernel_reductions: int
    reduction_kernel_perforations: int
    reduction_kernel_parallel_batches: int
    reduction_kernel_max_parallel_facets: int
    reduction_kernel_parallel_level_batches: int
    reduction_kernel_max_parallel_levels: int
    reduction_kernel_executor_workers: int
    reduction_kernel_facet_discovery_parallel_tasks: int
    reduction_kernel_essential_parallel_tasks: int
    reduction_kernel_aggregation_rounds: int
    reduction_kernel_aggregation_parallel_tasks: int
    reduction_kernel_cumulative_facet_seconds: float
    reduction_kernel_cumulative_essential_seconds: float
    reduction_kernel_cumulative_core_seconds: float
    reduction_kernel_cumulative_local_reduction_seconds: float
    reduction_kernel_cumulative_aggregation_seconds: float
    reduction_kernel_cumulative_merge_seconds: float


@dataclass(frozen=True)
class _Measurement:
    algorithm: str
    max_workers: int
    sequence: mp.MorseSequence
    sequence_seconds: float
    persistence_seconds: float
    total_seconds: float
    frame_metrics: dict[str, object]


PARALLEL_ALGORITHMS = {
    mp.PROCESS_LOWER_STARS_PARALLEL_SEQUENCE,
    mp.FLOODING_REDUCTION_KERNEL_PARALLEL_SEQUENCE,
}


def _metric_int(metrics: dict[str, object], name: str) -> int:
    return int(metrics.get(f"sequence_reduction_kernel_{name}", 0))


def _metric_seconds(metrics: dict[str, object], name: str) -> float:
    return 1.0e-9 * _metric_int(metrics, f"{name}_nanoseconds")


def make_injective_terrain(seed: int, grid_size: int) -> mp.FilteredComplex:
    """Create a connected triangulated terrain with strictly ordered vertices."""

    if grid_size < 2:
        raise ValueError("grid_size must be at least 2")
    rng = random.Random(seed)
    bumps = [
        (
            rng.random(),
            rng.random(),
            rng.uniform(-0.8, 0.8),
            rng.uniform(0.08, 0.22),
        )
        for _ in range(5)
    ]

    def vertex_id(row: int, col: int) -> int:
        return row * grid_size + col

    raw_values: dict[int, float] = {}
    for row in range(grid_size):
        x = row / float(grid_size - 1)
        for col in range(grid_size):
            y = col / float(grid_size - 1)
            value = (
                0.35 * x
                + 0.20 * y
                + 0.35 * math.sin(2.0 * math.pi * (x + 0.17 * seed))
                + 0.25 * math.cos(2.0 * math.pi * (1.7 * y - 0.11 * seed))
                + 0.15 * math.sin(2.0 * math.pi * (x + y))
            )
            for center_x, center_y, weight, sigma in bumps:
                squared_distance = (x - center_x) ** 2 + (y - center_y) ** 2
                value += weight * math.exp(-squared_distance / (2.0 * sigma * sigma))
            value += rng.uniform(-0.03, 0.03)
            raw_values[vertex_id(row, col)] = value

    ranked_vertices = sorted(
        raw_values, key=lambda vertex: (raw_values[vertex], vertex)
    )
    vertex_values = {vertex: float(rank) for rank, vertex in enumerate(ranked_vertices)}
    facets: list[tuple[int, int, int]] = []
    for row in range(grid_size - 1):
        for col in range(grid_size - 1):
            v00 = vertex_id(row, col)
            v10 = vertex_id(row + 1, col)
            v01 = vertex_id(row, col + 1)
            v11 = vertex_id(row + 1, col + 1)
            if (row + col + seed) % 2 == 0:
                facets.extend(((v00, v10, v11), (v00, v11, v01)))
            else:
                facets.extend(((v00, v10, v01), (v10, v11, v01)))
    return mp.FilteredComplex.from_lower_star(facets, vertex_values)


def make_injective_volume(seed: int, grid_size: int) -> mp.FilteredComplex:
    """Create a tetrahedral volume with a strict lower-star filtration."""

    if grid_size < 2:
        raise ValueError("grid_size must be at least 2")
    rng = random.Random(seed)
    bumps = [
        (
            rng.random(),
            rng.random(),
            rng.random(),
            rng.uniform(-0.8, 0.8),
            rng.uniform(0.10, 0.24),
        )
        for _ in range(7)
    ]

    def vertex_id(x: int, y: int, z: int) -> int:
        return (x * grid_size + y) * grid_size + z

    raw_values: dict[int, float] = {}
    for x_index in range(grid_size):
        x = x_index / float(grid_size - 1)
        for y_index in range(grid_size):
            y = y_index / float(grid_size - 1)
            for z_index in range(grid_size):
                z = z_index / float(grid_size - 1)
                value = (
                    0.25 * x
                    + 0.18 * y
                    + 0.12 * z
                    + 0.30 * math.sin(2.0 * math.pi * (x + 0.13 * seed))
                    + 0.22 * math.cos(2.0 * math.pi * (1.4 * y - 0.09 * seed))
                    + 0.18 * math.sin(2.0 * math.pi * (x + y + z))
                )
                for center_x, center_y, center_z, weight, sigma in bumps:
                    squared_distance = (
                        (x - center_x) ** 2 + (y - center_y) ** 2 + (z - center_z) ** 2
                    )
                    value += weight * math.exp(
                        -squared_distance / (2.0 * sigma * sigma)
                    )
                value += rng.uniform(-0.03, 0.03)
                raw_values[vertex_id(x_index, y_index, z_index)] = value

    ranked_vertices = sorted(
        raw_values, key=lambda vertex: (raw_values[vertex], vertex)
    )
    vertex_values = {vertex: float(rank) for rank, vertex in enumerate(ranked_vertices)}

    axes = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    facets: list[tuple[int, int, int, int]] = []
    for x in range(grid_size - 1):
        for y in range(grid_size - 1):
            for z in range(grid_size - 1):
                for order in itertools.permutations(range(3)):
                    coordinates = [x, y, z]
                    vertices = [vertex_id(*coordinates)]
                    for axis_index in order:
                        coordinates = [
                            coordinate + step
                            for coordinate, step in zip(
                                coordinates, axes[axis_index], strict=True
                            )
                        ]
                        vertices.append(vertex_id(*coordinates))
                    facets.append(tuple(vertices))
    return mp.FilteredComplex.from_lower_star(facets, vertex_values)


def _run_pipeline(
    complex_: mp.FilteredComplex,
    algorithm: str,
    max_workers: int | None,
) -> tuple[mp.MorseSequence, mp.PersistenceDiagram, float, float, float]:
    started = perf_counter()
    sequence = mp.compute_morse_sequence(
        complex_, algorithm=algorithm, max_workers=max_workers
    )
    sequence_finished = perf_counter()
    references = mp.compute_reference_map(complex_, sequence, algorithm=algorithm)
    diagram = mp.compute_morse_persistence(
        complex_, sequence, references, algorithm=algorithm
    )
    finished = perf_counter()
    return (
        sequence,
        diagram,
        sequence_finished - started,
        finished - sequence_finished,
        finished - started,
    )


def _critical_counts(
    complex_: mp.FilteredComplex,
    sequence: mp.MorseSequence,
) -> tuple[int, ...]:
    max_dimension = max(
        (complex_.dimension(simplex) for simplex in range(complex_.size)),
        default=-1,
    )
    counts = [0] * (max_dimension + 1)
    for simplex in sequence.critical_simplices:
        counts[complex_.dimension(simplex)] += 1
    return tuple(counts)


def _benchmark_complex(
    complex_: mp.FilteredComplex,
    *,
    family: str,
    name: str,
    seed: int,
    grid_size: int,
    num_vertices: int,
    algorithms: Iterable[str] = DEFAULT_ALGORITHMS,
    parallel_workers: int = 8,
    repeats: int = 3,
    warmups: int = 1,
) -> list[SimplicialStrategyBenchmarkRow]:
    if parallel_workers < 1:
        raise ValueError("parallel_workers must be positive")
    if repeats < 1:
        raise ValueError("repeats must be positive")
    if warmups < 0:
        raise ValueError("warmups must be non-negative")
    selected_algorithms = tuple(dict.fromkeys(algorithms))
    if mp.F_MAX_SEQUENCE not in selected_algorithms:
        raise ValueError("The F-Max baseline must be included")

    standard = mp.compute_standard_persistence(complex_)
    standard_finite = standard.finite_barcode()
    standard_essential = standard.essential_barcode()
    process_lower_stars = mp.compute_morse_sequence(
        complex_, algorithm=mp.PROCESS_LOWER_STARS_SEQUENCE
    )
    measurements: list[_Measurement] = []

    for algorithm in selected_algorithms:
        max_workers = parallel_workers if algorithm in PARALLEL_ALGORITHMS else None
        for _ in range(warmups):
            _run_pipeline(complex_, algorithm, max_workers)

        best: _Measurement | None = None
        for _ in range(repeats):
            sequence, diagram, sequence_seconds, persistence_seconds, total_seconds = (
                _run_pipeline(complex_, algorithm, max_workers)
            )
            if (
                diagram.finite_barcode() != standard_finite
                or diagram.essential_barcode() != standard_essential
            ):
                raise AssertionError(f"{algorithm} barcode differs from standard")
            measurement = _Measurement(
                algorithm=algorithm,
                max_workers=parallel_workers if max_workers is not None else 1,
                sequence=sequence,
                sequence_seconds=sequence_seconds,
                persistence_seconds=persistence_seconds,
                total_seconds=total_seconds,
                frame_metrics={},
            )
            if best is None or measurement.total_seconds < best.total_seconds:
                best = measurement
        if best is None:
            raise RuntimeError("Benchmark did not run")
        profile = mp.profile_morse_reference_frame(
            complex_, algorithm=algorithm, max_workers=max_workers
        )
        measurements.append(
            _Measurement(
                algorithm=best.algorithm,
                max_workers=best.max_workers,
                sequence=best.sequence,
                sequence_seconds=best.sequence_seconds,
                persistence_seconds=best.persistence_seconds,
                total_seconds=best.total_seconds,
                frame_metrics=profile.frame_metrics,
            )
        )

    f_max = next(
        measurement
        for measurement in measurements
        if measurement.algorithm == mp.F_MAX_SEQUENCE
    )
    f_max_critical_count = len(f_max.sequence.critical_simplices)
    max_dimension = max(
        (complex_.dimension(simplex) for simplex in range(complex_.size)),
        default=-1,
    )
    rows: list[SimplicialStrategyBenchmarkRow] = []
    for measurement in measurements:
        critical_count = len(measurement.sequence.critical_simplices)
        eliminated = complex_.size - critical_count
        kernel_algorithm = measurement.algorithm in {
            mp.FLOODING_REDUCTION_KERNEL_SEQUENCE,
            mp.FLOODING_REDUCTION_KERNEL_PARALLEL_SEQUENCE,
        }
        kernel_metrics = measurement.frame_metrics
        rows.append(
            SimplicialStrategyBenchmarkRow(
                family=family,
                name=name,
                seed=seed,
                grid_size=grid_size,
                algorithm=measurement.algorithm,
                max_workers=measurement.max_workers,
                repeats=repeats,
                warmups=warmups,
                cpp_backend=complex_.cpp_backend_active(),
                num_simplices=complex_.size,
                num_vertices=num_vertices,
                num_levels=complex_.num_levels,
                max_dimension=max_dimension,
                num_critical_simplices=critical_count,
                critical_simplices_by_dimension=_critical_counts(
                    complex_, measurement.sequence
                ),
                num_regular_pairs=eliminated // 2,
                critical_ratio=float(critical_count) / complex_.size,
                critical_count_delta_vs_f_max=(critical_count - f_max_critical_count),
                critical_count_ratio_vs_f_max=(
                    float(critical_count) / f_max_critical_count
                    if f_max_critical_count
                    else math.inf
                ),
                sequence_matches_process_lower_stars=(
                    measurement.sequence.steps == process_lower_stars.steps
                ),
                barcode_matches_standard=True,
                sequence_seconds=measurement.sequence_seconds,
                persistence_seconds=measurement.persistence_seconds,
                total_seconds=measurement.total_seconds,
                sequence_seconds_per_eliminated_simplex=(
                    measurement.sequence_seconds / eliminated
                    if eliminated
                    else math.inf
                ),
                sequence_speedup_vs_f_max=(
                    f_max.sequence_seconds / measurement.sequence_seconds
                    if measurement.sequence_seconds
                    else math.inf
                ),
                total_speedup_vs_f_max=(
                    f_max.total_seconds / measurement.total_seconds
                    if measurement.total_seconds
                    else math.inf
                ),
                reduction_kernel_metrics_available=(
                    kernel_algorithm and complex_.cpp_backend_active()
                ),
                reduction_kernel_levels=_metric_int(kernel_metrics, "levels"),
                reduction_kernel_rounds=_metric_int(kernel_metrics, "rounds"),
                reduction_kernel_facet_kernels=_metric_int(
                    kernel_metrics, "facet_kernels"
                ),
                reduction_kernel_reductions=_metric_int(kernel_metrics, "reductions"),
                reduction_kernel_perforations=_metric_int(
                    kernel_metrics, "perforations"
                ),
                reduction_kernel_parallel_batches=_metric_int(
                    kernel_metrics, "parallel_batches"
                ),
                reduction_kernel_max_parallel_facets=_metric_int(
                    kernel_metrics, "max_parallel_facets"
                ),
                reduction_kernel_parallel_level_batches=_metric_int(
                    kernel_metrics, "parallel_level_batches"
                ),
                reduction_kernel_max_parallel_levels=_metric_int(
                    kernel_metrics, "max_parallel_levels"
                ),
                reduction_kernel_executor_workers=_metric_int(
                    kernel_metrics, "executor_workers"
                ),
                reduction_kernel_facet_discovery_parallel_tasks=_metric_int(
                    kernel_metrics, "facet_discovery_parallel_tasks"
                ),
                reduction_kernel_essential_parallel_tasks=_metric_int(
                    kernel_metrics, "essential_parallel_tasks"
                ),
                reduction_kernel_aggregation_rounds=_metric_int(
                    kernel_metrics, "aggregation_rounds"
                ),
                reduction_kernel_aggregation_parallel_tasks=_metric_int(
                    kernel_metrics, "aggregation_parallel_tasks"
                ),
                reduction_kernel_cumulative_facet_seconds=_metric_seconds(
                    kernel_metrics, "facet"
                ),
                reduction_kernel_cumulative_essential_seconds=_metric_seconds(
                    kernel_metrics, "essential"
                ),
                reduction_kernel_cumulative_core_seconds=_metric_seconds(
                    kernel_metrics, "core"
                ),
                reduction_kernel_cumulative_local_reduction_seconds=_metric_seconds(
                    kernel_metrics, "local_reduction"
                ),
                reduction_kernel_cumulative_aggregation_seconds=_metric_seconds(
                    kernel_metrics, "aggregation"
                ),
                reduction_kernel_cumulative_merge_seconds=_metric_seconds(
                    kernel_metrics, "merge"
                ),
            )
        )
    return rows


def benchmark_terrain(
    seed: int,
    grid_size: int,
    *,
    algorithms: Iterable[str] = DEFAULT_ALGORITHMS,
    parallel_workers: int = 8,
    repeats: int = 3,
    warmups: int = 1,
) -> list[SimplicialStrategyBenchmarkRow]:
    return _benchmark_complex(
        make_injective_terrain(seed, grid_size),
        family="injective-terrain",
        name=f"injective-terrain-n{grid_size}-seed{seed}",
        seed=seed,
        grid_size=grid_size,
        num_vertices=grid_size**2,
        algorithms=algorithms,
        parallel_workers=parallel_workers,
        repeats=repeats,
        warmups=warmups,
    )


def benchmark_volume(
    seed: int,
    grid_size: int,
    *,
    algorithms: Iterable[str] = DEFAULT_ALGORITHMS,
    parallel_workers: int = 8,
    repeats: int = 3,
    warmups: int = 1,
) -> list[SimplicialStrategyBenchmarkRow]:
    return _benchmark_complex(
        make_injective_volume(seed, grid_size),
        family="injective-volume",
        name=f"injective-volume-n{grid_size}-seed{seed}",
        seed=seed,
        grid_size=grid_size,
        num_vertices=grid_size**3,
        algorithms=algorithms,
        parallel_workers=parallel_workers,
        repeats=repeats,
        warmups=warmups,
    )


def write_rows(
    rows: list[SimplicialStrategyBenchmarkRow],
    output: TextIO,
    output_format: str,
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
        for row in benchmark_terrain(
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
