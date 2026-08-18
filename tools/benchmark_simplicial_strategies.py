#!/usr/bin/env python3
"""Compare MorseFrames strategies on connected injective lower-star terrains."""

from __future__ import annotations

import argparse
import csv
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


@dataclass(frozen=True)
class _Measurement:
    algorithm: str
    max_workers: int
    sequence: mp.MorseSequence
    sequence_seconds: float
    persistence_seconds: float
    total_seconds: float


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
                value += weight * math.exp(
                    -squared_distance / (2.0 * sigma * sigma)
                )
            value += rng.uniform(-0.03, 0.03)
            raw_values[vertex_id(row, col)] = value

    ranked_vertices = sorted(raw_values, key=lambda vertex: (raw_values[vertex], vertex))
    vertex_values = {
        vertex: float(rank) for rank, vertex in enumerate(ranked_vertices)
    }
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


def benchmark_terrain(
    seed: int,
    grid_size: int,
    *,
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

    complex_ = make_injective_terrain(seed, grid_size)
    standard = mp.compute_standard_persistence(complex_)
    standard_finite = standard.finite_barcode()
    standard_essential = standard.essential_barcode()
    process_lower_stars = mp.compute_morse_sequence(
        complex_, algorithm=mp.PROCESS_LOWER_STARS_SEQUENCE
    )
    measurements: list[_Measurement] = []

    for algorithm in selected_algorithms:
        max_workers = (
            parallel_workers
            if algorithm == mp.PROCESS_LOWER_STARS_PARALLEL_SEQUENCE
            else None
        )
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
            )
            if best is None or measurement.total_seconds < best.total_seconds:
                best = measurement
        if best is None:
            raise RuntimeError("Benchmark did not run")
        measurements.append(best)

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
        rows.append(
            SimplicialStrategyBenchmarkRow(
                family="injective-terrain",
                name=f"injective-terrain-n{grid_size}-seed{seed}",
                seed=seed,
                grid_size=grid_size,
                algorithm=measurement.algorithm,
                max_workers=measurement.max_workers,
                repeats=repeats,
                warmups=warmups,
                cpp_backend=complex_.cpp_backend_active(),
                num_simplices=complex_.size,
                num_vertices=grid_size * grid_size,
                num_levels=complex_.num_levels,
                max_dimension=max_dimension,
                num_critical_simplices=critical_count,
                critical_simplices_by_dimension=_critical_counts(
                    complex_, measurement.sequence
                ),
                num_regular_pairs=eliminated // 2,
                critical_ratio=float(critical_count) / complex_.size,
                critical_count_delta_vs_f_max=(
                    critical_count - f_max_critical_count
                ),
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
            )
        )
    return rows


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
