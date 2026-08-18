#!/usr/bin/env python3
"""Benchmark sequential and parallel simplicial ProcessLowerStars."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from time import perf_counter
from typing import Iterable, TextIO


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

import morseframes as mp  # noqa: E402


@dataclass(frozen=True)
class LowerStarCase:
    name: str
    distribution: str
    fan_sizes: tuple[int, ...]
    complex: mp.FilteredComplex
    lower_star_sizes: tuple[int, ...]


@dataclass(frozen=True)
class ProcessLowerStarsBenchmarkRow:
    case: str
    distribution: str
    anchor_count: int
    min_anchor_fan_size: int
    max_anchor_fan_size: int
    total_anchor_fan_size: int
    algorithm: str
    max_workers: int
    repeats: int
    warmups: int
    cpp_backend: bool
    num_simplices: int
    num_vertices: int
    max_dimension: int
    lower_star_count: int
    max_lower_star_size: int
    mean_lower_star_size: float
    estimated_min_task_load: int
    estimated_max_task_load: int
    estimated_task_load_ratio: float
    num_critical_simplices: int
    critical_simplices_by_dimension: tuple[int, ...]
    num_regular_pairs: int
    critical_ratio: float
    eliminated_simplices: int
    exact_sequence_matches_sequential: bool
    sequence_seconds: float
    persistence_seconds: float
    total_seconds: float
    sequence_seconds_per_eliminated_simplex: float
    sequence_speedup_vs_sequential: float
    total_speedup_vs_sequential: float
    sequence_parallel_efficiency: float


def _validate_fan_sizes(fan_sizes: Iterable[int]) -> tuple[int, ...]:
    result = tuple(int(size) for size in fan_sizes)
    if not result:
        raise ValueError("At least one anchor is required.")
    if any(size < 1 for size in result):
        raise ValueError("Every fan size must be positive.")
    return result


def make_fan_case(name: str, distribution: str, fan_sizes: Iterable[int]) -> LowerStarCase:
    """Build disconnected triangle fans with an injective lower-star filtration."""

    sizes = _validate_fan_sizes(fan_sizes)
    facets: list[tuple[int, int, int]] = []
    base_vertices: list[int] = []
    next_vertex = 0
    bases_by_anchor: list[tuple[tuple[int, int], ...]] = []
    for fan_size in sizes:
        bases: list[tuple[int, int]] = []
        for _ in range(fan_size):
            left = next_vertex
            right = next_vertex + 1
            next_vertex += 2
            base_vertices.extend((left, right))
            bases.append((left, right))
        bases_by_anchor.append(tuple(bases))

    anchors = tuple(range(next_vertex, next_vertex + len(sizes)))
    for anchor, bases in zip(anchors, bases_by_anchor):
        facets.extend((left, right, anchor) for left, right in bases)

    vertex_values = {
        vertex: float(rank)
        for rank, vertex in enumerate((*base_vertices, *anchors))
    }
    complex_ = mp.FilteredComplex.from_lower_star(facets, vertex_values)
    return LowerStarCase(
        name=name,
        distribution=distribution,
        fan_sizes=sizes,
        complex=complex_,
        lower_star_sizes=lower_star_sizes(complex_),
    )


def make_balanced_case(anchor_count: int, fan_size: int) -> LowerStarCase:
    return make_fan_case(
        f"balanced-a{anchor_count}-f{fan_size}",
        "balanced",
        (fan_size,) * anchor_count,
    )


def make_skewed_case(
    anchor_count: int,
    light_fan_size: int,
    heavy_fan_size: int,
) -> LowerStarCase:
    if anchor_count < 1:
        raise ValueError("At least one anchor is required.")
    if heavy_fan_size < light_fan_size:
        raise ValueError("The heavy fan size must be at least the light fan size.")
    return make_fan_case(
        f"skewed-a{anchor_count}-f{light_fan_size}-{heavy_fan_size}",
        "skewed",
        (heavy_fan_size,) + (light_fan_size,) * (anchor_count - 1),
    )


def lower_star_sizes(complex_: mp.FilteredComplex) -> tuple[int, ...]:
    vertex_order = tuple(
        simplex
        for simplex in complex_.filtration_order
        if complex_.dimension(simplex) == 0
    )
    vertex_rank = {
        complex_.vertices(simplex)[0]: rank
        for rank, simplex in enumerate(vertex_order)
    }
    counts = [0] * len(vertex_order)
    for simplex in range(complex_.size):
        owner_rank = max(vertex_rank[vertex] for vertex in complex_.vertices(simplex))
        counts[owner_rank] += 1
    return tuple(counts)


def estimated_task_loads(star_sizes: tuple[int, ...], max_workers: int) -> tuple[int, ...]:
    if max_workers < 1:
        raise ValueError("max_workers must be positive.")
    task_count = min(max_workers, len(star_sizes))
    loads = [0] * task_count
    ordered = sorted(enumerate(star_sizes), key=lambda item: (-item[1], item[0]))
    for _, star_size in ordered:
        task = min(range(task_count), key=lambda index: (loads[index], index))
        loads[task] += star_size
    return tuple(loads)


def _critical_counts_by_dimension(
    complex_: mp.FilteredComplex,
    critical_simplices: tuple[int, ...],
) -> tuple[int, ...]:
    max_dimension = max(
        (complex_.dimension(simplex) for simplex in range(complex_.size)),
        default=-1,
    )
    counts = [0] * (max_dimension + 1)
    for simplex in critical_simplices:
        counts[complex_.dimension(simplex)] += 1
    return tuple(counts)


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


def benchmark_case(
    case: LowerStarCase,
    *,
    workers: Iterable[int] = (1, 2, 4),
    repeats: int = 3,
    warmups: int = 1,
) -> list[ProcessLowerStarsBenchmarkRow]:
    if repeats < 1:
        raise ValueError("repeats must be positive.")
    if warmups < 0:
        raise ValueError("warmups must be non-negative.")
    worker_counts = tuple(dict.fromkeys(int(worker) for worker in workers))
    if not worker_counts or any(worker < 1 for worker in worker_counts):
        raise ValueError("Worker counts must be positive.")

    complex_ = case.complex
    standard = mp.compute_standard_persistence(complex_)
    standard_finite = standard.finite_barcode()
    standard_essential = standard.essential_barcode()
    baseline = mp.compute_morse_sequence(
        complex_, algorithm=mp.PROCESS_LOWER_STARS_SEQUENCE
    )

    configurations = [(mp.PROCESS_LOWER_STARS_SEQUENCE, 1, None)]
    configurations.extend(
        (mp.PROCESS_LOWER_STARS_PARALLEL_SEQUENCE, worker, worker)
        for worker in worker_counts
    )
    rows: list[ProcessLowerStarsBenchmarkRow] = []
    baseline_sequence_seconds: float | None = None
    baseline_total_seconds: float | None = None
    max_dimension = max(
        (complex_.dimension(simplex) for simplex in range(complex_.size)),
        default=-1,
    )

    for algorithm, reported_workers, max_workers in configurations:
        for _ in range(warmups):
            _run_pipeline(complex_, algorithm, max_workers)

        best: tuple[mp.MorseSequence, float, float, float] | None = None
        for _ in range(repeats):
            sequence, diagram, sequence_seconds, persistence_seconds, total_seconds = (
                _run_pipeline(complex_, algorithm, max_workers)
            )
            if sequence.steps != baseline.steps:
                raise AssertionError(
                    f"{algorithm} with {reported_workers} workers differs from sequential"
                )
            if (
                diagram.finite_barcode() != standard_finite
                or diagram.essential_barcode() != standard_essential
            ):
                raise AssertionError("Morse and standard persistence barcodes differ.")
            if best is None or total_seconds < best[3]:
                best = (
                    sequence,
                    sequence_seconds,
                    persistence_seconds,
                    total_seconds,
                )

        if best is None:
            raise RuntimeError("Benchmark did not run.")
        sequence, sequence_seconds, persistence_seconds, total_seconds = best
        if algorithm == mp.PROCESS_LOWER_STARS_SEQUENCE:
            baseline_sequence_seconds = sequence_seconds
            baseline_total_seconds = total_seconds
        if baseline_sequence_seconds is None or baseline_total_seconds is None:
            raise RuntimeError("The sequential baseline must run first.")
        sequence_speedup = (
            baseline_sequence_seconds / sequence_seconds
            if sequence_seconds
            else math.inf
        )
        total_speedup = (
            baseline_total_seconds / total_seconds if total_seconds else math.inf
        )
        critical_count = len(sequence.critical_simplices)
        eliminated = complex_.size - critical_count
        loads = estimated_task_loads(case.lower_star_sizes, reported_workers)
        rows.append(
            ProcessLowerStarsBenchmarkRow(
                case=case.name,
                distribution=case.distribution,
                anchor_count=len(case.fan_sizes),
                min_anchor_fan_size=min(case.fan_sizes),
                max_anchor_fan_size=max(case.fan_sizes),
                total_anchor_fan_size=sum(case.fan_sizes),
                algorithm=algorithm,
                max_workers=reported_workers,
                repeats=repeats,
                warmups=warmups,
                cpp_backend=complex_.cpp_backend_active(),
                num_simplices=complex_.size,
                num_vertices=len(case.lower_star_sizes),
                max_dimension=max_dimension,
                lower_star_count=len(case.lower_star_sizes),
                max_lower_star_size=max(case.lower_star_sizes),
                mean_lower_star_size=(
                    float(sum(case.lower_star_sizes)) / len(case.lower_star_sizes)
                ),
                estimated_min_task_load=min(loads),
                estimated_max_task_load=max(loads),
                estimated_task_load_ratio=(
                    float(max(loads)) / min(loads) if min(loads) else math.inf
                ),
                num_critical_simplices=critical_count,
                critical_simplices_by_dimension=_critical_counts_by_dimension(
                    complex_, sequence.critical_simplices
                ),
                num_regular_pairs=(complex_.size - critical_count) // 2,
                critical_ratio=float(critical_count) / complex_.size,
                eliminated_simplices=eliminated,
                exact_sequence_matches_sequential=True,
                sequence_seconds=sequence_seconds,
                persistence_seconds=persistence_seconds,
                total_seconds=total_seconds,
                sequence_seconds_per_eliminated_simplex=(
                    sequence_seconds / eliminated if eliminated else math.inf
                ),
                sequence_speedup_vs_sequential=sequence_speedup,
                total_speedup_vs_sequential=total_speedup,
                sequence_parallel_efficiency=(
                    sequence_speedup / reported_workers
                ),
            )
        )
    return rows


def write_rows(
    rows: list[ProcessLowerStarsBenchmarkRow],
    output: TextIO,
    output_format: str,
) -> None:
    if not rows:
        raise ValueError("No benchmark rows were generated.")
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
    parser.add_argument("--anchors", type=int, default=16)
    parser.add_argument("--balanced-fan", type=int, default=8)
    parser.add_argument("--light-fan", type=int, default=2)
    parser.add_argument("--heavy-fan", type=int, default=98)
    parser.add_argument("--workers", type=int, nargs="+", default=(1, 2, 4, 8))
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--format", choices=("csv", "json"), default="csv")
    parser.add_argument("--output", type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    cases = (
        make_balanced_case(args.anchors, args.balanced_fan),
        make_skewed_case(args.anchors, args.light_fan, args.heavy_fan),
    )
    rows = [
        row
        for case in cases
        for row in benchmark_case(
            case,
            workers=args.workers,
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
