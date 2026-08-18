#!/usr/bin/env python3
"""Profile parallel phase costs on injective tetrahedral volumes."""

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

import benchmark_simplicial_strategies as strategy_benchmark  # noqa: E402
import morseframes as mp  # noqa: E402


STRATEGIES = {
    "process-lower-stars": mp.PROCESS_LOWER_STARS_PARALLEL_SEQUENCE,
    "reduction-kernel": mp.FLOODING_REDUCTION_KERNEL_PARALLEL_SEQUENCE,
}


@dataclass(frozen=True)
class TetrahedralPhaseProfileRow:
    family: str
    name: str
    seed: int
    grid_size: int
    strategy: str
    algorithm: str
    max_workers: int
    repeats: int
    num_simplices: int
    num_levels: int
    num_critical_simplices: int
    profile_seconds: float
    sequence_total_seconds: float
    sequence_core_seconds: float
    callback_seconds: float
    callback_share: float
    process_lower_stars_builder_init_seconds: float
    process_lower_stars_builder_init_share: float
    process_lower_stars_setup_seconds: float
    process_lower_stars_setup_share: float
    process_lower_stars_local_wall_seconds: float
    process_lower_stars_local_wall_share: float
    process_lower_stars_replay_seconds: float
    process_lower_stars_replay_share: float
    process_lower_stars_cumulative_task_seconds: float
    process_lower_stars_task_parallelism: float
    process_lower_stars_task_time_imbalance: float
    process_lower_stars_task_load_imbalance: float
    process_lower_stars_parallel_tasks: int
    reduction_kernel_max_parallel_levels: int


def _seconds(metrics: dict[str, object], name: str) -> float:
    return 1.0e-9 * int(metrics.get(name, 0))


def _ratio(numerator: float, denominator: float) -> float:
    return numerator / denominator if denominator else math.nan


def benchmark_profile(
    seed: int,
    grid_size: int,
    *,
    workers: tuple[int, ...] = (1, 8),
    strategies: tuple[str, ...] = tuple(STRATEGIES),
    repeats: int = 5,
) -> list[TetrahedralPhaseProfileRow]:
    selected_workers = tuple(dict.fromkeys(workers))
    selected_strategies = tuple(dict.fromkeys(strategies))
    if not selected_workers or any(worker < 1 for worker in selected_workers):
        raise ValueError("workers must contain positive integers")
    if repeats < 1:
        raise ValueError("repeats must be positive")
    unknown = set(selected_strategies) - STRATEGIES.keys()
    if unknown:
        raise ValueError(f"Unknown strategies: {', '.join(sorted(unknown))}")

    complex_ = strategy_benchmark.make_injective_volume(seed, grid_size)
    rows: list[TetrahedralPhaseProfileRow] = []
    for strategy in selected_strategies:
        algorithm = STRATEGIES[strategy]
        for worker_count in selected_workers:
            profiles = [
                mp.profile_morse_reference_frame(
                    complex_, algorithm=algorithm, max_workers=worker_count
                )
                for _ in range(repeats)
            ]
            profile = min(profiles, key=lambda candidate: candidate.profile_seconds)
            metrics = profile.frame_metrics
            sequence_total = _seconds(metrics, "sequence_total_nanoseconds")
            sequence_core = _seconds(metrics, "sequence_core_nanoseconds")
            callback = sum(
                _seconds(metrics, name)
                for name in (
                    "reference_update_nanoseconds",
                    "reduction_plan_nanoseconds",
                    "release_cleanup_nanoseconds",
                )
            )
            builder_init = _seconds(
                metrics, "sequence_process_lower_stars_builder_init_nanoseconds"
            )
            setup = _seconds(metrics, "sequence_process_lower_stars_setup_nanoseconds")
            local_wall = _seconds(
                metrics, "sequence_process_lower_stars_local_wall_nanoseconds"
            )
            replay = _seconds(
                metrics, "sequence_process_lower_stars_replay_nanoseconds"
            )
            cumulative_task = _seconds(
                metrics,
                "sequence_process_lower_stars_cumulative_task_nanoseconds",
            )
            min_task = _seconds(
                metrics, "sequence_process_lower_stars_min_task_nanoseconds"
            )
            max_task = _seconds(
                metrics, "sequence_process_lower_stars_max_task_nanoseconds"
            )
            min_load = int(metrics.get("sequence_process_lower_stars_min_task_load", 0))
            max_load = int(metrics.get("sequence_process_lower_stars_max_task_load", 0))
            is_process_lower_stars = strategy == "process-lower-stars"
            rows.append(
                TetrahedralPhaseProfileRow(
                    family="injective-volume",
                    name=f"injective-volume-n{grid_size}-seed{seed}",
                    seed=seed,
                    grid_size=grid_size,
                    strategy=strategy,
                    algorithm=algorithm,
                    max_workers=worker_count,
                    repeats=repeats,
                    num_simplices=profile.num_simplices,
                    num_levels=profile.num_levels,
                    num_critical_simplices=profile.num_critical_simplices,
                    profile_seconds=profile.profile_seconds,
                    sequence_total_seconds=sequence_total,
                    sequence_core_seconds=sequence_core,
                    callback_seconds=callback,
                    callback_share=_ratio(callback, sequence_total),
                    process_lower_stars_builder_init_seconds=builder_init,
                    process_lower_stars_builder_init_share=(
                        _ratio(builder_init, sequence_total)
                        if is_process_lower_stars
                        else math.nan
                    ),
                    process_lower_stars_setup_seconds=setup,
                    process_lower_stars_setup_share=(
                        _ratio(setup, sequence_total)
                        if is_process_lower_stars
                        else math.nan
                    ),
                    process_lower_stars_local_wall_seconds=local_wall,
                    process_lower_stars_local_wall_share=(
                        _ratio(local_wall, sequence_total)
                        if is_process_lower_stars
                        else math.nan
                    ),
                    process_lower_stars_replay_seconds=replay,
                    process_lower_stars_replay_share=(
                        _ratio(replay, sequence_total)
                        if is_process_lower_stars
                        else math.nan
                    ),
                    process_lower_stars_cumulative_task_seconds=cumulative_task,
                    process_lower_stars_task_parallelism=(
                        _ratio(cumulative_task, local_wall)
                        if is_process_lower_stars
                        else math.nan
                    ),
                    process_lower_stars_task_time_imbalance=(
                        _ratio(max_task, min_task)
                        if is_process_lower_stars
                        else math.nan
                    ),
                    process_lower_stars_task_load_imbalance=(
                        _ratio(float(max_load), float(min_load))
                        if is_process_lower_stars
                        else math.nan
                    ),
                    process_lower_stars_parallel_tasks=int(
                        metrics.get("sequence_process_lower_stars_parallel_tasks", 0)
                    ),
                    reduction_kernel_max_parallel_levels=int(
                        metrics.get("sequence_reduction_kernel_max_parallel_levels", 0)
                    ),
                )
            )
    return rows


def write_rows(
    rows: list[TetrahedralPhaseProfileRow], output: TextIO, output_format: str
) -> None:
    if not rows:
        raise ValueError("No profile rows were generated")
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
    parser.add_argument("--workers", type=int, nargs="+", default=(1, 8))
    parser.add_argument(
        "--strategies", nargs="+", choices=tuple(STRATEGIES), default=tuple(STRATEGIES)
    )
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--format", choices=("csv", "json"), default="csv")
    parser.add_argument("--output", type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    rows = [
        row
        for grid_size in args.sizes
        for seed in args.seeds
        for row in benchmark_profile(
            seed,
            grid_size,
            workers=tuple(args.workers),
            strategies=tuple(args.strategies),
            repeats=args.repeats,
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
