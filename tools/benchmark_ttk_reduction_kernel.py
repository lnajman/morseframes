#!/usr/bin/env python3
"""Directly compare TTK ProcessLowerStars with parallel ReductionKernel."""

from __future__ import annotations

import argparse
import csv
import sys
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "tools"))

import benchmark_simplicial_strategies as generators  # noqa: E402
import benchmark_ttk_process_lower_stars as ttk_support  # noqa: E402
import morseframes as mf  # noqa: E402


@dataclass(frozen=True)
class DirectComparisonRow:
    family: str
    dimension: int
    seed: int
    grid_size: int
    workers: int
    repeats: int
    warmups: int
    ttk_revision: str
    num_vertices: int
    num_simplices: int
    critical_simplices_by_dimension: tuple[int, ...]
    critical_counts_match: bool
    f_max_seconds: float
    reduction_kernel_seconds: float
    ttk_process_lower_stars_seconds: float
    reduction_kernel_ratio_vs_f_max: float
    ttk_ratio_vs_f_max: float
    reduction_kernel_ratio_vs_ttk: float
    reduction_kernel_speedup_vs_one: float
    ttk_speedup_vs_one: float


def benchmark_case(
    executable: Path,
    complex_: mf.FilteredComplex,
    family: str,
    dimension: int,
    seed: int,
    grid_size: int,
    workers: tuple[int, ...],
    repeats: int,
    warmups: int,
) -> list[DirectComparisonRow]:
    raw_rows: list[tuple[int, dict[str, object]]] = []
    with tempfile.TemporaryDirectory(prefix="morseframes-ttk-rk-") as directory:
        input_path = Path(directory) / "complex.txt"
        num_vertices, num_simplices = ttk_support.write_ttk_input(
            complex_, input_path
        )
        for worker_count in workers:
            result = ttk_support.run_ttk(
                executable, input_path, worker_count, repeats, warmups
            )
            if (
                int(result["dimension"]) != dimension
                or int(result["num_vertices"]) != num_vertices
                or int(result["num_simplices"]) != num_simplices
            ):
                raise AssertionError("TTK and MorseFrames constructed different complexes")
            raw_rows.append((worker_count, result))

    one_worker = next(
        (result for worker_count, result in raw_rows if worker_count == 1), None
    )
    if one_worker is None:
        raise ValueError("workers must include the one-worker baseline")
    one_reduction_kernel = float(one_worker["reduction_kernel_seconds"])
    one_ttk = float(one_worker["gradient_seconds"])

    rows: list[DirectComparisonRow] = []
    for worker_count, result in raw_rows:
        f_max_seconds = float(result["f_max_seconds"])
        reduction_kernel_seconds = float(result["reduction_kernel_seconds"])
        ttk_seconds = float(result["gradient_seconds"])
        ttk_critical = tuple(
            int(value) for value in result["critical_simplices_by_dimension"]
        )
        f_max_critical = tuple(
            int(value)
            for value in result["f_max_critical_simplices_by_dimension"]
        )
        reduction_kernel_critical = tuple(
            int(value)
            for value in result[
                "reduction_kernel_critical_simplices_by_dimension"
            ]
        )
        counts_match = ttk_critical == f_max_critical == reduction_kernel_critical
        if not counts_match:
            raise AssertionError(
                "Critical counts differ: "
                f"TTK={ttk_critical}, F-Max={f_max_critical}, "
                f"ReductionKernel={reduction_kernel_critical}"
            )
        rows.append(
            DirectComparisonRow(
                family=family,
                dimension=dimension,
                seed=seed,
                grid_size=grid_size,
                workers=worker_count,
                repeats=repeats,
                warmups=warmups,
                ttk_revision=ttk_support.TTK_REVISION,
                num_vertices=num_vertices,
                num_simplices=num_simplices,
                critical_simplices_by_dimension=ttk_critical,
                critical_counts_match=counts_match,
                f_max_seconds=f_max_seconds,
                reduction_kernel_seconds=reduction_kernel_seconds,
                ttk_process_lower_stars_seconds=ttk_seconds,
                reduction_kernel_ratio_vs_f_max=(
                    reduction_kernel_seconds / f_max_seconds
                ),
                ttk_ratio_vs_f_max=ttk_seconds / f_max_seconds,
                reduction_kernel_ratio_vs_ttk=(
                    reduction_kernel_seconds / ttk_seconds
                ),
                reduction_kernel_speedup_vs_one=(
                    one_reduction_kernel / reduction_kernel_seconds
                ),
                ttk_speedup_vs_one=one_ttk / ttk_seconds,
            )
        )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ttk-benchmark", type=Path, required=True)
    parser.add_argument(
        "--families",
        nargs="+",
        choices=("terrain", "volume"),
        default=("terrain", "volume"),
    )
    parser.add_argument("--terrain-sizes", type=int, nargs="+", default=(16, 32, 64))
    parser.add_argument("--volume-sizes", type=int, nargs="+", default=(4, 8, 12, 16))
    parser.add_argument("--seeds", type=int, nargs="+", default=(0, 1, 2))
    parser.add_argument("--workers", type=int, nargs="+", default=(1, 2, 4, 8))
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.ttk_benchmark.is_file():
        raise FileNotFoundError(args.ttk_benchmark)
    workers = tuple(dict.fromkeys(args.workers))
    if 1 not in workers or any(worker < 1 for worker in workers):
        raise ValueError("workers must contain positive integers including one")
    if args.repeats < 1 or args.warmups < 0:
        raise ValueError("repeats must be positive and warmups nonnegative")

    rows: list[DirectComparisonRow] = []
    if "terrain" in args.families:
        for grid_size in args.terrain_sizes:
            for seed in args.seeds:
                rows.extend(
                    benchmark_case(
                        args.ttk_benchmark,
                        generators.make_injective_terrain(seed, grid_size),
                        "injective-terrain",
                        2,
                        seed,
                        grid_size,
                        workers,
                        args.repeats,
                        args.warmups,
                    )
                )
    if "volume" in args.families:
        for grid_size in args.volume_sizes:
            for seed in args.seeds:
                rows.extend(
                    benchmark_case(
                        args.ttk_benchmark,
                        generators.make_injective_volume(seed, grid_size),
                        "injective-volume",
                        3,
                        seed,
                        grid_size,
                        workers,
                        args.repeats,
                        args.warmups,
                    )
                )

    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
    output = sys.stdout if args.output is None else args.output.open("w", newline="")
    try:
        writer = csv.DictWriter(output, fieldnames=list(asdict(rows[0])))
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)
    finally:
        if output is not sys.stdout:
            output.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
