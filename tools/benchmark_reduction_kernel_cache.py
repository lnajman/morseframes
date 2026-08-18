#!/usr/bin/env python3
"""Benchmark explicit ReductionKernel cache construction and reuse."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "tools"))

import benchmark_simplicial_strategies as fixtures  # noqa: E402
import morseframes as mf  # noqa: E402


@dataclass(frozen=True)
class CacheBenchmarkRow:
    family: str
    name: str
    seed: int
    grid_size: int
    num_simplices: int
    num_levels: int
    repeats: int
    warmups: int
    cache_build_seconds: float
    cache_entries: int
    cache_coboundary_entries: int
    cache_bytes: int
    uncached_gradient_seconds: float
    cached_gradient_seconds: float
    cached_speedup: float
    saving_seconds: float
    break_even_gradients: float
    sequence_matches_uncached: bool


def _profile(complex_: mf.FilteredComplex) -> tuple[float, tuple[object, ...]]:
    profile = mf.profile_morse_sequence(
        complex_, algorithm=mf.FLOODING_REDUCTION_KERNEL_SEQUENCE
    )
    sequence = mf.compute_morse_sequence(
        complex_, algorithm=mf.FLOODING_REDUCTION_KERNEL_SEQUENCE
    )
    return profile.construction_seconds, sequence.steps


def benchmark_case(
    family: str,
    grid_size: int,
    seed: int,
    *,
    repeats: int,
    warmups: int,
) -> CacheBenchmarkRow:
    maker = (
        fixtures.make_injective_terrain
        if family == "terrain"
        else fixtures.make_injective_volume
    )
    uncached = maker(seed, grid_size)
    cached = maker(seed, grid_size)
    if not uncached.cpp_backend_active() or not cached.cpp_backend_active():
        raise RuntimeError("The cache benchmark requires the native backend.")

    cache = cached.prepare_reduction_kernel_cache()
    for _ in range(warmups):
        _profile(uncached)
        _profile(cached)

    uncached_times: list[float] = []
    cached_times: list[float] = []
    uncached_steps: tuple[object, ...] = ()
    cached_steps: tuple[object, ...] = ()
    for repeat in range(repeats):
        if repeat % 2 == 0:
            uncached_time, uncached_steps = _profile(uncached)
            cached_time, cached_steps = _profile(cached)
        else:
            cached_time, cached_steps = _profile(cached)
            uncached_time, uncached_steps = _profile(uncached)
        uncached_times.append(uncached_time)
        cached_times.append(cached_time)

    uncached_median = statistics.median(uncached_times)
    cached_median = statistics.median(cached_times)
    saving = uncached_median - cached_median
    cache_build_seconds = cache["build_nanoseconds"] * 1.0e-9
    return CacheBenchmarkRow(
        family=f"injective-{family}",
        name=f"injective-{family}-n{grid_size}-seed{seed}",
        seed=seed,
        grid_size=grid_size,
        num_simplices=uncached.size,
        num_levels=uncached.num_levels,
        repeats=repeats,
        warmups=warmups,
        cache_build_seconds=cache_build_seconds,
        cache_entries=cache["entries"],
        cache_coboundary_entries=cache["coboundary_entries"],
        cache_bytes=cache["bytes"],
        uncached_gradient_seconds=uncached_median,
        cached_gradient_seconds=cached_median,
        cached_speedup=(uncached_median / cached_median),
        saving_seconds=saving,
        break_even_gradients=(
            cache_build_seconds / saving if saving > 0 else math.inf
        ),
        sequence_matches_uncached=(cached_steps == uncached_steps),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--terrain-sizes", nargs="+", type=int, default=[16, 32, 64])
    parser.add_argument("--volume-sizes", nargs="+", type=int, default=[4, 8, 12, 16])
    parser.add_argument("--seeds", nargs="+", type=int, default=[0, 1, 2])
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.repeats < 1 or args.warmups < 0:
        parser.error("repeats must be positive and warmups nonnegative")

    rows = [
        benchmark_case(family, size, seed, repeats=args.repeats, warmups=args.warmups)
        for family, sizes in (
            ("terrain", args.terrain_sizes),
            ("volume", args.volume_sizes),
        )
        for size in sizes
        for seed in args.seeds
    ]
    if not all(row.sequence_matches_uncached for row in rows):
        raise AssertionError("A cached gradient differs from its uncached gradient.")

    fieldnames = list(asdict(rows[0]))
    stream = args.output.open("w", newline="") if args.output else sys.stdout
    try:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)
    finally:
        if args.output:
            stream.close()


if __name__ == "__main__":
    main()
