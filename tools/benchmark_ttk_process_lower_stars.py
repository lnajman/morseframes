#!/usr/bin/env python3
"""Compare TTK and MorseFrames ProcessLowerStars on identical 2D/3D inputs."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "tools"))

import benchmark_simplicial_strategies as generators  # noqa: E402
import morseframes as mf  # noqa: E402

TTK_REVISION = "f4ffd1a1049d0ccf6e8f3eb4f7c096a6cc251ba0"


@dataclass(frozen=True)
class ComparisonRow:
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
    f_max_critical_simplices: int
    f_max_critical_by_dimension: tuple[int, ...]
    morseframes_critical_simplices: int
    morseframes_critical_by_dimension: tuple[int, ...]
    ttk_critical_simplices: int
    ttk_critical_by_dimension: tuple[int, ...]
    critical_counts_match: bool
    f_max_seconds: float
    morseframes_process_lower_stars_seconds: float
    ttk_process_lower_stars_seconds: float
    morseframes_ratio_vs_f_max: float
    ttk_ratio_vs_f_max: float
    morseframes_ratio_vs_ttk: float
    ttk_setup_seconds: float
    ttk_precondition_seconds: float


def best_profile(
    complex_: mf.FilteredComplex,
    algorithm: str,
    workers: int | None,
    repeats: int,
    warmups: int,
) -> mf.MorseSequenceProfile:
    for _ in range(warmups):
        mf.profile_morse_sequence(complex_, algorithm=algorithm, max_workers=workers)
    profiles = [
        mf.profile_morse_sequence(complex_, algorithm=algorithm, max_workers=workers)
        for _ in range(repeats)
    ]
    return min(profiles, key=lambda profile: profile.construction_seconds)


def write_ttk_input(complex_: mf.FilteredComplex, path: Path) -> tuple[int, int]:
    records = tuple(complex_.simplices())
    dimension = max(record.dimension for record in records)
    if dimension not in (2, 3):
        raise ValueError("The TTK comparison supports 2D and 3D complexes.")
    vertices = [record for record in records if record.dimension == 0]
    vertex_ids = sorted(record.vertices[0] for record in vertices)
    if vertex_ids != list(range(len(vertex_ids))):
        raise ValueError("TTK benchmark inputs require contiguous vertex identifiers.")
    values = [0.0] * len(vertices)
    for record in vertices:
        values[record.vertices[0]] = record.filtration
    facets = [record.vertices for record in records if record.dimension == dimension]

    with path.open("w", encoding="utf-8") as output:
        output.write(
            f"morseframes-ttk-v1 {dimension} {len(vertices)} {len(facets)}\n"
        )
        output.write(" ".join(format(value, ".17g") for value in values) + "\n")
        for facet in facets:
            output.write(" ".join(str(vertex) for vertex in facet) + "\n")
    return len(vertices), len(records)


def run_ttk(
    executable: Path,
    input_path: Path,
    workers: int,
    repeats: int,
    warmups: int,
) -> dict[str, object]:
    completed = subprocess.run(
        [
            str(executable),
            "--input",
            str(input_path),
            "--workers",
            str(workers),
            "--repeats",
            str(repeats),
            "--warmups",
            str(warmups),
        ],
        check=True,
        text=True,
        capture_output=True,
    )
    lines = [line for line in completed.stdout.splitlines() if line.startswith("{")]
    if not lines:
        raise RuntimeError(f"TTK benchmark produced no JSON row: {completed.stdout}")
    return json.loads(lines[-1])


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
) -> list[ComparisonRow]:
    if not mf.cpp_backend_active(complex_):
        raise RuntimeError("The comparison requires the MorseFrames native backend.")
    f_max = best_profile(
        complex_, mf.F_MAX_SEQUENCE, None, repeats=repeats, warmups=warmups
    )
    rows: list[ComparisonRow] = []
    with tempfile.TemporaryDirectory(prefix="morseframes-ttk-") as directory:
        input_path = Path(directory) / "complex.txt"
        num_vertices, num_simplices = write_ttk_input(complex_, input_path)
        for worker_count in workers:
            algorithm = (
                mf.PROCESS_LOWER_STARS_SEQUENCE
                if worker_count == 1
                else mf.PROCESS_LOWER_STARS_PARALLEL_SEQUENCE
            )
            morseframes = best_profile(
                complex_,
                algorithm,
                None if worker_count == 1 else worker_count,
                repeats=repeats,
                warmups=warmups,
            )
            ttk = run_ttk(
                executable, input_path, worker_count, repeats, warmups
            )
            if (
                int(ttk["dimension"]) != dimension
                or int(ttk["num_vertices"]) != num_vertices
                or int(ttk["num_simplices"]) != num_simplices
            ):
                raise AssertionError(
                    "TTK and MorseFrames did not construct the same complex: "
                    f"TTK=({ttk['dimension']}, {ttk['num_vertices']}, "
                    f"{ttk['num_simplices']}), MorseFrames=({dimension}, "
                    f"{num_vertices}, {num_simplices})"
                )
            ttk_critical = int(ttk["num_critical_simplices"])
            ttk_by_dimension = tuple(
                int(value) for value in ttk["critical_simplices_by_dimension"]
            )
            counts_match = (
                f_max.num_critical_simplices
                == morseframes.num_critical_simplices
                == ttk_critical
                and f_max.critical_simplices_by_dimension
                == morseframes.critical_simplices_by_dimension
                == ttk_by_dimension
            )
            if not counts_match:
                raise AssertionError(
                    "Critical counts differ: "
                    f"F-Max={f_max.critical_simplices_by_dimension}, "
                    f"MorseFrames={morseframes.critical_simplices_by_dimension}, "
                    f"TTK={ttk_by_dimension}"
                )
            ttk_seconds = float(ttk["gradient_seconds"])
            rows.append(
                ComparisonRow(
                    family=family,
                    dimension=dimension,
                    seed=seed,
                    grid_size=grid_size,
                    workers=worker_count,
                    repeats=repeats,
                    warmups=warmups,
                    ttk_revision=TTK_REVISION,
                    num_vertices=num_vertices,
                    num_simplices=num_simplices,
                    f_max_critical_simplices=f_max.num_critical_simplices,
                    f_max_critical_by_dimension=f_max.critical_simplices_by_dimension,
                    morseframes_critical_simplices=morseframes.num_critical_simplices,
                    morseframes_critical_by_dimension=(
                        morseframes.critical_simplices_by_dimension
                    ),
                    ttk_critical_simplices=ttk_critical,
                    ttk_critical_by_dimension=ttk_by_dimension,
                    critical_counts_match=counts_match,
                    f_max_seconds=f_max.construction_seconds,
                    morseframes_process_lower_stars_seconds=(
                        morseframes.construction_seconds
                    ),
                    ttk_process_lower_stars_seconds=ttk_seconds,
                    morseframes_ratio_vs_f_max=(
                        morseframes.construction_seconds / f_max.construction_seconds
                    ),
                    ttk_ratio_vs_f_max=ttk_seconds / f_max.construction_seconds,
                    morseframes_ratio_vs_ttk=(
                        morseframes.construction_seconds / ttk_seconds
                    ),
                    ttk_setup_seconds=float(ttk["setup_seconds"]),
                    ttk_precondition_seconds=float(ttk["precondition_seconds"]),
                )
            )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ttk-benchmark", type=Path, required=True)
    parser.add_argument(
        "--families", nargs="+", choices=("terrain", "volume"), default=("terrain", "volume")
    )
    parser.add_argument("--terrain-sizes", type=int, nargs="+", default=(16, 32, 64))
    parser.add_argument("--volume-sizes", type=int, nargs="+", default=(4, 8, 12, 16))
    parser.add_argument("--seeds", type=int, nargs="+", default=(0, 1, 2))
    parser.add_argument("--workers", type=int, nargs="+", default=(1, 2, 4, 8))
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.ttk_benchmark.is_file():
        raise FileNotFoundError(args.ttk_benchmark)
    if args.repeats < 1 or args.warmups < 0:
        raise ValueError("repeats must be positive and warmups nonnegative")
    workers = tuple(dict.fromkeys(args.workers))
    if not workers or any(worker < 1 for worker in workers):
        raise ValueError("workers must contain positive integers")

    rows: list[ComparisonRow] = []
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
