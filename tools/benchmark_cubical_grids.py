#!/usr/bin/env python3
"""Benchmark MorseFrames 2D cubical grids against GUDHI CubicalComplex."""

from __future__ import annotations

import argparse
import csv
import math
import random
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from time import perf_counter
from typing import Callable, Iterable, Sequence, TypeVar


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

import morseframes as mp  # noqa: E402


DEFAULT_FAMILIES = ("sinusoidal", "noisy-sinusoidal", "plateau", "random")
DEFAULT_SIZES = (16, 32)
DEFAULT_SEEDS = (0, 1)
DEFAULT_ALGORITHMS = (
    mp.F_MAX_SEQUENCE,
    mp.F_MIN_SEQUENCE,
    mp.SATURATED_SEQUENCE,
    mp.COREDUCTION_SEQUENCE,
)

T = TypeVar("T")
BarcodeSignature = tuple[
    tuple[tuple[int, float, float], ...],
    tuple[tuple[int, float], ...],
]


@dataclass(frozen=True)
class CubicalBenchmarkRow:
    family: str
    name: str
    seed: int
    width: int
    height: int
    num_cells: int
    num_levels: int
    algorithm: str
    critical_cells: int
    critical_ratio: float
    modulus: int
    repeats: int
    build_seconds: float
    sequence_seconds: float
    morse_reference_seconds: float
    morse_coreference_seconds: float
    standard_seconds: float
    gudhi_seconds: float | None
    morse_reference_total_seconds: float
    morse_coreference_total_seconds: float
    standard_total_seconds: float
    gudhi_over_morse_reference: float | None
    gudhi_over_morse_coreference: float | None
    gudhi_over_standard: float | None
    standard_over_morse_reference: float
    standard_over_morse_coreference: float
    finite_count: int
    essential_count: int
    matches_standard: bool
    matches_coreference: bool
    matches_gudhi: bool | None


def time_best(action: Callable[[], T], *, repeats: int) -> tuple[float, T]:
    best_seconds = float("inf")
    best_result: T | None = None
    for _ in range(repeats):
        started = perf_counter()
        result = action()
        elapsed = perf_counter() - started
        if elapsed < best_seconds:
            best_seconds = elapsed
            best_result = result
    if best_result is None:
        raise RuntimeError("No benchmark repetitions were run.")
    return best_seconds, best_result


def diagram_signature(diagram: mp.PersistenceDiagram) -> BarcodeSignature:
    return diagram.finite_barcode(include_zero=False), diagram.essential_barcode()


def _normalize(values: Sequence[float]) -> list[float]:
    lower = min(values)
    upper = max(values)
    if upper == lower:
        return [0.0 for _ in values]
    scale = upper - lower
    return [(value - lower) / scale for value in values]


def _base_sinusoidal(width: int, height: int, *, seed: int = 0) -> list[float]:
    rng = random.Random(seed)
    phase_x = rng.random() * 2.0 * math.pi
    phase_y = rng.random() * 2.0 * math.pi
    values: list[float] = []
    for y in range(height):
        yy = (2.0 * math.pi * y) / float(max(1, height - 1))
        for x in range(width):
            xx = (2.0 * math.pi * x) / float(max(1, width - 1))
            value = math.sin(xx) * math.sin(yy)
            value += 0.18 * math.cos(2.0 * xx + yy)
            value += 0.08 * math.sin(3.0 * xx + phase_x) * math.cos(yy + phase_y)
            values.append(value)
    return values


def make_vertex_values(
    family: str,
    *,
    width: int,
    height: int,
    seed: int,
    plateau_levels: int,
) -> tuple[list[float], str]:
    if width < 2 or height < 2:
        raise ValueError("width and height must be at least 2.")

    rng = random.Random(seed)
    if family == "sinusoidal":
        values = [round(value, 8) for value in _base_sinusoidal(width, height, seed=seed)]
    elif family == "noisy-sinusoidal":
        values = [
            round(value + 0.10 * rng.uniform(-1.0, 1.0), 8)
            for value in _base_sinusoidal(width, height, seed=seed)
        ]
    elif family == "plateau":
        if plateau_levels < 2:
            raise ValueError("plateau_levels must be at least 2.")
        normalized = _normalize(_base_sinusoidal(width, height, seed=seed))
        values = [
            round(round(value * float(plateau_levels - 1)) / float(plateau_levels - 1), 8)
            for value in normalized
        ]
    elif family == "random":
        values = [round(rng.random(), 8) for _ in range(width * height)]
    else:
        raise ValueError(f"Unknown cubical benchmark family {family!r}.")
    return values, f"{family}-g{width}x{height}-seed{seed}"


def gudhi_available() -> bool:
    try:
        import gudhi  # noqa: F401
    except ImportError:
        return False
    return True


def gudhi_cubical_barcode_from_values(
    width: int,
    height: int,
    values: Sequence[float],
    *,
    modulus: int,
) -> BarcodeSignature:
    try:
        import gudhi  # type: ignore[import-not-found]
    except ImportError as exc:
        raise RuntimeError("GUDHI is not importable in this Python environment.") from exc

    vertices = [
        [float(values[y * width + x]) for x in range(width)]
        for y in range(height)
    ]
    cubical = gudhi.CubicalComplex(vertices=vertices)
    persistence = cubical.persistence(homology_coeff_field=modulus, min_persistence=0)
    finite: list[tuple[int, float, float]] = []
    essential: list[tuple[int, float]] = []
    for dimension, (birth, death) in persistence:
        if death == math.inf:
            essential.append((int(dimension), float(birth)))
        elif birth < death:
            finite.append((int(dimension), float(birth), float(death)))
    return tuple(sorted(finite)), tuple(sorted(essential))


def _optional_ratio(numerator: float | None, denominator: float) -> float | None:
    if numerator is None or denominator == 0.0:
        return None
    return numerator / denominator


def benchmark_case(
    *,
    family: str,
    width: int,
    height: int,
    seed: int,
    algorithms: Sequence[str],
    modulus: int,
    repeats: int,
    plateau_levels: int,
    time_gudhi: bool,
) -> list[CubicalBenchmarkRow]:
    values, name = make_vertex_values(
        family,
        width=width,
        height=height,
        seed=seed,
        plateau_levels=plateau_levels,
    )
    build_seconds, grid = time_best(
        lambda: mp.CubicalGrid2DComplex.from_vertex_values(width, height, values),
        repeats=repeats,
    )

    standard_seconds, standard = time_best(
        lambda: mp.compute_standard_persistence(grid, modulus=modulus),
        repeats=repeats,
    )
    standard_signature = diagram_signature(standard)
    gudhi_seconds: float | None = None
    gudhi_signature: BarcodeSignature | None = None
    if time_gudhi:
        gudhi_seconds, gudhi_signature = time_best(
            lambda: gudhi_cubical_barcode_from_values(
                width,
                height,
                values,
                modulus=modulus,
            ),
            repeats=repeats,
        )
        if gudhi_signature != standard_signature:
            raise RuntimeError(f"GUDHI cubical barcode mismatch on {name}.")

    rows: list[CubicalBenchmarkRow] = []
    for algorithm in algorithms:
        sequence_seconds, sequence = time_best(
            lambda: mp.compute_morse_sequence(grid, algorithm=algorithm),
            repeats=repeats,
        )
        reference_seconds, reference = time_best(
            lambda: mp.compute_morse_persistence(
                grid,
                sequence,
                algorithm=algorithm,
                modulus=modulus,
            ),
            repeats=repeats,
        )
        coreference_seconds, coreference = time_best(
            lambda: mp.compute_morse_coreference_persistence(
                grid,
                sequence,
                algorithm=algorithm,
                modulus=modulus,
            ),
            repeats=repeats,
        )
        reference_signature = diagram_signature(reference)
        coreference_signature = diagram_signature(coreference)
        matches_standard = reference_signature == standard_signature
        matches_coreference = coreference_signature == standard_signature
        if not matches_standard:
            raise RuntimeError(f"Morse reference barcode mismatch on {name} with {algorithm}.")
        if not matches_coreference:
            raise RuntimeError(f"Morse coreference barcode mismatch on {name} with {algorithm}.")

        reference_total = build_seconds + sequence_seconds + reference_seconds
        coreference_total = build_seconds + sequence_seconds + coreference_seconds
        standard_total = build_seconds + standard_seconds
        critical_cells = len(sequence.critical_simplices)
        rows.append(
            CubicalBenchmarkRow(
                family=family,
                name=name,
                seed=seed,
                width=width,
                height=height,
                num_cells=grid.size,
                num_levels=grid.num_levels,
                algorithm=algorithm,
                critical_cells=critical_cells,
                critical_ratio=critical_cells / float(grid.size) if grid.size else 0.0,
                modulus=modulus,
                repeats=repeats,
                build_seconds=build_seconds,
                sequence_seconds=sequence_seconds,
                morse_reference_seconds=reference_seconds,
                morse_coreference_seconds=coreference_seconds,
                standard_seconds=standard_seconds,
                gudhi_seconds=gudhi_seconds,
                morse_reference_total_seconds=reference_total,
                morse_coreference_total_seconds=coreference_total,
                standard_total_seconds=standard_total,
                gudhi_over_morse_reference=_optional_ratio(gudhi_seconds, reference_total),
                gudhi_over_morse_coreference=_optional_ratio(gudhi_seconds, coreference_total),
                gudhi_over_standard=_optional_ratio(gudhi_seconds, standard_total),
                standard_over_morse_reference=standard_total / reference_total
                if reference_total
                else math.inf,
                standard_over_morse_coreference=standard_total / coreference_total
                if coreference_total
                else math.inf,
                finite_count=len(standard.finite_pairs),
                essential_count=len(standard.essential),
                matches_standard=matches_standard,
                matches_coreference=matches_coreference,
                matches_gudhi=(
                    None if gudhi_signature is None else gudhi_signature == standard_signature
                ),
            )
        )
    return rows


def run_benchmarks(
    *,
    families: Sequence[str] = DEFAULT_FAMILIES,
    sizes: Sequence[int] = DEFAULT_SIZES,
    seeds: Sequence[int] = DEFAULT_SEEDS,
    algorithms: Sequence[str] = DEFAULT_ALGORITHMS,
    modulus: int = 3,
    repeats: int = 3,
    plateau_levels: int = 8,
    time_gudhi: bool = True,
) -> list[CubicalBenchmarkRow]:
    if repeats < 1:
        raise ValueError("repeats must be positive.")
    if not mp.cpp_backend_available() or mp.CppCubicalGrid2DComplex is None:
        raise RuntimeError("The native C++ backend is required for cubical benchmarks.")
    if time_gudhi and not gudhi_available():
        raise RuntimeError("GUDHI is not importable; pass --skip-gudhi to omit it.")

    rows: list[CubicalBenchmarkRow] = []
    for family in families:
        for size in sizes:
            for seed in seeds:
                rows.extend(
                    benchmark_case(
                        family=family,
                        width=size,
                        height=size,
                        seed=seed,
                        algorithms=algorithms,
                        modulus=modulus,
                        repeats=repeats,
                        plateau_levels=plateau_levels,
                        time_gudhi=time_gudhi,
                    )
                )
    return rows


def _format_optional(value: float | None, *, digits: int = 3) -> str:
    if value is None:
        return "-"
    return f"{value:.{digits}f}"


def render_summary(rows: Iterable[CubicalBenchmarkRow]) -> str:
    lines = [
        "family size seed algorithm cells crit% ref(ms) coref(ms) std(ms) gudhi(ms) "
        "G/ref G/coref",
    ]
    for row in rows:
        lines.append(
            f"{row.family:17} "
            f"{row.width:4} "
            f"{row.seed:4} "
            f"{row.algorithm:20} "
            f"{row.num_cells:7} "
            f"{100.0 * row.critical_ratio:5.1f} "
            f"{1000.0 * row.morse_reference_total_seconds:7.2f} "
            f"{1000.0 * row.morse_coreference_total_seconds:9.2f} "
            f"{1000.0 * row.standard_total_seconds:7.2f} "
            f"{_format_optional(None if row.gudhi_seconds is None else 1000.0 * row.gudhi_seconds):>9} "
            f"{_format_optional(row.gudhi_over_morse_reference):>6} "
            f"{_format_optional(row.gudhi_over_morse_coreference):>7}"
        )
    return "\n".join(lines)


def write_csv(rows: Iterable[CubicalBenchmarkRow], output: Path | None) -> None:
    materialized = list(rows)
    if materialized:
        fieldnames = list(asdict(materialized[0]).keys())
    else:
        fieldnames = list(CubicalBenchmarkRow.__dataclass_fields__.keys())

    if output is None:
        handle = sys.stdout
        close = False
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        handle = output.open("w", newline="")
        close = True
    try:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in materialized:
            writer.writerow(asdict(row))
    finally:
        if close:
            handle.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--families",
        nargs="+",
        default=list(DEFAULT_FAMILIES),
        choices=DEFAULT_FAMILIES,
        help="cubical grid value families",
    )
    parser.add_argument("--sizes", nargs="+", type=int, default=list(DEFAULT_SIZES))
    parser.add_argument("--seeds", nargs="+", type=int, default=list(DEFAULT_SEEDS))
    parser.add_argument(
        "--algorithms",
        nargs="+",
        default=list(DEFAULT_ALGORITHMS),
        choices=mp.MORSE_SEQUENCE_ALGORITHMS,
        help="Morse sequence strategies",
    )
    parser.add_argument("--modulus", type=int, default=3, help="prime field modulus")
    parser.add_argument("--repeats", type=int, default=3, help="best-of repetitions")
    parser.add_argument("--plateau-levels", type=int, default=8)
    parser.add_argument(
        "--skip-gudhi",
        action="store_true",
        help="do not time GUDHI CubicalComplex",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--format", choices=("summary", "csv"), default="summary")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = run_benchmarks(
        families=tuple(args.families),
        sizes=tuple(args.sizes),
        seeds=tuple(args.seeds),
        algorithms=tuple(args.algorithms),
        modulus=args.modulus,
        repeats=args.repeats,
        plateau_levels=args.plateau_levels,
        time_gudhi=not args.skip_gudhi,
    )

    if args.format == "csv":
        write_csv(rows, args.output)
        return

    summary = render_summary(rows)
    if args.output is None:
        print(summary)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(summary)


if __name__ == "__main__":
    main()
