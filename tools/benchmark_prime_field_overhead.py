#!/usr/bin/env python3
"""Benchmark Morse persistence coefficient overhead on local synthetic cases."""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from time import perf_counter
from typing import Callable


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "tools"))

import morseframes as mp  # noqa: E402
from benchmark_persistence import DEFAULT_ROADMAP_CACHE, make_benchmark_complex  # noqa: E402


DEFAULT_FAMILIES = ("lower-star", "plateau", "rips")
DEFAULT_SIZES = (8, 12, 16)
DEFAULT_SEEDS = (0, 1)
DEFAULT_ALGORITHMS = (
    mp.SATURATED_SEQUENCE,
    mp.F_MAX_SEQUENCE,
    mp.COREDUCTION_SEQUENCE,
)
DEFAULT_PRIMES = (3, 5)


@dataclass(frozen=True)
class TimingResult:
    seconds: float
    finite_count: int
    essential_count: int
    signature: tuple[
        tuple[tuple[int, float, float], ...],
        tuple[tuple[int, float], ...],
    ]


@dataclass(frozen=True)
class PrimeFieldOverheadRow:
    family: str
    name: str
    seed: int
    size: int
    num_simplices: int
    num_levels: int
    algorithm: str
    critical_simplices: int
    modulus: int
    repeats: int
    sequence_seconds: float
    morse_seconds: float
    standard_seconds: float
    morse_over_z2: float
    standard_over_z2: float
    morse_speedup_vs_standard: float
    matches_standard: bool
    matches_z2: bool
    finite_count: int
    essential_count: int


def diagram_signature(diagram: mp.PersistenceDiagram) -> tuple[
    tuple[tuple[int, float, float], ...],
    tuple[tuple[int, float], ...],
]:
    return (diagram.finite_barcode(include_zero=False), diagram.essential_barcode())


def time_best(
    action: Callable[[], mp.PersistenceDiagram],
    *,
    repeats: int,
) -> TimingResult:
    best_seconds = float("inf")
    best_diagram: mp.PersistenceDiagram | None = None
    for _ in range(repeats):
        start = perf_counter()
        diagram = action()
        elapsed = perf_counter() - start
        if elapsed < best_seconds:
            best_seconds = elapsed
            best_diagram = diagram
    if best_diagram is None:
        raise RuntimeError("No benchmark repetitions were run.")
    return TimingResult(
        seconds=best_seconds,
        finite_count=len(best_diagram.finite_pairs),
        essential_count=len(best_diagram.essential),
        signature=diagram_signature(best_diagram),
    )


def time_sequence(
    complex_: mp.FilteredComplex,
    algorithm: str,
    *,
    repeats: int,
) -> tuple[float, mp.MorseSequence]:
    best_seconds = float("inf")
    best_sequence: mp.MorseSequence | None = None
    for _ in range(repeats):
        start = perf_counter()
        sequence = mp.compute_morse_sequence(complex_, algorithm=algorithm)
        elapsed = perf_counter() - start
        if elapsed < best_seconds:
            best_seconds = elapsed
            best_sequence = sequence
    if best_sequence is None:
        raise RuntimeError("No sequence benchmark repetitions were run.")
    return best_seconds, best_sequence


def benchmark_case(
    *,
    family: str,
    seed: int,
    size: int,
    algorithm: str,
    primes: tuple[int, ...],
    repeats: int,
    plateau_levels: int,
    roadmap_cache: Path,
    download_roadmap_data: bool,
) -> list[PrimeFieldOverheadRow]:
    complex_, name = make_benchmark_complex(
        family,
        seed,
        size,
        plateau_levels=plateau_levels,
        roadmap_cache=roadmap_cache,
        download_roadmap_data=download_roadmap_data,
    )
    sequence_seconds, sequence = time_sequence(complex_, algorithm, repeats=repeats)
    critical_simplices = len(sequence.critical_simplices)

    z2_morse = time_best(
        lambda: mp.compute_morse_persistence(complex_, sequence=sequence),
        repeats=repeats,
    )
    z2_standard = time_best(
        lambda: mp.compute_standard_persistence(complex_),
        repeats=repeats,
    )
    if z2_morse.signature != z2_standard.signature:
        raise RuntimeError(f"Z2 Morse/standard mismatch on {name} with {algorithm}.")

    rows = [
        PrimeFieldOverheadRow(
            family=family,
            name=name,
            seed=seed,
            size=size,
            num_simplices=complex_.size,
            num_levels=complex_.num_levels,
            algorithm=algorithm,
            critical_simplices=critical_simplices,
            modulus=2,
            repeats=repeats,
            sequence_seconds=sequence_seconds,
            morse_seconds=z2_morse.seconds,
            standard_seconds=z2_standard.seconds,
            morse_over_z2=1.0,
            standard_over_z2=1.0,
            morse_speedup_vs_standard=(
                z2_standard.seconds / z2_morse.seconds if z2_morse.seconds else float("inf")
            ),
            matches_standard=True,
            matches_z2=True,
            finite_count=z2_morse.finite_count,
            essential_count=z2_morse.essential_count,
        )
    ]

    for modulus in primes:
        morse = time_best(
            lambda modulus=modulus: mp.compute_morse_persistence(
                complex_, sequence=sequence, modulus=modulus
            ),
            repeats=repeats,
        )
        standard = time_best(
            lambda modulus=modulus: mp.compute_standard_persistence(complex_, modulus=modulus),
            repeats=repeats,
        )
        rows.append(
            PrimeFieldOverheadRow(
                family=family,
                name=name,
                seed=seed,
                size=size,
                num_simplices=complex_.size,
                num_levels=complex_.num_levels,
                algorithm=algorithm,
                critical_simplices=critical_simplices,
                modulus=modulus,
                repeats=repeats,
                sequence_seconds=sequence_seconds,
                morse_seconds=morse.seconds,
                standard_seconds=standard.seconds,
                morse_over_z2=morse.seconds / z2_morse.seconds if z2_morse.seconds else float("inf"),
                standard_over_z2=(
                    standard.seconds / z2_standard.seconds
                    if z2_standard.seconds
                    else float("inf")
                ),
                morse_speedup_vs_standard=(
                    standard.seconds / morse.seconds if morse.seconds else float("inf")
                ),
                matches_standard=morse.signature == standard.signature,
                matches_z2=morse.signature == z2_morse.signature,
                finite_count=morse.finite_count,
                essential_count=morse.essential_count,
            )
        )

    return rows


def summarize_markdown(rows: list[PrimeFieldOverheadRow]) -> str:
    lines = [
        "# Prime-field overhead benchmark",
        "",
        "Best-of-repeat timings. Morse timings exclude sequence construction; "
        "`sequence_seconds` is reported separately in the CSV.",
        "",
        "## Aggregate",
        "",
        "| modulus | cases | avg p-field / Z2 Morse | median p-field / Z2 Morse | "
        "avg p-field / Z2 standard | avg Morse speedup vs standard |",
        "|---:|---:|---:|---:|---:|---:|",
    ]

    for modulus in sorted({row.modulus for row in rows if row.modulus != 2}):
        modulus_rows = [row for row in rows if row.modulus == modulus]
        morse_over_z2 = [row.morse_over_z2 for row in modulus_rows]
        standard_over_z2 = [row.standard_over_z2 for row in modulus_rows]
        speedup = [row.morse_speedup_vs_standard for row in modulus_rows]
        lines.append(
            f"| F_{modulus} | {len(modulus_rows)} | "
            f"{statistics.fmean(morse_over_z2):.3f} | "
            f"{statistics.median(morse_over_z2):.3f} | "
            f"{statistics.fmean(standard_over_z2):.3f} | "
            f"{statistics.fmean(speedup):.3f}x |"
        )

    lines.extend([
        "",
        "## By Family And Strategy",
        "",
        "| group | cases | avg p-field / Z2 Morse | avg p-field / Z2 standard | "
        "avg Morse speedup vs standard |",
        "|---|---:|---:|---:|---:|",
    ])

    groups: dict[tuple[str, str, int], list[PrimeFieldOverheadRow]] = {}
    for row in rows:
        if row.modulus == 2:
            continue
        groups.setdefault((row.family, row.algorithm, row.modulus), []).append(row)

    for (family, algorithm, modulus), group_rows in sorted(groups.items()):
        avg_morse_over_z2 = statistics.fmean(row.morse_over_z2 for row in group_rows)
        avg_standard_over_z2 = statistics.fmean(row.standard_over_z2 for row in group_rows)
        avg_speedup = statistics.fmean(row.morse_speedup_vs_standard for row in group_rows)
        lines.append(
            f"| {family} / {algorithm} / F_{modulus} | {len(group_rows)} | "
            f"{avg_morse_over_z2:.3f} | {avg_standard_over_z2:.3f} | "
            f"{avg_speedup:.3f}x |"
        )

    mismatches = [
        row
        for row in rows
        if not row.matches_standard or (row.modulus != 2 and not row.matches_z2)
    ]
    lines.extend(["", f"Validation mismatches: {len(mismatches)}"])
    if mismatches:
        for row in mismatches[:10]:
            lines.append(
                f"- {row.name}, {row.algorithm}, F_{row.modulus}: "
                f"matches_standard={row.matches_standard}, matches_z2={row.matches_z2}"
            )
    return "\n".join(lines) + "\n"


def write_csv(rows: list[PrimeFieldOverheadRow], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--families", nargs="+", default=list(DEFAULT_FAMILIES))
    parser.add_argument("--sizes", nargs="+", type=int, default=list(DEFAULT_SIZES))
    parser.add_argument("--seeds", nargs="+", type=int, default=list(DEFAULT_SEEDS))
    parser.add_argument("--algorithms", nargs="+", default=list(DEFAULT_ALGORITHMS))
    parser.add_argument("--primes", nargs="+", type=int, default=list(DEFAULT_PRIMES))
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--plateau-levels", type=int, default=3)
    parser.add_argument("--roadmap-cache", type=Path, default=DEFAULT_ROADMAP_CACHE)
    parser.add_argument("--download-roadmap-data", action="store_true")
    parser.add_argument("--output-csv", type=Path)
    parser.add_argument("--output-md", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows: list[PrimeFieldOverheadRow] = []
    for family in args.families:
        for size in args.sizes:
            for seed in args.seeds:
                for algorithm in args.algorithms:
                    rows.extend(
                        benchmark_case(
                            family=family,
                            seed=seed,
                            size=size,
                            algorithm=algorithm,
                            primes=tuple(args.primes),
                            repeats=args.repeats,
                            plateau_levels=args.plateau_levels,
                            roadmap_cache=args.roadmap_cache,
                            download_roadmap_data=args.download_roadmap_data,
                        )
                    )

    if not rows:
        raise RuntimeError("No benchmark rows were produced.")

    if args.output_csv is not None:
        write_csv(rows, args.output_csv)
    summary = summarize_markdown(rows)
    if args.output_md is not None:
        args.output_md.parent.mkdir(parents=True, exist_ok=True)
        args.output_md.write_text(summary)
    sys.stdout.write(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
