#!/usr/bin/env python3
"""Render the unified gradient-only strategy comparison."""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT.parent / "gradient_strategy_benchmark.csv"
DEFAULT_FIGURE = ROOT / "docs" / "gradient_strategy_comparison.svg"
DEFAULT_TABLE = ROOT / "docs" / "gradient_strategy_comparison_table.tex"

FAMILY_ORDER = ("injective-terrain", "injective-volume")
FAMILY_LABELS = {
    "injective-terrain": "2D triangulated terrain",
    "injective-volume": "3D tetrahedral volume",
}
PARALLEL_ALGORITHMS = (
    "process-lower-stars-parallel",
    "flooding-reduction-kernel-parallel",
)
ALGORITHM_LABELS = {
    "f-max": "F-Max",
    "process-lower-stars": "ProcessLowerStars (sequential)",
    "process-lower-stars-parallel": "ProcessLowerStars",
    "flooding-reduction-kernel": "Reduction kernel (sequential)",
    "flooding-reduction-kernel-parallel": "Reduction kernel",
}
STYLES = {
    "process-lower-stars-parallel": ("#174A7E", "o", "-"),
    "flooding-reduction-kernel-parallel": ("#C76822", "s", "--"),
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as source:
        rows = list(csv.DictReader(source))
    if not rows:
        raise ValueError(f"No benchmark rows found in {path}")
    return rows


def _median(rows: list[dict[str, str]], column: str) -> float:
    return statistics.median(float(row[column]) for row in rows)


def render_figure(rows: list[dict[str, str]], output: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is required to render the comparison") from exc

    plt.rcParams["svg.hashsalt"] = "morseframes-gradient-strategies"
    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.9), sharex=True)
    for axis, family in zip(axes, FAMILY_ORDER, strict=True):
        family_rows = [row for row in rows if row["family"] == family]
        workers = sorted(
            {
                int(row["max_workers"])
                for row in family_rows
                if row["algorithm"] in PARALLEL_ALGORITHMS
            }
        )
        for algorithm in PARALLEL_ALGORITHMS:
            medians = []
            minima = []
            maxima = []
            for worker in workers:
                group = [
                    row
                    for row in family_rows
                    if row["algorithm"] == algorithm
                    and int(row["max_workers"]) == worker
                ]
                values = [
                    float(row["construction_time_ratio_vs_f_max"])
                    for row in group
                ]
                medians.append(statistics.median(values))
                minima.append(min(values))
                maxima.append(max(values))
            color, marker, linestyle = STYLES[algorithm]
            axis.errorbar(
                workers,
                medians,
                yerr=(
                    [median - minimum for median, minimum in zip(medians, minima)],
                    [maximum - median for median, maximum in zip(medians, maxima)],
                ),
                color=color,
                marker=marker,
                linestyle=linestyle,
                linewidth=1.8,
                capsize=3,
                label=ALGORITHM_LABELS[algorithm],
            )
        axis.axhline(
            1.0,
            color="#30343A",
            linestyle=(0, (3, 3)),
            linewidth=1.1,
            label="F-Max",
        )
        axis.set_title(FAMILY_LABELS[family], loc="left")
        axis.set_xlabel("Workers")
        axis.set_ylabel("Gradient time relative to F-Max")
        axis.set_xticks(workers)
        axis.grid(color="#D9DDE3", linewidth=0.8)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.legend(frameon=False, fontsize=8)
    case_count = len({row["name"] for row in rows})
    figure.suptitle(
        "Gradient construction: F-Max, ProcessLowerStars, and ReductionKernel",
        x=0.08,
        y=0.97,
        ha="left",
        fontsize=14,
    )
    figure.text(
        0.08,
        0.90,
        f"Median and min-max range across {case_count} cases; F-Max = 1",
        fontsize=9,
        color="#50555C",
    )
    figure.subplots_adjust(left=0.09, right=0.98, bottom=0.15, top=0.78, wspace=0.28)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(
        output,
        metadata={"Title": "Gradient strategy comparison", "Date": None},
    )
    plt.close(figure)
    if output.suffix.lower() == ".svg":
        output.write_text(
            "\n".join(line.rstrip() for line in output.read_text().splitlines())
            + "\n"
        )


def render_table(rows: list[dict[str, str]], output: Path) -> None:
    grouped: dict[tuple[str, str, int], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[
            (row["family"], row["algorithm"], int(row["max_workers"]))
        ].append(row)

    parallel_workers = sorted(
        {
            int(row["max_workers"])
            for row in rows
            if row["algorithm"] in PARALLEL_ALGORITHMS
        }
    )
    algorithm_order = (
        ("f-max", 1),
        ("process-lower-stars", 1),
        *(("process-lower-stars-parallel", worker) for worker in parallel_workers),
        ("flooding-reduction-kernel", 1),
        *(
            ("flooding-reduction-kernel-parallel", worker)
            for worker in parallel_workers
        ),
    )
    lines = [
        r"\begin{tabular}{llrrrrrr}",
        r"\toprule",
        (
            r"Data & Strategy & Workers & Critical/F-Max & Time/F-Max & "
            r"Speedup$_1$ & Efficiency & M simplex/s \\"
        ),
        r"\midrule",
    ]
    for family in FAMILY_ORDER:
        for algorithm, workers in algorithm_order:
            group = grouped.get((family, algorithm, workers))
            if not group:
                continue
            lines.append(
                " & ".join(
                    (
                        FAMILY_LABELS[family],
                        ALGORITHM_LABELS[algorithm],
                        str(workers),
                        f'{_median(group, "critical_count_ratio_vs_f_max"):.2f}',
                        f'{_median(group, "construction_time_ratio_vs_f_max"):.2f}',
                        f'{_median(group, "worker_speedup_vs_one"):.2f}',
                        f'{_median(group, "parallel_efficiency"):.2f}',
                        f'{1.0e-6 * _median(group, "simplices_per_second"):.2f}',
                    )
                )
                + r" \\"
            )
        lines.append(r"\addlinespace")
    lines.extend((r"\bottomrule", r"\end{tabular}"))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--figure-output", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument("--table-output", type=Path, default=DEFAULT_TABLE)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows = read_rows(args.input)
    render_figure(rows, args.figure_output)
    render_table(rows, args.table_output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
