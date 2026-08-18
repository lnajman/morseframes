#!/usr/bin/env python3
"""Render critical-quality and timing comparisons for simplicial strategies."""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT.parent / "simplicial_strategy_benchmark.csv"
DEFAULT_FIGURE = ROOT / "docs" / "simplicial_strategy_comparison.svg"
DEFAULT_TABLE = ROOT / "docs" / "simplicial_strategy_comparison_table.tex"

ALGORITHM_ORDER = (
    "process-lower-stars",
    "process-lower-stars-parallel",
    "f-max",
    "f-min",
    "same-level-reduction",
    "saturated",
)
ALGORITHM_LABELS = {
    "process-lower-stars": "ProcessLowerStars",
    "process-lower-stars-parallel": "ProcessLowerStars (8 workers)",
    "f-max": "F-Max",
    "f-min": "F-Min",
    "same-level-reduction": "Same-level reduction",
    "saturated": "Saturated",
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as source:
        rows = list(csv.DictReader(source))
    if not rows:
        raise ValueError(f"No benchmark rows found in {path}")
    return rows


def _summary(values: list[float]) -> tuple[float, float, float]:
    return min(values), statistics.median(values), max(values)


def render_figure(rows: list[dict[str, str]], output: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is required to render the comparison") from exc

    plt.rcParams["svg.hashsalt"] = "morseframes-simplicial-strategies"
    by_algorithm: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_algorithm[row["algorithm"]].append(row)

    algorithms = [algorithm for algorithm in ALGORITHM_ORDER if algorithm in by_algorithm]
    quality = {
        algorithm: _summary(
            [float(row["critical_count_ratio_vs_f_max"]) for row in by_algorithm[algorithm]]
        )
        for algorithm in algorithms
    }
    time_ratio = {
        algorithm: _summary(
            [1.0 / float(row["total_speedup_vs_f_max"]) for row in by_algorithm[algorithm]]
        )
        for algorithm in algorithms
    }

    figure, axes = plt.subplots(1, 2, figsize=(10.6, 4.8), sharey=True)
    y_positions = list(range(len(algorithms)))
    colors = [
        "#174A7E" if algorithm.startswith("process-lower-stars") else "#60656F"
        for algorithm in algorithms
    ]
    markers = ["o", "s", "D", "^", "v", "P"]
    for axis, summaries, title, xlabel in (
        (
            axes[0],
            quality,
            "Critical simplices relative to F-Max",
            "Critical-count ratio (lower is better)",
        ),
        (
            axes[1],
            time_ratio,
            "End-to-end time relative to F-Max",
            "Time ratio (lower is better)",
        ),
    ):
        panel_maximum = max(summary[2] for summary in summaries.values())
        for index, algorithm in enumerate(algorithms):
            minimum, median, maximum = summaries[algorithm]
            axis.errorbar(
                median,
                index,
                xerr=((median - minimum,), (maximum - median,)),
                fmt=markers[index % len(markers)],
                color=colors[index],
                ecolor=colors[index],
                elinewidth=1.5,
                capsize=3,
                markersize=6,
            )
            near_right_edge = median > 0.78 * panel_maximum
            axis.annotate(
                f"{median:.2f}x",
                (median, index),
                xytext=(-7 if near_right_edge else 7, 0),
                textcoords="offset points",
                ha="right" if near_right_edge else "left",
                va="center",
                fontsize=8,
                color="#30343A",
            )
        axis.axvline(1.0, color="#30343A", linestyle=(0, (2, 3)), linewidth=1.1)
        axis.set_title(title, loc="left")
        axis.set_xlabel(xlabel)
        axis.set_yticks(y_positions)
        axis.set_yticklabels([ALGORITHM_LABELS[algorithm] for algorithm in algorithms])
        axis.set_xlim(0, panel_maximum * 1.12)
        axis.grid(axis="x", color="#D9DDE3", linewidth=0.8)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.spines["left"].set_color("#60656F")
        axis.spines["bottom"].set_color("#60656F")
    axes[0].invert_yaxis()
    axes[1].tick_params(axis="y", labelleft=False)
    figure.suptitle(
        "MorseFrames strategies on injective simplicial terrains",
        x=0.08,
        y=0.97,
        ha="left",
        fontsize=15,
    )
    figure.text(
        0.08,
        0.91,
        "Median and min-max range across 9 cases (3 sizes x 3 seeds); F-Max = 1",
        ha="left",
        fontsize=9,
        color="#50555C",
    )
    figure.subplots_adjust(left=0.23, right=0.98, bottom=0.15, top=0.79, wspace=0.32)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(
        output,
        metadata={"Title": "Simplicial strategy comparison", "Date": None},
    )
    plt.close(figure)
    if output.suffix.lower() == ".svg":
        normalized = "\n".join(line.rstrip() for line in output.read_text().splitlines())
        output.write_text(normalized + "\n")


def render_table(rows: list[dict[str, str]], output: Path) -> None:
    grouped: dict[tuple[int, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(int(row["grid_size"]), row["algorithm"])].append(row)

    lines = [
        r"\begin{tabular}{rlrrrrrr}",
        r"\toprule",
        (
            r"Grid & Strategy & Simplices & Critical & Critical/F-Max & "
            r"Sequence (ms) & Total (ms) & Total speedup \\"
        ),
        r"\midrule",
    ]
    for grid_size in sorted({key[0] for key in grouped}):
        for algorithm in ALGORITHM_ORDER:
            group = grouped.get((grid_size, algorithm))
            if not group:
                continue
            median = statistics.median
            lines.append(
                " & ".join(
                    (
                        str(grid_size),
                        ALGORITHM_LABELS[algorithm],
                        f'{int(group[0]["num_simplices"]):,}',
                        f'{median(int(row["num_critical_simplices"]) for row in group):.0f}',
                        f'{median(float(row["critical_count_ratio_vs_f_max"]) for row in group):.2f}',
                        f'{1000.0 * median(float(row["sequence_seconds"]) for row in group):.2f}',
                        f'{1000.0 * median(float(row["total_seconds"]) for row in group):.2f}',
                        f'{median(float(row["total_speedup_vs_f_max"]) for row in group):.2f}',
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
