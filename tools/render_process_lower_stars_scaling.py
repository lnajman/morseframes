#!/usr/bin/env python3
"""Render the ProcessLowerStars worker-scaling figure and LaTeX table."""

from __future__ import annotations

import argparse
import ast
import csv
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT.parent / "process_lower_stars_scaling.csv"
DEFAULT_FIGURE = ROOT / "docs" / "process_lower_stars_scaling.svg"
DEFAULT_TABLE = ROOT / "docs" / "process_lower_stars_scaling_table.tex"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as source:
        rows = list(csv.DictReader(source))
    if not rows:
        raise ValueError(f"No benchmark rows found in {path}")
    return rows


def _parallel_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [
        row
        for row in rows
        if row["algorithm"] == "process-lower-stars-parallel"
    ]


def render_figure(rows: list[dict[str, str]], output: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is required to render the scaling figure") from exc

    plt.rcParams["svg.hashsalt"] = "morseframes-process-lower-stars"

    parallel = _parallel_rows(rows)
    distributions = ("balanced", "skewed")
    scales = sorted({int(row["num_simplices"]) for row in parallel})
    colors = ("#9CC4E4", "#4F81BD", "#174A7E", "#102F4C")
    markers = ("o", "s", "^", "D")
    line_styles = ("-", "--", "-.", ":")
    workers = sorted({int(row["max_workers"]) for row in parallel})
    max_observed_speedup = max(
        float(row["sequence_speedup_vs_sequential"]) for row in parallel
    )

    figure, axes = plt.subplots(1, 2, figsize=(10.2, 4.4), sharex=True, sharey=True)
    for axis, distribution in zip(axes, distributions):
        for index, scale in enumerate(scales):
            selected = sorted(
                (
                    row
                    for row in parallel
                    if row["distribution"] == distribution
                    and int(row["num_simplices"]) == scale
                ),
                key=lambda row: int(row["max_workers"]),
            )
            axis.plot(
                [int(row["max_workers"]) for row in selected],
                [float(row["sequence_speedup_vs_sequential"]) for row in selected],
                color=colors[index % len(colors)],
                marker=markers[index % len(markers)],
                linestyle=line_styles[index % len(line_styles)],
                linewidth=2.0,
                markersize=5.5,
                label=f"{scale:,} simplices",
            )
        axis.plot(
            workers,
            workers,
            color="#60656F",
            linestyle=(0, (2, 3)),
            linewidth=1.2,
            label="ideal scaling",
        )
        axis.set_title(f"{distribution.capitalize()} anchor workloads", loc="left")
        axis.set_xlabel("Worker budget")
        axis.set_xticks(workers)
        axis.set_xlim(min(workers), max(workers))
        axis.set_ylim(0, max(2.0, max_observed_speedup * 1.15))
        axis.grid(axis="y", color="#D9DDE3", linewidth=0.8)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.spines["left"].set_color("#60656F")
        axis.spines["bottom"].set_color("#60656F")
    axes[0].set_ylabel("Sequence speedup vs. sequential")
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.63, 0.84),
        ncol=4,
        fontsize=8,
    )
    figure.suptitle(
        "ProcessLowerStars sequence scaling",
        x=0.08,
        y=0.97,
        ha="left",
        fontsize=15,
    )
    figure.text(
        0.08,
        0.91,
        (
            "Matched simplex counts; best of repeated runs; exact sequences validated; "
            "ideal reference exceeds display range"
        ),
        ha="left",
        fontsize=9,
        color="#50555C",
    )
    figure.subplots_adjust(left=0.08, right=0.98, bottom=0.14, top=0.72, wspace=0.08)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(
        output,
        bbox_inches="tight",
        metadata={"Title": "ProcessLowerStars scaling", "Date": None},
    )
    plt.close(figure)


def _latex_escape(value: str) -> str:
    return value.replace("_", r"\_").replace("%", r"\%")


def render_table(rows: list[dict[str, str]], output: Path) -> None:
    parallel = _parallel_rows(rows)
    by_case: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in parallel:
        by_case[row["case"]].append(row)

    lines = [
        r"\begin{tabular}{lrrrrrrr}",
        r"\toprule",
        (
            r"Workload & Simplices & Max star & Load ratio & Critical (by dim.) "
            r"& Critical \% & Best seq. speedup & Best total speedup \\"
        ),
        r"\midrule",
    ]
    for case_rows in sorted(
        by_case.values(),
        key=lambda group: (int(group[0]["num_simplices"]), group[0]["distribution"]),
    ):
        max_worker_row = max(case_rows, key=lambda row: int(row["max_workers"]))
        best_sequence = max(float(row["sequence_speedup_vs_sequential"]) for row in case_rows)
        best_total = max(float(row["total_speedup_vs_sequential"]) for row in case_rows)
        critical = tuple(ast.literal_eval(max_worker_row["critical_simplices_by_dimension"]))
        critical_text = "(" + ",".join(str(value) for value in critical) + ")"
        lines.append(
            " & ".join(
                (
                    _latex_escape(max_worker_row["distribution"].capitalize()),
                    f'{int(max_worker_row["num_simplices"]):,}',
                    str(int(max_worker_row["max_lower_star_size"])),
                    f'{float(max_worker_row["estimated_task_load_ratio"]):.2f}',
                    critical_text,
                    f'{100.0 * float(max_worker_row["critical_ratio"]):.1f}',
                    f"{best_sequence:.2f}",
                    f"{best_total:.2f}",
                )
            )
            + r" \\"
        )
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
