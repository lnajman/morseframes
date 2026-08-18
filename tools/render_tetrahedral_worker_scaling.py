#!/usr/bin/env python3
"""Render tetrahedral worker-scaling benchmark artifacts."""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT.parent / "tetrahedral_worker_scaling.csv"
DEFAULT_FIGURE = ROOT / "docs" / "tetrahedral_worker_scaling.svg"
DEFAULT_TABLE = ROOT / "docs" / "tetrahedral_worker_scaling_table.tex"

STRATEGY_ORDER = ("process-lower-stars", "reduction-kernel")
STRATEGY_LABELS = {
    "process-lower-stars": "ProcessLowerStars",
    "reduction-kernel": "Reduction kernel",
}
STRATEGY_STYLES = {
    "process-lower-stars": ("#174A7E", "o", "-"),
    "reduction-kernel": ("#C76822", "s", "--"),
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
        raise RuntimeError("matplotlib is required to render the scaling plot") from exc

    plt.rcParams["svg.hashsalt"] = "morseframes-tetrahedral-worker-scaling"
    grouped: dict[tuple[str, int], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["strategy"], int(row["max_workers"]))].append(row)
    strategies = [
        strategy
        for strategy in STRATEGY_ORDER
        if any(key[0] == strategy for key in grouped)
    ]
    workers = sorted({key[1] for key in grouped})

    figure, axes = plt.subplots(1, 2, figsize=(10.8, 4.8), sharex=True)
    for axis, column, title in (
        (axes[0], "gradient_speedup_vs_one_worker", "Gradient speedup"),
        (axes[1], "simplices_per_second", "Gradient throughput"),
    ):
        for strategy in strategies:
            summaries = {
                worker: _summary(
                    [float(row[column]) for row in grouped[(strategy, worker)]]
                )
                for worker in workers
            }
            medians = [summaries[worker][1] for worker in workers]
            lower = [summaries[worker][1] - summaries[worker][0] for worker in workers]
            upper = [summaries[worker][2] - summaries[worker][1] for worker in workers]
            color, marker, linestyle = STRATEGY_STYLES[strategy]
            axis.errorbar(
                workers,
                medians,
                yerr=(lower, upper),
                color=color,
                marker=marker,
                linestyle=linestyle,
                linewidth=1.8,
                capsize=3,
                label=STRATEGY_LABELS[strategy],
            )
        if column == "gradient_speedup_vs_one_worker":
            axis.plot(
                workers,
                workers,
                color="#60656F",
                linestyle=(0, (3, 3)),
                linewidth=1.0,
                label="Ideal scaling",
            )
        axis.set_title(title, loc="left")
        axis.set_xlabel("Workers")
        axis.set_ylabel(
            "Speedup relative to one worker"
            if column == "gradient_speedup_vs_one_worker"
            else "Simplices per second"
        )
        axis.set_xticks(workers)
        if column == "gradient_speedup_vs_one_worker":
            axis.set_ylim(0, max(workers) * 1.08)
        axis.grid(color="#D9DDE3", linewidth=0.8)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
    axes[0].legend(frameon=False, loc="upper left", fontsize=8)
    axes[1].legend(frameon=False, loc="upper left", fontsize=8)

    case_count = len({row["name"] for row in rows})
    figure.suptitle(
        "Worker scaling on injective tetrahedral volumes",
        x=0.08,
        y=0.97,
        ha="left",
        fontsize=15,
    )
    figure.text(
        0.08,
        0.90,
        f"Median and min-max range across {case_count} cases; one worker = 1",
        ha="left",
        fontsize=9,
        color="#50555C",
    )
    figure.subplots_adjust(left=0.10, right=0.98, bottom=0.15, top=0.78, wspace=0.30)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(
        output,
        metadata={"Title": "Tetrahedral worker scaling", "Date": None},
    )
    plt.close(figure)
    if output.suffix.lower() == ".svg":
        normalized = "\n".join(
            line.rstrip() for line in output.read_text().splitlines()
        )
        output.write_text(normalized + "\n")


def render_table(rows: list[dict[str, str]], output: Path) -> None:
    grouped: dict[tuple[int, str, int], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[
            (int(row["grid_size"]), row["strategy"], int(row["max_workers"]))
        ].append(row)

    lines = [
        r"\begin{tabular}{rlrrrrrr}",
        r"\toprule",
        (
            r"Grid & Strategy & Workers & Simplices & Gradient (ms) & "
            r"Speedup & Efficiency & M simplex/s \\"
        ),
        r"\midrule",
    ]
    median = statistics.median
    for grid_size in sorted({key[0] for key in grouped}):
        for strategy in STRATEGY_ORDER:
            for workers in sorted(
                key[2] for key in grouped if key[:2] == (grid_size, strategy)
            ):
                group = grouped[(grid_size, strategy, workers)]
                lines.append(
                    " & ".join(
                        (
                            str(grid_size),
                            STRATEGY_LABELS[strategy],
                            str(workers),
                            f'{int(group[0]["num_simplices"]):,}',
                            f'{1000.0 * median(float(row["gradient_seconds"]) for row in group):.2f}',
                            f'{median(float(row["gradient_speedup_vs_one_worker"]) for row in group):.2f}',
                            f'{median(float(row["parallel_efficiency"]) for row in group):.2f}',
                            f'{1.0e-6 * median(float(row["simplices_per_second"]) for row in group):.2f}',
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
