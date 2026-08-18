#!/usr/bin/env python3
"""Render reduction-kernel worker-scaling benchmark artifacts."""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT.parent / "reduction_kernel_scaling.csv"
DEFAULT_FIGURE = ROOT / "docs" / "reduction_kernel_scaling.svg"
DEFAULT_TABLE = ROOT / "docs" / "reduction_kernel_scaling_table.tex"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as source:
        rows = list(csv.DictReader(source))
    if not rows:
        raise ValueError(f"No benchmark rows found in {path}")
    return rows


def _summary(values: list[float]) -> tuple[float, float, float]:
    return min(values), statistics.median(values), max(values)


def _median_column(rows: list[dict[str, str]], column: str) -> float:
    return statistics.median(float(row[column]) for row in rows)


def render_figure(rows: list[dict[str, str]], output: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is required to render the scaling plot") from exc

    plt.rcParams["svg.hashsalt"] = "morseframes-reduction-kernel-scaling"
    by_workers: dict[int, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_workers[int(row["max_workers"])].append(row)
    workers = sorted(by_workers)

    sequence = {
        worker: _summary(
            [float(row["sequence_speedup_vs_one_worker"]) for row in group]
        )
        for worker, group in by_workers.items()
    }
    total = {
        worker: _summary(
            [float(row["total_speedup_vs_one_worker"]) for row in group]
        )
        for worker, group in by_workers.items()
    }

    figure, axes = plt.subplots(1, 2, figsize=(10.8, 4.6), sharex=True)
    for axis, summaries, title, ideal in (
        (axes[0], sequence, "Kernel construction speedup", True),
        (axes[1], total, "End-to-end speedup", False),
    ):
        medians = [summaries[worker][1] for worker in workers]
        lower = [
            summaries[worker][1] - summaries[worker][0] for worker in workers
        ]
        upper = [
            summaries[worker][2] - summaries[worker][1] for worker in workers
        ]
        axis.errorbar(
            workers,
            medians,
            yerr=(lower, upper),
            color="#B4473D",
            marker="o",
            linewidth=1.8,
            capsize=3,
        )
        if ideal:
            axis.plot(
                workers,
                workers,
                color="#60656F",
                linestyle=(0, (3, 3)),
                linewidth=1.0,
                label="Ideal construction scaling",
            )
            axis.legend(frameon=False, loc="upper left", fontsize=8)
        else:
            axis.axhline(
                1.0, color="#60656F", linestyle=(0, (3, 3)), linewidth=1.0
            )
        axis.set_title(title, loc="left")
        axis.set_xlabel("Workers")
        axis.set_ylabel("Speedup relative to one worker")
        axis.set_xticks(workers)
        maximum = max(
            max(summary[2] for summary in summaries.values()),
            max(workers) if ideal else 1,
        )
        axis.set_ylim(0, maximum * 1.08)
        axis.grid(color="#D9DDE3", linewidth=0.8)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    case_count = len({row["name"] for row in rows})
    figure.suptitle(
        "Reduction-kernel scaling after coarse level scheduling",
        x=0.08,
        y=0.97,
        ha="left",
        fontsize=15,
    )
    figure.text(
        0.08,
        0.90,
        f"Median and min-max range across {case_count} injective terrains",
        ha="left",
        fontsize=9,
        color="#50555C",
    )
    figure.subplots_adjust(left=0.10, right=0.98, bottom=0.15, top=0.78, wspace=0.30)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(
        output,
        metadata={"Title": "Reduction-kernel worker scaling", "Date": None},
    )
    plt.close(figure)
    if output.suffix.lower() == ".svg":
        normalized = "\n".join(line.rstrip() for line in output.read_text().splitlines())
        output.write_text(normalized + "\n")


def render_table(rows: list[dict[str, str]], output: Path) -> None:
    grouped: dict[tuple[int, int], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(int(row["grid_size"]), int(row["max_workers"]))].append(row)

    lines = [
        r"\begin{tabular}{rrrrrrrr}",
        r"\toprule",
        (
            r"Grid & Workers & Simplices & Sequence (ms) & Sequence speedup & "
            r"Total (ms) & Total speedup & Level batches \\"
        ),
        r"\midrule",
    ]
    median = statistics.median
    for grid_size in sorted({key[0] for key in grouped}):
        for workers in sorted(key[1] for key in grouped if key[0] == grid_size):
            group = grouped[(grid_size, workers)]
            lines.append(
                " & ".join(
                    (
                        str(grid_size),
                        str(workers),
                        f'{int(group[0]["num_simplices"]):,}',
                        f'{1000.0 * median(float(row["sequence_seconds"]) for row in group):.2f}',
                        f'{_median_column(group, "sequence_speedup_vs_one_worker"):.2f}',
                        f'{1000.0 * median(float(row["total_seconds"]) for row in group):.2f}',
                        f'{_median_column(group, "total_speedup_vs_one_worker"):.2f}',
                        f'{median(int(row["parallel_level_batches"]) for row in group):.0f}',
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
