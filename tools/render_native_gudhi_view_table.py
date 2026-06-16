#!/usr/bin/env python3
"""Render the native GUDHI Simplex_tree benchmark table from CSV rows."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "docs" / "native_gudhi_view_quick.csv"
DEFAULT_OUTPUT = ROOT / "docs" / "native_gudhi_view_quick_table.tex"
DEFAULT_CAPTION_TITLE = r"Native \texttt{Gudhi::Simplex\_tree} quick benchmark."
DEFAULT_LABEL = "tab:native-gudhi-view-quick"
DEFAULT_STRATEGIES = (
    "f-max",
    "f-min",
    "plateau-greedy",
    "same-level-reduction",
)
STRATEGY_LABELS = {
    "f-max": r"\texttt{f-max}",
    "f-min": r"\texttt{f-min}",
    "plateau-greedy": r"\texttt{plateau-greedy}",
    "same-level-reduction": r"\slred",
}


@dataclass(frozen=True)
class RatioSummary:
    median: float
    q1: float
    q3: float
    samples: int


@dataclass(frozen=True)
class NativeGudhiSummary:
    name: str
    strategy: str
    simplices: int
    direct_full_ms: float
    compact_full_ms: float
    gudhi_persistence_ms: float
    direct_reducer_ms: float
    gudhi_over_direct: RatioSummary
    gudhi_over_reducer: RatioSummary


def _float(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    if value == "":
        raise ValueError(f"Missing required numeric column {key!r}.")
    return float(value)


def _int(row: dict[str, str], key: str) -> int:
    value = row.get(key, "")
    if value == "":
        raise ValueError(f"Missing required integer column {key!r}.")
    return int(value)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _quantile(values: list[float], probability: float) -> float:
    if not values:
        raise ValueError("Cannot compute a quantile of an empty list.")
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def _ratio_summary(values: Iterable[float]) -> RatioSummary:
    materialized = list(values)
    return RatioSummary(
        median=statistics.median(materialized),
        q1=_quantile(materialized, 0.25),
        q3=_quantile(materialized, 0.75),
        samples=len(materialized),
    )


def summarize_rows(
    rows: Iterable[dict[str, str]],
    *,
    strategies: Iterable[str] = DEFAULT_STRATEGIES,
    statistic: str = "mean",
) -> list[NativeGudhiSummary]:
    strategy_order = tuple(strategies)
    strategy_rank = {strategy: index for index, strategy in enumerate(strategy_order)}
    if statistic == "mean":
        summarize = statistics.fmean
    elif statistic == "median":
        summarize = statistics.median
    else:
        raise ValueError(f"Unknown statistic {statistic!r}.")
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        strategy = row.get("strategy", "")
        if strategy in strategy_rank:
            grouped[(row["name"], strategy)].append(row)

    summaries: list[NativeGudhiSummary] = []
    for key in sorted(grouped, key=lambda item: (item[0], strategy_rank[item[1]])):
        group = grouped[key]
        summaries.append(
            NativeGudhiSummary(
                name=key[0],
                strategy=key[1],
                simplices=_int(group[0], "simplices"),
                direct_full_ms=summarize(
                    _float(row, "direct_full_ms") for row in group
                ),
                compact_full_ms=summarize(
                    _float(row, "compact_full_ms") for row in group
                ),
                gudhi_persistence_ms=summarize(
                    _float(row, "gudhi_persistence_ms") for row in group
                ),
                direct_reducer_ms=summarize(
                    _float(row, "direct_reducer_ms") for row in group
                ),
                gudhi_over_direct=_ratio_summary(
                    _float(row, "gudhi_persistence_ms")
                    / _float(row, "direct_full_ms")
                    for row in group
                ),
                gudhi_over_reducer=_ratio_summary(
                    _float(row, "gudhi_persistence_ms")
                    / _float(row, "direct_reducer_ms")
                    for row in group
                ),
            )
        )
    return summaries


def _texttt(value: str) -> str:
    return r"\texttt{" + value.replace("_", r"\_") + "}"


def _format_ratio(summary: RatioSummary) -> str:
    if summary.samples <= 1:
        return f"{summary.median:.2f}"
    return (
        r"\ratioiqr{"
        f"{summary.median:.2f}"
        r"}{"
        f"{summary.q1:.2f}"
        r"}{"
        f"{summary.q3:.2f}"
        r"}"
    )


def render_latex_table(
    rows: list[NativeGudhiSummary],
    *,
    caption_title: str = DEFAULT_CAPTION_TITLE,
    label: str = DEFAULT_LABEL,
) -> str:
    lines = [
        "% Generated by morseframes/tools/render_native_gudhi_view_table.py.",
        "% Re-run that script instead of editing this table by hand.",
        r"\providecommand{\ratioiqr}[3]{\shortstack[c]{#1\\[-0.25ex]{\tiny #2--#3}}}",
        r"\begin{table}[ht]",
        r"\centering",
        r"\scriptsize",
        r"\setlength{\tabcolsep}{2.4pt}",
        r"\begin{tabular}{llrrrrcc}",
        r"\toprule",
        (
            r"Case & Strategy & Simp. & Direct (ms) & Import (ms) & GUDHI (ms) "
            r"& GUDHI/Direct & GUDHI/Reducer \\"
        ),
        r"\midrule",
    ]
    previous_name: str | None = None
    for row in rows:
        if previous_name is not None and row.name != previous_name:
            lines.append(r"\midrule")
        previous_name = row.name
        lines.append(
            f"{_texttt(row.name)} & "
            f"{STRATEGY_LABELS.get(row.strategy, _texttt(row.strategy))} & "
            f"{row.simplices:,} & "
            f"{row.direct_full_ms:.2f} & "
            f"{row.compact_full_ms:.2f} & "
            f"{row.gudhi_persistence_ms:.2f} & "
            f"{_format_ratio(row.gudhi_over_direct)} & "
            f"{_format_ratio(row.gudhi_over_reducer)} \\\\"
        )
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            (
                rf"\caption{{{caption_title}  "
                r"\texttt{Direct} uses the read-only \texttt{Simplex\_tree} view "
                r"and our Morse reducer; \texttt{Import} first copies the complex "
                r"into the compact owning structure.  \texttt{GUDHI} is GUDHI's "
                r"in-process persistent cohomology on the same \texttt{Simplex\_tree}.  "
                r"\texttt{GUDHI/Direct}$<1$ means GUDHI is faster end-to-end; "
                r"\texttt{GUDHI/Reducer}$>1$ means our reducer alone is faster than "
                r"GUDHI persistence, excluding view construction and Morse input "
                r"construction.  Ratio columns show the median with $Q_1$--$Q_3$ "
                r"underneath when repeat-level rows are available.}"
            ),
            rf"\label{{{label}}}",
            r"\end{table}",
            "",
        ]
    )
    return "\n".join(lines)


def render_summary(rows: list[NativeGudhiSummary]) -> str:
    direct_ratios = [row.gudhi_over_direct.median for row in rows]
    reducer_ratios = [
        row.gudhi_over_reducer.median
        for row in rows
        if row.strategy in {"f-max", "f-min", "same-level-reduction"}
    ]
    return (
        f"Median GUDHI/Direct range: {min(direct_ratios):.2f}--{max(direct_ratios):.2f}\n"
        f"Main-strategy median GUDHI/Reducer range: "
        f"{min(reducer_ratios):.2f}--{max(reducer_ratios):.2f}\n"
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render the native GUDHI benchmark LaTeX table."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--strategies", nargs="+", default=list(DEFAULT_STRATEGIES))
    parser.add_argument("--caption-title", default=DEFAULT_CAPTION_TITLE)
    parser.add_argument("--label", default=DEFAULT_LABEL)
    parser.add_argument("--statistic", choices=("mean", "median"), default="mean")
    parser.add_argument("--stdout", action="store_true")
    parser.add_argument("--summary", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    rows = summarize_rows(
        read_rows(args.input), strategies=args.strategies, statistic=args.statistic
    )
    if not rows:
        raise ValueError("No native GUDHI benchmark rows matched the requested filters.")
    table = render_latex_table(
        rows, caption_title=args.caption_title, label=args.label
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(table)
    if args.stdout:
        sys.stdout.write(table)
    if args.summary:
        sys.stdout.write(render_summary(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
