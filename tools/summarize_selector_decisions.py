#!/usr/bin/env python3
"""Summarize fair selector-validation CSVs into a default-policy decision table."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, TextIO


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = ROOT.parent / "work"
DEFAULT_INPUTS = (
    ("report", DEFAULT_OUTPUT_DIR / "profile_metric_fair_validation_comparison.csv"),
    ("extended-holdout", DEFAULT_OUTPUT_DIR / "profile_metric_extended_holdout_comparison.csv"),
    ("plateau-holdout", DEFAULT_OUTPUT_DIR / "profile_metric_plateau_holdout_comparison.csv"),
    ("terrain-holdout", DEFAULT_OUTPUT_DIR / "profile_metric_terrain_holdout_comparison.csv"),
    ("image-holdout", DEFAULT_OUTPUT_DIR / "profile_metric_image_holdout_comparison.csv"),
    ("image-stress", DEFAULT_OUTPUT_DIR / "profile_metric_image_stress_comparison.csv"),
)
DEFAULT_TABLE_OUTPUT = DEFAULT_OUTPUT_DIR / "profile_selector_decision_summary.txt"
DEFAULT_CSV_OUTPUT = DEFAULT_OUTPUT_DIR / "profile_selector_decision_summary.csv"
DEFAULT_JSON_OUTPUT = DEFAULT_OUTPUT_DIR / "profile_selector_decision_summary.json"
DEFAULT_LATEX_OUTPUT = ROOT / "docs" / "profile_selector_decision_summary_table.tex"
DEFAULT_PROSE_OUTPUT = ROOT / "docs" / "profile_selector_decision_summary_prose.tex"

METRIC_DISPLAY_NAMES = {
    "profile_total_work": "PTW",
    "critical_ratio": "CR",
    "adaptive": "Adapt",
    "adaptive_structured": "Adapt-S",
    "estimated_reducer_work": "ERW",
    "boundary_annotation_work": "BAW",
    "initial_annotation_size": "IAS",
    "profile_seconds": "Psec",
    "working_set_size": "Wset",
}

SPLIT_DISPLAY_NAMES = {
    "base": "Base",
    "base-holdout": "Base hold-out",
    "cam": "CAM",
    "cam-holdout": "CAM hold-out",
    "plateau-holdout": "Plateau hold-out",
    "terrain-holdout": "Terrain hold-out",
    "image-holdout": "Image hold-out",
    "image-stress": "Image stress",
    "roadmap": "Roadmap",
    "roadmap-holdout": "Roadmap hold-out",
}

SOURCE_DISPLAY_NAMES = {
    "report": "Report",
    "extended-holdout": "Extended",
    "plateau-holdout": "Plateau",
    "terrain-holdout": "Terrain",
    "image-holdout": "Image",
    "image-stress": "Image+",
}


@dataclass(frozen=True)
class MetricAggregate:
    source: str
    split: str
    metric: str
    cases: int
    match_percent: float
    unavailable_count: int
    avg_penalty_percent: float
    max_penalty_percent: float
    avg_selected_rank: float


@dataclass(frozen=True)
class SelectorDecisionRow:
    source: str
    split: str
    cases: int
    current_metric: str
    best_metrics: str
    current_avg_penalty_percent: float
    best_avg_penalty_percent: float
    gap_percent_points: float
    current_match_percent: float
    best_match_percent: float
    current_avg_selected_rank: float
    best_avg_selected_rank: float
    current_max_penalty_percent: float
    best_max_penalty_percent: float
    decision: str


def _optional_float(row: dict[str, str], key: str, default: float = 0.0) -> float:
    value = row.get(key, "")
    if value == "":
        return default
    return float(value)


def _optional_int(row: dict[str, str], key: str, default: int = 0) -> int:
    value = row.get(key, "")
    if value == "":
        return default
    return int(float(value))


def parse_input_spec(spec: str) -> tuple[str, Path]:
    if "=" in spec:
        label, path = spec.split("=", 1)
        label = label.strip()
        if not label:
            raise ValueError(f"Invalid input spec {spec!r}: empty label")
        return label, Path(path)
    path = Path(spec)
    return path.stem, path


def parse_metric_label(label: str, default_split: str) -> tuple[str, str]:
    if ":" in label:
        split, metric = label.split(":", 1)
        return split, metric
    return default_split, label


def read_metric_aggregates(
    inputs: Iterable[tuple[str, Path]],
    *,
    family: str = "all",
) -> list[MetricAggregate]:
    aggregates: list[MetricAggregate] = []
    for source, path in inputs:
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                if row.get("family") != family:
                    continue
                split, metric = parse_metric_label(row["label"], source)
                aggregates.append(
                    MetricAggregate(
                        source=source,
                        split=split,
                        metric=metric,
                        cases=_optional_int(row, "cases"),
                        match_percent=_optional_float(row, "match_percent"),
                        unavailable_count=_optional_int(row, "unavailable_count"),
                        avg_penalty_percent=_optional_float(row, "avg_penalty_percent"),
                        max_penalty_percent=_optional_float(row, "max_penalty_percent"),
                        avg_selected_rank=_optional_float(row, "avg_selected_rank"),
                    )
                )
    return aggregates


def _decision(
    *,
    current: MetricAggregate | None,
    best: MetricAggregate,
    best_metrics: Iterable[str],
    tie_tolerance: float,
    promotion_margin: float,
) -> str:
    if current is None:
        return "missing-current"
    if current.metric in set(best_metrics):
        return "keep-current"
    gap = current.avg_penalty_percent - best.avg_penalty_percent
    if gap <= tie_tolerance:
        return "keep-current"
    if gap < promotion_margin:
        return "watch-best"
    return "switch-candidate"


def summarize_decisions(
    inputs: Iterable[tuple[str, Path]],
    *,
    current_metric: str = "adaptive",
    family: str = "all",
    tie_tolerance: float = 0.005,
    promotion_margin: float = 0.5,
) -> list[SelectorDecisionRow]:
    aggregates = read_metric_aggregates(inputs, family=family)
    grouped: dict[tuple[str, str], list[MetricAggregate]] = defaultdict(list)
    for aggregate in aggregates:
        grouped[(aggregate.source, aggregate.split)].append(aggregate)

    decisions: list[SelectorDecisionRow] = []
    for (source, split), group in grouped.items():
        best_penalty = min(row.avg_penalty_percent for row in group)
        best_group = [
            row for row in group
            if row.avg_penalty_percent <= best_penalty + tie_tolerance
        ]
        best_group.sort(key=lambda row: (row.avg_penalty_percent, row.metric))
        best = best_group[0]
        current = next((row for row in group if row.metric == current_metric), None)
        if current is None:
            current_penalty = 0.0
            current_match = 0.0
            current_rank = 0.0
            current_max = 0.0
            cases = best.cases
        else:
            current_penalty = current.avg_penalty_percent
            current_match = current.match_percent
            current_rank = current.avg_selected_rank
            current_max = current.max_penalty_percent
            cases = current.cases

        best_metrics = "/".join(row.metric for row in best_group)
        decisions.append(
            SelectorDecisionRow(
                source=source,
                split=split,
                cases=cases,
                current_metric=current_metric,
                best_metrics=best_metrics,
                current_avg_penalty_percent=current_penalty,
                best_avg_penalty_percent=best.avg_penalty_percent,
                gap_percent_points=max(0.0, current_penalty - best.avg_penalty_percent),
                current_match_percent=current_match,
                best_match_percent=max(row.match_percent for row in best_group),
                current_avg_selected_rank=current_rank,
                best_avg_selected_rank=min(row.avg_selected_rank for row in best_group),
                current_max_penalty_percent=current_max,
                best_max_penalty_percent=min(row.max_penalty_percent for row in best_group),
                decision=_decision(
                    current=current,
                    best=best,
                    best_metrics=[row.metric for row in best_group],
                    tie_tolerance=tie_tolerance,
                    promotion_margin=promotion_margin,
                ),
            )
        )
    return decisions


def _display_metric(metric: str) -> str:
    return METRIC_DISPLAY_NAMES.get(metric, metric.replace("_", "-"))


def _display_metric_list(metrics: str) -> str:
    return "/".join(_display_metric(metric) for metric in metrics.split("/"))


def _latex_metric(metric: str) -> str:
    return r"\texttt{" + metric.replace("_", "-") + "}"


def write_table(rows: list[SelectorDecisionRow], output: TextIO) -> None:
    headers = [
        "source",
        "split",
        "cases",
        "current",
        "best",
        "cur_pen%",
        "best_pen%",
        "gap_pp",
        "decision",
    ]
    output.write(
        f"{headers[0]:>17} {headers[1]:>17} {headers[2]:>7} "
        f"{headers[3]:>9} {headers[4]:>16} "
        + " ".join(f"{header:>10}" for header in headers[5:8])
        + f" {headers[8]:>17}\n"
    )
    for row in rows:
        output.write(
            f"{SOURCE_DISPLAY_NAMES.get(row.source, row.source):>17} "
            f"{row.split:>17} "
            f"{row.cases:7d} "
            f"{_display_metric(row.current_metric):>9} "
            f"{_display_metric_list(row.best_metrics):>16} "
            f"{row.current_avg_penalty_percent:10.2f} "
            f"{row.best_avg_penalty_percent:10.2f} "
            f"{row.gap_percent_points:10.2f} "
            f"{row.decision:>17}\n"
        )


def write_csv(rows: list[SelectorDecisionRow], output: TextIO) -> None:
    writer = csv.DictWriter(output, fieldnames=list(asdict(rows[0]).keys()))
    writer.writeheader()
    for row in rows:
        writer.writerow(asdict(row))


def write_json(rows: list[SelectorDecisionRow], output: TextIO) -> None:
    json.dump([asdict(row) for row in rows], output, indent=2)
    output.write("\n")


def _latex_escape(value: str) -> str:
    return (
        value
        .replace("\\", r"\textbackslash{}")
        .replace("_", r"\_")
        .replace("%", r"\%")
        .replace("&", r"\&")
        .replace("#", r"\#")
    )


def render_latex_table(rows: Iterable[SelectorDecisionRow]) -> str:
    row_list = list(rows)
    current_metric = row_list[0].current_metric if row_list else "adaptive"
    lines = [
        "% Generated by morseframes/tools/summarize_selector_decisions.py.",
        "% Re-run that script instead of editing this table by hand.",
        r"\begin{table}[ht]",
        r"\centering",
        r"\scriptsize",
        r"\setlength{\tabcolsep}{3pt}",
        r"\begin{tabular}{@{}llrlrrl@{}}",
        r"\toprule",
        r"Source & Split & Cases & Best metric(s) & Cur./best pen. \% & Gap pp & Decision \\",
        r"\midrule",
    ]
    for row in row_list:
        source = _latex_escape(SOURCE_DISPLAY_NAMES.get(row.source, row.source))
        split = _latex_escape(SPLIT_DISPLAY_NAMES.get(row.split, row.split))
        best_metrics = _latex_escape(_display_metric_list(row.best_metrics))
        decision = _latex_escape(row.decision.replace("keep-current", "keep").replace("-", " "))
        lines.append(
            f"{source} & "
            f"{split} & "
            f"{row.cases:d} & "
            f"{best_metrics} & "
            f"{row.current_avg_penalty_percent:.2f}/{row.best_avg_penalty_percent:.2f} & "
            f"{row.gap_percent_points:.2f} & "
            f"{decision} " + r"\\"
        )
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            r"\caption{Selector-decision audit over the aggregate rows of the fair",
            r"profile-selection validation summaries.  The current metric is",
            f"{_latex_metric(current_metric)}.  A zero gap means that the current selector is tied with",
            r"the best observed average penalty for that validation scope.}",
            r"\label{tab:selector-decision-summary}",
            r"\end{table}",
            "",
        ]
    )
    return "\n".join(lines)


def render_prose(rows: Iterable[SelectorDecisionRow]) -> str:
    row_list = list(rows)
    current_metric = row_list[0].current_metric if row_list else "adaptive"
    keep_count = sum(1 for row in row_list if row.decision == "keep-current")
    watch_rows = [row for row in row_list if row.decision == "watch-best"]
    switch_rows = [row for row in row_list if row.decision == "switch-candidate"]
    lines = [
        "% Generated by morseframes/tools/summarize_selector_decisions.py.",
        "% Re-run that script instead of editing this paragraph by hand.",
    ]
    if keep_count == len(row_list):
        if current_metric == "adaptive_structured":
            tail = (
                "The structured-grid branch removes the previous image-grid "
                "switch-candidate gap while preserving the previous tied scopes."
            )
        else:
            tail = (
                "In each scope, the current selector is tied with the best observed "
                "average penalty."
            )
        lines.append(
            f"The selector-decision audit keeps the current {_latex_metric(current_metric)} metric in "
            f"all {len(row_list)} aggregate validation scopes.  {tail}"
        )
    else:
        lines.append(
            f"The selector-decision audit keeps the current {_latex_metric(current_metric)} metric in "
            f"{keep_count} of {len(row_list)} aggregate validation scopes."
        )
        if watch_rows:
            watched = ", ".join(
                f"{row.split} ({_display_metric_list(row.best_metrics)}, "
                f"{row.gap_percent_points:.2f} pp)"
                for row in watch_rows
            )
            lines.append(f"The following scopes have close alternatives to watch: {watched}.")
        if switch_rows:
            switched = ", ".join(
                f"{row.split} ({_display_metric_list(row.best_metrics)}, "
                f"{row.gap_percent_points:.2f} pp)"
                for row in switch_rows
            )
            lines.append(f"The following scopes exceed the promotion margin: {switched}.")
    lines.append("")
    return "\n".join(lines)


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def write_outputs(
    rows: list[SelectorDecisionRow],
    *,
    table_output: Path | None,
    csv_output: Path | None,
    json_output: Path | None,
    latex_output: Path | None,
    prose_output: Path | None,
) -> None:
    if not rows:
        raise ValueError("No selector-decision rows were generated.")
    if table_output is not None:
        table_output.parent.mkdir(parents=True, exist_ok=True)
        with table_output.open("w") as handle:
            write_table(rows, handle)
    if csv_output is not None:
        csv_output.parent.mkdir(parents=True, exist_ok=True)
        with csv_output.open("w", newline="") as handle:
            write_csv(rows, handle)
    if json_output is not None:
        json_output.parent.mkdir(parents=True, exist_ok=True)
        with json_output.open("w") as handle:
            write_json(rows, handle)
    if latex_output is not None:
        _write_text(latex_output, render_latex_table(rows))
    if prose_output is not None:
        _write_text(prose_output, render_prose(rows))


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read fair profile-selection comparison CSVs and summarize whether the "
            "current adaptive selector should be kept, watched, or replaced."
        )
    )
    parser.add_argument(
        "inputs",
        nargs="*",
        help=(
            "Comparison CSV inputs as PATH or LABEL=PATH. Defaults to the report, "
            "extended-holdout, and plateau-holdout comparison CSVs in work/."
        ),
    )
    parser.add_argument("--current-metric", default="adaptive")
    parser.add_argument(
        "--family",
        default="all",
        help="Comparison-summary family/group row to audit. Defaults to aggregate 'all'.",
    )
    parser.add_argument(
        "--tie-tolerance",
        type=float,
        default=0.005,
        help="Penalty-percent tolerance for treating metrics as tied.",
    )
    parser.add_argument(
        "--promotion-margin",
        type=float,
        default=0.5,
        help="Penalty-percent-point gap needed before recommending a switch.",
    )
    parser.add_argument("--table-output", type=Path, default=DEFAULT_TABLE_OUTPUT)
    parser.add_argument("--csv-output", type=Path, default=DEFAULT_CSV_OUTPUT)
    parser.add_argument("--json-output", type=Path, default=None)
    parser.add_argument("--latex-output", type=Path, default=DEFAULT_LATEX_OUTPUT)
    parser.add_argument("--prose-output", type=Path, default=DEFAULT_PROSE_OUTPUT)
    parser.add_argument("--skip-table", action="store_true")
    parser.add_argument("--skip-csv", action="store_true")
    parser.add_argument("--skip-latex", action="store_true")
    parser.add_argument("--skip-prose", action="store_true")
    parser.add_argument(
        "--stdout",
        action="store_true",
        help="Print the text table to stdout in addition to writing requested outputs.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    inputs = (
        [parse_input_spec(spec) for spec in args.inputs]
        if args.inputs
        else list(DEFAULT_INPUTS)
    )
    rows = summarize_decisions(
        inputs,
        current_metric=args.current_metric,
        family=args.family,
        tie_tolerance=args.tie_tolerance,
        promotion_margin=args.promotion_margin,
    )
    write_outputs(
        rows,
        table_output=None if args.skip_table else args.table_output,
        csv_output=None if args.skip_csv else args.csv_output,
        json_output=args.json_output,
        latex_output=None if args.skip_latex else args.latex_output,
        prose_output=None if args.skip_prose else args.prose_output,
    )
    if args.stdout:
        write_table(rows, sys.stdout)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
