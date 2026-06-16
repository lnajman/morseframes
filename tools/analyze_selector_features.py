#!/usr/bin/env python3
"""Diagnose profile-selector margins against cheap complex descriptors."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, TextIO


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = ROOT.parent / "work"
DEFAULT_INPUTS = (
    ("base", DEFAULT_OUTPUT_DIR / "profile_vs_measured_base_validation_multi_metric.csv"),
    ("cam", DEFAULT_OUTPUT_DIR / "profile_vs_measured_cam_validation_multi_metric.csv"),
    ("roadmap", DEFAULT_OUTPUT_DIR / "profile_vs_measured_roadmap_validation_multi_metric.csv"),
    ("base-holdout", DEFAULT_OUTPUT_DIR / "profile_vs_measured_base_extended_holdout_multi_metric.csv"),
    ("cam-holdout", DEFAULT_OUTPUT_DIR / "profile_vs_measured_cam_extended_holdout_multi_metric.csv"),
    ("roadmap-holdout", DEFAULT_OUTPUT_DIR / "profile_vs_measured_roadmap_extended_holdout_multi_metric.csv"),
    ("plateau-holdout", DEFAULT_OUTPUT_DIR / "profile_vs_measured_plateau_holdout_multi_metric.csv"),
    ("terrain-holdout", DEFAULT_OUTPUT_DIR / "profile_vs_measured_terrain_holdout_multi_metric.csv"),
    ("image-holdout", DEFAULT_OUTPUT_DIR / "profile_vs_measured_image_holdout_multi_metric.csv"),
    ("image-stress", DEFAULT_OUTPUT_DIR / "profile_vs_measured_image_stress_multi_metric.csv"),
)
DEFAULT_TABLE_OUTPUT = DEFAULT_OUTPUT_DIR / "selector_feature_diagnostic.txt"
DEFAULT_CSV_OUTPUT = DEFAULT_OUTPUT_DIR / "selector_feature_diagnostic.csv"
DEFAULT_JSON_OUTPUT = DEFAULT_OUTPUT_DIR / "selector_feature_diagnostic.json"
DEFAULT_LATEX_OUTPUT = ROOT / "docs" / "selector_feature_diagnostic_table.tex"
DEFAULT_PROSE_OUTPUT = ROOT / "docs" / "selector_feature_diagnostic_prose.tex"

SPLIT_DISPLAY_NAMES = {
    "base": "Base",
    "base-holdout": "Base hold-out",
    "cam": "CAM",
    "cam-holdout": "CAM hold-out",
    "roadmap": "Roadmap",
    "roadmap-holdout": "Roadmap hold-out",
    "plateau-holdout": "Plateau hold-out",
    "terrain-holdout": "Terrain hold-out",
    "image-holdout": "Image hold-out",
    "image-stress": "Image stress",
    "all": "All",
}

FEATURE_COLUMNS = (
    "num_simplices",
    "shape_simplex_vertex_ratio",
    "shape_largest_level_size",
    "shape_largest_level_ratio",
    "shape_level_concentration",
    "shape_level_entropy_ratio",
)


@dataclass(frozen=True)
class SelectorCase:
    split: str
    family: str
    name: str
    seed: str
    size: str
    critical_metric: str
    time_metric: str
    critical_algorithm: str
    time_algorithm: str
    measured_best_algorithm: str
    critical_penalty_percent: float
    time_penalty_percent: float
    margin_percent_points: float
    winner: str
    num_simplices: float | None
    shape_simplex_vertex_ratio: float | None
    shape_largest_level_size: float | None
    shape_largest_level_ratio: float | None
    shape_level_concentration: float | None
    shape_level_entropy_ratio: float | None


@dataclass(frozen=True)
class SelectorFeatureSummary:
    split: str
    family: str
    cases: int
    time_wins: int
    critical_wins: int
    ties: int
    critical_avg_penalty_percent: float
    time_avg_penalty_percent: float
    avg_margin_percent_points: float
    max_time_gain_percent_points: float
    max_critical_gain_percent_points: float
    descriptor_coverage_percent: float
    avg_num_simplices: float | None
    avg_simplex_vertex_ratio: float | None
    avg_largest_level_size: float | None
    avg_largest_level_ratio: float | None
    avg_level_concentration: float | None
    avg_level_entropy_ratio: float | None


def _optional_float(row: dict[str, str], key: str) -> float | None:
    value = row.get(key, "")
    if value == "":
        return None
    return float(value)


def _required_float(row: dict[str, str], key: str) -> float:
    value = _optional_float(row, key)
    if value is None:
        raise ValueError(f"Missing required numeric column {key!r}.")
    return value


def _case_key(row: dict[str, str]) -> tuple[str, str, str, str]:
    return (
        row.get("family", ""),
        row.get("name", ""),
        row.get("seed", ""),
        row.get("size", ""),
    )


def parse_input_spec(spec: str) -> tuple[str, Path]:
    if "=" in spec:
        label, path = spec.split("=", 1)
        label = label.strip()
        if not label:
            raise ValueError(f"Invalid input spec {spec!r}: empty label.")
        return label, Path(path)
    path = Path(spec)
    return path.stem, path


def read_profile_rows(inputs: Iterable[tuple[str, Path]]) -> list[tuple[str, dict[str, str]]]:
    rows: list[tuple[str, dict[str, str]]] = []
    for split, path in inputs:
        with path.open(newline="") as handle:
            rows.extend((split, row) for row in csv.DictReader(handle))
    return rows


def build_selector_cases(
    inputs: Iterable[tuple[str, Path]],
    *,
    critical_metric: str = "critical_ratio",
    time_metrics: Iterable[str] = ("profile_total_work", "profile_seconds"),
    tie_tolerance: float = 0.005,
) -> list[SelectorCase]:
    time_metric_set = set(time_metrics)
    grouped: dict[tuple[str, str, str, str, str], dict[str, dict[str, str]]] = defaultdict(dict)
    for split, row in read_profile_rows(inputs):
        metric = row.get("profile_selection_metric", "")
        if metric == critical_metric or metric in time_metric_set:
            grouped[(split,) + _case_key(row)][metric] = row

    cases: list[SelectorCase] = []
    for (split, family, name, seed, size), rows_by_metric in sorted(grouped.items()):
        critical_row = rows_by_metric.get(critical_metric)
        time_rows = [
            row for metric, row in rows_by_metric.items()
            if metric in time_metric_set
        ]
        if critical_row is None or not time_rows:
            continue
        time_row = min(
            time_rows,
            key=lambda row: (
                _required_float(row, "profile_penalty_percent"),
                row.get("profile_selection_metric", ""),
            ),
        )
        critical_penalty = _required_float(critical_row, "profile_penalty_percent")
        time_penalty = _required_float(time_row, "profile_penalty_percent")
        margin = critical_penalty - time_penalty
        if margin > tie_tolerance:
            winner = "time"
        elif margin < -tie_tolerance:
            winner = "critical"
        else:
            winner = "tie"

        features = {
            column: _optional_float(critical_row, column)
            for column in FEATURE_COLUMNS
        }
        cases.append(
            SelectorCase(
                split=split,
                family=family,
                name=name,
                seed=seed,
                size=size,
                critical_metric=critical_metric,
                time_metric=time_row.get("profile_selection_metric", ""),
                critical_algorithm=critical_row.get("profile_selected_algorithm", ""),
                time_algorithm=time_row.get("profile_selected_algorithm", ""),
                measured_best_algorithm=critical_row.get("measured_best_algorithm", ""),
                critical_penalty_percent=critical_penalty,
                time_penalty_percent=time_penalty,
                margin_percent_points=margin,
                winner=winner,
                num_simplices=features["num_simplices"],
                shape_simplex_vertex_ratio=features["shape_simplex_vertex_ratio"],
                shape_largest_level_size=features["shape_largest_level_size"],
                shape_largest_level_ratio=features["shape_largest_level_ratio"],
                shape_level_concentration=features["shape_level_concentration"],
                shape_level_entropy_ratio=features["shape_level_entropy_ratio"],
            )
        )
    return cases


def _average(values: Iterable[float | None]) -> float | None:
    materialized = [value for value in values if value is not None]
    if not materialized:
        return None
    return sum(materialized) / float(len(materialized))


def _summarize_group(split: str, family: str, cases: list[SelectorCase]) -> SelectorFeatureSummary:
    if not cases:
        raise ValueError("Cannot summarize an empty selector-case group.")
    descriptor_rows = sum(
        1
        for case in cases
        if case.shape_level_concentration is not None
    )
    return SelectorFeatureSummary(
        split=split,
        family=family,
        cases=len(cases),
        time_wins=sum(1 for case in cases if case.winner == "time"),
        critical_wins=sum(1 for case in cases if case.winner == "critical"),
        ties=sum(1 for case in cases if case.winner == "tie"),
        critical_avg_penalty_percent=(
            sum(case.critical_penalty_percent for case in cases) / float(len(cases))
        ),
        time_avg_penalty_percent=(
            sum(case.time_penalty_percent for case in cases) / float(len(cases))
        ),
        avg_margin_percent_points=(
            sum(case.margin_percent_points for case in cases) / float(len(cases))
        ),
        max_time_gain_percent_points=max(case.margin_percent_points for case in cases),
        max_critical_gain_percent_points=(
            -min(case.margin_percent_points for case in cases)
        ),
        descriptor_coverage_percent=100.0 * descriptor_rows / float(len(cases)),
        avg_num_simplices=_average(case.num_simplices for case in cases),
        avg_simplex_vertex_ratio=_average(
            case.shape_simplex_vertex_ratio for case in cases
        ),
        avg_largest_level_size=_average(case.shape_largest_level_size for case in cases),
        avg_largest_level_ratio=_average(case.shape_largest_level_ratio for case in cases),
        avg_level_concentration=_average(
            case.shape_level_concentration for case in cases
        ),
        avg_level_entropy_ratio=_average(
            case.shape_level_entropy_ratio for case in cases
        ),
    )


def summarize_selector_cases(
    cases: Iterable[SelectorCase],
    *,
    include_all: bool = True,
) -> list[SelectorFeatureSummary]:
    case_list = list(cases)
    grouped: dict[tuple[str, str], list[SelectorCase]] = defaultdict(list)
    for case in case_list:
        grouped[(case.split, case.family)].append(case)
    summaries = [
        _summarize_group(split, family, group)
        for (split, family), group in sorted(grouped.items())
    ]
    if include_all and case_list:
        summaries.append(_summarize_group("all", "all", case_list))
    return summaries


def _format_optional(value: float | None, digits: int = 2) -> str:
    if value is None or math.isnan(value):
        return ""
    return f"{value:.{digits}f}"


def write_summary_table(rows: list[SelectorFeatureSummary], output: TextIO) -> None:
    headers = [
        "split",
        "family",
        "cases",
        "time",
        "crit",
        "tie",
        "crit%",
        "time%",
        "margin",
        "desc%",
        "simp",
        "s/v",
        "levelconc",
    ]
    output.write(
        f"{headers[0]:>16} {headers[1]:>14} "
        + " ".join(f"{header:>9}" for header in headers[2:])
        + "\n"
    )
    for row in rows:
        output.write(
            f"{row.split:>16} "
            f"{row.family:>14} "
            f"{row.cases:9d} "
            f"{row.time_wins:9d} "
            f"{row.critical_wins:9d} "
            f"{row.ties:9d} "
            f"{row.critical_avg_penalty_percent:9.2f} "
            f"{row.time_avg_penalty_percent:9.2f} "
            f"{row.avg_margin_percent_points:9.2f} "
            f"{row.descriptor_coverage_percent:9.1f} "
            f"{_format_optional(row.avg_num_simplices, 1):>9} "
            f"{_format_optional(row.avg_simplex_vertex_ratio, 2):>9} "
            f"{_format_optional(row.avg_level_concentration, 3):>9}\n"
        )


def write_summary_csv(rows: list[SelectorFeatureSummary], output: TextIO) -> None:
    writer = csv.DictWriter(output, fieldnames=list(asdict(rows[0]).keys()))
    writer.writeheader()
    for row in rows:
        writer.writerow(asdict(row))


def write_summary_json(rows: list[SelectorFeatureSummary], output: TextIO) -> None:
    json.dump([asdict(row) for row in rows], output, indent=2)
    output.write("\n")


def _latex_escape(value: str) -> str:
    return value.replace("_", r"\_")


def _latex_number(value: float | None, digits: int = 2) -> str:
    if value is None:
        return "--"
    return f"{value:.{digits}f}"


def render_latex_table(rows: list[SelectorFeatureSummary]) -> str:
    selected = [
        row
        for row in rows
        if row.split in {
            "plateau-holdout",
            "terrain-holdout",
            "image-holdout",
            "image-stress",
            "all",
        }
    ]
    lines = [
        "% Generated by morseframes/tools/analyze_selector_features.py.",
        "% Re-run that script instead of editing this table by hand.",
        r"\begin{table}[ht]",
        r"\centering",
        r"\small",
        r"\begin{tabular}{llrrrrrr}",
        r"\toprule",
        r"Split & Family & Cases & TW wins & CR wins & Ties & CR pen. & TW pen. \\",
        r"\midrule",
    ]
    for row in selected:
        lines.append(
            f"{_latex_escape(SPLIT_DISPLAY_NAMES.get(row.split, row.split))} & "
            f"{_latex_escape(row.family)} & "
            f"{row.cases} & "
            f"{row.time_wins} & "
            f"{row.critical_wins} & "
            f"{row.ties} & "
            f"{row.critical_avg_penalty_percent:.2f} & "
            f"{row.time_avg_penalty_percent:.2f} \\\\"
        )
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            (
                r"\caption{Selector-feature diagnostic comparing requested "
                r"\texttt{critical-ratio} against the best requested total-work/time "
                r"metric on each validation case.  Positive evidence for the structured "
                r"branch appears when the TW penalty is lower than the CR penalty.}"
            ),
            r"\label{tab:selector-feature-diagnostic}",
            r"\end{table}",
            "",
        ]
    )
    return "\n".join(lines)


def render_prose(rows: list[SelectorFeatureSummary]) -> str:
    lookup = {(row.split, row.family): row for row in rows}
    all_row = lookup.get(("all", "all"))
    image_holdout = lookup.get(("image-holdout", "image-grid"))
    image_stress = lookup.get(("image-stress", "image-grid"))
    terrain = lookup.get(("terrain-holdout", "terrain"))
    lines = [
        "% Generated by morseframes/tools/analyze_selector_features.py.",
        "% Re-run that script instead of editing this paragraph by hand.",
    ]
    if all_row is not None:
        lines.append(
            "The selector-feature diagnostic compares requested "
            r"\texttt{critical-ratio} with the best requested total-work/time metric "
            f"on {all_row.cases} validation cases.  The aggregate average penalties are "
            f"${all_row.critical_avg_penalty_percent:.2f}\\%$ for critical-ratio and "
            f"${all_row.time_avg_penalty_percent:.2f}\\%$ for total-work/time."
        )
    if image_holdout is not None and image_stress is not None:
        lines.append(
            "For image-grid filtrations, the total-work/time branch is best on the small "
            f"hold-out by {image_holdout.avg_margin_percent_points:.2f} percentage points "
            f"and on the larger stress split by {image_stress.avg_margin_percent_points:.2f} "
            "percentage points.  This supports keeping the structured branch as a "
            "diagnostic candidate while requiring more validation before promoting it."
        )
        if (
            image_holdout.avg_simplex_vertex_ratio is not None
            and image_stress.avg_simplex_vertex_ratio is not None
        ):
            lines.append(
                "The regenerated image rows now expose descriptor coverage for this "
                f"comparison: the average simplex-to-vertex ratio rises from "
                f"{image_holdout.avg_simplex_vertex_ratio:.2f} on the hold-out to "
                f"{image_stress.avg_simplex_vertex_ratio:.2f} on the stress split, "
                "which gives the next calibration pass a concrete scale variable to test."
            )
    if terrain is not None:
        lines.append(
            "The terrain hold-out remains a useful control because it has grid structure "
            f"but does not show the same large total-work/time margin "
            f"({terrain.avg_margin_percent_points:.2f} percentage points)."
        )
    lines.append("")
    return "\n".join(lines)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def write_table(path: Path, rows: list[SelectorFeatureSummary]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        write_summary_table(rows, handle)


def write_csv(path: Path, rows: list[SelectorFeatureSummary]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        write_summary_csv(rows, handle)


def write_json(path: Path, rows: list[SelectorFeatureSummary]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        write_summary_json(rows, handle)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare critical-ratio profile selection with total-work/time selection "
            "and summarize the margins against cheap shape descriptors."
        )
    )
    parser.add_argument(
        "inputs",
        nargs="*",
        help=(
            "Optional CSV inputs as PATH or LABEL=PATH. Defaults to the validation "
            "artifacts in work/."
        ),
    )
    parser.add_argument("--critical-metric", default="critical_ratio")
    parser.add_argument(
        "--time-metrics",
        nargs="+",
        default=["profile_total_work", "profile_seconds"],
    )
    parser.add_argument("--tie-tolerance", type=float, default=0.005)
    parser.add_argument("--table-output", type=Path, default=DEFAULT_TABLE_OUTPUT)
    parser.add_argument("--csv-output", type=Path, default=DEFAULT_CSV_OUTPUT)
    parser.add_argument("--json-output", type=Path, default=DEFAULT_JSON_OUTPUT)
    parser.add_argument("--latex-output", type=Path, default=DEFAULT_LATEX_OUTPUT)
    parser.add_argument("--prose-output", type=Path, default=DEFAULT_PROSE_OUTPUT)
    parser.add_argument("--stdout", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    inputs = (
        tuple(parse_input_spec(spec) for spec in args.inputs)
        if args.inputs
        else DEFAULT_INPUTS
    )
    cases = build_selector_cases(
        inputs,
        critical_metric=args.critical_metric,
        time_metrics=args.time_metrics,
        tie_tolerance=args.tie_tolerance,
    )
    summaries = summarize_selector_cases(cases)
    write_table(args.table_output, summaries)
    write_csv(args.csv_output, summaries)
    write_json(args.json_output, summaries)
    write_text(args.latex_output, render_latex_table(summaries))
    write_text(args.prose_output, render_prose(summaries))
    if args.stdout:
        write_summary_table(summaries, sys.stdout)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
