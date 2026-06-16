#!/usr/bin/env python3
"""Regenerate the fair multi-metric profile-validation artifacts."""

from __future__ import annotations

import argparse
import csv
import os
import shlex
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
DEFAULT_OUTPUT_DIR = ROOT.parent / "work"
DEFAULT_PROFILE_METRICS = (
    "profile_total_work",
    "critical_ratio",
    "adaptive",
    "adaptive_structured",
)
DEFAULT_TABLE_OUTPUT = ROOT / "docs" / "profile_metric_fair_validation_table.tex"
DEFAULT_PROSE_OUTPUT = ROOT / "docs" / "profile_metric_fair_validation_prose.tex"
DEFAULT_MANIFEST_OUTPUT = ROOT / "docs" / "fair_profile_validation_manifest.md"
VALIDATION_PRESETS = (
    "report",
    "extended-holdout",
    "plateau-holdout",
    "terrain-holdout",
    "image-holdout",
    "image-stress",
)
SPLIT_DISPLAY_NAMES = {
    "base": "Base",
    "base-holdout": "Base hold-out",
    "cam": r"CAM $S^4$",
    "cam-holdout": r"CAM $S^4$ hold-out",
    "plateau-holdout": "Plateau hold-out",
    "roadmap": "Roadmap",
    "roadmap-holdout": "Roadmap hold-out",
    "terrain-holdout": "Terrain hold-out",
    "image-holdout": "Image hold-out",
    "image-stress": "Image stress",
}
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
METRIC_LATEX_NAMES = {
    "profile_total_work": r"\texttt{profile-total-work}",
    "critical_ratio": r"\texttt{critical-ratio}",
    "adaptive": r"\texttt{adaptive}",
    "adaptive_structured": r"\texttt{adaptive-structured}",
    "estimated_reducer_work": r"\texttt{estimated-reducer-work}",
    "boundary_annotation_work": r"\texttt{boundary-annotation-work}",
    "initial_annotation_size": r"\texttt{initial-annotation-size}",
    "profile_seconds": r"\texttt{profile-seconds}",
    "working_set_size": r"\texttt{working-set-size}",
}
SMALL_NUMBER_WORDS = {
    0: "zero",
    1: "one",
    2: "two",
    3: "three",
    4: "four",
    5: "five",
    6: "six",
    7: "seven",
    8: "eight",
    9: "nine",
    10: "ten",
}


@dataclass(frozen=True)
class ValidationSplit:
    label: str
    families: tuple[str, ...]
    sizes: tuple[int, ...]
    seeds: tuple[int, ...]
    output_name: str


@dataclass(frozen=True)
class FairMetricTableRow:
    split: str
    requested_metric: str
    effective_metric: str
    cases: int
    match_percent: float
    unavailable_count: int
    avg_penalty_percent: float
    max_penalty_percent: float
    avg_selected_rank: float


def default_validation_splits() -> tuple[ValidationSplit, ...]:
    return validation_splits("report")


def validation_splits(preset: str) -> tuple[ValidationSplit, ...]:
    normalized = preset.lower().replace("_", "-")
    if normalized == "report":
        return (
            ValidationSplit(
                label="base",
                families=("lower-star", "plateau", "rips"),
                sizes=(16, 24, 32),
                seeds=(10, 11, 12),
                output_name="profile_vs_measured_base_validation_multi_metric.csv",
            ),
            ValidationSplit(
                label="cam",
                families=("cam-s4-rips",),
                sizes=(16, 24, 32),
                seeds=(10, 11, 12),
                output_name="profile_vs_measured_cam_validation_multi_metric.csv",
            ),
            ValidationSplit(
                label="roadmap",
                families=("roadmap-rips",),
                sizes=(50, 100, 104),
                seeds=(0,),
                output_name="profile_vs_measured_roadmap_validation_multi_metric.csv",
            ),
        )
    if normalized == "extended-holdout":
        return (
            ValidationSplit(
                label="base-holdout",
                families=("lower-star", "plateau", "rips"),
                sizes=(20, 28, 36, 44),
                seeds=(20, 21, 22, 23, 24),
                output_name="profile_vs_measured_base_extended_holdout_multi_metric.csv",
            ),
            ValidationSplit(
                label="cam-holdout",
                families=("cam-s4-rips",),
                sizes=(20, 28, 36),
                seeds=(20, 21, 22),
                output_name="profile_vs_measured_cam_extended_holdout_multi_metric.csv",
            ),
            ValidationSplit(
                label="roadmap-holdout",
                families=("roadmap-rips",),
                sizes=(50, 100, 104),
                seeds=(0,),
                output_name="profile_vs_measured_roadmap_extended_holdout_multi_metric.csv",
            ),
        )
    if normalized == "plateau-holdout":
        return (
            ValidationSplit(
                label="plateau-holdout",
                families=("plateau",),
                sizes=(12, 16, 20, 24, 28, 32, 36, 40, 44, 48),
                seeds=(30, 31, 32, 33, 34, 35, 36, 37, 38, 39),
                output_name="profile_vs_measured_plateau_holdout_multi_metric.csv",
            ),
        )
    if normalized == "terrain-holdout":
        return (
            ValidationSplit(
                label="terrain-holdout",
                families=("terrain",),
                sizes=(6, 8, 10, 12),
                seeds=(40, 41, 42, 43, 44),
                output_name="profile_vs_measured_terrain_holdout_multi_metric.csv",
            ),
        )
    if normalized == "image-holdout":
        return (
            ValidationSplit(
                label="image-holdout",
                families=("image-grid",),
                sizes=(6, 8, 10, 12),
                seeds=(50, 51, 52, 53, 54),
                output_name="profile_vs_measured_image_holdout_multi_metric.csv",
            ),
        )
    if normalized == "image-stress":
        return (
            ValidationSplit(
                label="image-stress",
                families=("image-grid",),
                sizes=(14, 18, 22, 26),
                seeds=(60, 61, 62, 63, 64),
                output_name="profile_vs_measured_image_stress_multi_metric.csv",
            ),
        )
    supported = ", ".join(VALIDATION_PRESETS)
    raise ValueError(f"Unknown validation preset {preset!r}. Supported: {supported}.")


def default_profile_metrics(validation_preset: str) -> tuple[str, ...]:
    normalized = validation_preset.lower().replace("_", "-")
    if normalized in {
        "plateau-holdout",
        "terrain-holdout",
        "image-holdout",
        "image-stress",
    }:
        return (
            "profile_total_work",
            "critical_ratio",
            "adaptive",
            "adaptive_structured",
            "estimated_reducer_work",
            "boundary_annotation_work",
            "initial_annotation_size",
            "profile_seconds",
            "working_set_size",
        )
    return DEFAULT_PROFILE_METRICS


def artifact_slug(validation_preset: str) -> str:
    normalized = validation_preset.lower().replace("_", "-")
    if normalized == "report":
        return "fair_validation"
    return normalized.replace("-", "_")


def default_table_output(validation_preset: str) -> Path:
    if validation_preset == "report":
        return DEFAULT_TABLE_OUTPUT
    return ROOT / "docs" / f"profile_metric_{artifact_slug(validation_preset)}_table.tex"


def default_prose_output(validation_preset: str) -> Path:
    if validation_preset == "report":
        return DEFAULT_PROSE_OUTPUT
    return ROOT / "docs" / f"profile_metric_{artifact_slug(validation_preset)}_prose.tex"


def default_manifest_output(validation_preset: str) -> Path:
    if validation_preset == "report":
        return DEFAULT_MANIFEST_OUTPUT
    return ROOT / "docs" / f"profile_metric_{artifact_slug(validation_preset)}_manifest.md"


def resolve_generated_outputs(args: argparse.Namespace) -> None:
    if args.table_output is None:
        args.table_output = default_table_output(args.validation_preset)
    if args.prose_output is None:
        args.prose_output = default_prose_output(args.validation_preset)
    if args.manifest_output is None:
        args.manifest_output = default_manifest_output(args.validation_preset)


def summary_output_name(args: argparse.Namespace, *, effective: bool, output_format: str) -> str:
    suffix = "effective_comparison" if effective else "comparison"
    return f"profile_metric_{artifact_slug(args.validation_preset)}_{suffix}.{output_format}"


def _strings(values: Iterable[int | str]) -> list[str]:
    return [str(value) for value in values]


def build_benchmark_command(
    *,
    python_executable: str,
    split: ValidationSplit,
    output_dir: Path,
    repeats: int,
    profile_repeats: int,
    sequence_algorithm: str,
    profile_candidate_gate: str,
    profile_selection_metrics: Iterable[str],
    roadmap_cache: Path,
    time_gudhi_cam: bool,
    download_roadmap_data: bool,
) -> list[str]:
    command = [
        python_executable,
        str(TOOLS / "benchmark_persistence.py"),
        "--mode",
        "profile-vs-measured",
        "--families",
        *split.families,
        "--sizes",
        *_strings(split.sizes),
        "--seeds",
        *_strings(split.seeds),
        "--repeats",
        str(repeats),
        "--profile-repeats",
        str(profile_repeats),
        "--sequence-algorithm",
        sequence_algorithm,
        "--profile-selection-metrics",
        *profile_selection_metrics,
        "--profile-candidate-gate",
        profile_candidate_gate,
        "--roadmap-cache",
        str(roadmap_cache),
        "--format",
        "csv",
        "--output",
        str(output_dir / split.output_name),
    ]
    if time_gudhi_cam:
        command.append("--time-gudhi-cam")
    if download_roadmap_data:
        command.append("--download-roadmap-data")
    return command


def build_summary_command(
    *,
    python_executable: str,
    splits: Iterable[ValidationSplit],
    output_dir: Path,
    label_column: str,
    output_format: str,
    output_name: str,
) -> list[str]:
    inputs = [
        f"{split.label}={output_dir / split.output_name}"
        for split in splits
    ]
    return [
        python_executable,
        str(TOOLS / "compare_profile_gates.py"),
        *inputs,
        "--label-column",
        label_column,
        "--format",
        output_format,
        "--output",
        str(output_dir / output_name),
    ]


def build_all_commands(args: argparse.Namespace) -> list[list[str]]:
    output_dir = args.output_dir
    splits = validation_splits(args.validation_preset)
    commands: list[list[str]] = []
    if not args.summaries_only:
        commands.extend(
            build_benchmark_command(
                python_executable=args.python_executable,
                split=split,
                output_dir=output_dir,
                repeats=args.repeats,
                profile_repeats=args.profile_repeats,
                sequence_algorithm=args.sequence_algorithm,
                profile_candidate_gate=args.profile_candidate_gate,
                profile_selection_metrics=args.profile_selection_metrics,
                roadmap_cache=args.roadmap_cache,
                time_gudhi_cam=not args.skip_gudhi_cam,
                download_roadmap_data=args.download_roadmap_data,
            )
            for split in splits
        )

    commands.extend(
        [
            build_summary_command(
                python_executable=args.python_executable,
                splits=splits,
                output_dir=output_dir,
                label_column="profile_selection_metric",
                output_format="table",
                output_name=summary_output_name(
                    args,
                    effective=False,
                    output_format="txt",
                ),
            ),
            build_summary_command(
                python_executable=args.python_executable,
                splits=splits,
                output_dir=output_dir,
                label_column="profile_selection_metric",
                output_format="csv",
                output_name=summary_output_name(
                    args,
                    effective=False,
                    output_format="csv",
                ),
            ),
            build_summary_command(
                python_executable=args.python_executable,
                splits=splits,
                output_dir=output_dir,
                label_column="effective_profile_selection_metric",
                output_format="table",
                output_name=summary_output_name(
                    args,
                    effective=True,
                    output_format="txt",
                ),
            ),
            build_summary_command(
                python_executable=args.python_executable,
                splits=splits,
                output_dir=output_dir,
                label_column="effective_profile_selection_metric",
                output_format="csv",
                output_name=summary_output_name(
                    args,
                    effective=True,
                    output_format="csv",
                ),
            ),
        ]
    )
    return commands


def format_command(command: Iterable[str]) -> str:
    return " ".join(shlex.quote(part) for part in command)


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


def _optional_bool(row: dict[str, str], key: str) -> bool:
    return row.get(key, "").strip().lower() == "true"


def _display_metric(metric: str) -> str:
    normalized = metric.lower().replace("-", "_")
    if normalized in METRIC_DISPLAY_NAMES:
        return METRIC_DISPLAY_NAMES[normalized]
    return metric.replace("_", r"\_")


def _display_effective_metric(metrics: Iterable[str]) -> str:
    values = sorted({_display_metric(metric) for metric in metrics})
    return "/".join(values)


def summarize_fair_metric_rows(
    *,
    splits: Iterable[ValidationSplit],
    output_dir: Path,
    metric_order: Iterable[str],
) -> list[FairMetricTableRow]:
    ordered_metrics = list(metric_order)
    summaries: list[FairMetricTableRow] = []
    for split in splits:
        path = output_dir / split.output_name
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            grouped[row["profile_selection_metric"]].append(row)

        metric_names = [
            metric for metric in ordered_metrics
            if metric in grouped
        ]
        metric_names.extend(
            metric for metric in sorted(grouped)
            if metric not in metric_names
        )
        for requested_metric in metric_names:
            group = grouped[requested_metric]
            cases = len(group)
            if cases == 0:
                continue
            summaries.append(
                FairMetricTableRow(
                    split=split.label,
                    requested_metric=requested_metric,
                    effective_metric=_display_effective_metric(
                        row["effective_profile_selection_metric"] for row in group
                    ),
                    cases=cases,
                    match_percent=(
                        100.0
                        * sum(1 for row in group if _optional_bool(row, "profile_matches_measured"))
                        / cases
                    ),
                    unavailable_count=sum(
                        1 for row in group
                        if _optional_int(row, "measured_best_profile_rank") == 0
                    ),
                    avg_penalty_percent=(
                        sum(_optional_float(row, "profile_penalty_percent") for row in group)
                        / cases
                    ),
                    max_penalty_percent=max(
                        _optional_float(row, "profile_penalty_percent") for row in group
                    ),
                    avg_selected_rank=(
                        sum(_optional_float(row, "profile_selected_measured_rank") for row in group)
                        / cases
                    ),
                )
            )
    return summaries


def render_fair_metric_table(rows: Iterable[FairMetricTableRow]) -> str:
    lines = [
        "% Generated by morseframes/tools/run_fair_profile_validation.py.",
        "% Re-run that script instead of editing this table by hand.",
        r"\begin{table}[ht]",
        r"\centering",
        r"\small",
        r"\begin{tabular}{llrrrrrr}",
        r"\toprule",
        r"Split & Metric (eff.) & Cases & Match \% & Unavail. & Avg. pen. \% & Max pen. \% & Rank \\",
        r"\midrule",
    ]
    previous_split: str | None = None
    for row in rows:
        if previous_split is not None and previous_split != row.split:
            lines.append(r"\midrule")
        metric = f"{_display_metric(row.requested_metric)} ({row.effective_metric})"
        lines.append(
            f"{SPLIT_DISPLAY_NAMES.get(row.split, row.split)} & "
            f"{metric} & "
            f"{row.cases:d} & "
            f"{row.match_percent:.1f} & "
            f"{row.unavailable_count:d} & "
            f"{row.avg_penalty_percent:.2f} & "
            f"{row.max_penalty_percent:.2f} & "
            f"{row.avg_selected_rank:.2f} " + r"\\"
        )
        previous_split = row.split

    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            r"\caption{Fair multi-metric validation of profile-selection scores with the current",
            r"\texttt{family-aware} candidate gate.  PTW denotes \texttt{profile-total-work}, CR denotes",
            r"\texttt{critical-ratio}, Adapt denotes \texttt{adaptive}, and Adapt-S denotes",
            r"\texttt{adaptive-structured}; the metric in parentheses is",
            r"the effective metric after resolving the adaptive policy.  Within each split, all",
            r"requested metrics are evaluated against the same profiles and the same measured portfolio",
            r"timings.  ``Unavail.'' counts cases where the measured fastest strategy was not in the",
            r"profiled candidate set.  ``Rank'' is the average measured-time rank of the",
            r"profile-selected strategy.}",
            r"\label{tab:fair-profile-metrics}",
            r"\end{table}",
            "",
        ]
    )
    return "\n".join(lines)


def _metric_latex(metric: str) -> str:
    normalized = metric.lower().replace("-", "_")
    if normalized in METRIC_LATEX_NAMES:
        return METRIC_LATEX_NAMES[normalized]
    return r"\texttt{" + metric.replace("_", "-") + "}"


def _percent(value: float) -> str:
    return f"${value:.2f}\\%$"


def _same_value(*values: float) -> bool:
    if not values:
        return True
    first = values[0]
    return all(abs(value - first) < 0.005 for value in values[1:])


def _case_count(value: int) -> str:
    return SMALL_NUMBER_WORDS.get(value, str(value))


def _split_role(split: str) -> str:
    if split.startswith("base"):
        return "base"
    if split.startswith("cam"):
        return "cam"
    if split.startswith("plateau"):
        return "plateau"
    if split.startswith("terrain"):
        return "terrain"
    if split.startswith("image"):
        return "image"
    if split.startswith("roadmap"):
        return "roadmap"
    return split


def _metric_penalty_list(rows: Iterable[FairMetricTableRow]) -> str:
    return ", ".join(
        f"{_metric_latex(row.requested_metric)} {_percent(row.avg_penalty_percent)}"
        for row in rows
    )


def render_fair_metric_prose(rows: Iterable[FairMetricTableRow]) -> str:
    row_list = list(rows)
    lookup = {
        (_split_role(row.split), row.requested_metric): row
        for row in row_list
    }
    roles = {role for role, _metric in lookup}
    lines = [
        "% Generated by morseframes/tools/run_fair_profile_validation.py.",
        "% Re-run that script instead of editing this paragraph by hand.",
    ]

    base_ptw = lookup.get(("base", "profile_total_work"))
    base_cr = lookup.get(("base", "critical_ratio"))
    base_adaptive = lookup.get(("base", "adaptive"))
    if base_ptw is not None and base_cr is not None and base_adaptive is not None:
        if _same_value(base_cr.avg_penalty_percent, base_adaptive.avg_penalty_percent):
            target = (
                f"{_percent(base_cr.avg_penalty_percent)} for "
                f"{_metric_latex('critical_ratio')}/{_metric_latex('adaptive')}"
            )
        else:
            target = (
                f"{_percent(base_cr.avg_penalty_percent)} for "
                f"{_metric_latex('critical_ratio')} and "
                f"{_percent(base_adaptive.avg_penalty_percent)} for "
                f"{_metric_latex('adaptive')}"
            )
        lines.append(
            "On this validation split with five timing repeats, the controlled comparison "
            f"reduces the base-family aggregate average penalty from "
            f"{_percent(base_ptw.avg_penalty_percent)} for "
            f"{_metric_latex('profile_total_work')} to {target}."
        )

    cam_rows = [
        lookup[key]
        for key in (
            ("cam", "profile_total_work"),
            ("cam", "critical_ratio"),
            ("cam", "adaptive"),
        )
        if key in lookup
    ]
    if len(cam_rows) == 3 and _same_value(
        *(row.avg_penalty_percent for row in cam_rows)
    ):
        lines.append(
            "For CAM-style Rips all three metrics select the same strategies on this "
            f"split, with {_percent(cam_rows[0].avg_penalty_percent)} average penalty."
        )
    elif cam_rows:
        lines.append(
            "For CAM-style Rips, the average penalties are "
            f"{_metric_penalty_list(cam_rows)}."
        )

    roadmap_ptw = lookup.get(("roadmap", "profile_total_work"))
    roadmap_cr = lookup.get(("roadmap", "critical_ratio"))
    roadmap_adaptive = lookup.get(("roadmap", "adaptive"))
    if roadmap_ptw is not None and roadmap_cr is not None and roadmap_adaptive is not None:
        if (
            _same_value(roadmap_ptw.avg_penalty_percent, 0.0)
            and _same_value(roadmap_adaptive.avg_penalty_percent, 0.0)
            and _same_value(roadmap_ptw.match_percent, 100.0)
            and _same_value(roadmap_adaptive.match_percent, 100.0)
        ):
            lines.append(
                f"For Roadmap Rips, {_metric_latex('profile_total_work')} and "
                f"{_metric_latex('adaptive')} both select the measured best in all "
                f"{_case_count(roadmap_ptw.cases)} cases, while fixed "
                f"{_metric_latex('critical_ratio')} has "
                f"{_percent(roadmap_cr.avg_penalty_percent)} average penalty."
            )
        else:
            lines.append(
                "For Roadmap Rips, the average penalties are "
                f"{_metric_penalty_list([roadmap_ptw, roadmap_cr, roadmap_adaptive])}."
            )

    plateau_rows = [
        row
        for (role, _metric), row in lookup.items()
        if role == "plateau"
    ]
    if plateau_rows:
        ranked = sorted(plateau_rows, key=lambda row: row.avg_penalty_percent)
        best = ranked[0]
        lines.append(
            "For the plateau-focused hold-out, the best average penalty is "
            f"{_percent(best.avg_penalty_percent)} from "
            f"{_metric_latex(best.requested_metric)}; the full ranking is "
            f"{_metric_penalty_list(ranked)}."
        )

    terrain_rows = [
        row
        for (role, _metric), row in lookup.items()
        if role == "terrain"
    ]
    if terrain_rows:
        ranked = sorted(terrain_rows, key=lambda row: row.avg_penalty_percent)
        best = ranked[0]
        lines.append(
            "For the terrain hold-out, the best average penalty is "
            f"{_percent(best.avg_penalty_percent)} from "
            f"{_metric_latex(best.requested_metric)}; the full ranking is "
            f"{_metric_penalty_list(ranked)}."
        )

    image_rows = [
        row
        for (role, _metric), row in lookup.items()
        if role == "image"
    ]
    if image_rows:
        ranked = sorted(image_rows, key=lambda row: row.avg_penalty_percent)
        best = ranked[0]
        split_names = sorted({row.split for row in image_rows})
        if len(split_names) == 1 and split_names[0] == "image-holdout":
            split_label = "image hold-out"
        elif len(split_names) == 1 and split_names[0] == "image-stress":
            split_label = "image stress split"
        elif len(split_names) == 1:
            split_label = SPLIT_DISPLAY_NAMES.get(split_names[0], split_names[0]).lower()
        else:
            split_label = "image-grid splits"
        lines.append(
            f"For the {split_label}, the best average penalty is "
            f"{_percent(best.avg_penalty_percent)} from "
            f"{_metric_latex(best.requested_metric)}; the full ranking is "
            f"{_metric_penalty_list(ranked)}."
        )

    if roles <= {"plateau"}:
        lines.append(
            "This plateau-focused diagnostic should be used to decide whether plateau "
            "needs a different effective metric from the rest of the base families."
        )
    elif roles <= {"terrain"}:
        lines.append(
            "This terrain-focused diagnostic tests whether the adaptive selector trained "
            "on random complexes transfers to a structured grid filtration with plateaus."
        )
    elif roles <= {"image"}:
        lines.append(
            "This image-grid diagnostic tests whether stronger blocky plateaus favor a "
            "different cheap selector from the current adaptive policy."
        )
    else:
        lines.append(
            "Thus adaptive scoring with the current \\texttt{family-aware} gate remains the "
            "most promising experimental selector, and future selector comparisons should use "
            "the multi-metric mode rather than independent timing runs."
        )
    lines.append("")
    return "\n".join(lines)


def write_fair_metric_table(
    *,
    splits: Iterable[ValidationSplit],
    output_dir: Path,
    metric_order: Iterable[str],
    table_output: Path,
) -> None:
    table_output.parent.mkdir(parents=True, exist_ok=True)
    table_rows = summarize_fair_metric_rows(
        splits=splits,
        output_dir=output_dir,
        metric_order=metric_order,
    )
    table_output.write_text(render_fair_metric_table(table_rows))


def write_fair_metric_prose(
    *,
    splits: Iterable[ValidationSplit],
    output_dir: Path,
    metric_order: Iterable[str],
    prose_output: Path,
) -> None:
    prose_output.parent.mkdir(parents=True, exist_ok=True)
    table_rows = summarize_fair_metric_rows(
        splits=splits,
        output_dir=output_dir,
        metric_order=metric_order,
    )
    prose_output.write_text(render_fair_metric_prose(table_rows))


def render_experiment_manifest(
    *,
    args: argparse.Namespace,
    splits: Iterable[ValidationSplit],
    commands: Iterable[list[str]],
    full_validation_commands: Iterable[list[str]] | None = None,
) -> str:
    split_list = list(splits)
    command_list = list(commands)
    full_command_list = (
        list(full_validation_commands)
        if full_validation_commands is not None
        else command_list
    )
    lines = [
        "# Fair Profile Validation Manifest",
        "",
        "Generated by `morseframes/tools/run_fair_profile_validation.py`.",
        "This file records the reproducible validation recipe used by the experiments note.",
        "",
        "## Validation Parameters",
        "",
        f"- Output directory: `{args.output_dir}`",
        f"- Roadmap cache: `{args.roadmap_cache}`",
        f"- Python executable: `{args.python_executable}`",
        f"- Inherited PYTHONPATH: `{os.environ.get('PYTHONPATH', '')}`",
        f"- Validation preset: `{args.validation_preset}`",
        f"- Timing repeats: `{args.repeats}`",
        f"- Profile repeats: `{args.profile_repeats}`",
        f"- Sequence algorithm setting: `{args.sequence_algorithm}`",
        f"- Profile candidate gate: `{args.profile_candidate_gate}`",
        f"- Profile selection metrics: `{', '.join(args.profile_selection_metrics)}`",
        f"- GUDHI/CAM timing: `{'disabled' if args.skip_gudhi_cam else 'enabled'}`",
        f"- Roadmap downloading: `{'enabled' if args.download_roadmap_data else 'disabled'}`",
        f"- Summaries-only mode: `{'enabled' if args.summaries_only else 'disabled'}`",
        "",
        "## Validation Splits",
        "",
        "| Split | Families | Sizes | Seeds | CSV artifact |",
        "| --- | --- | --- | --- | --- |",
    ]
    for split in split_list:
        lines.append(
            f"| `{split.label}` | "
            f"`{', '.join(split.families)}` | "
            f"`{', '.join(str(size) for size in split.sizes)}` | "
            f"`{', '.join(str(seed) for seed in split.seeds)}` | "
            f"`{args.output_dir / split.output_name}` |"
        )

    lines.extend(
        [
            "",
            "## Derived Artifacts",
            "",
            f"- Requested-metric summary table: `{args.output_dir / summary_output_name(args, effective=False, output_format='txt')}`",
            f"- Requested-metric summary CSV: `{args.output_dir / summary_output_name(args, effective=False, output_format='csv')}`",
            f"- Effective-metric summary table: `{args.output_dir / summary_output_name(args, effective=True, output_format='txt')}`",
            f"- Effective-metric summary CSV: `{args.output_dir / summary_output_name(args, effective=True, output_format='csv')}`",
            f"- LaTeX table fragment: `{args.table_output}`",
            f"- LaTeX prose fragment: `{args.prose_output}`",
            f"- This manifest: `{args.manifest_output}`",
            "",
            "## Environment Assumptions",
            "",
            "- Run from the repository workspace root, or pass absolute paths for outputs and caches.",
            "- The Python path must expose `morseframes/python`.",
            "- To include GUDHI/CAM timings, the Python path must also expose the local GUDHI build.",
            "- Roadmap datasets must already exist in the Roadmap cache, unless `--download-roadmap-data` is used.",
            "- Timing comparisons should be run sequentially on an otherwise quiet machine; the driver intentionally executes the benchmark splits in order.",
            "",
            "## Full Validation Recipe",
            "",
            "```bash",
        ]
    )
    lines.extend(format_command(command) for command in full_command_list)
    lines.extend(
        [
            "```",
            "",
        ]
    )
    if command_list != full_command_list:
        lines.extend(
            [
                "## Commands Executed by This Invocation",
                "",
                "This manifest was regenerated in `--summaries-only` mode, so the commands below "
                "refresh summaries and generated report fragments from existing CSV files.",
                "",
                "```bash",
            ]
        )
        lines.extend(format_command(command) for command in command_list)
        lines.extend(["```", ""])
    lines.extend(
        [
            "Use `--dry-run` to print the command recipe without executing timings.",
            "",
        ]
    )
    return "\n".join(lines)


def write_experiment_manifest(
    *,
    args: argparse.Namespace,
    splits: Iterable[ValidationSplit],
    commands: Iterable[list[str]],
    full_validation_commands: Iterable[list[str]] | None = None,
    manifest_output: Path,
) -> None:
    manifest_output.parent.mkdir(parents=True, exist_ok=True)
    manifest_output.write_text(
        render_experiment_manifest(
            args=args,
            splits=splits,
            commands=commands,
            full_validation_commands=full_validation_commands,
        )
    )


def run_commands(commands: Iterable[list[str]], *, dry_run: bool = False) -> None:
    for command in commands:
        print(f"$ {format_command(command)}", flush=True)
        if not dry_run:
            subprocess.run(command, check=True)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Regenerate the fair profile-selection validation CSVs and comparison "
            "summaries used by the experiments note."
        )
    )
    parser.add_argument(
        "--python-executable",
        default=sys.executable,
        help="Python executable used to run the benchmark and comparison tools.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for generated CSV and summary files.",
    )
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--profile-repeats", type=int, default=5)
    parser.add_argument("--sequence-algorithm", default="portfolio")
    parser.add_argument("--profile-candidate-gate", default="family-aware")
    parser.add_argument(
        "--validation-preset",
        choices=list(VALIDATION_PRESETS),
        default="report",
        help=(
            "Validation split recipe. 'report' reproduces the experiments-note table; "
            "'extended-holdout', 'plateau-holdout', 'terrain-holdout', and "
            "'image-holdout' write separate diagnostic hold-out artifacts."
        ),
    )
    parser.add_argument(
        "--profile-selection-metrics",
        nargs="+",
        help=(
            "Profile metrics evaluated fairly against the same measured timings. "
            "Defaults depend on --validation-preset."
        ),
    )
    parser.add_argument(
        "--roadmap-cache",
        type=Path,
        default=ROOT.parent / "work" / "roadmap-data",
    )
    parser.add_argument(
        "--skip-gudhi-cam",
        action="store_true",
        help="Do not request GUDHI/CAM timing during the validation runs.",
    )
    parser.add_argument(
        "--download-roadmap-data",
        action="store_true",
        help="Download missing Roadmap datasets before the Roadmap validation split.",
    )
    parser.add_argument(
        "--summaries-only",
        action="store_true",
        help="Skip benchmark runs and regenerate comparison summaries from existing CSVs.",
    )
    parser.add_argument(
        "--table-output",
        type=Path,
        default=None,
        help="LaTeX table fragment generated from the fair validation CSVs.",
    )
    parser.add_argument(
        "--skip-table",
        action="store_true",
        help="Do not regenerate the LaTeX table fragment.",
    )
    parser.add_argument(
        "--prose-output",
        type=Path,
        default=None,
        help="LaTeX prose fragment generated from the fair validation CSVs.",
    )
    parser.add_argument(
        "--skip-prose",
        action="store_true",
        help="Do not regenerate the LaTeX prose fragment.",
    )
    parser.add_argument(
        "--manifest-output",
        type=Path,
        default=None,
        help="Markdown manifest recording validation parameters and artifact paths.",
    )
    parser.add_argument(
        "--skip-manifest",
        action="store_true",
        help="Do not regenerate the experiment manifest.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing them.",
    )
    args = parser.parse_args(argv)
    if args.profile_selection_metrics is None:
        args.profile_selection_metrics = list(default_profile_metrics(args.validation_preset))
    resolve_generated_outputs(args)
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    splits = validation_splits(args.validation_preset)
    commands = build_all_commands(args)
    full_args = argparse.Namespace(**vars(args))
    full_args.summaries_only = False
    full_validation_commands = build_all_commands(full_args)
    run_commands(commands, dry_run=args.dry_run)
    if not args.skip_table:
        if args.dry_run:
            print(f"# would write LaTeX table fragment: {args.table_output}", flush=True)
        else:
            write_fair_metric_table(
                splits=splits,
                output_dir=args.output_dir,
                metric_order=args.profile_selection_metrics,
                table_output=args.table_output,
            )
    if not args.skip_prose:
        if args.dry_run:
            print(f"# would write LaTeX prose fragment: {args.prose_output}", flush=True)
        else:
            write_fair_metric_prose(
                splits=splits,
                output_dir=args.output_dir,
                metric_order=args.profile_selection_metrics,
                prose_output=args.prose_output,
            )
    if not args.skip_manifest:
        if args.dry_run:
            print(f"# would write experiment manifest: {args.manifest_output}", flush=True)
        else:
            write_experiment_manifest(
                args=args,
                splits=splits,
                commands=commands,
                full_validation_commands=full_validation_commands,
                manifest_output=args.manifest_output,
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
