#!/usr/bin/env python3
"""Compare profile-vs-measured CSV outputs for candidate-gate experiments."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, TextIO

from calibrate_profile_gate import shape_bin


@dataclass(frozen=True)
class GateComparisonRow:
    label: str
    family: str
    cases: int
    avg_profile_candidates: float
    match_percent: float
    unavailable_count: int
    unavailable_percent: float
    avg_penalty_percent: float
    max_penalty_percent: float
    avg_selected_rank: float
    profile_overhead_percent: float


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


def read_profile_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def parse_input_spec(spec: str) -> tuple[str, Path]:
    if "=" in spec:
        label, path = spec.split("=", 1)
        label = label.strip()
        if not label:
            raise ValueError(f"Invalid input spec {spec!r}: empty label")
        return label, Path(path)
    path = Path(spec)
    return path.stem, path


def group_key(row: dict[str, str], group_by: str) -> str:
    normalized = group_by.lower().replace("_", "-")
    family = row.get("family", "unknown-family")
    size = row.get("size", "unknown-size")
    if normalized == "family":
        return family
    if normalized == "size":
        return f"n={size}"
    if normalized == "family-size":
        return f"{family}:n={size}"
    if normalized == "shape-bin":
        return shape_bin(row)
    if normalized == "family-shape":
        return f"{family}:{shape_bin(row)}"
    raise ValueError(
        f"Unknown grouping {group_by!r}. "
        "Supported: family, size, family-size, shape-bin, family-shape."
    )


def _summarize_group(label: str, group: str, rows: list[dict[str, str]]) -> GateComparisonRow:
    cases = len(rows)
    if cases == 0:
        raise ValueError("Cannot summarize an empty row group.")

    matches = sum(1 for row in rows if _optional_bool(row, "profile_matches_measured"))
    unavailable = sum(1 for row in rows if _optional_int(row, "measured_best_profile_rank") == 0)
    total_profile_seconds = sum(_optional_float(row, "total_profile_seconds") for row in rows)
    total_measured_seconds = sum(
        _optional_float(row, "total_measured_morse_seconds") for row in rows
    )
    return GateComparisonRow(
        label=label,
        family=group,
        cases=cases,
        avg_profile_candidates=(
            sum(_optional_float(row, "profile_candidate_count") for row in rows) / cases
        ),
        match_percent=100.0 * matches / cases,
        unavailable_count=unavailable,
        unavailable_percent=100.0 * unavailable / cases,
        avg_penalty_percent=(
            sum(_optional_float(row, "profile_penalty_percent") for row in rows) / cases
        ),
        max_penalty_percent=max(
            _optional_float(row, "profile_penalty_percent") for row in rows
        ),
        avg_selected_rank=(
            sum(_optional_float(row, "profile_selected_measured_rank") for row in rows) / cases
        ),
        profile_overhead_percent=(
            100.0 * total_profile_seconds / total_measured_seconds
            if total_measured_seconds > 0.0
            else 0.0
        ),
    )


def summarize_gate_rows(
    label: str,
    rows: Iterable[dict[str, str]],
    *,
    include_all: bool = True,
    group_by: str = "family",
) -> list[GateComparisonRow]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    all_rows = list(rows)
    for row in all_rows:
        grouped[group_key(row, group_by)].append(row)

    summaries = [
        _summarize_group(label, family, group)
        for family, group in sorted(grouped.items())
    ]
    if include_all and all_rows:
        summaries.append(_summarize_group(label, "all", all_rows))
    return summaries


def compare_gate_outputs(
    inputs: Iterable[tuple[str, Path]],
    *,
    include_all: bool = True,
    group_by: str = "family",
    label_column: str | None = None,
) -> list[GateComparisonRow]:
    input_list = list(inputs)
    rows: list[GateComparisonRow] = []
    for label, path in input_list:
        profile_rows = read_profile_rows(path)
        if label_column is None:
            rows.extend(
                summarize_gate_rows(
                    label,
                    profile_rows,
                    include_all=include_all,
                    group_by=group_by,
                )
            )
            continue

        labeled_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in profile_rows:
            labeled_rows[row.get(label_column, "")].append(row)
        for column_value, group_rows in sorted(labeled_rows.items()):
            if not column_value:
                column_value = f"missing-{label_column}"
            row_label = column_value if len(input_list) == 1 else f"{label}:{column_value}"
            rows.extend(
                summarize_gate_rows(
                    row_label,
                    group_rows,
                    include_all=include_all,
                    group_by=group_by,
                )
            )
    return rows


def write_table(rows: list[GateComparisonRow], output: TextIO) -> None:
    headers = [
        "label",
        "group",
        "cases",
        "cands",
        "match%",
        "unavail",
        "unav%",
        "avg_pen%",
        "max_pen%",
        "rank",
        "prof_over%",
    ]
    output.write(
        f"{headers[0]:>18} {headers[1]:>14} "
        + " ".join(f"{header:>10}" for header in headers[2:])
        + "\n"
    )
    for row in rows:
        output.write(
            f"{row.label:>18} "
            f"{row.family:>14} "
            f"{row.cases:10d} "
            f"{row.avg_profile_candidates:10.2f} "
            f"{row.match_percent:10.1f} "
            f"{row.unavailable_count:10d} "
            f"{row.unavailable_percent:10.1f} "
            f"{row.avg_penalty_percent:10.2f} "
            f"{row.max_penalty_percent:10.2f} "
            f"{row.avg_selected_rank:10.2f} "
            f"{row.profile_overhead_percent:10.1f}\n"
        )


def write_csv(rows: list[GateComparisonRow], output: TextIO) -> None:
    writer = csv.DictWriter(output, fieldnames=list(asdict(rows[0]).keys()))
    writer.writeheader()
    for row in rows:
        writer.writerow(asdict(row))


def write_json(rows: list[GateComparisonRow], output: TextIO) -> None:
    json.dump([asdict(row) for row in rows], output, indent=2)
    output.write("\n")


def write_rows(
    rows: list[GateComparisonRow],
    *,
    output_format: str,
    output_path: Path | None,
) -> None:
    if not rows:
        raise ValueError("No profile-gate comparison rows were generated.")

    handle = None
    if output_path is None:
        output = sys.stdout
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        handle = output_path.open("w", newline="")
        output = handle

    try:
        if output_format == "table":
            write_table(rows, output)
        elif output_format == "csv":
            write_csv(rows, output)
        elif output_format == "json":
            write_json(rows, output)
        else:
            raise ValueError(f"Unknown output format: {output_format}")
    finally:
        if handle is not None:
            handle.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "inputs",
        nargs="+",
        help="CSV inputs as PATH or LABEL=PATH.",
    )
    parser.add_argument(
        "--no-all",
        action="store_true",
        help="Do not add an aggregate 'all' row for each input.",
    )
    parser.add_argument(
        "--group-by",
        default="family",
        choices=["family", "size", "family-size", "shape-bin", "family-shape"],
        help="Grouping used for per-input summaries.",
    )
    parser.add_argument(
        "--label-column",
        help=(
            "Use values from this CSV column as comparison labels, useful when a single "
            "CSV contains several selector metrics."
        ),
    )
    parser.add_argument("--format", choices=["table", "csv", "json"], default="table")
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows = compare_gate_outputs(
        [parse_input_spec(spec) for spec in args.inputs],
        include_all=not args.no_all,
        group_by=args.group_by,
        label_column=args.label_column,
    )
    write_rows(rows, output_format=args.format, output_path=args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
