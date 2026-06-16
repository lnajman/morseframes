#!/usr/bin/env python3
"""Render stage timing tables from the native GUDHI Simplex_tree benchmark CSV."""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "docs" / "native_gudhi_stage_profile_quick.csv"
DEFAULT_TABLE_OUTPUT = ROOT / "docs" / "native_gudhi_stage_profile_quick_table.tex"
DEFAULT_PROSE_OUTPUT = ROOT / "docs" / "native_gudhi_stage_profile_quick_prose.tex"
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
    "same-level-reduction": r"\texttt{same-level}",
}
STAGE_PROSE_LABELS = {
    "boundary": "boundary lookup",
    "coboundary": "coboundary construction",
    "extract": "simplex extraction",
    "frame-extra": "reference/reduction-input overhead",
    "order": "filtration ordering",
    "sequence": "sequence construction",
}


@dataclass(frozen=True)
class StageSummary:
    name: str
    strategy: str
    simplices: int
    view_build_ms: float
    view_extract_ms: float
    view_boundary_ms: float
    view_coboundary_ms: float
    view_order_ms: float
    sequence_ms: float
    frame_extra_ms: float
    reducer_ms: float
    direct_full_ms: float
    gudhi_ms: float

    @property
    def gudhi_over_direct(self) -> float:
        return self.gudhi_ms / self.direct_full_ms

    @property
    def gudhi_over_reducer(self) -> float:
        return self.gudhi_ms / self.reducer_ms

    @property
    def largest_pre_reducer_stage(self) -> str:
        stages = {
            "extract": self.view_extract_ms,
            "boundary": self.view_boundary_ms,
            "coboundary": self.view_coboundary_ms,
            "order": self.view_order_ms,
            "sequence": self.sequence_ms,
            "frame-extra": self.frame_extra_ms,
        }
        return max(stages, key=stages.__getitem__)

    @property
    def largest_pre_reducer_ms(self) -> float:
        return max(
            self.view_extract_ms,
            self.view_boundary_ms,
            self.view_coboundary_ms,
            self.view_order_ms,
            self.sequence_ms,
            self.frame_extra_ms,
        )

    def share(self, value: float) -> float:
        if self.direct_full_ms == 0.0:
            return 0.0
        return 100.0 * value / self.direct_full_ms


def _float(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    if value == "":
        raise ValueError(f"Missing required numeric column {key!r}.")
    return float(value)


def _float_or(row: dict[str, str], key: str, default: float) -> float:
    value = row.get(key, "")
    if value == "":
        return default
    return float(value)


def _int(row: dict[str, str], key: str) -> int:
    value = row.get(key, "")
    if value == "":
        raise ValueError(f"Missing required integer column {key!r}.")
    return int(value)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def summarize_rows(
    rows: Iterable[dict[str, str]],
    *,
    strategies: Iterable[str] = DEFAULT_STRATEGIES,
) -> list[StageSummary]:
    strategy_order = tuple(strategies)
    strategy_rank = {strategy: index for index, strategy in enumerate(strategy_order)}
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        strategy = row.get("strategy", "")
        if strategy in strategy_rank:
            grouped[(row["name"], strategy)].append(row)

    summaries: list[StageSummary] = []
    for key in sorted(grouped, key=lambda item: (item[0], strategy_rank[item[1]])):
        group = grouped[key]
        summaries.append(
            StageSummary(
                name=key[0],
                strategy=key[1],
                simplices=_int(group[0], "simplices"),
                view_build_ms=statistics.fmean(
                    _float(row, "view_build_ms") for row in group
                ),
                view_extract_ms=statistics.fmean(
                    _float_or(row, "view_extract_ms", _float(row, "view_build_ms"))
                    for row in group
                ),
                view_boundary_ms=statistics.fmean(
                    _float_or(row, "view_boundary_ms", 0.0) for row in group
                ),
                view_coboundary_ms=statistics.fmean(
                    _float_or(row, "view_coboundary_ms", 0.0) for row in group
                ),
                view_order_ms=statistics.fmean(
                    _float_or(row, "view_order_ms", 0.0) for row in group
                ),
                sequence_ms=statistics.fmean(
                    _float(row, "view_sequence_ms") for row in group
                ),
                frame_extra_ms=statistics.fmean(
                    _float(row, "view_frame_extra_ms") for row in group
                ),
                reducer_ms=statistics.fmean(
                    _float(row, "direct_reducer_ms") for row in group
                ),
                direct_full_ms=statistics.fmean(
                    _float(row, "direct_full_ms") for row in group
                ),
                gudhi_ms=statistics.fmean(
                    _float(row, "gudhi_persistence_ms") for row in group
                ),
            )
        )
    return summaries


def _texttt(value: str) -> str:
    return r"\texttt{" + value.replace("_", r"\_") + "}"


def _stage_count_summary(counts: dict[str, int]) -> str:
    parts = [
        f"{STAGE_PROSE_LABELS.get(stage, stage)} in {count} row"
        + ("" if count == 1 else "s")
        for stage, count in sorted(counts.items())
    ]
    if len(parts) == 1:
        return parts[0]
    return ", ".join(parts[:-1]) + f", and {parts[-1]}"


def render_latex_table(rows: list[StageSummary]) -> str:
    lines = [
        "% Generated by morseframes/tools/render_native_gudhi_stage_profile.py.",
        "% Re-run that script instead of editing this table by hand.",
        r"\begin{table}[ht]",
        r"\centering",
        r"\scriptsize",
        r"\setlength{\tabcolsep}{1.8pt}",
        r"\begin{tabular}{llrrrrrrrrl}",
        r"\toprule",
        (
            r"Case & Strategy & Extract & Bound. & Cobd. & Order & Seq. "
            r"& Frame & Reducer & Direct & Largest \\"
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
            f"{row.view_extract_ms:.2f} & "
            f"{row.view_boundary_ms:.2f} & "
            f"{row.view_coboundary_ms:.2f} & "
            f"{row.view_order_ms:.2f} & "
            f"{row.sequence_ms:.2f} & "
            f"{row.frame_extra_ms:.2f} & "
            f"{row.reducer_ms:.2f} & "
            f"{row.direct_full_ms:.2f} & "
            f"{_texttt(row.largest_pre_reducer_stage)} \\\\"
        )
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            (
                r"\caption{Native GUDHI plugin stage profile in milliseconds.  "
                r"The first four timing columns split construction of the read-only "
                r"\texttt{Simplex\_tree} view into simplex extraction, boundary lookup, "
                r"coboundary construction, and filtration ordering.  "
                r"\texttt{Seq.} times sequence construction alone.  "
                r"\texttt{Frame extra} is the fused reduction-input time minus "
                r"sequence-only time, so it estimates reference/reduction-input overhead.}"
            ),
            r"\label{tab:native-gudhi-stage-profile-quick}",
            r"\end{table}",
            "",
        ]
    )
    return "\n".join(lines)


def render_latex_prose(rows: list[StageSummary]) -> str:
    if not rows:
        raise ValueError("Cannot summarize an empty stage profile.")

    extract_share = statistics.fmean(row.share(row.view_extract_ms) for row in rows)
    boundary_share = statistics.fmean(row.share(row.view_boundary_ms) for row in rows)
    coboundary_share = statistics.fmean(row.share(row.view_coboundary_ms) for row in rows)
    order_share = statistics.fmean(row.share(row.view_order_ms) for row in rows)
    sequence_share = statistics.fmean(row.share(row.sequence_ms) for row in rows)
    frame_extra_share = statistics.fmean(row.share(row.frame_extra_ms) for row in rows)
    reducer_share = statistics.fmean(row.share(row.reducer_ms) for row in rows)
    largest_counts: dict[str, int] = defaultdict(int)
    for row in rows:
        largest_counts[row.largest_pre_reducer_stage] += 1
    largest_summary = _stage_count_summary(largest_counts)
    direct_ratios = [row.gudhi_over_direct for row in rows]
    reducer_ratios = [
        row.gudhi_over_reducer
        for row in rows
        if row.strategy in {"f-max", "f-min", "same-level-reduction"}
    ]

    return (
        "% Generated by morseframes/tools/render_native_gudhi_stage_profile.py.\n"
        "% Re-run that script instead of editing this paragraph by hand.\n"
        "In this quick native profile, the largest pre-reducer stage is "
        f"{largest_summary}, out of {len(rows)} summarized rows.  Averaged over the "
        "reported rows, the direct-path time is split approximately as "
        f"{extract_share:.1f}\\% simplex extraction, {boundary_share:.1f}\\% boundary "
        f"lookup, {coboundary_share:.1f}\\% coboundary construction, "
        f"{order_share:.1f}\\% filtration ordering, {sequence_share:.1f}\\% "
        f"sequence-only construction, {frame_extra_share:.1f}\\% reference/reduction-input "
        f"overhead, and {reducer_share:.1f}\\% reducer time.  The GUDHI/direct "
        "end-to-end ratio "
        f"ranges from {min(direct_ratios):.2f} to {max(direct_ratios):.2f}, while the "
        "main-strategy GUDHI/reducer ratio ranges from "
        f"{min(reducer_ratios):.2f} to {max(reducer_ratios):.2f}.  This split now shows "
        "a more balanced cost distribution: reducer time is no longer the single "
        "dominant component, while extraction and sequence construction remain the "
        "largest non-reducer costs.\n"
    )


def render_stdout_summary(rows: list[StageSummary]) -> str:
    largest_counts: dict[str, int] = defaultdict(int)
    for row in rows:
        largest_counts[row.largest_pre_reducer_stage] += 1
    return (
        "Largest pre-reducer stage counts: "
        + ", ".join(f"{stage}={count}" for stage, count in sorted(largest_counts.items()))
        + "\n"
        + render_latex_prose(rows)
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--table-output", type=Path, default=DEFAULT_TABLE_OUTPUT)
    parser.add_argument("--prose-output", type=Path, default=DEFAULT_PROSE_OUTPUT)
    parser.add_argument("--strategies", nargs="+", default=list(DEFAULT_STRATEGIES))
    parser.add_argument("--stdout", action="store_true")
    parser.add_argument("--summary", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    rows = summarize_rows(read_rows(args.input), strategies=args.strategies)
    if not rows:
        raise ValueError("No native GUDHI benchmark rows matched the requested filters.")

    table = render_latex_table(rows)
    prose = render_latex_prose(rows)
    args.table_output.parent.mkdir(parents=True, exist_ok=True)
    args.table_output.write_text(table)
    args.prose_output.parent.mkdir(parents=True, exist_ok=True)
    args.prose_output.write_text(prose)
    if args.stdout:
        sys.stdout.write(table)
    if args.summary:
        sys.stdout.write(render_stdout_summary(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
