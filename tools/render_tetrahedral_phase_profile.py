#!/usr/bin/env python3
"""Render tetrahedral parallel phase-profile tables."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT.parent / "tetrahedral_phase_profile.csv"
DEFAULT_TABLE = ROOT / "docs" / "tetrahedral_phase_profile_table.tex"

STRATEGY_ORDER = ("process-lower-stars", "reduction-kernel")
STRATEGY_LABELS = {
    "process-lower-stars": "ProcessLowerStars",
    "reduction-kernel": "Reduction kernel",
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as source:
        rows = list(csv.DictReader(source))
    if not rows:
        raise ValueError(f"No profile rows found in {path}")
    return rows


def _median(rows: list[dict[str, str]], column: str) -> float:
    return statistics.median(float(row[column]) for row in rows)


def _percent(value: float) -> str:
    return "--" if math.isnan(value) else f"{100.0 * value:.1f}"


def _ratio(value: float) -> str:
    return "--" if math.isnan(value) else f"{value:.2f}"


def render_table(rows: list[dict[str, str]], output: Path) -> None:
    grouped: dict[tuple[int, str, int], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[
            (int(row["grid_size"]), row["strategy"], int(row["max_workers"]))
        ].append(row)

    lines = [
        r"\begin{tabular}{rlrrrrrrrrr}",
        r"\toprule",
        (
            r"Grid & Strategy & Workers & Gradient (ms) & Build (\%) & M simplex/s & "
            r"Builder (\%) & Setup (\%) & Local (\%) & Replay (\%) & Task parallelism \\"
        ),
        r"\midrule",
    ]
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
                            f'{1000.0 * _median(group, "construction_seconds"):.2f}',
                            _percent(
                                _median(group, "sequence_build_seconds")
                                / _median(group, "construction_seconds")
                            ),
                            f'{1.0e-6 * _median(group, "simplices_per_second"):.2f}',
                            _percent(
                                _median(
                                    group,
                                    "process_lower_stars_builder_init_share",
                                )
                            ),
                            _percent(_median(group, "process_lower_stars_setup_share")),
                            _percent(
                                _median(group, "process_lower_stars_local_wall_share")
                            ),
                            _percent(
                                _median(group, "process_lower_stars_replay_share")
                            ),
                            _ratio(
                                _median(group, "process_lower_stars_task_parallelism")
                            ),
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
    parser.add_argument("--table-output", type=Path, default=DEFAULT_TABLE)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    render_table(read_rows(args.input), args.table_output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
