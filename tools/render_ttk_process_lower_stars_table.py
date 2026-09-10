#!/usr/bin/env python3
"""Render aggregate TTK ProcessLowerStars comparison ratios as LaTeX."""

from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=ROOT.parent / "ttk_process_lower_stars.csv"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "docs" / "ttk_process_lower_stars_table.tex",
    )
    return parser.parse_args()


def median(rows: list[dict[str, str]], field: str) -> float:
    return statistics.median(float(row[field]) for row in rows)


def main() -> int:
    args = parse_args()
    with args.input.open(newline="") as input_file:
        rows = list(csv.DictReader(input_file))
    if not rows:
        raise ValueError("No TTK comparison rows found")
    if any(row["critical_counts_match"] != "True" for row in rows):
        raise ValueError("Cannot render a table containing critical-count mismatches")

    labels = {
        "injective-terrain": "2D terrain",
        "injective-volume": "3D volume",
    }
    lines = [
        r"\begin{tabular}{llrrrr}",
        r"\toprule",
        r"Input & Workers & MF/F-Max & TTK/F-Max & MF/TTK & TTK+prep./F-Max \\",
        r"\midrule",
    ]
    for family in ("injective-terrain", "injective-volume"):
        family_rows = [row for row in rows if row["family"] == family]
        for workers in sorted({int(row["workers"]) for row in family_rows}):
            selected = [
                row for row in family_rows if int(row["workers"]) == workers
            ]
            ttk_initial = statistics.median(
                (
                    float(row["ttk_setup_seconds"])
                    + float(row["ttk_precondition_seconds"])
                    + float(row["ttk_process_lower_stars_seconds"])
                )
                / float(row["f_max_seconds"])
                for row in selected
            )
            lines.append(
                f"{labels[family]} & {workers} & "
                f"{median(selected, 'morseframes_ratio_vs_f_max'):.2f} & "
                f"{median(selected, 'ttk_ratio_vs_f_max'):.2f} & "
                f"{median(selected, 'morseframes_ratio_vs_ttk'):.2f} & "
                f"{ttk_initial:.2f} \\\\"
            )
        if family == "injective-terrain":
            lines.append(r"\midrule")
    lines.extend([r"\bottomrule", r"\end{tabular}", ""])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
