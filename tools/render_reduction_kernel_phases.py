#!/usr/bin/env python3
"""Render the largest tested case in each RK phase-profile workload family."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import statistics


def largest_case_rows(data):
    if "completed_utc" not in data:
        raise ValueError("Cannot publish an incomplete benchmark")
    if not data["cases"] or not all(c["all_sequences_match"] for c in data["cases"]):
        raise ValueError("Missing cases or failed exact sequence checks")
    rows = []
    for mode in ("lower-star", "plateau"):
        for family in ("terrain", "volume"):
            candidates = [c for c in data["cases"] if c["mode"] == mode and c["family"] == family]
            if not candidates:
                continue
            size = max(c["size"] for c in candidates)
            cases = [c for c in candidates if c["size"] == size]
            for case in cases:
                if not {1, 8}.issubset({m["workers"] for m in case["measurements"]}):
                    raise ValueError("Table requires one- and eight-worker measurements")
            one = [next(m["summary"] for m in c["measurements"] if m["workers"] == 1) for c in cases]
            eight = [next(m["summary"] for m in c["measurements"] if m["workers"] == 8) for c in cases]
            row = dict(mode=mode, family=family, size=size, cases=len(cases),
                       simplices=cases[0]["num_simplices"],
                       one_ms=1e3 * statistics.median(m["performance_seconds"]["median"] for m in one),
                       eight_ms=1e3 * statistics.median(m["performance_seconds"]["median"] for m in eight),
                       speedup=statistics.median(m["speedup_vs_one"] for m in eight))
            for name in ("builder", "setup", "level_wall", "replay"):
                row[name + "_percent"] = 100 * statistics.median(
                    m["coarse"]["shares"][name + "_seconds"]["median"] for m in eight)
            row["other_percent"] = 100 * statistics.median(m["coarse"]["other_share"]["median"] for m in eight)
            rows.append(row)
    return rows


def render_table(data):
    lines = [r"\begin{tabular}{lrrrrrr}", r"\hline",
             r"Workload & $n$ & 1 worker (ms) & 8 workers (ms) & Speedup & Builder (\%) & Levels (\%) \\",
             r"\hline"]
    for row in largest_case_rows(data):
        label = ("Lower-star" if row["mode"] == "lower-star" else "Plateau") + (
            " 2D" if row["family"] == "terrain" else " 3D")
        lines.append(f"{label} & {row['size']} & {row['one_ms']:.2f} & {row['eight_ms']:.2f} & "
                     f"{row['speedup']:.3g} & {row['builder_percent']:.1f} & {row['level_wall_percent']:.1f} " + r"\\")
    lines.extend([r"\hline", r"\end{tabular}",
                  "% Time and speedup: uninstrumented; phase shares: separate coarse runs at 8 workers.",
                  "% Medians over case medians at the largest measured size in each family."])
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--table-output", type=Path, required=True)
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    args.table_output.write_text(render_table(data))
    print(json.dumps(largest_case_rows(data), indent=2))


if __name__ == "__main__":
    main()
