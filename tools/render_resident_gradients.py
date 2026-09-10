#!/usr/bin/env python3
"""Render resident-input total and phase tables, never mixing legacy timings."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import statistics

from benchmark_resident_gradients import ALGORITHMS, summarize

LABELS = {"f_max": "F-Max", "reduction_kernel": "RK", "ttk": "TTK"}
PHASE_ORDER = {
    "f_max": ("representation_and_filtration", "builder_setup", "gradient"),
    "reduction_kernel": ("representation_and_filtration", "builder_setup", "gradient"),
    "ttk": ("native_object_init", "vertex_order", "representation_setup",
            "connectivity_precondition", "gradient"),
}


def validated_groups(data):
    if data.get("schema") != "resident-gradient-study-v1" or not data.get("completed_utc"):
        raise ValueError("Expected a completed resident-input benchmark, not legacy kernel timings")
    if not data.get("cases"):
        raise ValueError("Empty benchmark")
    args = data["arguments"]
    expected_cases = {(family, size, seed)
                      for family in ("terrain", "volume")
                      for size in args[family + "_sizes"] for seed in args["seeds"]}
    actual_cases = [(c["family"], c["size"], c["seed"]) for c in data["cases"]]
    if len(set(actual_cases)) != len(actual_cases) or set(actual_cases) != expected_cases:
        raise ValueError("Missing or duplicate input cases")
    groups = {}
    for case in data["cases"]:
        measurements = case["measurements"]
        if not measurements or len({m["workers"] for m in measurements}) != len(measurements):
            raise ValueError("Missing or duplicate worker measurements")
        if {m["workers"] for m in measurements} != set(data["arguments"]["workers"]):
            raise ValueError("Incomplete worker configurations")
        for measurement in measurements:
            if measurement["raw"]["workers"] != measurement["workers"]:
                raise ValueError("Worker metadata mismatch")
            # Recompute from raw evidence instead of trusting stored summaries.
            summary = summarize(measurement["raw"], data["arguments"]["repeats"],
                                data["arguments"]["diagnostics"])
            key = (case["family"], case["size"], measurement["workers"])
            groups.setdefault(key, []).append(summary)
    return groups


def render_tables(data):
    groups = validated_groups(data)
    totals = [r"\begin{tabular}{llrrrrr}", r"\hline",
              r"Input & $n$ & Workers & F-Max (ms) & RK (ms) & TTK (ms) & RK/TTK \\",
              r"\hline"]
    phases = [r"\begin{tabular}{lrrllr}", r"\hline",
              r"Input & $n$ & Workers & Algorithm & Phase & Time (ms) \\", r"\hline"]
    for (family, size, workers), summaries in sorted(groups.items()):
        label = "2D" if family == "terrain" else "3D"
        times = [1e3 * statistics.median(s[a]["total_seconds"]["median"] for s in summaries)
                 for a in ALGORITHMS]
        ratio = statistics.median(s["paired_ratios"]["reduction_kernel/ttk"]["median"]
                                  for s in summaries)
        totals.append(f"{label} & {size} & {workers} & " +
                      " & ".join(f"{v:.3f}" for v in times) + f" & {ratio:.3f} " + r"\\")
        # Keep the phase fragment compact: largest case per family, 1/8 workers.
        if workers not in (1, 8) or size != max(k[1] for k in groups if k[0] == family):
            continue
        for algorithm in ALGORITHMS:
            total = 1e3 * statistics.median(
                s[algorithm]["diagnostic_total_seconds"]["median"] for s in summaries)
            prefix = f"{label} & {size} & {workers} & {LABELS[algorithm]}"
            phases.append(f"{prefix} & Total (diagnostic) & {total:.3f} " + r"\\")
            for phase in PHASE_ORDER[algorithm]:
                value = 1e3 * statistics.median(
                    s[algorithm]["phases_seconds"][phase]["median"] for s in summaries)
                phases.append(f"{prefix} & {phase.replace('_', ' ')} & {value:.3f} " + r"\\")
    totals.extend([r"\hline", r"\end{tabular}",
                   "% Uninstrumented totals from common resident arrays, all native preparation included.",
                   "% Medians over case medians; ratios use paired per-repetition ratios before aggregation."])
    phases.extend([r"\hline", r"\end{tabular}",
                   "% Separate diagnostic runs, largest input per family at 1/8 workers.",
                   "% Phases partition each raw diagnostic total; medians and rounded values need not add.",
                   "% Nested RK/F-Max details remain in raw JSON; do not add them to parent phases.",
                   "% TTK gradient includes lower-star construction and matching together."])
    return "\n".join(totals) + "\n", "\n".join(phases) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--table-output", type=Path, required=True)
    parser.add_argument("--phases-output", type=Path, required=True)
    args = parser.parse_args()
    totals, phases = render_tables(json.loads(args.input.read_text()))
    args.table_output.write_text(totals)
    args.phases_output.write_text(phases)


if __name__ == "__main__":
    main()
