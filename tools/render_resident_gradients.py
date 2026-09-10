#!/usr/bin/env python3
"""Render construction-separated comparisons, retaining historical v1 totals."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import statistics

from benchmark_resident_gradients import ALGORITHMS, ALGORITHM_PHASES, CONSTRUCTION_PHASES, summarize

LABELS = {"f_max": "F-Max", "reduction_kernel": "RK", "ttk": "TTK"}
PHASE_ORDER = {
    "f_max": ("representation_and_filtration", "builder_setup", "gradient"),
    "reduction_kernel": ("representation_and_filtration", "builder_setup", "gradient"),
    "ttk": ("native_object_init", "vertex_order", "representation_setup",
            "connectivity_precondition", "gradient"),
}


def validated_groups(data):
    if data.get("schema") not in ("resident-gradient-study-v1", "resident-gradient-study-v2") or not data.get("completed_utc"):
        raise ValueError("Expected a completed resident-input benchmark, not legacy kernel timings")
    split = data["schema"] == "resident-gradient-study-v2"
    if split and (data.get("construction_phases") != {a: sorted(CONSTRUCTION_PHASES[a]) for a in ALGORITHMS}
                  or data.get("algorithm_phases") != {a: sorted(ALGORITHM_PHASES[a]) for a in ALGORITHMS}):
        raise ValueError("Missing or changed construction/algorithm boundary")
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
            expected_schema = "resident-gradient-v2" if split else "resident-gradient-v1"
            if measurement["raw"].get("schema") != expected_schema:
                raise ValueError("Mixed timing schemas")
            if measurement["raw"]["workers"] != measurement["workers"]:
                raise ValueError("Worker metadata mismatch")
            # Recompute from raw evidence instead of trusting stored summaries.
            summary = summarize(measurement["raw"], data["arguments"]["repeats"],
                                data["arguments"]["diagnostics"])
            key = (case["family"], case["size"], measurement["workers"])
            groups.setdefault(key, []).append(summary)
    return groups


def render_tables(data, comparison=None):
    groups = validated_groups(data)
    split = data["schema"] == "resident-gradient-study-v2"
    comparison = comparison or ("algorithm" if split else "total")
    if comparison not in ("algorithm", "total") or (comparison == "algorithm" and not split):
        raise ValueError("Algorithm comparisons require v2 non-profiled phase samples")
    metric = "algorithm_seconds" if comparison == "algorithm" else "total_seconds"
    ratios = "algorithm_paired_ratios" if comparison == "algorithm" else "paired_ratios"
    header = r"Input & $n$ & Workers & F-Max (ms) & RK (ms) & TTK (ms) & RK/TTK \\"
    if split:
        label = "algorithm" if comparison == "algorithm" else "total"
        header = (f"Input & $n$ & Workers & F-Max {label} (ms) & RK {label} (ms) & "
                  f"TTK {label} (ms) & RK/TTK " + r"\\")
    totals = [r"\begin{tabular}{llrrrrr}", r"\hline",
              header, r"\hline"]
    phases = [r"\begin{tabular}{lrrllr}", r"\hline",
              r"Input & $n$ & Workers & Algorithm & Phase & Time (ms) \\", r"\hline"]
    for (family, size, workers), summaries in sorted(groups.items()):
        label = "2D" if family == "terrain" else "3D"
        times = [1e3 * statistics.median(s[a][metric]["median"] for s in summaries)
                 for a in ALGORITHMS]
        ratio = statistics.median(s[ratios]["reduction_kernel/ttk"]["median"]
                                  for s in summaries)
        totals.append(f"{label} & {size} & {workers} & " +
                      " & ".join(f"{v:.3f}" for v in times) + f" & {ratio:.3f} " + r"\\")
        # Keep the phase fragment compact: largest case per family, 1/8 workers.
        if workers not in (1, 8) or size != max(k[1] for k in groups if k[0] == family):
            continue
        for algorithm in ALGORITHMS:
            total = 1e3 * statistics.median(
                s[algorithm]["total_seconds" if split else "diagnostic_total_seconds"]["median"] for s in summaries)
            prefix = f"{label} & {size} & {workers} & {LABELS[algorithm]}"
            total_label = "Total (boundary clocks)" if split else "Total (diagnostic)"
            phases.append(f"{prefix} & {total_label} & {total:.3f} " + r"\\")
            for phase in PHASE_ORDER[algorithm]:
                value = 1e3 * statistics.median(
                    s[algorithm]["performance_phases_seconds" if split else "phases_seconds"][phase]["median"] for s in summaries)
                phases.append(f"{prefix} & {phase.replace('_', ' ')} & {value:.3f} " + r"\\")
    scope = ("% Algorithm time excludes native construction; includes MF builder and TTK ordering/lower stars."
             if comparison == "algorithm" else
             "% Resident-to-gradient totals include all native preparation.")
    totals.extend([r"\hline", r"\end{tabular}",
                   scope if split else "% Uninstrumented totals from common resident arrays, all native preparation included.",
                   "% Medians over case medians; ratios use paired per-repetition ratios before aggregation."])
    phases.extend([r"\hline", r"\end{tabular}",
                   ("% Non-profiled phase-boundary samples, largest input per family at 1/8 workers." if split else
                    "% Separate diagnostic runs, largest input per family at 1/8 workers."),
                   ("% Phases partition each raw performance total; medians and rounded values need not add." if split else
                    "% Phases partition each raw diagnostic total; medians and rounded values need not add."),
                   "% Nested RK/F-Max details remain in raw JSON; do not add them to parent phases.",
                   "% TTK gradient includes lower-star construction and matching together."])
    return "\n".join(totals) + "\n", "\n".join(phases) + "\n"


def render_construction(data):
    groups = validated_groups(data)
    if data["schema"] != "resident-gradient-study-v2":
        raise ValueError("Construction reporting requires v2 non-profiled phase samples")
    rows = [r"\begin{tabular}{llrrrrr}", r"\hline",
            r"Input & $n$ & Workers & F-Max construction (ms) & RK construction (ms) & TTK construction (ms) & Shared loading (ms) \\", r"\hline"]
    for (family, size, workers), summaries in sorted(groups.items()):
        times = [1e3 * statistics.median(s[a]["construction_seconds"]["median"] for s in summaries)
                 for a in ALGORITHMS]
        loading = 1e3 * statistics.median(
            m["raw"]["input_loading_seconds"] for c in data["cases"]
            if c["family"] == family and c["size"] == size
            for m in c["measurements"] if m["workers"] == workers)
        label = "2D" if family == "terrain" else "3D"
        rows.append(f"{label} & {size} & {workers} & " +
                    " & ".join(f"{v:.3f}" for v in times + [loading]) + " " + r"\\")
    rows.extend([r"\hline", r"\end{tabular}",
                 "% Native construction only; excluded from the algorithm comparison.",
                 "% MF: representation/filtration. TTK: object/representation/connectivity; ordering stays in algorithm time.",
                 "% Shared loading: one parse per native invocation, medians across seeds, not per algorithm or cold-cache I/O.",
                 "% Per-sample component sums are exact; separately aggregated medians need not add."])
    return "\n".join(rows) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--table-output", type=Path, required=True)
    parser.add_argument("--phases-output", type=Path, required=True)
    parser.add_argument("--construction-output", type=Path)
    parser.add_argument("--comparison", choices=("algorithm", "total"))
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    totals, phases = render_tables(data, args.comparison)
    construction = render_construction(data) if args.construction_output else None
    args.table_output.write_text(totals)
    args.phases_output.write_text(phases)
    if args.construction_output:
        args.construction_output.write_text(construction)


if __name__ == "__main__":
    main()
