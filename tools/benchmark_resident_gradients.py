#!/usr/bin/env python3
"""Compare fresh RK, F-Max, TTK and optional PLS gradients from resident arrays.

Report native construction separately and compare the remaining algorithm
work. Fresh resident-to-gradient totals remain available. File parsing is
measured once per native invocation, outside all per-algorithm times. Internal
phase diagnostics use separate runs and never enter performance ratios.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import random
import statistics
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from benchmark_reduction_kernel_ab import header_digest
from benchmark_ttk_process_lower_stars import TTK_REVISION, write_ttk_input
import benchmark_simplicial_strategies as generators

ALGORITHMS = ("f_max", "reduction_kernel", "ttk")
PLS_ALGORITHMS = ALGORITHMS + ("process_lower_stars",)
OUTER_PHASES = {
    "f_max": {"representation_and_filtration", "builder_setup", "gradient"},
    "reduction_kernel": {"representation_and_filtration", "builder_setup", "gradient"},
    "process_lower_stars": {"representation_and_filtration", "builder_setup", "gradient"},
    "ttk": {"native_object_init", "vertex_order", "representation_setup",
            "connectivity_precondition", "gradient"},
}
CONSTRUCTION_PHASES = {
    "f_max": {"representation_and_filtration"},
    "reduction_kernel": {"representation_and_filtration"},
    "process_lower_stars": {"representation_and_filtration"},
    "ttk": {"native_object_init", "representation_setup", "connectivity_precondition"},
}
ALGORITHM_PHASES = {a: OUTER_PHASES[a] - CONSTRUCTION_PHASES[a] for a in PLS_ALGORITHMS}
DETAIL_PHASES = {
    "f_max": {"workspace_init", "candidate_seeding", "candidate_selection",
              "emission_and_updates", "callbacks"},
    "reduction_kernel": {"workspace_and_pool", "level_processing", "replay"},
    "process_lower_stars": {"lower_star_setup", "local_processing", "replay"},
    "ttk": set(),
}
ORDERS = ((0, 1, 2), (2, 1, 0), (1, 2, 0), (0, 2, 1), (2, 0, 1), (1, 0, 2))
PLS_ORDERS = ((0, 1, 3, 2), (1, 2, 0, 3), (2, 3, 1, 0), (3, 0, 2, 1))


def algorithm_names(schema):
    return PLS_ALGORITHMS if schema in ("resident-gradient-v3", "resident-gradient-study-v3") else ALGORITHMS


def validate_pls_orders(result, repeats):
    actual = result.get("performance_orders", [])
    if repeats % 4 or len(actual) != repeats:
        raise ValueError("PLS comparison requires complete four-method order cycles")
    for offset in range(4):
        expected = [[PLS_ALGORITHMS[a] for a in PLS_ORDERS[(i + offset) % 4]]
                    for i in range(repeats)]
        if actual == expected:
            return
    raise ValueError("Incorrect or unbalanced four-method execution order")


def paired_distribution(numerator, denominator, bootstrap=False):
    ratios = [x / y for x, y in zip(numerator, denominator, strict=True)]
    result = distribution(ratios, positive=True)
    if bootstrap:
        rng = random.Random(0)
        boot = sorted(statistics.median(rng.choices(ratios, k=len(ratios)))
                      for _ in range(2000))
        result["paired_repetition_bootstrap_95_interval"] = [boot[49], boot[1949]]
    return result


def distribution(values, *, positive=False):
    if not values or any(not math.isfinite(v) or v < 0 or (positive and v == 0)
                         for v in values):
        raise ValueError("Expected finite nonnegative samples (positive for elapsed totals)")
    q = (statistics.quantiles(values, n=4, method="inclusive")
         if len(values) > 1 else [values[0]] * 3)
    return dict(median=statistics.median(values), q1=q[0], q3=q[2],
                min=min(values), max=max(values))


def validate_outer(outer, total, algorithm):
    if set(outer) != OUTER_PHASES[algorithm]:
        raise ValueError("Missing or unexpected phase measurements")
    distribution(list(outer.values()))
    if not math.isclose(sum(outer.values()), total, rel_tol=1e-9, abs_tol=1e-12):
        raise ValueError("Outer phases do not account for the total")
    if outer["gradient"] <= 0:
        raise ValueError("Gradient duration must be positive")


def performance_components(result, algorithm):
    """Partition each non-profiled sample, never subtract aggregate medians."""
    if result.get("schema") not in ("resident-gradient-v2", "resident-gradient-v3"):
        raise ValueError("Construction-separated comparisons require v2/v3 performance phases")
    raw = result["algorithms"][algorithm]
    totals, rows = raw["performance_seconds"], raw.get("performance_phases_seconds", [])
    if len(rows) != len(totals):
        raise ValueError("Incomplete performance phases")
    components = {k: [] for k in ("construction_seconds", "algorithm_seconds", "gradient_seconds")}
    for total, outer in zip(totals, rows, strict=True):
        validate_outer(outer, total, algorithm)
        components["construction_seconds"].append(sum(outer[k] for k in CONSTRUCTION_PHASES[algorithm]))
        components["algorithm_seconds"].append(sum(outer[k] for k in ALGORITHM_PHASES[algorithm]))
        components["gradient_seconds"].append(outer["gradient"])
    return components


def summarize(result, repeats, diagnostics):
    extended = result.get("schema") == "resident-gradient-v3"
    split = result.get("schema") in ("resident-gradient-v2", "resident-gradient-v3")
    if result.get("schema") not in ("resident-gradient-v1", "resident-gradient-v2", "resident-gradient-v3"):
        raise ValueError("Wrong timing schema: rebuild the resident-gradient executable")
    algorithms = algorithm_names(result["schema"])
    if extended:
        validate_pls_orders(result, repeats)
        if type(result.get("euler_characteristic")) is not int:
            raise ValueError("Missing Euler characteristic")
    if split:
        distribution([result.get("input_loading_seconds", -1)], positive=True)
    if result.get("ttk_revision") != TTK_REVISION:
        raise ValueError("Native benchmark was not built against the pinned TTK revision")
    if result.get("exact_reference_checks") is not True:
        raise ValueError("Missing successful exact reference checks")
    if result.get("dimension") not in (1, 2, 3) or result.get("simplices", 0) <= 0:
        raise ValueError("Invalid complex metadata")
    if type(result.get("critical_counts_match")) is not bool:
        raise ValueError("Missing critical-count agreement flag")
    if set(result.get("algorithms", {})) != set(algorithms):
        raise ValueError("Missing algorithm measurements")
    summary = {}
    for algorithm in algorithms:
        raw = result["algorithms"][algorithm]
        counts = raw.get("critical_counts", [])
        if (len(counts) != result["dimension"] + 1 or
                any(type(c) is not int or c < 0 for c in counts) or
                sum(counts) > result["simplices"] or
                (result["simplices"] - sum(counts)) % 2):
            raise ValueError("Invalid critical-simplex counts")
        if extended and sum((-1) ** dim * n for dim, n in enumerate(counts)) != result["euler_characteristic"]:
            raise ValueError("Critical-cell Euler characteristic differs")
        performance = raw["performance_seconds"]
        rows = raw["diagnostics"]
        if len(performance) != repeats or len(rows) != diagnostics:
            raise ValueError("Incomplete timing repetitions")
        entry = {"total_seconds": distribution(performance, positive=True)}
        if split:
            entry.update({k: distribution(v, positive=True)
                          for k, v in performance_components(result, algorithm).items()})
            entry["performance_phases_seconds"] = {
                k: distribution([r[k] for r in raw["performance_phases_seconds"]])
                for k in OUTER_PHASES[algorithm]}
        totals, phases, details, shares = [], {}, {}, {}
        for row in rows:
            total = row["total_seconds"]
            distribution([total], positive=True)
            outer = row["phases_seconds"]
            inner = row["gradient_details_seconds"]
            if set(inner) != DETAIL_PHASES[algorithm]:
                raise ValueError("Missing or unexpected phase measurements")
            validate_outer(outer, total, algorithm)
            if inner:
                distribution(list(inner.values()))
                remaining = outer["gradient"] - sum(inner.values())
                if remaining < -max(1e-12, outer["gradient"] * 1e-9):
                    raise ValueError("Nested phases exceed gradient time (double counting)")
                # Residual is explicit, not silently attributed to local work.
                for name, value in {**inner, "unattributed": max(0.0, remaining)}.items():
                    details.setdefault(name, []).append(value)
            totals.append(total)
            for name, value in outer.items():
                phases.setdefault(name, []).append(value)
                shares.setdefault(name, []).append(value / total)
        entry["diagnostic_total_seconds"] = distribution(totals, positive=True)
        entry["phases_seconds"] = {k: distribution(v) for k, v in phases.items()}
        entry["phase_shares"] = {k: distribution(v) for k, v in shares.items()}
        entry["gradient_details_seconds"] = (
            {k: distribution(v) for k, v in details.items()} if details else None)
        summary[algorithm] = entry
    counts = [result["algorithms"][a]["critical_counts"] for a in algorithms]
    if result.get("critical_counts_match") != all(c == counts[0] for c in counts):
        raise ValueError("Critical-count agreement flag is inconsistent")
    # Different algorithms may legitimately create different critical counts.
    # Exact consistency with each algorithm's own reference is mandatory.
    summary["paired_ratios"] = {}
    if split:
        summary["algorithm_paired_ratios"] = {}
    comparisons = [("reduction_kernel", "ttk"), ("reduction_kernel", "f_max"), ("ttk", "f_max")]
    if extended:
        comparisons.extend([("process_lower_stars", "ttk"), ("process_lower_stars", "f_max"),
                            ("reduction_kernel", "process_lower_stars")])
    for numerator, denominator in comparisons:
        a = result["algorithms"][numerator]["performance_seconds"]
        b = result["algorithms"][denominator]["performance_seconds"]
        summary["paired_ratios"][f"{numerator}/{denominator}"] = paired_distribution(a, b, extended)
        if split:
            a = performance_components(result, numerator)["algorithm_seconds"]
            b = performance_components(result, denominator)["algorithm_seconds"]
            summary["algorithm_paired_ratios"][f"{numerator}/{denominator}"] = paired_distribution(a, b, extended)
    return summary


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", type=Path, required=True)
    parser.add_argument("--terrain-sizes", type=int, nargs="*", default=[16, 64])
    parser.add_argument("--volume-sizes", type=int, nargs="*", default=[8, 16])
    parser.add_argument("--seeds", type=int, nargs="+", default=[0, 2])
    parser.add_argument("--workers", type=int, nargs="+", default=[1, 2, 4, 8])
    parser.add_argument("--repeats", type=int, default=12)
    parser.add_argument("--diagnostics", type=int)
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--include-pls", action="store_true", help="Add native PLS using balanced four-method orders")
    parser.add_argument("--order-offset", type=int, default=0)
    parser.add_argument("--reverse-inputs", action="store_true")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    cycle = 4 if args.include_pls else 6
    if args.diagnostics is None:
        args.diagnostics = cycle
    if (not args.terrain_sizes and not args.volume_sizes) or any(
            n < 2 for n in args.terrain_sizes + args.volume_sizes):
        parser.error("Choose at least one input size, each >= 2")
    if (1 not in args.workers or min(args.workers) < 1 or
            len(set(args.workers)) != len(args.workers)):
        parser.error("Workers must be unique, positive, and include one")
    if (len(set(args.seeds)) != len(args.seeds) or
            len(set(args.terrain_sizes)) != len(args.terrain_sizes) or
            len(set(args.volume_sizes)) != len(args.volume_sizes)):
        parser.error("Sizes and seeds must be unique")
    if (args.repeats < cycle or args.repeats % cycle or args.diagnostics < cycle or
            args.diagnostics % cycle or args.warmups < 0 or not 0 <= args.order_offset < cycle):
        parser.error(f"Repetitions/diagnostics must be positive multiples of {cycle}; order offset must be 0..{cycle - 1}")
    if args.output.exists():
        parser.error("Output exists; choose a fresh path to preserve evidence")
    return args


def command(*args):
    return subprocess.check_output(args, cwd=ROOT, text=True).strip()


def run_native(executable, path, workers, repeats, diagnostics, warmups,
               include_pls=False, order_offset=0):
    completed = subprocess.run(
        [str(executable.resolve()), "--input", str(path), "--workers", str(workers),
         "--repeats", str(repeats), "--diagnostics", str(diagnostics),
         "--warmups", str(warmups)] + (["--include-pls", "1"] if include_pls else []) +
        (["--order-offset", str(order_offset)] if order_offset else []),
        text=True, capture_output=True, check=True)
    lines = [line for line in completed.stdout.splitlines() if line.startswith("{")]
    if len(lines) != 1:
        raise ValueError("Expected one native JSON result")
    return json.loads(lines[0])


def main():
    args = parse_args()
    algorithms = PLS_ALGORITHMS if args.include_pls else ALGORITHMS
    orders = PLS_ORDERS if args.include_pls else ORDERS
    if not args.benchmark.is_file():
        raise FileNotFoundError(args.benchmark)
    args.input_dir.mkdir(parents=True, exist_ok=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "schema": "resident-gradient-study-v3" if args.include_pls else "resident-gradient-study-v2",
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(),
        "binary_format": command("file", str(args.benchmark.resolve())),
        "binary_sha256": hashlib.sha256(args.benchmark.read_bytes()).hexdigest(),
        "source_revision_at_run": command("git", "rev-parse", "HEAD"),
        "source_dirty_at_run": bool(command("git", "status", "--porcelain")),
        "headers_sha256_at_run": header_digest(ROOT / "include"),
        "driver_sha256_at_run": hashlib.sha256(
            (ROOT / "benchmarks/ttk/benchmark_resident_gradient.cpp").read_bytes()).hexdigest(),
        "runner_sha256_at_run": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "requested_ttk_revision": TTK_REVISION,
        "omp_wait_policy": os.environ.get("OMP_WAIT_POLICY", "runtime-default"),
        "arguments": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        "timing_scope": "Common resident vertex values and maximal-cell arrays to a ready native gradient; all native preparation included",
        "primary_comparison": "algorithm_seconds, excluding separately reported native construction; includes MorseFrames builder setup and TTK vertex ordering/lower stars",
        "construction_phases": {a: sorted(CONSTRUCTION_PHASES[a]) for a in algorithms},
        "algorithm_phases": {a: sorted(ALGORITHM_PHASES[a]) for a in algorithms},
        "performance_clock_note": "Outer phase-boundary clocks only; no internal profiling; components partition every total",
        "loading_note": "input_loading_seconds is one shared file read/parse/validation per native invocation, not a per-algorithm cost or cold-cache disk benchmark",
        "excluded": ["input generation and file I/O", "reference validation", "post-readiness teardown", "persistence"],
        "diagnostic_note": "Separate runs; nested gradient details are not additive to outer phases; medians need not add; TTK lower stars and matching are combined",
        "algorithm_order_cycle": [[algorithms[i] for i in order] for order in orders],
        "interval_note": "v3 intervals resample paired repetition ratios within one session; "
                         "not independent-session/machine uncertainty or a causal comparison to older runs",
        "cases": [],
    }
    inputs = []
    # Finish all input generation before any performance measurements.
    for family, sizes in (("terrain", args.terrain_sizes), ("volume", args.volume_sizes)):
        for size in sizes:
            for seed in args.seeds:
                path = args.input_dir / f"{family}-n{size}-seed{seed}.txt"
                if not path.exists():
                    generator = (generators.make_injective_terrain if family == "terrain"
                                 else generators.make_injective_volume)
                    write_ttk_input(generator(seed, size), path)
                inputs.append((family, size, seed, path))
    if args.reverse_inputs:
        inputs.reverse()
    for index, (family, size, seed, path) in enumerate(inputs):
        order = list(args.workers)
        if not args.include_pls:
            shift = (index // 2) % len(order)
            order = order[shift:] + order[:shift]
            if index % 2:
                order.reverse()
        case = dict(family=family, size=size, seed=seed,
                    input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                    worker_order=order, measurements=[])
        for workers in order:
            print(f"Resident {family} n={size} seed={seed} workers={workers}", flush=True)
            raw = run_native(args.benchmark, path, workers, args.repeats,
                             args.diagnostics, args.warmups, args.include_pls, args.order_offset)
            if raw.get("schema") != ("resident-gradient-v3" if args.include_pls else "resident-gradient-v2"):
                raise ValueError("Rebuild the native driver to collect performance phase splits")
            if args.include_pls and raw.get("performance_orders") != [
                    [algorithms[a] for a in orders[(i + args.order_offset) % len(orders)]]
                    for i in range(args.repeats)]:
                raise ValueError("Native execution order differs from requested schedule")
            expected_dim = 2 if family == "terrain" else 3
            if (raw["workers"] != workers or raw["dimension"] != expected_dim or
                    raw["vertices"] != size ** expected_dim):
                raise ValueError("Native metadata differs from requested case")
            case["measurements"].append(dict(
                workers=workers, raw=raw,
                summary=summarize(raw, args.repeats, args.diagnostics)))
        baseline = next(m for m in case["measurements"] if m["workers"] == 1)
        for measurement in case["measurements"]:
            if (measurement["raw"]["simplices"] != baseline["raw"]["simplices"] or
                    measurement["raw"]["vertices"] != baseline["raw"]["vertices"]):
                raise ValueError("Topology differs across worker configurations")
            for algorithm in algorithms:
                if (measurement["raw"]["algorithms"][algorithm]["critical_counts"] !=
                        baseline["raw"]["algorithms"][algorithm]["critical_counts"]):
                    raise ValueError("Critical counts changed across worker counts")
                measurement["summary"][algorithm]["speedup_vs_one"] = (
                    baseline["summary"][algorithm]["total_seconds"]["median"] /
                    measurement["summary"][algorithm]["total_seconds"]["median"])
                measurement["summary"][algorithm]["algorithm_speedup_vs_one"] = (
                    baseline["summary"][algorithm]["algorithm_seconds"]["median"] /
                    measurement["summary"][algorithm]["algorithm_seconds"]["median"])
        data["cases"].append(case)
        args.output.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    data["completed_utc"] = datetime.now(timezone.utc).isoformat()
    args.output.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    print(f"Saved {len(data['cases'])} inputs to {args.output}", flush=True)


if __name__ == "__main__":
    main()
