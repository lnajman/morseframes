#!/usr/bin/env python3
"""Compare fresh RK, F-Max and TTK gradients from common resident mesh arrays.

Native preparation is inside timing. Input generation, I/O, reference checks,
and destruction after gradient readiness are outside. Phase diagnostics use
separate runs; their samples never enter the headline performance ratios.
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
import statistics
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from benchmark_reduction_kernel_ab import header_digest
from benchmark_ttk_process_lower_stars import TTK_REVISION, write_ttk_input
import benchmark_simplicial_strategies as generators

ALGORITHMS = ("f_max", "reduction_kernel", "ttk")
OUTER_PHASES = {
    "f_max": {"representation_and_filtration", "builder_setup", "gradient"},
    "reduction_kernel": {"representation_and_filtration", "builder_setup", "gradient"},
    "ttk": {"native_object_init", "vertex_order", "representation_setup",
            "connectivity_precondition", "gradient"},
}
DETAIL_PHASES = {
    "f_max": {"workspace_init", "candidate_seeding", "candidate_selection",
              "emission_and_updates", "callbacks"},
    "reduction_kernel": {"workspace_and_pool", "level_processing", "replay"},
    "ttk": set(),
}
ORDERS = ((0, 1, 2), (2, 1, 0), (1, 2, 0), (0, 2, 1), (2, 0, 1), (1, 0, 2))


def distribution(values, *, positive=False):
    if not values or any(not math.isfinite(v) or v < 0 or (positive and v == 0)
                         for v in values):
        raise ValueError("Expected finite nonnegative samples (positive for elapsed totals)")
    q = (statistics.quantiles(values, n=4, method="inclusive")
         if len(values) > 1 else [values[0]] * 3)
    return dict(median=statistics.median(values), q1=q[0], q3=q[2],
                min=min(values), max=max(values))


def summarize(result, repeats, diagnostics):
    if result.get("schema") != "resident-gradient-v1":
        raise ValueError("Wrong timing schema: rebuild the resident-gradient executable")
    if result.get("ttk_revision") != TTK_REVISION:
        raise ValueError("Native benchmark was not built against the pinned TTK revision")
    if result.get("exact_reference_checks") is not True:
        raise ValueError("Missing successful exact reference checks")
    if result.get("dimension") not in (1, 2, 3) or result.get("simplices", 0) <= 0:
        raise ValueError("Invalid complex metadata")
    if type(result.get("critical_counts_match")) is not bool:
        raise ValueError("Missing critical-count agreement flag")
    if set(result.get("algorithms", {})) != set(ALGORITHMS):
        raise ValueError("Missing algorithm measurements")
    summary = {}
    for algorithm in ALGORITHMS:
        raw = result["algorithms"][algorithm]
        counts = raw.get("critical_counts", [])
        if (len(counts) != result["dimension"] + 1 or
                any(type(c) is not int or c < 0 for c in counts) or
                sum(counts) > result["simplices"] or
                (result["simplices"] - sum(counts)) % 2):
            raise ValueError("Invalid critical-simplex counts")
        performance = raw["performance_seconds"]
        rows = raw["diagnostics"]
        if len(performance) != repeats or len(rows) != diagnostics:
            raise ValueError("Incomplete timing repetitions")
        entry = {"total_seconds": distribution(performance, positive=True)}
        totals, phases, details, shares = [], {}, {}, {}
        for row in rows:
            total = row["total_seconds"]
            distribution([total], positive=True)
            outer = row["phases_seconds"]
            inner = row["gradient_details_seconds"]
            if set(outer) != OUTER_PHASES[algorithm] or set(inner) != DETAIL_PHASES[algorithm]:
                raise ValueError("Missing or unexpected phase measurements")
            distribution(list(outer.values()))
            if not math.isclose(sum(outer.values()), total, rel_tol=1e-9, abs_tol=1e-12):
                raise ValueError("Outer phases do not account for the diagnostic total")
            if outer["gradient"] <= 0:
                raise ValueError("Gradient duration must be positive")
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
    counts = [result["algorithms"][a]["critical_counts"] for a in ALGORITHMS]
    if result.get("critical_counts_match") != (counts[0] == counts[1] == counts[2]):
        raise ValueError("Critical-count agreement flag is inconsistent")
    # Different algorithms may legitimately create different critical counts.
    # Exact consistency with each algorithm's own reference is mandatory.
    summary["paired_ratios"] = {}
    for numerator, denominator in (("reduction_kernel", "ttk"),
                                   ("reduction_kernel", "f_max"), ("ttk", "f_max")):
        a = result["algorithms"][numerator]["performance_seconds"]
        b = result["algorithms"][denominator]["performance_seconds"]
        summary["paired_ratios"][f"{numerator}/{denominator}"] = distribution(
            [x / y for x, y in zip(a, b, strict=True)], positive=True)
    return summary


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", type=Path, required=True)
    parser.add_argument("--terrain-sizes", type=int, nargs="*", default=[16, 64])
    parser.add_argument("--volume-sizes", type=int, nargs="*", default=[8, 16])
    parser.add_argument("--seeds", type=int, nargs="+", default=[0, 2])
    parser.add_argument("--workers", type=int, nargs="+", default=[1, 2, 4, 8])
    parser.add_argument("--repeats", type=int, default=12)
    parser.add_argument("--diagnostics", type=int, default=6)
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
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
    if (args.repeats < 6 or args.repeats % 6 or args.diagnostics < 6 or
            args.diagnostics % 6 or args.warmups < 0):
        parser.error("Repetitions and diagnostics must be positive multiples of six")
    if args.output.exists():
        parser.error("Output exists; choose a fresh path to preserve evidence")
    return args


def command(*args):
    return subprocess.check_output(args, cwd=ROOT, text=True).strip()


def run_native(executable, path, workers, repeats, diagnostics, warmups):
    completed = subprocess.run(
        [str(executable.resolve()), "--input", str(path), "--workers", str(workers),
         "--repeats", str(repeats), "--diagnostics", str(diagnostics),
         "--warmups", str(warmups)], text=True, capture_output=True, check=True)
    lines = [line for line in completed.stdout.splitlines() if line.startswith("{")]
    if len(lines) != 1:
        raise ValueError("Expected one native JSON result")
    return json.loads(lines[0])


def main():
    args = parse_args()
    if not args.benchmark.is_file():
        raise FileNotFoundError(args.benchmark)
    args.input_dir.mkdir(parents=True, exist_ok=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "schema": "resident-gradient-study-v1",
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(),
        "binary_format": command("file", str(args.benchmark.resolve())),
        "binary_sha256": hashlib.sha256(args.benchmark.read_bytes()).hexdigest(),
        "source_revision_at_run": command("git", "rev-parse", "HEAD"),
        "source_dirty_at_run": bool(command("git", "status", "--porcelain")),
        "headers_sha256_at_run": header_digest(ROOT / "include"),
        "driver_sha256_at_run": hashlib.sha256(
            (ROOT / "benchmarks/ttk/benchmark_resident_gradient.cpp").read_bytes()).hexdigest(),
        "requested_ttk_revision": TTK_REVISION,
        "omp_wait_policy": os.environ.get("OMP_WAIT_POLICY", "runtime-default"),
        "arguments": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        "timing_scope": "Common resident vertex values and maximal-cell arrays to a ready native gradient; all native preparation included",
        "excluded": ["input generation and file I/O", "reference validation", "post-readiness teardown", "persistence"],
        "diagnostic_note": "Separate runs; nested gradient details are not additive to outer phases; medians need not add; TTK lower stars and matching are combined",
        "algorithm_order_cycle": [[ALGORITHMS[i] for i in order] for order in ORDERS],
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
    for index, (family, size, seed, path) in enumerate(inputs):
        order = list(args.workers)
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
                             args.diagnostics, args.warmups)
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
            for algorithm in ALGORITHMS:
                if (measurement["raw"]["algorithms"][algorithm]["critical_counts"] !=
                        baseline["raw"]["algorithms"][algorithm]["critical_counts"]):
                    raise ValueError("Critical counts changed across worker counts")
                measurement["summary"][algorithm]["speedup_vs_one"] = (
                    baseline["summary"][algorithm]["total_seconds"]["median"] /
                    measurement["summary"][algorithm]["total_seconds"]["median"])
        data["cases"].append(case)
        args.output.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    data["completed_utc"] = datetime.now(timezone.utc).isoformat()
    args.output.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    print(f"Saved {len(data['cases'])} inputs to {args.output}", flush=True)


if __name__ == "__main__":
    main()
