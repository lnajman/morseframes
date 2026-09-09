#!/usr/bin/env python3
"""Repeated native RK timings and separate coarse/detailed diagnostic runs.

The resident native worker checks every result against its validated sequential
gradient. Performance samples never enable instrumentation. Coarse runs retain
metrics-free facet kernels; detailed runs explain work but can perturb it.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shlex
import shutil
import statistics
import subprocess
import tempfile
from datetime import datetime, timezone

from benchmark_reduction_kernel_ab import (
    ROOT, Worker, command_output, header_digest, snapshot_headers,
)


def distribution(values):
    if not values or any(not math.isfinite(v) or v < 0 for v in values):
        raise ValueError("Expected nonempty, finite, nonnegative measurements")
    quartiles = (statistics.quantiles(values, n=4, method="inclusive")
                 if len(values) > 1 else [values[0]] * 3)
    return {"median": statistics.median(values), "q1": quartiles[0],
            "q3": quartiles[2], "min": min(values), "max": max(values)}


def summarize(samples):
    performance = [v for block in samples["run"] for v in block]
    if not performance or min(performance) <= 0:
        raise ValueError("Performance timings must be positive")
    summary = {"performance_seconds": distribution(performance)}
    for mode in ("coarse", "detailed"):
        rows = [dict(row) for block in samples[mode] for row in block]
        if not rows or any(set(row) != set(rows[0]) for row in rows):
            raise ValueError("Missing or inconsistent diagnostic samples")
        summary[mode] = {key: distribution([row[key] for row in rows])
                         for key in rows[0]}
        # These are disjoint outer wall-clock phases from the SAME run.
        for row in rows:
            if row["total_seconds"] <= 0:
                raise ValueError("Diagnostic total must be positive")
            row["accounted_seconds"] = sum(row[k] for k in (
                "builder_seconds", "setup_seconds", "level_wall_seconds", "replay_seconds"))
            if row["accounted_seconds"] > row["total_seconds"] * 1.001:
                raise ValueError("Outer phases exceed their measured total")
        summary[mode]["shares"] = {
            key: distribution([row[key] / row["total_seconds"] for row in rows])
            for key in ("builder_seconds", "setup_seconds", "level_wall_seconds", "replay_seconds")
        }
        summary[mode]["other_share"] = distribution([
            max(0, 1 - row["accounted_seconds"] / row["total_seconds"]) for row in rows])
        summary[mode]["overhead_ratio"] = (
            summary[mode]["total_seconds"]["median"] / summary["performance_seconds"]["median"])
        activity = [row["cumulative_level_task_seconds"] / row["level_wall_seconds"]
                    for row in rows if row["max_parallel_levels"] > 1 and row["level_wall_seconds"] > 0]
        summary[mode]["effective_level_task_parallelism"] = distribution(activity) if activity else None
    return summary


def worker_order(workers, block):
    order = list(workers)
    shift = (block // 2) % len(order)
    order = order[shift:] + order[:shift]
    return order[::-1] if block % 2 else order


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--terrain-sizes", type=int, nargs="*", default=[16, 64])
    parser.add_argument("--volume-sizes", type=int, nargs="*", default=[16, 32])
    parser.add_argument("--plateau-terrain-sizes", type=int, nargs="*", default=[16, 32])
    parser.add_argument("--plateau-volume-sizes", type=int, nargs="*", default=[4, 8, 12])
    parser.add_argument("--seeds", type=int, nargs="+", default=[0, 2])
    parser.add_argument("--workers", type=int, nargs="+", default=[1, 2, 4, 8])
    parser.add_argument("--blocks", type=int, default=8)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--coarse-blocks", type=int, default=8)
    parser.add_argument("--detailed-blocks", type=int, default=3)
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--compiler", default="clang++")
    parser.add_argument("--cxx-flags", default="-std=c++17 -O3 -DNDEBUG -pthread")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if (min(args.workers) < 1 or 1 not in args.workers or len(set(args.workers)) != len(args.workers)
            or min(args.blocks, args.repeats, args.coarse_blocks, args.detailed_blocks, args.warmups) < 1):
        parser.error("Require unique positive workers including 1, and positive repeat counts")
    if args.output.exists():
        parser.error("Output exists; use a new path to preserve evidence")
    import benchmark_simplicial_strategies as generators
    from benchmark_ttk_process_lower_stars import write_ttk_input
    cases = []
    args.input_dir.mkdir(parents=True, exist_ok=True)
    for family, sizes, mode, seeds in (
        ("terrain", args.terrain_sizes, "lower-star", args.seeds),
        ("volume", args.volume_sizes, "lower-star", args.seeds),
        ("terrain", args.plateau_terrain_sizes, "plateau", [args.seeds[0]]),
        ("volume", args.plateau_volume_sizes, "plateau", [args.seeds[0]]),
    ):
        for size in sizes:
            for seed in seeds:
                path = args.input_dir / f"{family}-n{size}-seed{seed}.txt"
                if not path.exists():
                    print(f"Generating {path.name}", flush=True)
                    generator = getattr(generators, f"make_injective_{family}")
                    complex_ = generator(seed, size)
                    write_ttk_input(complex_, path)
                    del complex_
                    gc.collect()
                cases.append({"family": family, "size": size, "seed": seed,
                              "mode": mode, "input": str(path.resolve()),
                              "input_sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    if not cases:
        parser.error("Select at least one case")
    source = ROOT / "benchmarks" / "benchmark_reduction_kernel_ab.cpp"
    compiler = shutil.which(args.compiler)
    if compiler is None:
        raise FileNotFoundError(args.compiler)
    flags = shlex.split(args.cxx_flags) + ["-DMORSEFRAMES_RK_PHASE_PROFILE"]
    result = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(), "cpu_count": os.cpu_count(),
        "compiler": command_output(compiler, "--version"), "flags": flags,
        "driver_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "orchestrator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "arguments": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        "timing_scope": "Fresh builder, workspace, pool, gradient and replay; excludes topology, "
                        "validation, protocol I/O and returned sequence destruction. No persistence.",
        "diagnostic_note": "Coarse and detailed samples are separate runs, not performance timings. "
                           "Detailed local times accumulate over workers and nest inside facet execution; "
                           "never add them to outer wall-clock phases. Other includes destruction/pool teardown.",
        "cases": [],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="morseframes-phases-") as scratch_name:
        scratch = Path(scratch_name)
        snapshot = scratch / "snapshot"
        result["revision"] = snapshot_headers("WORKTREE", snapshot)
        result["headers_sha256"] = header_digest(snapshot / "include")
        binary = scratch / "worker"
        print("Building native phase worker", flush=True)
        subprocess.run([compiler, *flags, "-I", str(snapshot / "include"),
                        str(source), "-o", str(binary)], check=True)
        result["binary_format"] = command_output("file", str(binary))
        result["binary_sha256"] = hashlib.sha256(binary.read_bytes()).hexdigest()
        # All generation and compilation have finished before measurements start.
        for case in cases:
            print(f"Measuring {case['mode']} {case['family']} n={case['size']} seed={case['seed']}", flush=True)
            worker = Worker(binary, Path(case["input"]), scratch / "reference.sequence", case["mode"])
            try:
                case.update(worker.metadata)
                case["sequence_sha256"] = hashlib.sha256((scratch / "reference.sequence").read_bytes()).hexdigest()
                samples = {count: {mode: [] for mode in ("run", "coarse", "detailed")}
                           for count in args.workers}
                orders = {}
                for mode, blocks, repeats in (("run", args.blocks, args.repeats),
                        ("coarse", args.coarse_blocks, 1), ("detailed", args.detailed_blocks, 1)):
                    for count in args.workers:
                        worker.process.stdin.write(f"{mode} {count} {args.warmups}\n")
                        worker.process.stdin.flush()
                        worker.read()
                    orders[mode] = []
                    for block in range(blocks):
                        order = worker_order(args.workers, block)
                        orders[mode].append(order)
                        for count in order:
                            worker.process.stdin.write(f"{mode} {count} {repeats}\n")
                            worker.process.stdin.flush()
                            sample = worker.read()
                            if len(sample) != repeats:
                                raise ValueError("Incomplete sample block")
                            samples[count][mode].append(sample)
                case["orders"] = orders
                case["measurements"] = [
                    {"workers": count, "summary": summarize(samples[count]), "samples": samples[count]}
                    for count in args.workers]
                baseline = next(m for m in case["measurements"] if m["workers"] == 1)
                for measurement in case["measurements"]:
                    measurement["summary"]["speedup_vs_one"] = (
                        baseline["summary"]["performance_seconds"]["median"] /
                        measurement["summary"]["performance_seconds"]["median"])
                case["all_sequences_match"] = True
                result["cases"].append(case)
                args.output.write_text(json.dumps(result, indent=2) + "\n")
            finally:
                worker.close()
    result["completed_utc"] = datetime.now(timezone.utc).isoformat()
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"Saved {len(result['cases'])} cases to {args.output}", flush=True)


if __name__ == "__main__":
    main()
