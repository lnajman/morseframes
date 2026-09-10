#!/usr/bin/env python3
"""Prepared-complex F-Max/PLS/RK comparison and exact old/new PLS checks.

Construction is separate. Fresh builders, lower-star partitioning, priorities,
worker pools, local work and replay remain inside gradient timing. No persistence.
"""
from __future__ import annotations

import argparse
from contextlib import ExitStack
from datetime import datetime, timezone
import filecmp
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import platform
import random
import shlex
import shutil
import subprocess
import tempfile
from pls_phase_profile import validate as validate_pls_profile

from benchmark_reduction_kernel_ab import (
    ROOT, Worker, command_output, header_digest, snapshot_headers, summarize,
)

ALGORITHMS = ("f_max", "process_lower_stars", "reduction_kernel")
PHASES = ("builder_seconds", "kernel_seconds", "algorithm_seconds")
VERSIONS = ("baseline", "candidate")


def grid_input(dimension: int, size: int, seed: int) -> str:
    """Freudenthal triangulation, with shuffled injective vertex ranks.

    Limit face-generation attempts before enumerating factorially many cells.
    These small synthetic grids test dimension, not a representative data corpus.
    """
    if not 1 <= dimension <= 7 or size < 2:
        raise ValueError("grid dimension must be 1..7 and side length >= 2")
    cells_count = (size - 1) ** dimension * math.factorial(dimension)
    if cells_count * (2 ** (dimension + 1) - 1) > 3_000_000:
        raise ValueError("grid exceeds the bounded face-generation budget")
    values = list(range(size ** dimension))
    random.Random(seed).shuffle(values)
    strides = [size ** i for i in range(dimension)]
    lines = [f"morseframes-ttk-v1 {dimension} {len(values)} {cells_count}",
             " ".join(map(str, values))]
    for coordinates in itertools.product(range(size - 1), repeat=dimension):
        base = sum(x * stride for x, stride in zip(coordinates, strides))
        for permutation in itertools.permutations(range(dimension)):
            vertices = [base]
            for axis in permutation:
                vertices.append(vertices[-1] + strides[axis])
            lines.append(" ".join(map(str, vertices)))
    return "\n".join(lines) + "\n"


def check_runs(runs, repeats):
    if len(runs) != repeats:
        raise AssertionError("Incomplete repetitions")
    for run in runs:
        if set(run) != set(ALGORITHMS):
            raise AssertionError("Missing algorithms")
        for times in run.values():
            if set(times) != set(PHASES) or any(
                not math.isfinite(t) or t <= 0 for t in times.values()
            ):
                raise AssertionError("Invalid timings")
            if not math.isclose(times["builder_seconds"] + times["kernel_seconds"],
                                times["algorithm_seconds"], rel_tol=1e-10):
                raise AssertionError("Phases do not sum to algorithm time")


def summaries(samples):
    result = {}
    for algorithm in ALGORITHMS:
        result[algorithm] = {}
        for phase in PHASES:
            blocks = [{v: [r[algorithm][phase] for r in b[v]] for v in VERSIONS}
                      for b in samples]
            result[algorithm][phase] = summarize(blocks)
    return result


def comparisons(samples):
    result = {}
    for numerator, denominator in [("process_lower_stars", "f_max"),
                                   ("reduction_kernel", "process_lower_stars"),
                                   ("reduction_kernel", "f_max")]:
        blocks = [{"baseline": [r[denominator]["algorithm_seconds"] for r in b["candidate"]],
                   "candidate": [r[numerator]["algorithm_seconds"] for r in b["candidate"]]}
                  for b in samples]
        result[f"{numerator}/{denominator}"] = summarize(blocks)
    return result


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--baseline", required=True)
    p.add_argument("--candidate", default="WORKTREE")
    p.add_argument("--inputs", type=Path, nargs="*", default=[])
    p.add_argument("--grids", nargs="*", default=[], help="DIMENSION:SIDE, e.g. 4:4 5:3 6:2 7:2")
    p.add_argument("--seeds", type=int, nargs="+", default=[0, 2])
    p.add_argument("--input-dir", type=Path)
    p.add_argument("--reverse-inputs", action="store_true")
    p.add_argument("--workers", type=int, nargs="+", default=[1, 8])
    p.add_argument("--blocks", type=int, default=8)
    p.add_argument("--repeats", type=int, default=3)
    p.add_argument("--warmups", type=int, default=2)
    p.add_argument("--profiles", type=int, default=3)
    p.add_argument("--memory-repeats", type=int, default=3)
    p.add_argument("--compiler", default="clang++")
    p.add_argument("--cxx-flags", default="-std=c++17 -O3 -DNDEBUG -pthread")
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args(argv)
    if (a.blocks < 4 or a.blocks % 2 or a.repeats < 1 or a.blocks * a.repeats % 6
            or min(a.workers) < 1 or a.warmups < 1 or a.profiles < 1 or a.memory_repeats < 0):
        p.error("positive counts, even blocks >= 4, and total repetitions divisible by 6 required")
    if len(set(a.workers)) != len(a.workers) or len(set(a.seeds)) != len(a.seeds):
        p.error("duplicate worker counts or seeds")
    if not a.inputs and not a.grids:
        p.error("supply --inputs or --grids")
    if a.output.exists():
        p.error("output exists; preserve previous evidence")
    return a


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    args = parse_args()
    source = ROOT / "benchmarks/benchmark_simplicial_gradients.cpp"
    compiler = shutil.which(args.compiler)
    if compiler is None:
        raise FileNotFoundError(args.compiler)
    flags = shlex.split(args.cxx_flags)
    result = {
        "schema": "simplicial-gradient-ab-v1", "completed": False,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(), "architecture": platform.machine(),
        "cpu_count": os.cpu_count(), "compiler": command_output(compiler, "--version"),
        "flags": flags, "source_status": command_output("git", "status", "--porcelain"),
        "headers_patch": command_output("git", "diff", "HEAD", "--", "include"),
        "driver_sha256": digest(source), "runner_sha256": digest(Path(__file__)),
        "helper_sha256": digest(ROOT / "tools/benchmark_reduction_kernel_ab.py"),
        "profile_helper_sha256": digest(ROOT / "benchmarks/pls_profile.hpp"),
        "settings": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()
                     if k != "inputs"},
        "timing_scope": "Prepared common native complex; fresh builder plus full gradient, "
                        "including function-specific setup, pools and replay. Loading and native "
                        "construction reported separately. No persistence or output teardown.",
        "memory_scope": "Fresh-process peak RSS after construction / gradient; includes input, "
                        "runtime and constructor high-water mark, not isolated live heap.",
        "interval_scope": "Paired-block bootstrap within this session, not across machines/sessions.",
        "builds": {}, "inputs": [], "cases": [],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="morseframes-simplicial-") as directory:
        scratch = Path(directory)
        inputs = [p.resolve(strict=True) for p in args.inputs]
        input_dir = (args.input_dir or scratch / "inputs").resolve()
        input_dir.mkdir(parents=True, exist_ok=True)
        for specification in args.grids:
            dimension, side = map(int, specification.split(":"))
            for seed in args.seeds:
                path = input_dir / f"grid-d{dimension}-n{side}-seed{seed}.txt"
                contents = grid_input(dimension, side, seed)
                if path.exists():
                    if path.read_text() != contents:
                        raise ValueError(f"Conflicting existing input: {path}")
                else:
                    path.write_text(contents)
                inputs.append(path)
        if len(set(inputs)) != len(inputs):
            raise ValueError("Duplicate input paths")
        if args.reverse_inputs:
            inputs.reverse()
        binaries = {}
        for version, ref in zip(VERSIONS, (args.baseline, args.candidate)):
            build = scratch / version
            revision = snapshot_headers(ref, build)
            binary = build / "worker"
            print(f"Building {version}: {revision}", flush=True)
            subprocess.run([compiler, *flags, "-I", str(build / "include"), str(source),
                            "-o", str(binary)], check=True)
            binaries[version] = binary
            result["builds"][version] = {"revision": revision,
                "headers_sha256": header_digest(build / "include"),
                "binary_sha256": digest(binary), "format": command_output("file", str(binary))}
        # No compilation, input generation or tests during the measurement loop.
        for case_index, path in enumerate(inputs):
            with ExitStack() as stack:
                workers = {}
                for v in VERSIONS:
                    workers[v] = Worker(binaries[v], path, scratch / f"{v}.dump")
                    stack.callback(workers[v].close)
                metadata = {v: w.metadata for v, w in workers.items()}
                if metadata["baseline"]["identity"] != metadata["candidate"]["identity"]:
                    raise AssertionError("Complex or sequential gradients changed")
                if not filecmp.cmp(scratch / "baseline.dump", scratch / "candidate.dump", shallow=False):
                    raise AssertionError("Complete ordered complex/gradient dumps differ")
                identity = metadata["candidate"]["identity"]
                for algorithm in identity["algorithms"].values():
                    counts = algorithm["critical_counts"]
                    if sum((-1) ** d * n for d, n in enumerate(counts)) != identity["euler"]:
                        raise AssertionError("Critical-cell Euler characteristic mismatch")
                result["inputs"].append({"path": str(path), "sha256": digest(path),
                    "metadata": metadata, "exact_dumps_match": True})
                for worker_index, count in enumerate(args.workers):
                    for w in workers.values():
                        check_runs(w.run(count, args.warmups), args.warmups)
                    samples = []
                    for block in range(args.blocks):
                        order = list(VERSIONS)
                        if (block + case_index + worker_index) % 2:
                            order.reverse()
                        sample = {"order": order}
                        for v in order:
                            sample[v] = workers[v].run(count, args.repeats)
                            check_runs(sample[v], args.repeats)
                        samples.append(sample)
                    profiles = {v: [] for v in VERSIONS}
                    for repeat in range(args.profiles):
                        for v in VERSIONS[::(-1 if repeat % 2 else 1)]:
                            w = workers[v]
                            w.process.stdin.write(f"profile {count}\n"); w.process.stdin.flush()
                            profile = w.read()
                            if profile["stars"] != identity["vertices"]:
                                raise AssertionError("Incorrect profile star count")
                            if "pls_profile_seconds" in profile:
                                profile["validated_fine_seconds"] = validate_pls_profile(
                                    profile["pls_profile_seconds"], profile["setup_seconds"],
                                    profile["local_wall_seconds"], profile["replay_seconds"],
                                    profile["algorithm_seconds"] - profile["builder_seconds"])
                            profiles[v].append(profile)
                    memory = {v: {a: [] for a in ALGORITHMS} for v in VERSIONS}
                    for repeat in range(args.memory_repeats):
                        for v in VERSIONS[::(-1 if repeat % 2 else 1)]:
                            for a, algorithm in enumerate(ALGORITHMS):
                                memory[v][algorithm].append(json.loads(subprocess.check_output(
                                    [str(binaries[v]), str(path), "--memory", str(a), str(count)], text=True)))
                    case = {"input_index": case_index, "workers": count, "samples": samples,
                            "profiles": profiles, "memory": memory,
                            "summary": summaries(samples), "comparison": comparisons(samples)}
                    result["cases"].append(case)
                    args.output.write_text(json.dumps(result, indent=2) + "\n")
                    ratio = case["summary"]["process_lower_stars"]["algorithm_seconds"]["median_paired_ratio"]
                    print(f"{path.name}, workers={count}: PLS new/old {ratio:.3f}", flush=True)
        if len(result["cases"]) != len(inputs) * len(args.workers):
            raise AssertionError("Incomplete case coverage")
        result["completed"] = True
        result["finished_utc"] = datetime.now(timezone.utc).isoformat()
        args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
