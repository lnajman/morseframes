#!/usr/bin/env python3
"""Controlled construction A/B from resident arrays; separate clocks and memory.

The same driver is compiled against both header snapshots. No builds or input
generation overlap measurements. Diagnostic clocks are never performance data.
"""
from __future__ import annotations

import argparse
from contextlib import ExitStack
from datetime import datetime, timezone
import filecmp
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import tempfile

from benchmark_reduction_kernel_ab import (
    ROOT, Worker, command_output, header_digest, snapshot_headers, summarize,
)

METRICS = ("construction_seconds", "enumeration_and_insertion_seconds",
           "finalize_seconds", "fmax_seconds", "rk_seconds")


def summarize_samples(samples):
    if not samples or len(samples) % 2:
        raise ValueError("An even number of paired blocks is required")
    for block in samples:
        if sorted(block["order"]) != ["baseline", "candidate"]:
            raise ValueError("Invalid version order")
        if not block["baseline"] or len(block["baseline"]) != len(block["candidate"]):
            raise ValueError("Incomplete paired block")
        for name in ("baseline", "candidate"):
            for sample in block[name]:
                if set(sample) != set(METRICS) or any(
                    not math.isfinite(x) or x <= 0 for x in sample.values()
                ):
                    raise ValueError("Invalid performance sample")
                if sample["construction_seconds"] <= 0 or not math.isclose(
                    sample["construction_seconds"],
                    sample["enumeration_and_insertion_seconds"] + sample["finalize_seconds"],
                    rel_tol=1e-10, abs_tol=1e-12,
                ):
                    raise ValueError("Construction phases do not partition the total")
    if sum(b["order"][0] == "baseline" for b in samples) * 2 != len(samples):
        raise ValueError("Version order must be balanced")
    return {metric: summarize([
        {name: [sample[metric] for sample in block[name]]
         for name in ("baseline", "candidate")} for block in samples
    ]) for metric in METRICS}


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--baseline", required=True)
    p.add_argument("--candidate", default="WORKTREE")
    p.add_argument("--inputs", type=Path, nargs="+", required=True)
    p.add_argument("--workers", type=int, nargs="+", default=[1, 8])
    p.add_argument("--blocks", type=int, default=8)
    p.add_argument("--repeats", type=int, default=2)
    p.add_argument("--warmups", type=int, default=2)
    p.add_argument("--diagnostics", type=int, default=3)
    p.add_argument("--compiler", default="clang++")
    p.add_argument("--cxx-flags", default="-std=c++17 -O3 -DNDEBUG -pthread")
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args(argv)
    if args.blocks < 4 or args.blocks % 2 or min(
        args.repeats, args.warmups, args.diagnostics, *args.workers
    ) < 1:
        p.error("blocks even and >= 4; other counts positive")
    if len(set(args.inputs)) != len(args.inputs) or len(set(args.workers)) != len(args.workers):
        p.error("duplicate inputs/workers")
    if args.output.exists():
        p.error("output exists; preserve previous evidence")
    if not all(path.is_file() for path in args.inputs):
        p.error("all inputs must already exist")
    return args


def main():
    args = parse_args()
    compiler = shutil.which(args.compiler)
    if compiler is None:
        raise FileNotFoundError(args.compiler)
    source = ROOT / "benchmarks/benchmark_complex_construction.cpp"
    result = dict(
        schema="complex-construction-ab-v1",
        started_utc=datetime.now(timezone.utc).isoformat(),
        platform=platform.platform(), architecture=platform.machine(),
        compiler=command_output(compiler, "--version"), flags=shlex.split(args.cxx_flags),
        cpu_count=os.cpu_count(), source_revision=command_output("git", "rev-parse", "HEAD"),
        source_status=command_output("git", "status", "--porcelain"),
        driver_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        worktree_header_patch=command_output("git", "diff", "HEAD", "--", "include"),
        arguments={k: ([str(x) for x in v] if k == "inputs" else str(v) if isinstance(v, Path) else v)
                   for k, v in vars(args).items()},
        timing_scope="Resident arrays to finalized complex; gradient timers separately include fresh builders. "
                     "Loading, validation, protocol, output destruction excluded. F-Max/RK order alternates. "
                     "Diagnostic insertion clocks perturb enumeration; their residual includes clock overhead.",
        memory_scope="Fresh-process peak RSS after construction, before gradients/validation; "
                     "includes runtime, resident input, allocator and constructor temporaries. Not live heap size.",
        interval_note="Paired-block bootstrap describes this session, not independent-session uncertainty.",
        builds={}, cases=[],
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="morseframes-construction-") as directory:
        scratch = Path(directory)
        binaries = {}
        for name, ref in (("baseline", args.baseline), ("candidate", args.candidate)):
            snapshot = scratch / name
            revision = snapshot_headers(ref, snapshot)
            binary = snapshot / "worker"
            print(f"Building {name}: {revision}", flush=True)
            subprocess.run([compiler, *result["flags"], "-I", str(snapshot / "include"),
                            str(source), "-o", str(binary)], check=True)
            result["builds"][name] = dict(
                revision=revision, headers_sha256=header_digest(snapshot / "include"),
                binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                binary_format=command_output("file", str(binary)),
            )
            binaries[name] = binary
        for index, input_path in enumerate(args.inputs):
            case = dict(input_path=str(input_path.resolve()),
                        input_sha256=hashlib.sha256(input_path.read_bytes()).hexdigest(),
                        memory={}, diagnostics={}, measurements=[])
            # Isolate memory samples from both reference construction and gradients.
            for name, binary in binaries.items():
                case["memory"][name] = [json.loads(subprocess.check_output(
                    [str(binary), str(input_path), "--memory"], text=True))
                    for _ in range(args.diagnostics)]
            with ExitStack() as stack:
                workers = {}
                for name, binary in binaries.items():
                    workers[name] = Worker(binary, input_path, scratch / (name + ".dump"))
                    stack.callback(workers[name].close)
                if workers["baseline"].metadata != workers["candidate"].metadata or not filecmp.cmp(
                    scratch / "baseline.dump", scratch / "candidate.dump", shallow=False
                ):
                    raise AssertionError("Exact complex/gradient comparison failed")
                case.update(workers["baseline"].metadata)
                for entries in case["memory"].values():
                    for entry in entries:
                        if entry["simplices"] != case["simplices"] or not (
                            0 <= entry["input_peak_bytes"] <= entry["construction_peak_bytes"]
                        ):
                            raise AssertionError("Invalid isolated memory sample")
                case["exact_complex_and_reference_gradients_match"] = True
                for count in (args.workers if index % 2 == 0 else list(reversed(args.workers))):
                    for worker in workers.values():
                        worker.run(count, args.warmups)
                    samples = []
                    for block in range(args.blocks):
                        order = ["baseline", "candidate"]
                        if (block + index + args.workers.index(count)) % 2:
                            order.reverse()
                        sample = {"order": order}
                        for name in order:
                            sample[name] = workers[name].run(count, args.repeats)
                        samples.append(sample)
                    summary = summarize_samples(samples)
                    case["measurements"].append(dict(workers=count, samples=samples, summary=summary))
                    print(f"{input_path.stem} workers={count}: construction ratio "
                          f"{summary['construction_seconds']['median_paired_ratio']:.3f}", flush=True)
                for name, worker in workers.items():
                    case["diagnostics"][name] = []
                    for _ in range(args.diagnostics):
                        worker.process.stdin.write("profile\n")
                        worker.process.stdin.flush()
                        diagnostic = worker.read()
                        if diagnostic["complex_fingerprint"] != case["complex_fingerprint"]:
                            raise AssertionError("Diagnostic complex differs")
                        if any(not math.isfinite(v) or v < 0 for k, v in diagnostic.items()
                               if k.endswith("seconds")):
                            raise AssertionError("Invalid diagnostic duration")
                        if diagnostic["simplices"] != case["simplices"] or diagnostic[
                            "insertion_attempts"
                        ] != case["insertion_attempts"]:
                            raise AssertionError("Diagnostic counters differ")
                        case["diagnostics"][name].append(diagnostic)
            result["cases"].append(case)
            args.output.write_text(json.dumps(result, indent=2) + "\n")
    result["completed_utc"] = datetime.now(timezone.utc).isoformat()
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
