#!/usr/bin/env python3
"""Alternate isolated native RK builds, retaining raw times and exact checks.

Use git revisions for immutable comparisons, or WORKTREE for current headers.
No test/build process is launched while measurements are running. Each worker
keeps its topology resident and blocks on stdin while the other is measured.
"""
from __future__ import annotations

import argparse
from contextlib import ExitStack
import filecmp
import gc
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import shutil
import statistics
import subprocess
import tempfile
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]


def command_output(*args: str) -> str:
    return subprocess.check_output(args, cwd=ROOT, text=True).strip()


def snapshot_headers(ref: str, destination: Path) -> str:
    destination.mkdir()
    if ref == "WORKTREE":
        shutil.copytree(ROOT / "include", destination / "include")
        return "WORKTREE at " + command_output("git", "rev-parse", "HEAD")
    revision = command_output("git", "rev-parse", "--verify", ref + "^{commit}")
    archive = subprocess.run(
        ["git", "archive", revision, "include"], cwd=ROOT,
        stdout=subprocess.PIPE, check=True,
    )
    subprocess.run(["tar", "-x", "-C", str(destination)],
                   input=archive.stdout, check=True)
    return revision


def header_digest(path: Path) -> str:
    digest = hashlib.sha256()
    for file in sorted(path.rglob("*.hpp")):
        digest.update(str(file.relative_to(path)).encode())
        digest.update(b"\0")
        digest.update(file.read_bytes())
    return digest.hexdigest()


class Worker:
    def __init__(self, binary: Path, input_path: Path, dump: Path,
                 input_mode: str | None = None):
        self.process = subprocess.Popen(
            [str(binary), str(input_path), str(dump)] +
            ([] if input_mode is None else [input_mode]), text=True,
            stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        )
        try:
            self.metadata = self.read()
        except BaseException:
            self.close()
            raise

    def read(self):
        line = self.process.stdout.readline()
        if not line:
            raise RuntimeError("Native benchmark worker exited unexpectedly")
        return json.loads(line)

    def run(self, workers: int, repeats: int) -> list[float]:
        self.process.stdin.write(f"run {workers} {repeats}\n")
        self.process.stdin.flush()
        return self.read()

    def close(self):
        try:
            if self.process.poll() is None:
                self.process.stdin.write("quit\n")
                self.process.stdin.flush()
                self.process.wait(timeout=10)
        except (BrokenPipeError, subprocess.TimeoutExpired):
            self.process.terminate()
            self.process.wait(timeout=10)
        finally:
            self.process.stdin.close()
            self.process.stdout.close()


def quartiles(values: list[float]) -> list[float]:
    return statistics.quantiles(values, n=4, method="inclusive")[::2]


def summarize(samples: list[dict]) -> dict:
    by_version = {
        name: [t for block in samples for t in block[name]]
        for name in ("baseline", "candidate")
    }
    ratios = [statistics.median(b["candidate"]) / statistics.median(b["baseline"])
              for b in samples]
    rng = random.Random(0)
    boot = sorted(statistics.median(rng.choices(ratios, k=len(ratios)))
                  for _ in range(2000))
    return {
        "median_seconds": {k: statistics.median(v) for k, v in by_version.items()},
        "iqr_seconds": {k: quartiles(v) for k, v in by_version.items()},
        "median_paired_ratio": statistics.median(ratios),
        "paired_ratio_iqr": quartiles(ratios),
        "paired_ratio_bootstrap_95_interval": [boot[49], boot[1949]],
        "faster_blocks": sum(r < 1 for r in ratios),
        "blocks": len(ratios),
    }


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--candidate", default="WORKTREE")
    parser.add_argument("--sizes", type=int, nargs="+", default=[16, 24, 32])
    parser.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    parser.add_argument("--workers", type=int, nargs="+", default=[1, 2, 4, 8])
    parser.add_argument("--blocks", type=int, default=12)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--compiler", default="clang++")
    parser.add_argument("--cxx-flags", default="-std=c++17 -O3 -DNDEBUG -pthread")
    parser.add_argument("--input-dir", type=Path,
                        help="Reuse/save deterministic input complexes here")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if (min(args.sizes) < 2 or min(args.workers) < 1 or args.blocks < 4 or
            args.blocks % 2 or args.repeats < 1 or args.warmups < 1):
        parser.error("sizes >= 2, workers/repeats/warmups >= 1; blocks even and >= 4")
    if args.output.exists():
        parser.error("output exists; use a new path to preserve previous evidence")
    return args


def main():
    import shlex
    import benchmark_simplicial_strategies as generators
    from benchmark_ttk_process_lower_stars import write_ttk_input

    args = parse_args()
    compiler = shutil.which(args.compiler)
    if compiler is None:
        raise FileNotFoundError(args.compiler)
    flags = shlex.split(args.cxx_flags)
    source = ROOT / "benchmarks" / "benchmark_reduction_kernel_ab.cpp"
    result = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(), "python_architecture": platform.machine(),
        "compiler": command_output(compiler, "--version"), "flags": flags,
        "cpu_count": os.cpu_count(),
        "benchmark_source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "timing_scope": "Fresh builder + gradient, including pool setup and replay; "
                        "excludes topology, validation, protocol and output destruction",
        "blocks": args.blocks, "repeats_per_block": args.repeats,
        "warmups_per_worker_count": args.warmups,
        "interval_note": "Paired-block bootstrap; describes this session only, "
                         "not independent-machine or independent-session uncertainty",
        "builds": {}, "cases": [],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="morseframes-ab-") as scratch:
        scratch = Path(scratch)
        binaries = {}
        for name, ref in [("baseline", args.baseline), ("candidate", args.candidate)]:
            directory = scratch / name
            revision = snapshot_headers(ref, directory)
            binary = directory / "worker"
            print(f"Building {name}: {revision}", flush=True)
            subprocess.run([compiler, *flags, "-I", str(directory / "include"),
                            str(source), "-o", str(binary)], check=True)
            result["builds"][name] = {
                "revision": revision, "headers_sha256": header_digest(directory / "include"),
                "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
                "binary_format": command_output("file", str(binary)),
            }
            binaries[name] = binary
        inputs = args.input_dir or scratch / "inputs"
        inputs.mkdir(parents=True, exist_ok=True)
        for size in args.sizes:
            for seed in args.seeds:
                input_path = inputs / f"volume-n{size}-seed{seed}.txt"
                if not input_path.exists():
                    print(f"Generating volume n={size}, seed={seed}", flush=True)
                    complex_ = generators.make_injective_volume(seed, size)
                    write_ttk_input(complex_, input_path)
                    del complex_
                    gc.collect()
                input_digest = hashlib.sha256(input_path.read_bytes()).hexdigest()
                with ExitStack() as stack:
                    workers = {}
                    for name in binaries:
                        workers[name] = Worker(binaries[name], input_path,
                                               scratch / (name + ".sequence"))
                        stack.callback(workers[name].close)
                    if workers["baseline"].metadata != workers["candidate"].metadata:
                        raise AssertionError("Complex/critical counts differ between builds")
                    if not filecmp.cmp(scratch / "baseline.sequence",
                                       scratch / "candidate.sequence", shallow=False):
                        raise AssertionError("Exact reference sequences differ between builds")
                    for count in args.workers:
                        for worker in workers.values():
                            worker.run(count, args.warmups)
                        samples = []
                        for block in range(args.blocks):
                            order = ["baseline", "candidate"]
                            if (block + seed + args.workers.index(count)) % 2:
                                order.reverse()
                            sample = {"order": order}
                            for name in order:
                                sample[name] = workers[name].run(count, args.repeats)
                            samples.append(sample)
                        case = {
                            "grid_size": size, "seed": seed, "workers": count,
                            "input_sha256": input_digest,
                            **workers["baseline"].metadata,
                            "exact_sequences_match": True,
                            "samples": samples, "summary": summarize(samples),
                        }
                        result["cases"].append(case)
                        args.output.write_text(json.dumps(result, indent=2) + "\n")
                        summary = case["summary"]
                        print(f"n={size} seed={seed} workers={count}: candidate/base "
                              f"{summary['median_paired_ratio']:.3f}, "
                              f"IQR {summary['paired_ratio_iqr']}, "
                              f"exact sequence OK", flush=True)
    result["finished_utc"] = datetime.now(timezone.utc).isoformat()
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
