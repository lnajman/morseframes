#!/usr/bin/env python3
"""Bounded diagnostic-only PLS pass on existing inputs and a frozen executable."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import statistics
import tempfile

from benchmark_reduction_kernel_ab import Worker, ROOT, command_output, header_digest
from pls_phase_profile import validate


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--inputs", type=Path, nargs="+", required=True)
    p.add_argument("--workers", type=int, nargs="+", default=[1, 8])
    p.add_argument("--repeats", type=int, default=6)
    p.add_argument("--warmups", type=int, default=2)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    if min(a.workers) < 1 or a.repeats < 1 or a.warmups < 1 or a.output.exists():
        p.error("positive counts and a fresh output path are required")
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    data = {"schema": "pls-fine-profile-v1", "started_utc": datetime.now(timezone.utc).isoformat(),
            "source_revision": command_output("git", "rev-parse", "HEAD"),
            "source_status": command_output("git", "status", "--porcelain"),
            "headers_sha256": header_digest(ROOT / "include"), "binary_sha256": digest(a.binary),
            "driver_sha256": digest(ROOT / "benchmarks/benchmark_simplicial_gradients.cpp"),
            "profile_helper_sha256": digest(ROOT / "benchmarks/pls_profile.hpp"),
            "settings": {k: str(v) if isinstance(v, Path) else v for k,v in vars(a).items() if k != "inputs"},
            "scope": "Diagnostic only; setup/local/cleanup children overlap their parents, not one another.",
            "cases": []}
    with a.output.open("x") as output, tempfile.TemporaryDirectory(prefix="pls-profile-") as directory:
        for path in a.inputs:
            w = Worker(a.binary.resolve(), path.resolve(), Path(directory)/"reference.dump")
            try:
                for count in a.workers:
                    rows = []
                    for i in range(a.warmups + a.repeats):
                        w.process.stdin.write(f"profile {count}\n"); w.process.stdin.flush()
                        row = w.read()
                        row["validated_fine_seconds"] = validate(row["pls_profile_seconds"], row["setup_seconds"],
                            row["local_wall_seconds"], row["replay_seconds"], row["algorithm_seconds"] - row["builder_seconds"])
                        if i >= a.warmups:
                            rows.append(row)
                    case = {"input": str(path.resolve()), "input_sha256": digest(path), "workers": count,
                            "metadata": w.metadata, "profiles": rows}
                    data["cases"].append(case)
                    print(path.name, count, {k: round(1000*statistics.median(r["validated_fine_seconds"][k] for r in rows),3)
                                            for k in rows[0]["validated_fine_seconds"]}, flush=True)
            finally:
                w.close()
        data["completed_utc"] = datetime.now(timezone.utc).isoformat()
        json.dump(data, output, indent=2); output.write("\n")


if __name__ == "__main__":
    main()
