#!/usr/bin/env python3
"""Run the smoke checks for the Morse persistence discussion package."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def run(
    command: list[str],
    *,
    cwd: Path = ROOT,
    env: dict[str, str] | None = None,
    capture_to: Path | None = None,
) -> None:
    print("$ " + " ".join(command), flush=True)
    if capture_to is None:
        subprocess.run(command, cwd=cwd, env=env, check=True)
        return

    completed = subprocess.run(
        command,
        cwd=cwd,
        env=env,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    capture_to.parent.mkdir(parents=True, exist_ok=True)
    capture_to.write_text(completed.stdout)
    print(f"wrote {rel(capture_to)}", flush=True)


def python_env() -> dict[str, str]:
    env = os.environ.copy()
    python_path = str(ROOT / "python")
    existing = env.get("PYTHONPATH")
    env["PYTHONPATH"] = python_path if not existing else f"{python_path}{os.pathsep}{existing}"
    return env


def build_dir_from_arg(value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return ROOT / path


def configure_if_needed(build_dir: Path, force: bool) -> None:
    cache = build_dir / "CMakeCache.txt"
    if force or not cache.exists():
        run(["cmake", "-S", str(ROOT), "-B", str(build_dir)])


def run_native_gudhi_benchmark(args: argparse.Namespace, build_dir: Path) -> None:
    executable = build_dir / "morse_benchmark_gudhi_view"
    if not executable.exists():
        raise FileNotFoundError(
            f"Native GUDHI benchmark executable not found: {executable}"
        )

    mode_flag = f"--{args.native_mode}"
    command = [str(executable), mode_flag, "--repeats", str(args.native_repeats)]
    for strategy in args.native_strategy:
        command.extend(["--strategy", strategy])

    output_prefix = f"reproduction_native_gudhi_{args.native_mode}"
    if args.native_strategy:
        strategy_suffix = "_".join(strategy.replace("-", "_") for strategy in args.native_strategy)
        output_prefix += f"_{strategy_suffix}"

    csv_path = DOCS / f"{output_prefix}.csv"
    run(command, capture_to=csv_path)

    table_path = DOCS / f"{output_prefix}_table.tex"
    run(
        [
            sys.executable,
            str(ROOT / "tools" / "render_native_gudhi_view_table.py"),
            "--input",
            str(csv_path),
            "--output",
            str(table_path),
            "--summary",
        ]
    )

    stage_table_path = DOCS / f"{output_prefix}_stage_table.tex"
    stage_prose_path = DOCS / f"{output_prefix}_stage_prose.tex"
    run(
        [
            sys.executable,
            str(ROOT / "tools" / "render_native_gudhi_stage_profile.py"),
            "--input",
            str(csv_path),
            "--table-output",
            str(stage_table_path),
            "--prose-output",
            str(stage_prose_path),
            "--summary",
        ]
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--build-dir",
        default="build",
        help="CMake build directory, relative to the project root unless absolute",
    )
    parser.add_argument(
        "--configure",
        action="store_true",
        help="rerun CMake configuration even if CMakeCache.txt exists",
    )
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--skip-tests", action="store_true")
    parser.add_argument("--skip-python-tutorial", action="store_true")
    parser.add_argument(
        "--with-native-gudhi-benchmark",
        action="store_true",
        help="rerun the native GUDHI benchmark into docs/reproduction_native_gudhi_*",
    )
    parser.add_argument(
        "--native-mode",
        choices=("quick", "large"),
        default="quick",
        help="native GUDHI benchmark mode",
    )
    parser.add_argument(
        "--native-repeats",
        type=int,
        default=3,
        help="number of benchmark repeats when --with-native-gudhi-benchmark is used",
    )
    parser.add_argument(
        "--native-strategy",
        action="append",
        default=[],
        help="optional strategy filter; can be passed more than once",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    build_dir = build_dir_from_arg(args.build_dir)

    configure_if_needed(build_dir, args.configure)

    if not args.skip_build:
        parallel = str(os.cpu_count() or 4)
        run(["cmake", "--build", str(build_dir), "--parallel", parallel])

    if not args.skip_tests:
        run(["ctest", "--test-dir", str(build_dir), "--output-on-failure"])

    if not args.skip_python_tutorial:
        run(
            [
                sys.executable,
                str(ROOT / "python" / "examples" / "prime_field_tutorial.py"),
                "--modulus",
                "3",
                "--algorithm",
                "f-max",
            ],
            env=python_env(),
        )

    if args.with_native_gudhi_benchmark:
        run_native_gudhi_benchmark(args, build_dir)

    viewer_pdf = DOCS / "experiments_morse_persistence_viewer.pdf"
    main_pdf = DOCS / "experiments_morse_persistence.pdf"
    if viewer_pdf.exists() and not main_pdf.exists():
        shutil.copyfile(viewer_pdf, main_pdf)

    print("\nDiscussion package checks completed.")
    print(f"Read first: {rel(DOCS / 'discussion_package.md')}")
    print(f"Viewer PDF: {rel(viewer_pdf)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
