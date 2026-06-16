#!/usr/bin/env python3
"""Compare built-in Morse examples with GUDHI when GUDHI is installed."""

from __future__ import annotations

import argparse
import math
import subprocess
from pathlib import Path


EXAMPLES = {
    "one_vertex": [((0,), 0.0)],
    "edge_later": [((0,), 0.0), ((1,), 0.0), ((0, 1), 1.0)],
    "edge_same": [((0,), 0.0), ((1,), 0.0), ((0, 1), 0.0)],
    "triangle_boundary": [
        ((0,), 0.0),
        ((1,), 0.0),
        ((2,), 0.0),
        ((0, 1), 1.0),
        ((0, 2), 1.0),
        ((1, 2), 1.0),
    ],
    "filled_triangle": [
        ((0,), 0.0),
        ((1,), 0.0),
        ((2,), 0.0),
        ((0, 1), 1.0),
        ((0, 2), 1.0),
        ((1, 2), 1.0),
        ((0, 1, 2), 2.0),
    ],
}


def morse_barcode(executable: Path, example: str) -> list[tuple[str, int, float, float]]:
    output = subprocess.check_output([str(executable), example], text=True)
    rows = []
    for line in output.splitlines():
        kind, dim, birth, death = line.split("\t")
        rows.append((kind, int(dim), float(birth), math.inf if death == "inf" else float(death)))
    return sorted(rows)


def gudhi_barcode(example: str) -> list[tuple[str, int, float, float]]:
    try:
        import gudhi
    except ModuleNotFoundError as exc:
        raise SystemExit("GUDHI is not installed. Install the Python package `gudhi` to use this tool.") from exc

    st = gudhi.SimplexTree()
    for simplex, filtration in EXAMPLES[example]:
        st.insert(simplex, filtration=filtration)
    st.make_filtration_non_decreasing()
    st.persistence(persistence_dim_max=True)

    rows = []
    max_dim = max(len(simplex) - 1 for simplex, _ in EXAMPLES[example])
    for dim in range(max_dim + 1):
        for birth, death in st.persistence_intervals_in_dimension(dim):
            if math.isinf(death):
                rows.append(("essential", dim, float(birth), math.inf))
            elif birth < death:
                rows.append(("finite", dim, float(birth), float(death)))
    return sorted(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--exe",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "build" / "morse_example_barcode",
    )
    args = parser.parse_args()

    for example in EXAMPLES:
        morse = morse_barcode(args.exe, example)
        gudhi = gudhi_barcode(example)
        if morse != gudhi:
            print(f"{example}: mismatch")
            print(f"  morse: {morse}")
            print(f"  gudhi: {gudhi}")
            return 1
        print(f"{example}: ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
