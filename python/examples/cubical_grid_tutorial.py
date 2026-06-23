#!/usr/bin/env python3
"""Tutorial: Morse persistence on a 2D cubical grid."""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from time import perf_counter

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import morseframes as mf  # noqa: E402


def vertex_field(width: int, height: int) -> list[float]:
    """Return row-major vertex values with plateaus and multiple extrema."""

    if width < 2 or height < 2:
        raise ValueError("width and height must be at least 2.")

    values: list[float] = []
    for y in range(height):
        yy = (2.0 * math.pi * y) / float(height - 1)
        for x in range(width):
            xx = (2.0 * math.pi * x) / float(width - 1)
            value = math.sin(xx) * math.sin(yy)
            value += 0.12 * math.cos(2.0 * xx + yy)
            values.append(round(value, 6))
    return values


def time_best(action, *, repeats: int):
    best_seconds = float("inf")
    best_result = None
    for _ in range(repeats):
        start = perf_counter()
        result = action()
        seconds = perf_counter() - start
        if seconds < best_seconds:
            best_seconds = seconds
            best_result = result
    if best_result is None:
        raise RuntimeError("No repetitions were run.")
    return best_seconds, best_result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--width", type=int, default=8, help="number of grid vertices in x")
    parser.add_argument("--height", type=int, default=8, help="number of grid vertices in y")
    parser.add_argument(
        "--algorithm",
        default=mf.F_MAX_SEQUENCE,
        choices=mf.MORSE_SEQUENCE_ALGORITHMS,
        help="Morse sequence strategy",
    )
    parser.add_argument(
        "--modulus",
        type=int,
        default=3,
        help="prime field modulus; use 2 for ordinary Z2",
    )
    parser.add_argument("--repeats", type=int, default=3, help="best-of timing repetitions")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.repeats < 1:
        raise ValueError("repeats must be positive.")
    if not mf.cpp_backend_available() or mf.CppCubicalGrid2DComplex is None:
        raise RuntimeError("The native C++ backend is required for cubical grids.")

    values = vertex_field(args.width, args.height)
    grid = mf.CubicalGrid2DComplex.from_vertex_values(args.width, args.height, values)

    sequence_seconds, sequence = time_best(
        lambda: mf.compute_morse_sequence(grid, algorithm=args.algorithm),
        repeats=args.repeats,
    )
    morse_seconds, morse = time_best(
        lambda: mf.compute_morse_persistence(
            grid,
            sequence,
            algorithm=args.algorithm,
            modulus=args.modulus,
        ),
        repeats=args.repeats,
    )
    coreference_seconds, coreference = time_best(
        lambda: mf.compute_morse_coreference_persistence(
            grid,
            sequence,
            algorithm=args.algorithm,
            modulus=args.modulus,
        ),
        repeats=args.repeats,
    )
    standard_seconds, standard = time_best(
        lambda: mf.compute_standard_persistence(grid, modulus=args.modulus),
        repeats=args.repeats,
    )

    assert morse.finite_barcode() == standard.finite_barcode()
    assert morse.essential_barcode() == standard.essential_barcode()
    assert coreference.finite_barcode() == standard.finite_barcode()
    assert coreference.essential_barcode() == standard.essential_barcode()

    square = grid.square(0, 0)
    print(f"MorseFrames {mf.__version__}")
    print(f"Grid: {args.width} x {args.height} vertices")
    print(f"Cells: {grid.size}")
    print(f"Levels: {grid.num_levels}")
    print(f"Algorithm: {sequence.algorithm}")
    print(f"Field: F_{args.modulus}")
    print(f"Critical cells: {len(sequence.critical_simplices)}")
    print(f"First square boundary: {grid.boundary(square)}")
    print(
        "First square F_p signs: "
        f"{tuple(grid.boundary_coefficient(square, i, args.modulus) for i in range(4))}"
    )
    print("\nTimings, best of repeats")
    print(f"  sequence:    {sequence_seconds:.6f}s")
    print(f"  Morse ref:   {morse_seconds:.6f}s")
    print(f"  Morse coref: {coreference_seconds:.6f}s")
    print(f"  standard:    {standard_seconds:.6f}s")
    if morse_seconds:
        print(f"  standard / Morse ref: {standard_seconds / morse_seconds:.3f}x")
    print("\nBarcodes")
    print(f"  finite:    {morse.finite_barcode()}")
    print(f"  essential: {morse.essential_barcode()}")
    print("\nMorse reference, coreference, and ordinary persistence agree.")


if __name__ == "__main__":
    main()
