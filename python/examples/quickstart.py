"""Minimal MorseFrames Python quickstart."""

from __future__ import annotations

import morseframes as mf


def main() -> None:
    complex_ = mf.FilteredComplex.from_simplices(
        [
            ([0], 0.0),
            ([1], 0.0),
            ([2], 0.0),
            ([0, 1], 1.0),
            ([1, 2], 1.0),
            ([0, 2], 2.0),
            ([0, 1, 2], 3.0),
        ]
    )

    sequence = mf.compute_morse_sequence(complex_, algorithm="f-max")
    references = mf.compute_reference_map(complex_, sequence)
    diagram = mf.compute_morse_persistence(complex_, sequence, references)
    diagram_f3 = mf.compute_morse_persistence(complex_, algorithm="f-max", modulus=3)

    print(f"MorseFrames {mf.__version__}")
    print(f"simplices: {complex_.size}")
    print(f"critical simplices: {len(sequence.critical_simplices)}")
    print(f"F2 finite barcode: {diagram.finite_barcode()}")
    print(f"F3 finite barcode: {diagram_f3.finite_barcode()}")


if __name__ == "__main__":
    main()

