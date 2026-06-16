#!/usr/bin/env python3
"""Tutorial: Morse persistence over a prime field F_p."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import morseframes as mp  # noqa: E402


def build_plateau_complex() -> mp.FilteredComplex:
    """Return a tiny lower-star complex with a genuine same-level plateau."""

    return mp.FilteredComplex.from_lower_star(
        [(0, 1, 2), (0, 2, 3)],
        {0: 1.0, 1: 0.0, 2: 1.0, 3: 0.0},
    )


def field_annotation_as_simplices(
    complex_: mp.FilteredComplex,
    sequence: mp.MorseSequence,
    annotation: mp.FieldAnnotation,
) -> tuple[tuple[mp.Simplex, int], ...]:
    """Replace critical ids in an F_p annotation by the critical simplices."""

    return tuple(
        (complex_.vertices(sequence.critical_simplices[critical_id]), coefficient)
        for critical_id, coefficient in annotation
    )


def print_nonzero_annotations(
    title: str,
    complex_: mp.FilteredComplex,
    sequence: mp.MorseSequence,
    annotations: tuple[mp.FieldAnnotation, ...],
    *,
    limit: int = 8,
) -> None:
    print(f"\n{title}")
    shown = 0
    for simplex_id, annotation in enumerate(annotations):
        if not annotation:
            continue
        simplex = complex_.vertices(simplex_id)
        as_simplices = field_annotation_as_simplices(complex_, sequence, annotation)
        print(f"  {simplex}: {as_simplices}")
        shown += 1
        if shown == limit:
            remaining = sum(1 for item in annotations if item) - shown
            if remaining > 0:
                print(f"  ... {remaining} more nonzero annotations")
            break


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--modulus",
        type=int,
        default=3,
        help="prime field modulus p; composite moduli are rejected",
    )
    parser.add_argument(
        "--algorithm",
        default=mp.SAME_LEVEL_REDUCTION_SEQUENCE,
        choices=mp.MORSE_SEQUENCE_ALGORITHMS,
        help="Morse sequence strategy",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    complex_ = build_plateau_complex()
    sequence = mp.compute_morse_sequence(complex_, algorithm=args.algorithm)

    references = mp.compute_reference_map_modp(
        complex_,
        sequence,
        modulus=args.modulus,
    )
    coreferences = mp.compute_coreference_map_modp(
        complex_,
        sequence,
        modulus=args.modulus,
    )

    morse = mp.compute_morse_persistence_modp(
        complex_,
        sequence=sequence,
        references=references,
        modulus=args.modulus,
    )
    dual_morse = mp.compute_morse_coreference_persistence_modp(
        complex_,
        sequence=sequence,
        coreferences=coreferences,
        modulus=args.modulus,
    )
    standard = mp.compute_standard_persistence(complex_, modulus=args.modulus)

    assert morse.finite_barcode() == standard.finite_barcode()
    assert morse.essential_barcode() == standard.essential_barcode()
    assert dual_morse.finite_barcode() == standard.finite_barcode()
    assert dual_morse.essential_barcode() == standard.essential_barcode()

    print(f"Prime field: F_{args.modulus}")
    print(f"Algorithm: {sequence.algorithm}")
    print(f"Simplices: {len(complex_)}")
    print(f"Critical simplices: {sequence.critical_simplices_as_simplices(complex_)}")

    print("\nFiltration")
    for simplex, value in complex_.filtration_list():
        print(f"  {simplex}: {value}")

    print_nonzero_annotations("Reference map over F_p", complex_, sequence, references)
    print_nonzero_annotations("Coreference map over F_p", complex_, sequence, coreferences)

    print("\nOff-diagonal barcodes")
    print(f"  Morse finite:       {morse.finite_barcode()}")
    print(f"  Morse essential:    {morse.essential_barcode()}")
    print(f"  Standard finite:    {standard.finite_barcode()}")
    print(f"  Standard essential: {standard.essential_barcode()}")
    print("\nMorse and standard persistence agree over the selected prime field.")


if __name__ == "__main__":
    main()
