# Python prime-field tutorial

This tutorial shows how to use the Python API over a prime field `F_p`.
The current implementation supports prime moduli such as `2`, `3`, `5`, and
`7`. Composite rings such as `Z_4` or `Z_6` are rejected because the barcode
reducers assume field coefficients.

The same public objects are used as in the `Z2` path:

- `FilteredComplex` stores the filtered simplicial complex.
- `MorseSequence` stores critical simplices and regular pairs.
- `compute_reference_map_modp` computes the Morse reference map over `F_p`.
- `compute_coreference_map_modp` computes the dual coreference map over `F_p`.
- `compute_morse_persistence_modp` reduces the Morse complex from references.
- `compute_morse_coreference_persistence_modp` reduces it from coreferences.
- `compute_standard_persistence(..., modulus=p)` is the full-complex oracle.

## Run the example

From the workspace root:

```sh
python3 morseframes/python/examples/prime_field_tutorial.py --modulus 3
```

To try another Morse sequence strategy:

```sh
python3 morseframes/python/examples/prime_field_tutorial.py \
  --modulus 5 \
  --algorithm f-max
```

The example constructs a small lower-star filtered complex with two triangles
and a genuine plateau:

```python
import morseframes as mp

complex_ = mp.FilteredComplex.from_lower_star(
    [(0, 1, 2), (0, 2, 3)],
    {0: 1.0, 1: 0.0, 2: 1.0, 3: 0.0},
)
```

A Morse sequence is independent of the coefficient field:

```python
sequence = mp.compute_morse_sequence(
    complex_,
    algorithm=mp.SAME_LEVEL_REDUCTION_SEQUENCE,
)
```

The coefficient field enters when the reference and coreference maps are built.
For odd primes, annotations carry signed coefficients modulo `p`:

```python
references = mp.compute_reference_map_modp(
    complex_,
    sequence,
    modulus=3,
)

coreferences = mp.compute_coreference_map_modp(
    complex_,
    sequence,
    modulus=3,
)
```

Those maps can be passed directly to the Morse reducers:

```python
morse = mp.compute_morse_persistence_modp(
    complex_,
    sequence=sequence,
    references=references,
    modulus=3,
)

dual_morse = mp.compute_morse_coreference_persistence_modp(
    complex_,
    sequence=sequence,
    coreferences=coreferences,
    modulus=3,
)
```

The standard reducer uses the same `modulus` keyword:

```python
standard = mp.compute_standard_persistence(complex_, modulus=3)

assert morse.finite_barcode() == standard.finite_barcode()
assert morse.essential_barcode() == standard.essential_barcode()
assert dual_morse.finite_barcode() == standard.finite_barcode()
assert dual_morse.essential_barcode() == standard.essential_barcode()
```

By default, `finite_barcode()` omits zero-length intervals. This is the right
comparison for persistence; zero-length plateau cancellations can differ between
the full complex and the Morse-reduced complex. Use
`finite_barcode(include_zero=True)` only when inspecting implementation-level
tie-breaking.

## Working with annotations

The `Z2` reference map is a tuple of critical ids. The `F_p` map is a tuple of
`(critical_id, coefficient)` pairs:

```python
for simplex_id, annotation in enumerate(references):
    if annotation:
        simplex = complex_.vertices(simplex_id)
        print(simplex, annotation)
```

To display the critical ids as simplices, index through
`sequence.critical_simplices`:

```python
def annotation_as_simplices(annotation):
    return tuple(
        (complex_.vertices(sequence.critical_simplices[critical_id]), coefficient)
        for critical_id, coefficient in annotation
    )
```

Coefficients are stored in the range `1..p-1`; for example, the coefficient `2`
in `F_3` represents `-1`.

## Convenience entry points

If you do not need to inspect the maps, the generic persistence functions are
usually enough:

```python
morse = mp.compute_morse_persistence(
    complex_,
    algorithm="f-max",
    modulus=3,
)

dual_morse = mp.compute_morse_coreference_persistence(
    complex_,
    algorithm="f-max",
    modulus=3,
)
```

The explicit `_modp` names are useful when you want code to make the coefficient
field visible:

```python
standard = mp.compute_standard_persistence_modp(complex_, 3)
```

## Current scope

This is prime-field persistence, not general `Z_n` persistence. In particular:

- `modulus=2` gives the usual `Z2` behavior.
- Odd prime moduli use oriented boundary and coboundary coefficients.
- Composite moduli raise `ValueError`.
- The Python prototype exposes this functionality now; the current GUDHI-facing
  patch is intentionally kept at `Z2` for a smaller first contribution.
