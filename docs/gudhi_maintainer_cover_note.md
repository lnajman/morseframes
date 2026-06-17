# GUDHI Maintainer Cover Note

This is a short cover note for opening a GUDHI maintainer discussion or draft
pull request. It deliberately describes the contribution in GUDHI terms and
does not rely on an external Morse persistence manuscript.

## Proposed Title

Add a header-only Morse persistence module for `Simplex_tree`

## Summary

This patch proposes a first C++ `Morse_persistence` module for GUDHI. The
module computes ordinary persistence over `Z/2Z` directly from a filtered
`Gudhi::Simplex_tree<>`.

The intended use case is a filtered simplicial complex with plateaus. Equal
filtration values are treated as genuine levels: the module constructs
same-level Morse pairings inside each level, builds the associated reference
map, and reduces the resulting annotations to obtain the persistence diagram.
This avoids forcing a lower-star or simplex-wise refinement before the Morse
sequence is constructed.

The module is meant as a complementary reduction path, not a replacement for
GUDHI's existing persistent cohomology implementation.

## What Is Included

The first patch is intentionally narrow:

- a header-only module under `src/Morse_persistence/include/gudhi`;
- one public umbrella header, `<gudhi/Morse_persistence.h>`;
- direct input from `Gudhi::Simplex_tree<>`;
- ordinary persistence over `Z2`;
- fused Morse sequence and reference-map construction;
- local simplex ids that can be mapped back to `Simplex_tree::Simplex_handle`;
- four public sequence strategies:
  `SAME_LEVEL_REDUCTION`, `F_MAX`, `F_MIN`, and `PLATEAU_GREEDY`;
- one small example;
- one C++ test file covering plateau and non-plateau cases;
- one Doxygen intro page.

The main entry point is:

```cpp
auto result = Gudhi::morse_persistence::compute_morse_persistence(
    simplex_tree,
    Gudhi::morse_persistence::Morse_sequence_strategy::F_MAX);
```

The result exposes the Morse sequence, finite intervals, essential intervals,
basic metrics, and a handle map back to the input `Simplex_tree`.

## What Is Not Included

The first patch does not include:

- Python bindings;
- benchmark scripts;
- prime-field or composite-ring coefficient APIs;
- cubical complexes;
- experimental flooding strategies as stable API;
- a new public complex data structure.

Those are better discussed after the first C++ `Simplex_tree` API is reviewed.

## Validation

The staging patch was checked against a shallow clean checkout of:

```text
https://github.com/GUDHI/gudhi-devel.git
HEAD 3a7d79c
```

The following checks passed:

```sh
cmake -S <gudhi-checkout> -B <build> \
  -DWITH_GUDHI_PYTHON=OFF \
  -DWITH_GUDHI_GUDHUI=OFF \
  -DWITH_GUDHI_EXAMPLE=ON \
  -DWITH_GUDHI_TEST=ON
cmake --build <build> --target Morse_persistence_example_from_simplex_tree
cmake --build <build> --target Morse_persistence_test_simplex_tree
ctest --test-dir <build> -R Morse_persistence --output-on-failure
cmake --build <build> --target doxygen
```

The C++ test matrix checks the public API on:

- a single vertex;
- a single increasing edge;
- a filled triangle entirely on one plateau;
- a filled triangle with a later tail edge;
- two components joined late;
- a one-cycle killed by a later triangle.

For each case, it checks the public strategies listed above, compares the
finite and essential intervals with expected `Z2` persistence, and verifies
that simplex ids can be mapped back to non-null `Simplex_tree` handles.

## Questions For Maintainers

The most useful first feedback would be:

- Is `Gudhi::morse_persistence` an acceptable namespace and module boundary?
- Is a header-only module acceptable for this first version?
- Should the result expose local simplex ids plus a handle map, or should
  intervals expose `Simplex_tree::Simplex_handle` directly?
- Should all four strategies be public initially, or should the first API expose
  only one strategy plus a simpler baseline?
- Should the test be converted to GUDHI's Boost-test helper style, or is the
  current standalone executable plus `add_test` acceptable?
- What extra Doxygen/API documentation would be expected before a non-draft pull
  request?

## Suggested Review Framing

The first review should focus on API shape, integration style, and correctness
checks. Performance claims can be discussed later once maintainers agree that
the module boundary and result types are suitable for GUDHI.
