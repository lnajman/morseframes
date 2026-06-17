# GUDHI Upstream Skeleton Branch

This branch is a staging branch for shaping the MorseFrames GUDHI adapter into
a small upstream-style C++ contribution.

It is not intended to be opened as a GUDHI pull request as-is. The purpose is to
keep the upstream-facing surface narrow while preserving the broader prototype
on `main`.

## First Patch Scope

The first GUDHI-facing patch should include:

- a header-only `Morse_persistence` module;
- direct input from `Gudhi::Simplex_tree<>`;
- same-level Morse sequence construction on plateaus;
- fused Morse sequence and reference-map construction;
- ordinary persistence over `Z2`;
- local simplex ids mapped back to `Simplex_tree::Simplex_handle`;
- a small example;
- a maintainer-style C++ test matrix.

It should not include:

- Python bindings;
- benchmark machinery;
- prime-field or composite-ring coefficient APIs;
- cubical complexes;
- experimental flooding strategies as stable API.

## Public Strategy Set

The initial public strategy set should stay conservative:

```text
SAME_LEVEL_REDUCTION
F_MAX
F_MIN
PLATEAU_GREEDY
```

The broader MorseFrames prototype can keep saturated and flooding strategies for
experiments, but those should not complicate the first GUDHI discussion.

## Current Prototype Files

The current wrapper layer lives in:

```text
include/gudhi/Morse_persistence.h
include/gudhi/Morse_persistence/complex_view.h
include/gudhi/Morse_persistence/diagram.h
include/gudhi/Morse_persistence/morse_sequence.h
include/gudhi/Morse_persistence/persistence_reducer.h
include/gudhi/Morse_persistence/reference_map.h
include/gudhi/Morse_persistence/strategy.h
```

The upstream-shaped example is:

```text
examples/example_morse_persistence_from_simplex_tree.cpp
```

The current adapter test matrix is:

```text
tests/test_gudhi_simplex_tree_view.cpp
```

## Build Check

From the repository root, with GUDHI and Boost headers available:

```sh
cmake -S . -B build-gudhi-skeleton \
  -DMORSEFRAMES_GUDHI_INCLUDE_DIR=/path/to/gudhi/include \
  -DMORSEFRAMES_BOOST_INCLUDE_DIR=/path/to/boost/include

cmake --build build-gudhi-skeleton
ctest --test-dir build-gudhi-skeleton --output-on-failure
```

The expected optional targets are:

```text
morseframes_gudhi_simplex_tree_view_tests
morseframes_gudhi_simplex_tree_example
morseframes_gudhi_style_simplex_tree_example
morseframes_benchmark_gudhi_view
```

## Next Mechanical Step

The actual upstream patch should mechanically move the current prototype kernel
headers used by the wrapper into:

```text
include/gudhi/Morse_persistence/internal/
```

and rewrite internal includes from `morseframes/...` to
`gudhi/Morse_persistence/internal/...`.

That move should be done after the public API and test matrix are agreed upon,
because it is mostly namespace/include churn rather than algorithmic work.
