# GUDHI Upstream Staging Bundle

This directory contains copy-ready GUDHI-side files for the first
`Morse_persistence` module discussion. It is intentionally smaller than the
full MorseFrames prototype: the goal is to give maintainers a focused C++ patch
that can be reviewed in GUDHI terms.

The staged API is header-only, accepts a `Gudhi::Simplex_tree<>`, computes
same-level Morse sequences on plateaus, builds the reference map while building
the sequence, and computes ordinary persistence over `Z2`.

## Scope

The first GUDHI-facing bundle includes:

- `src/Morse_persistence/include/gudhi/Morse_persistence.h`;
- public headers under
  `src/Morse_persistence/include/gudhi/Morse_persistence/`;
- private implementation headers under
  `src/Morse_persistence/include/gudhi/Morse_persistence/internal/`;
- one Doxygen intro page under `src/Morse_persistence/doc/`;
- one `Simplex_tree` example;
- one maintainer-style C++ test matrix.

The public strategy set is:

```text
SAME_LEVEL_REDUCTION
F_MAX
F_MIN
PLATEAU_GREEDY
```

The first bundle does not include Python bindings, benchmark machinery,
prime-field or composite-ring coefficient APIs, cubical complexes, or
experimental flooding strategies as stable API.

## Copy Targets

The copy targets are listed in `MANIFEST.md`. In short:

- copy the headers from this branch to
  `src/Morse_persistence/include/gudhi/Morse_persistence.h` and
  `src/Morse_persistence/include/gudhi/Morse_persistence/`;
- copy `examples/example_morse_persistence_from_simplex_tree.cpp` to
  `src/Morse_persistence/example/example_morse_persistence_from_simplex_tree.cpp`;
- copy `upstream/gudhi/src/Morse_persistence/example/CMakeLists.txt` to
  `src/Morse_persistence/example/CMakeLists.txt`;
- copy `upstream/gudhi/src/Morse_persistence/test/` to
  `src/Morse_persistence/test/`;
- copy `upstream/gudhi/src/Morse_persistence/doc/` to
  `src/Morse_persistence/doc/`;
- add `add_gudhi_module(Morse_persistence)` to the top-level GUDHI module list.

## Clean-Checkout Validation

The bundle was checked on June 17, 2026 against a shallow clean checkout of
the public GUDHI development repository:

```text
https://github.com/GUDHI/gudhi-devel.git
HEAD 3a7d79c
```

The following commands passed in the disposable GUDHI tree:

```sh
cmake -S <gudhi-checkout> -B <build> \
  -DWITH_GUDHI_PYTHON=OFF \
  -DWITH_GUDHI_GUDHUI=OFF \
  -DWITH_GUDHI_EXAMPLE=ON \
  -DWITH_GUDHI_TEST=ON
cmake --build <build> --target Morse_persistence_test_simplex_tree
cmake --build <build> --target Morse_persistence_example_from_simplex_tree
ctest --test-dir <build> -R Morse_persistence --output-on-failure
```

The final `ctest` run passed both:

```text
Morse_persistence_example_from_simplex_tree
Morse_persistence_test_simplex_tree
```

This directory is therefore a review bundle and placement guide. It is not a
submitted pull request.
