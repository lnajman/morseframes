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

- `include/gudhi/Morse_persistence.h`;
- public headers under `include/gudhi/Morse_persistence/`;
- private implementation headers under
  `include/gudhi/Morse_persistence/internal/`;
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

- copy the headers from this branch to `include/gudhi/Morse_persistence.h` and
  `include/gudhi/Morse_persistence/`;
- copy `examples/example_morse_persistence_from_simplex_tree.cpp` to
  `example/Morse_persistence/example_morse_persistence_from_simplex_tree.cpp`;
- copy `upstream/gudhi/example/Morse_persistence/CMakeLists.txt` to
  `example/Morse_persistence/CMakeLists.txt`;
- copy `upstream/gudhi/test/Morse_persistence/` to
  `test/Morse_persistence/`;
- add `add_gudhi_module(Morse_persistence)` to the top-level GUDHI module list.

## Dry-Run Validation

The bundle was checked on June 17, 2026 against a disposable copy of the local
GUDHI 3.12.0 source tree:

```text
/Users/laurentnajman/Documents/A trier au propre/PortableDisqueE/src/gudhi.3.12.0/
```

The local source tree was not a git checkout and already contained older
`Morse_persistence` artifacts, so the dry run first excluded those directories
and then copied in only the files from this staging branch.

The following commands passed in the disposable GUDHI tree:

```sh
cmake --build <build> --target Morse_persistence_test_simplex_tree
cmake --build <build-with-examples> --target Morse_persistence_example_from_simplex_tree
ctest --test-dir <build-with-examples> -R Morse_persistence --output-on-failure
```

The final `ctest` run passed both:

```text
Morse_persistence_example_from_simplex_tree
Morse_persistence_test_simplex_tree
```

This directory is therefore a review bundle and placement guide. It is not a
`git format-patch` generated from an upstream GUDHI repository.
