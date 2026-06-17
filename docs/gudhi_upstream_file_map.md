# GUDHI Upstream File Map

This note records the mechanical placement of the current GUDHI-facing
MorseFrames branch in a GUDHI development checkout. It is a staging checklist,
not yet a pull request.

The current branch is:

```text
gudhi-upstream-skeleton
```

The clean checkout used for validation was:

```text
https://github.com/GUDHI/gudhi-devel.git
HEAD 3a7d79c
```

## Header Files

Copy the public umbrella header to:

```text
src/Morse_persistence/include/gudhi/Morse_persistence.h
```

Copy the public module headers to:

```text
src/Morse_persistence/include/gudhi/Morse_persistence/complex_view.h
src/Morse_persistence/include/gudhi/Morse_persistence/diagram.h
src/Morse_persistence/include/gudhi/Morse_persistence/morse_sequence.h
src/Morse_persistence/include/gudhi/Morse_persistence/persistence_reducer.h
src/Morse_persistence/include/gudhi/Morse_persistence/reference_map.h
src/Morse_persistence/include/gudhi/Morse_persistence/strategy.h
```

Copy the trimmed internal implementation snapshot to:

```text
src/Morse_persistence/include/gudhi/Morse_persistence/internal/annotation.h
src/Morse_persistence/include/gudhi/Morse_persistence/internal/complex_view.h
src/Morse_persistence/include/gudhi/Morse_persistence/internal/field_annotation_store.h
src/Morse_persistence/include/gudhi/Morse_persistence/internal/field_arithmetic.h
src/Morse_persistence/include/gudhi/Morse_persistence/internal/filtered_complex.h
src/Morse_persistence/include/gudhi/Morse_persistence/internal/inverse_annotation_store.h
src/Morse_persistence/include/gudhi/Morse_persistence/internal/morse_reference_api.h
src/Morse_persistence/include/gudhi/Morse_persistence/internal/morse_sequence.h
src/Morse_persistence/include/gudhi/Morse_persistence/internal/reference_persistence.h
src/Morse_persistence/include/gudhi/Morse_persistence/internal/simplex_tree_builder.h
src/Morse_persistence/include/gudhi/Morse_persistence/internal/simplex_tree_morse.h
src/Morse_persistence/include/gudhi/Morse_persistence/internal/working_sets.h
```

The copied internal headers are intentionally private implementation details.
The public surface should stay under `Gudhi::morse_persistence`; the internal
code lives under `Gudhi::morse_persistence::internal`.

## Documentation

Copy the module intro page:

```text
upstream/gudhi/src/Morse_persistence/doc/Intro_morse_persistence.h
```

to:

```text
src/Morse_persistence/doc/Intro_morse_persistence.h
```

## Example

Copy the upstream-shaped example:

```text
examples/example_morse_persistence_from_simplex_tree.cpp
```

to:

```text
src/Morse_persistence/example/example_morse_persistence_from_simplex_tree.cpp
```

The matching copy-ready CMake stub is:

```text
upstream/gudhi/src/Morse_persistence/example/CMakeLists.txt
```

Place it at:

```text
src/Morse_persistence/example/CMakeLists.txt
```

## Tests

The copy-ready GUDHI test candidate is:

```text
upstream/gudhi/src/Morse_persistence/test/test_morse_persistence_simplex_tree.cpp
```

Place it at:

```text
src/Morse_persistence/test/test_morse_persistence_simplex_tree.cpp
```

The matching copy-ready CMake stub is:

```text
upstream/gudhi/src/Morse_persistence/test/CMakeLists.txt
```

Place it at:

```text
src/Morse_persistence/test/CMakeLists.txt
```

The candidate file covers:

- one vertex;
- one increasing edge;
- one filled triangle entirely on one plateau;
- a filled triangle with a later tail edge;
- two components joined late;
- a one-cycle killed by a later triangle.

For each case, it tests all public strategies:

```text
SAME_LEVEL_REDUCTION
F_MAX
F_MIN
PLATEAU_GREEDY
```

The expected checks include:

- finite off-diagonal intervals;
- essential intervals;
- dimensions of finite intervals;
- ability to map diagram simplex ids back to `Simplex_tree::Simplex_handle`;
- rejection of non-public strategy names such as `saturated`.

The older MorseFrames integration test remains:

```text
tests/test_gudhi_simplex_tree_view.cpp
```

It intentionally keeps additional prototype-level oracle checks through
`morseframes::...`; those checks should not be copied into GUDHI.

## Module Registration

In a clean upstream GUDHI checkout, add:

```cmake
add_gudhi_module(Morse_persistence)
```

to the module list in the top-level `CMakeLists.txt`, near the other
header-only topology modules.

## Local Build Check

From a GUDHI build directory, the intended smoke check is:

```sh
cmake --build . --target Morse_persistence_example_from_simplex_tree
cmake --build . --target Morse_persistence_test_simplex_tree
ctest -R Morse_persistence --output-on-failure
```

For the MorseFrames staging branch itself, use:

```sh
cmake -S . -B build-gudhi-skeleton \
  -DMORSEFRAMES_GUDHI_INCLUDE_DIR=/path/to/gudhi/src/Simplex_tree/include \
  -DMORSEFRAMES_BOOST_INCLUDE_DIR=/path/to/boost/include

cmake --build build-gudhi-skeleton
ctest --test-dir build-gudhi-skeleton --output-on-failure
```

## Not In The First Patch

Do not include these in the first GUDHI discussion patch:

- Python bindings;
- benchmark scripts and generated benchmark tables;
- prime-field or composite-ring coefficient APIs;
- cubical complexes;
- flooding and saturated experimental strategies.
