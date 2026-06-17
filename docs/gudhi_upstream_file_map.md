# GUDHI Upstream File Map

This note records the mechanical placement of the current GUDHI-facing
MorseFrames branch in a GUDHI checkout. It is a staging checklist, not yet a
claim that the patch is ready for submission.

The current branch is:

```text
gudhi-upstream-skeleton
```

## Header Files

Copy the public umbrella header to:

```text
include/gudhi/Morse_persistence.h
```

Copy the public module headers to:

```text
include/gudhi/Morse_persistence/complex_view.h
include/gudhi/Morse_persistence/diagram.h
include/gudhi/Morse_persistence/morse_sequence.h
include/gudhi/Morse_persistence/persistence_reducer.h
include/gudhi/Morse_persistence/reference_map.h
include/gudhi/Morse_persistence/strategy.h
```

Copy the trimmed internal implementation snapshot to:

```text
include/gudhi/Morse_persistence/internal/annotation.h
include/gudhi/Morse_persistence/internal/complex_view.h
include/gudhi/Morse_persistence/internal/field_annotation_store.h
include/gudhi/Morse_persistence/internal/field_arithmetic.h
include/gudhi/Morse_persistence/internal/filtered_complex.h
include/gudhi/Morse_persistence/internal/inverse_annotation_store.h
include/gudhi/Morse_persistence/internal/morse_reference_api.h
include/gudhi/Morse_persistence/internal/morse_sequence.h
include/gudhi/Morse_persistence/internal/reference_persistence.h
include/gudhi/Morse_persistence/internal/simplex_tree_builder.h
include/gudhi/Morse_persistence/internal/simplex_tree_morse.h
include/gudhi/Morse_persistence/internal/working_sets.h
```

The copied internal headers are intentionally private implementation details.
The public surface should stay under `Gudhi::morse_persistence`; the internal
code lives under `Gudhi::morse_persistence::internal`.

## Example

Copy the upstream-shaped example:

```text
examples/example_morse_persistence_from_simplex_tree.cpp
```

to:

```text
example/Morse_persistence/example_morse_persistence_from_simplex_tree.cpp
```

A minimal GUDHI `example/Morse_persistence/CMakeLists.txt` is:

```cmake
add_executable_with_targets(
  Morse_persistence_example_from_simplex_tree
  example_morse_persistence_from_simplex_tree.cpp
  TBB::tbb)

add_test(
  NAME Morse_persistence_example_from_simplex_tree
  COMMAND $<TARGET_FILE:Morse_persistence_example_from_simplex_tree>)
```

## Tests

The current MorseFrames test:

```text
tests/test_gudhi_simplex_tree_view.cpp
```

contains two kinds of checks:

- public API checks through `Gudhi::morse_persistence`;
- prototype-level oracle checks through `morseframes::...`.

For a GUDHI patch, keep only the public API checks and expected barcode checks,
then place the result at:

```text
test/Morse_persistence/test_morse_persistence_simplex_tree.cpp
```

A minimal GUDHI `test/Morse_persistence/CMakeLists.txt` is:

```cmake
add_executable_with_targets(
  Morse_persistence_test_simplex_tree
  test_morse_persistence_simplex_tree.cpp
  TBB::tbb)

add_test(
  NAME Morse_persistence_test_simplex_tree
  COMMAND $<TARGET_FILE:Morse_persistence_test_simplex_tree>)
```

The upstream test should cover at least:

- one vertex;
- one increasing edge;
- one filled triangle entirely on one plateau;
- a filled triangle with a later tail edge;
- two components joined late;
- a one-cycle killed by a later triangle.

For each case, test all public strategies:

```text
SAME_LEVEL_REDUCTION
F_MAX
F_MIN
PLATEAU_GREEDY
```

The expected checks should include:

- finite off-diagonal intervals;
- essential intervals;
- dimensions of finite intervals;
- ability to map diagram simplex ids back to `Simplex_tree::Simplex_handle`;
- rejection of non-public strategy names such as `saturated`.

## Module Registration

In a clean upstream GUDHI checkout, add:

```cmake
add_gudhi_module(Morse_persistence)
```

to the module list in the top-level `CMakeLists.txt`, near the other
header-only topology modules.

The local GUDHI checkout used during MorseFrames development already contains a
`Morse_persistence` module entry and matching `example/` and `test/`
directories. A clean checkout should not assume those are present.

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
  -DMORSEFRAMES_GUDHI_INCLUDE_DIR=/path/to/gudhi/include \
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
