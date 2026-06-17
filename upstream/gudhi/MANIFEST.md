# GUDHI Morse Persistence Manifest

This manifest lists the files that belong to the first upstream-style
`Morse_persistence` bundle and their intended destinations in a clean GUDHI
development checkout.

## Headers

Copy:

```text
include/gudhi/Morse_persistence.h
include/gudhi/Morse_persistence/complex_view.h
include/gudhi/Morse_persistence/diagram.h
include/gudhi/Morse_persistence/morse_sequence.h
include/gudhi/Morse_persistence/persistence_reducer.h
include/gudhi/Morse_persistence/reference_map.h
include/gudhi/Morse_persistence/strategy.h
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

to:

```text
src/Morse_persistence/include/gudhi/Morse_persistence.h
src/Morse_persistence/include/gudhi/Morse_persistence/complex_view.h
src/Morse_persistence/include/gudhi/Morse_persistence/diagram.h
src/Morse_persistence/include/gudhi/Morse_persistence/morse_sequence.h
src/Morse_persistence/include/gudhi/Morse_persistence/persistence_reducer.h
src/Morse_persistence/include/gudhi/Morse_persistence/reference_map.h
src/Morse_persistence/include/gudhi/Morse_persistence/strategy.h
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

The public namespace is:

```cpp
Gudhi::morse_persistence
```

The implementation namespace is:

```cpp
Gudhi::morse_persistence::internal
```

## Documentation

Copy:

```text
upstream/gudhi/src/Morse_persistence/doc/Intro_morse_persistence.h
```

to:

```text
src/Morse_persistence/doc/Intro_morse_persistence.h
```

## Example

Copy:

```text
examples/example_morse_persistence_from_simplex_tree.cpp
```

to:

```text
src/Morse_persistence/example/example_morse_persistence_from_simplex_tree.cpp
```

Copy:

```text
upstream/gudhi/src/Morse_persistence/example/CMakeLists.txt
```

to:

```text
src/Morse_persistence/example/CMakeLists.txt
```

## Test

Copy:

```text
upstream/gudhi/src/Morse_persistence/test/CMakeLists.txt
upstream/gudhi/src/Morse_persistence/test/test_morse_persistence_simplex_tree.cpp
```

to:

```text
src/Morse_persistence/test/CMakeLists.txt
src/Morse_persistence/test/test_morse_persistence_simplex_tree.cpp
```

The test is intentionally pure GUDHI-facing: it includes only GUDHI headers and
does not depend on `morseframes/...`.

## Module Registration

Add the module to the top-level GUDHI `CMakeLists.txt`:

```cmake
add_gudhi_module(Morse_persistence)
```

This should be placed near the other header-only topology modules.

## Smoke Commands

From a GUDHI build directory configured with examples and tests enabled:

```sh
cmake --build . --target Morse_persistence_test_simplex_tree
cmake --build . --target Morse_persistence_example_from_simplex_tree
ctest -R Morse_persistence --output-on-failure
```

The corresponding MorseFrames staging check is:

```sh
cmake -S . -B build-gudhi-skeleton \
  -DMORSEFRAMES_GUDHI_INCLUDE_DIR=/path/to/gudhi/src/Simplex_tree/include \
  -DMORSEFRAMES_BOOST_INCLUDE_DIR=/path/to/boost/include
cmake --build build-gudhi-skeleton
ctest --test-dir build-gudhi-skeleton --output-on-failure
```
