# GUDHI upstream patch map

This document maps the current Morse persistence prototype to a small first
GUDHI patch. It is not a patch itself; it is the checklist we can use before
copying code into a GUDHI branch.

The GUDHI source inspected for the current clean-checkout map is:

```text
https://github.com/GUDHI/gudhi-devel.git
HEAD 3a7d79c
```

## First patch objective

The first upstream patch should add a header-only C++ module that computes
Morse persistence directly on `Gudhi::Simplex_tree<>`:

- same-level Morse sequence construction;
- fused Morse sequence and reference-map construction;
- Morse persistence over `Z2` as the first public coefficient scope;
- local simplex ids mapped back to `Simplex_tree::Simplex_handle`;
- a tiny example and a small maintainer-style test matrix.

The first patch should not include Python bindings, benchmark machinery,
experimental flooding strategies as stable API, or cubical complexes.
The prototype already has prime-field experiments, but exposing coefficient
selection should be a later discussion rather than part of this first API.

## GUDHI module hook

GUDHI modules are enabled in the root `CMakeLists.txt` with `add_gudhi_module`.
The first upstream patch would add:

```cmake
add_gudhi_module(Morse_persistence)
```

near the existing module list in:

```text
CMakeLists.txt
```

The GUDHI development module machinery then looks for optional subdirectories
inside `src/<module>/`, such as:

```text
src/Morse_persistence/example/CMakeLists.txt
src/Morse_persistence/test/CMakeLists.txt
src/Morse_persistence/utilities/CMakeLists.txt
```

The `test` subdirectory is included when `WITH_GUDHI_TEST` and
`Boost::unit_test_framework` are available.

## Header files to move

The prototype currently exposes the GUDHI-shaped wrapper under:

```text
morseframes/include/gudhi/Morse_persistence.h
morseframes/include/gudhi/Morse_persistence/complex_view.h
morseframes/include/gudhi/Morse_persistence/diagram.h
morseframes/include/gudhi/Morse_persistence/morse_sequence.h
morseframes/include/gudhi/Morse_persistence/persistence_reducer.h
morseframes/include/gudhi/Morse_persistence/reference_map.h
morseframes/include/gudhi/Morse_persistence/strategy.h
```

In a GUDHI branch, these should become:

```text
src/Morse_persistence/include/gudhi/Morse_persistence.h
src/Morse_persistence/include/gudhi/Morse_persistence/complex_view.h
src/Morse_persistence/include/gudhi/Morse_persistence/diagram.h
src/Morse_persistence/include/gudhi/Morse_persistence/morse_sequence.h
src/Morse_persistence/include/gudhi/Morse_persistence/persistence_reducer.h
src/Morse_persistence/include/gudhi/Morse_persistence/reference_map.h
src/Morse_persistence/include/gudhi/Morse_persistence/strategy.h
```

The implementation currently depends on prototype kernel headers under:

```text
morseframes/include/morseframes/annotation.hpp
morseframes/include/morseframes/complex_view.hpp
morseframes/include/morseframes/debug_checks.hpp
morseframes/include/morseframes/filtered_complex.hpp
morseframes/include/morseframes/instrumentation.hpp
morseframes/include/morseframes/inverse_annotation_store.hpp
morseframes/include/morseframes/morse_reference_api.hpp
morseframes/include/morseframes/morse_sequence.hpp
morseframes/include/morseframes/reference_persistence.hpp
morseframes/include/morseframes/simplex_tree_builder.hpp
morseframes/include/morseframes/simplex_tree_morse.hpp
morseframes/include/morseframes/working_sets.hpp
```

For the first GUDHI branch, these should not remain under the local
`morseframes/` prototype include root. A clean layout would be:

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

The public wrappers should include only public GUDHI paths. Internal headers
should include other internal headers with `gudhi/Morse_persistence/internal/...`
paths. This is a mechanical namespace/include rewrite, not an algorithmic
change.

## Namespace map

The current prototype kernel uses:

```cpp
namespace morseframes { ... }
```

For GUDHI, the public API already uses:

```cpp
namespace Gudhi {
namespace morse_persistence {
  ...
}
}
```

Recommended upstream structure:

```cpp
namespace Gudhi {
namespace morse_persistence {
  // Public API.

namespace internal {
  // Implementation details.
}

}
}
```

This keeps the public namespace consistent with the wrapper already used in the
prototype and makes it clear which types are stable API.

## Public API to keep in the first patch

Keep the public API small:

```cpp
#include <gudhi/Morse_persistence.h>

auto result = Gudhi::morse_persistence::compute_morse_persistence(
    simplex_tree,
    Gudhi::morse_persistence::Morse_sequence_strategy::F_MAX);
```

Stable public types and functions:

```text
Morse_sequence_strategy
Simplex_id
Vertex_id
Level_id
Persistence_pair
Essential_interval
Persistence_diagram
Morse_step_type
Morse_step
Morse_sequence
Simplex_tree_view<SimplexTree>
Simplex_tree_morse_result<SimplexTree>
Simplex_tree_morse_frame<SimplexTree>
make_simplex_tree_view(simplex_tree)
compute_morse_sequence_and_reference_map(...)
compute_morse_persistence(...)
compute_morse_persistence_with_metrics(...)
off_diagonal_pairs(diagram)
```

Expose these strategies as the stable first set:

```text
SAME_LEVEL_REDUCTION
F_MAX
F_MIN
PLATEAU_GREEDY
```

Keep these internal for experiments and omit them from the first public
GUDHI-facing API:

```text
EXPERIMENTAL_SATURATED
EXPERIMENTAL_FLOODING_MAX
EXPERIMENTAL_FLOODING_MIN
EXPERIMENTAL_FLOODING_MINMAX
EXPERIMENTAL_FLOODING_MAXMIN
```

The wrapper enum now follows this conservative choice. The internal kernel still
keeps these strategies so the experimental report and benchmarks can continue to
compare them while the strategy taxonomy is settled.

## Documentation files

Add a Doxygen intro page:

```text
src/Morse_persistence/doc/Intro_morse_persistence.h
```

The minimal first page states:

- the module computes ordinary persistence via same-level Morse sequences and
  reference maps;
- it works directly on `Gudhi::Simplex_tree`;
- plateaus are handled directly, without lower-star refinement;
- the first public coefficient field is `Z2`;
- simplex handles in the result are valid only while the input Simplex-tree is
  alive and structurally unchanged.

The first staging patch does not update the generated
`src/common/doc/examples.h` file manually.

## Example files

The prototype candidate example is:

```text
morseframes/examples/example_morse_persistence_from_simplex_tree.cpp
```

In GUDHI it should move to:

```text
src/Morse_persistence/example/example_morse_persistence_from_simplex_tree.cpp
```

Add a module example CMake file:

```text
src/Morse_persistence/example/CMakeLists.txt
```

with:

```cmake
add_executable_with_targets(
  Morse_persistence_example_from_simplex_tree
  example_morse_persistence_from_simplex_tree.cpp
  TBB::tbb)

add_test(
  NAME Morse_persistence_example_from_simplex_tree
  COMMAND $<TARGET_FILE:Morse_persistence_example_from_simplex_tree>)
```

If maintainers prefer no `TBB::tbb` dependency for a header-only example, this
can be reduced to the same pattern used by simpler modules, but most GUDHI
examples use `add_executable_with_targets(... TBB::tbb)`.

## Test files

The prototype maintainer-style test is currently embedded in:

```text
morseframes/tests/test_gudhi_simplex_tree_view.cpp
```

For a GUDHI branch, split or rename it as:

```text
src/Morse_persistence/test/test_morse_persistence_simplex_tree.cpp
src/Morse_persistence/test/CMakeLists.txt
```

The current skeleton test uses a small runtime `CHECK(...)` helper rather than
plain `assert(...)`, because GUDHI release builds define `NDEBUG`. The CMake
file can be kept as an executable plus `add_test`, or converted to the GUDHI
Boost-test helper if maintainers prefer that convention:

```cmake
include(GUDHI_boost_test)

add_executable(
  Morse_persistence_test_simplex_tree
  test_morse_persistence_simplex_tree.cpp)

gudhi_add_boost_test(Morse_persistence_test_simplex_tree)
```

If the exact GUDHI test archive used for submission follows a different
module-test convention, keep the same test body but adapt the CMake wrapper to
that branch.

The first test matrix should include:

```text
single vertex
single increasing edge
filled triangle on one plateau
triangle with tail
two components joined late
1-cycle killed by later triangle
direct Simplex_tree view versus compact import
all stable public strategies on the same complex
```

Each strategy test should check:

```text
every simplex appears exactly once
regular pairs are same-level face/coface pairs
reference recurrence is valid
off-diagonal intervals match ordinary Z2 persistence
essential intervals match ordinary Z2 persistence
local ids map back to non-null Simplex_tree handles
```

For plateau examples, do not make zero-length interval multiplicity part of the
public test contract. Valid reductions may represent diagonal intervals
differently.

under a new `Morse_persistence` section, unless GUDHI's example-list generator
is run instead.

## CMake summary

Root file:

```text
CMakeLists.txt
```

Add:

```cmake
add_gudhi_module(Morse_persistence)
```

Example file:

```text
src/Morse_persistence/example/CMakeLists.txt
```

Test file:

```text
src/Morse_persistence/test/CMakeLists.txt
```

No utility executable is needed in the first patch. Benchmarks should stay out
of the first upstream patch.

## What to leave out of the first patch

Leave these in the prototype for now:

- Python bindings;
- benchmark scripts and CSV/report generation;
- flooding strategies as stable public API;
- prime-field coefficient selection as public API;
- composite `Z_n` coefficients;
- cubical-complex support;
- terrain/image benchmark datasets;
- adaptive strategy selection;
- the research experiment report.

This keeps the first GUDHI discussion about one thing: whether the
`Simplex_tree` Morse persistence API and implementation shape are acceptable.

## Pre-submission checklist

Before copying into a GUDHI branch:

1. Rewrite `morseframes::` implementation namespace into
   `Gudhi::morse_persistence::internal`.
2. Rewrite prototype include paths from `morseframes/...` to
   `gudhi/Morse_persistence/internal/...`.
3. Decide whether experimental strategies are omitted or left behind an
   explicit experimental name.
4. Keep the public coefficient scope at `Z2`; leave prime-field and composite
   coefficient choices outside the first patch.
5. Keep the runtime `CHECK(...)` helper or convert it into Boost unit tests if
   the target GUDHI branch expects Boost tests.
6. Build with:

```sh
cmake -S . -B build -DWITH_GUDHI_EXAMPLE=ON -DWITH_GUDHI_TEST=ON
cmake --build build --target Morse_persistence_example_from_simplex_tree
ctest --test-dir build -R Morse_persistence
```

7. Run the example manually and check that it prints one finite off-diagonal
   `H0` interval and one essential `H0` interval on the plateau-plus-tail
   complex.
