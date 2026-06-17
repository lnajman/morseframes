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

The public wrapper now reaches the kernel only through:

```text
include/gudhi/Morse_persistence/internal/
```

Those internal files are a copied, GUDHI-namespaced snapshot of the prototype
headers needed by the wrapper. They live in
`Gudhi::morse_persistence::internal` and include each other through
`gudhi/Morse_persistence/internal/...`. The copied snapshot is trimmed to the
public strategy set: `SAME_LEVEL_REDUCTION`, `F_MAX`, `F_MIN`, and
`PLATEAU_GREEDY`.

The upstream-shaped example is:

```text
examples/example_morse_persistence_from_simplex_tree.cpp
```

The current adapter test matrix is:

```text
tests/test_gudhi_simplex_tree_view.cpp
```

Detailed follow-up notes:

```text
docs/gudhi_upstream_file_map.md
docs/gudhi_maintainer_scope_note.md
upstream/gudhi/README.md
upstream/gudhi/MANIFEST.md
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
morseframes_gudhi_upstream_simplex_tree_tests
morseframes_gudhi_simplex_tree_example
morseframes_gudhi_style_simplex_tree_example
morseframes_benchmark_gudhi_view
```

## Current Test Shape

The adapter test matrix exercises the public `Gudhi::morse_persistence` API for
all four public strategies on the maintainer-style examples, then keeps a small
set of prototype-level checks for invariants that are easier to inspect through
the lower-level complex view.

## Pure GUDHI Test Candidate

The pure GUDHI test candidate, extracted from the public-API subset of
`tests/test_gudhi_simplex_tree_view.cpp`, lives at:

```text
upstream/gudhi/test/Morse_persistence/test_morse_persistence_simplex_tree.cpp
```

## Next Mechanical Step

The dry run on June 17, 2026 used a disposable copy of the local GUDHI 3.12.0
source tree, with pre-existing Morse_persistence artifacts excluded before
copying the staging files in.

The following passed:

```sh
ctest -R Morse_persistence --output-on-failure
```

with both:

```text
Morse_persistence_example_from_simplex_tree
Morse_persistence_test_simplex_tree
```

The compact review bundle is now recorded in:

```text
upstream/gudhi/README.md
upstream/gudhi/MANIFEST.md
```

The next mechanical step is to decide whether to generate an actual patch file
from a clean GUDHI git checkout, or keep this as a copy-ready staging bundle
until the maintainer-facing scope is agreed.
