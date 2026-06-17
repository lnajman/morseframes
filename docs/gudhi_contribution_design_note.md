# GUDHI contribution design note

This note sketches how the current Morse persistence prototype could be shaped
into a small GUDHI-compatible C++ contribution. It is intentionally focused on
the software interface: the goal is to describe the API surface, the internal
adapter, the assumptions on `Gudhi::Simplex_tree`, and the tests that should be
in place before discussing integration upstream.

For the concrete file-by-file patch plan against a GUDHI source tree, see
`gudhi_upstream_patch_map.md`.

## Scope

The first candidate contribution should compute, directly from a filtered
`Gudhi::Simplex_tree<>`:

- a same-level Morse sequence;
- the associated reference map, preferably fused with sequence construction;
- a reduced Morse persistence diagram over `Z2`;
- optional counters for timing and structural diagnostics.

This should not require lower-star refinement. Plateaus are part of the input
filtration and are handled directly by same-level pairings.

The prototype kernel and Python layer now also support prime fields `F_p` for
experiments. That should remain outside the first GUDHI-facing patch. A `Z2`
public API matches the existing ordinary-persistence level we want to compare
against first, and a coefficient parameter can be promoted later if maintainers
want it. Composite `Z_n` coefficients are intentionally out of scope here.

The first version should stay focused on simplicial complexes. Cubical
complexes and alternative Morse sequence algorithms should remain possible
extensions, but they should not complicate the initial Simplex-tree API.

## Proposed public layout

The prototype now provides one umbrella header and a small internal
subdirectory:

```text
include/gudhi/Morse_persistence.h
include/gudhi/Morse_persistence/complex_view.h
include/gudhi/Morse_persistence/morse_sequence.h
include/gudhi/Morse_persistence/reference_map.h
include/gudhi/Morse_persistence/persistence_reducer.h
include/gudhi/Morse_persistence/strategy.h
include/gudhi/Morse_persistence/diagram.h
```

These public wrapper headers are currently backed by the prototype kernel files:

```text
include/morseframes/morse_reference_api.hpp
include/morseframes/simplex_tree_morse.hpp
include/morseframes/simplex_tree_builder.hpp
include/morseframes/morse_sequence.hpp
include/morseframes/reference_persistence.hpp
include/morseframes/annotation.hpp
```

For an upstream discussion, `include/morseframes/...` should remain a local
prototype namespace only. The GUDHI-facing names now live under
`include/gudhi/...`.

## Namespace and naming

The public wrapper namespace is:

```cpp
namespace Gudhi::morse_persistence {
  // Public API.
}
```

This is consistent with nested GUDHI components such as
`Gudhi::persistence_matrix`, while keeping the module separate from
`Gudhi::Simplex_tree` itself.

The wrapper layer currently uses the following GUDHI-shaped names:

```text
morseframes::MorseSequenceStrategy  -> Gudhi::morse_persistence::Morse_sequence_strategy
morseframes::MorseSequence          -> Gudhi::morse_persistence::Morse_sequence
morseframes::SimplexTreeComplexView -> Gudhi::morse_persistence::Simplex_tree_view
morseframes::SimplexTreeMorseReferenceResult
                                    -> Gudhi::morse_persistence::Simplex_tree_morse_result
compute_simplex_tree_morse_reference_persistence
                                    -> compute_morse_persistence
build_morse_reference_frame         -> compute_morse_sequence_and_reference_map
```

The public API should avoid prototype-specific terms where possible. In
particular, "reference-map persistence" is precise for the paper, but the
GUDHI entry point can simply be `compute_morse_persistence`, with detailed
documentation explaining that the implementation uses the reference map.

## Public API sketch

A minimal first API could look like this:

```cpp
#include <gudhi/Simplex_tree.h>
#include <gudhi/Morse_persistence.h>

Gudhi::Simplex_tree<> st;
// Insert simplices, assign filtrations, then initialize the filtration cache.
st.initialize_filtration();

auto result = Gudhi::morse_persistence::compute_morse_persistence(
    st,
    Gudhi::morse_persistence::Morse_sequence_strategy::F_MAX);

for (const auto& interval : result.finite_intervals()) {
  auto birth_simplex = result.simplex_tree_handle(interval.birth_simplex());
  auto death_simplex = result.simplex_tree_handle(interval.death_simplex());
}
```

The result should own the temporary view so that local simplex ids used in the
Morse sequence and diagram can still be mapped back to `Simplex_tree` handles.
The input simplex tree must outlive the result.

A slightly more explicit API should also be available for experiments:

```cpp
auto frame = Gudhi::morse_persistence::compute_morse_sequence_and_reference_map(
    st,
    Gudhi::morse_persistence::Morse_sequence_strategy::F_MAX);

auto diagram = Gudhi::morse_persistence::compute_morse_persistence(st, frame);
```

This keeps the interface open to new sequence builders, including flooding and
future algorithms, without changing the reducer.

## Internal data structure

The implementation should not add another owning simplex-tree data structure
for the main path. The current prototype uses a `SimplexTreeComplexView`: a
read-only adapter over `Gudhi::Simplex_tree<>`.

The view builds the pieces needed by the Morse algorithms:

- contiguous local simplex ids in `[0, size())`;
- maps from local ids to `Simplex_tree::Simplex_handle`;
- sorted vertex tuples for deterministic tie-breaking;
- boundary and coboundary lists in local ids;
- filtration values, level ids, and simplices grouped by level.

This is the useful middle ground. The algorithm gets dense ids and direct
boundary/coboundary access, while the user keeps the original `Simplex_tree`.
After this view is built, a separate trie-like structure is unlikely to help
the Morse reduction itself, because the hot operations are over local ids,
annotations, inverse annotation lists, and same-level boundary/coboundary
queries.

The compact owning complex used in the prototype should remain useful for
testing, Python bindings, file import, and synthetic benchmarks, but it should
not be the first public GUDHI entry point.

## Required Simplex_tree assumptions

The first version should document these assumptions explicitly:

- The input is a valid filtered simplicial complex stored in
  `Gudhi::Simplex_tree`.
- Filtration values are monotone on faces.
- The filtration cache has been initialized, or the API calls
  `initialize_filtration()` on a local copy when mutation is acceptable.
- The input tree is not modified while the Morse result is used.
- Simplex handles returned by the result are only valid as long as the input
  tree remains alive and structurally unchanged.
- The first public GUDHI implementation works over `Z2`.
- The first implementation targets ordinary persistence, not extended
  persistence.
- Equal filtration values are allowed and are treated as genuine plateaus.

For deterministic behavior, tie-breaking should use filtration level,
dimension, and vertex tuple, not internal Simplex-tree handle order.

## Strategy set for the first version

The first GUDHI-facing implementation should expose only the strategies we can
explain and test clearly:

```text
SAME_LEVEL_REDUCTION
F_MAX
F_MIN
PLATEAU_GREEDY
```

The prototype also contains saturated and flooding variants. They are useful for
experiments, but they can be kept behind an experimental option until the names
and relation to the papers are settled.

The reducer should accept a precomputed `Morse_sequence` and reference map.
This is important: it lets future sequence algorithms be added without
rewriting persistence.

## Tiny maintainer test matrix

The first upstream-style test set should be small and precise:

| Test | Purpose |
| --- | --- |
| single vertex | essential `H0` interval and no finite interval |
| single edge with increasing filtration | one finite `H0` interval and one essential `H0` interval |
| filled triangle with all simplices on one plateau | direct plateau handling, no lower-star refinement |
| triangle plus tail with several levels | comparison against standard persistence |
| disconnected components with one late edge | multiple `H0` births and one merge |
| one 1-cycle killed by a later 2-simplex | finite `H1` interval |
| direct Simplex-tree view vs compact import | identical barcode and identical strategy signature |
| all public strategies on the same complex | valid sequence, valid references, same persistence barcode |

For each strategy, tests should check:

- every simplex appears exactly once, either as critical or in one regular pair;
- each regular pair is a same-level face/coface pair;
- the reference recurrence holds;
- the off-diagonal persistence barcode and essential intervals agree with
  GUDHI's ordinary persistence or with an independent standard `Z2` reducer;
- local simplex ids can be mapped back to non-null Simplex-tree handles.

For plateau examples, zero-length intervals should not be the main public
contract. Different valid reductions may represent diagonal intervals
differently, while the off-diagonal intervals and essential classes should be
stable.

## Benchmark shape

The benchmark should not try to prove absolute superiority in a first PR. It
should answer narrower engineering questions:

- cost of building the `Simplex_tree` view;
- cost of computing the Morse sequence and reference map;
- number of critical simplices;
- reducer time after Morse compression;
- comparison with ordinary persistence on the same `Simplex_tree` input.

The current native benchmark already separates direct view construction,
compact import, sequence/reference construction, and reducer time. That shape is
good for deciding whether the adapter is worth integrating.

## Current prototype status

The prototype already has the key ingredients:

- generic `ComplexView` templates in `morse/morse_reference_api.hpp`;
- a direct `Gudhi::Simplex_tree<>` adapter in `morse/simplex_tree_morse.hpp`;
- GUDHI-shaped wrapper headers under `include/gudhi/Morse_persistence`;
- deterministic strategy signatures across direct view and compact import;
- a minimal example in `examples/gudhi_simplex_tree_morse.cpp`;
- a candidate upstream-style example in
  `examples/example_morse_persistence_from_simplex_tree.cpp`;
- C++ tests in `tests/test_gudhi_simplex_tree_view.cpp`, including the small
  maintainer matrix above;
- a native benchmark in `benchmarks/benchmark_gudhi_view.cpp`;
- prototype prime-field support and coefficient-overhead benchmarks, kept out
  of the first GUDHI-facing API.

The remaining work before a real GUDHI patch is mostly upstream polishing:

1. Decide whether the wrapper names match GUDHI maintainers' preferred style.
2. Expand Doxygen comments where maintainers expect concept-level documentation.
3. Keep prototype-only experimental strategies out of the public API.
4. Move the candidate example into GUDHI's `example/` layout if the module is
   accepted.
5. Port the maintainer matrix into GUDHI's test layout when the module layout is
   fixed.
6. Decide after API review whether prime-field coefficients should become a
   public follow-up; keep composite `Z_n` out of the first patch.
