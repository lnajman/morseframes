# GUDHI Maintainer Scope Note

This is a draft note for discussing the GUDHI-facing Morse persistence
prototype with maintainers. It avoids relying on the unreleased Morse
persistence manuscript.

## Proposed Scope

The proposed contribution is a small header-only C++ module:

```cpp
#include <gudhi/Morse_persistence.h>
```

The first entry point computes ordinary persistence over `Z2` directly from a
filtered `Gudhi::Simplex_tree<>`:

```cpp
auto result = Gudhi::morse_persistence::compute_morse_persistence(
    simplex_tree,
    Gudhi::morse_persistence::Morse_sequence_strategy::F_MAX);
```

The result stores:

- the Morse sequence;
- finite persistence pairs;
- essential intervals;
- basic construction and reduction metrics;
- a view that maps local simplex ids back to `Simplex_tree::Simplex_handle`.

## What Is New

The module constructs a same-level Morse sequence on the filtered complex and
uses the associated reference map as the input to a persistence reducer.

The important point for GUDHI discussion is plateau handling. Equal filtration
values are treated as genuine filtration levels. The algorithm can pair
simplexes within a plateau and does not need to first impose a lower-star or
simplex-wise refinement of that plateau as part of the Morse sequence
construction.

The first public strategy set is deliberately small:

```text
SAME_LEVEL_REDUCTION
F_MAX
F_MIN
PLATEAU_GREEDY
```

The broader MorseFrames prototype contains other experimental strategies, but
they are not part of the proposed first GUDHI API.

## What Is Not Included

The first patch should not include:

- Python bindings;
- benchmarks;
- coefficient fields beyond `Z2`;
- cubical complexes;
- flooding or saturated experimental strategies;
- a new public complex data structure.

Those can be discussed later if the C++ `Simplex_tree` API is considered useful.

## Relation To Existing GUDHI Persistence

This module is not intended to replace `Persistent_cohomology` or the existing
`Simplex_tree` persistence path. It is a complementary reduction path that first
builds a Morse frame and then reduces annotations indexed by critical
simplexes.

The current validation approach is conservative:

- compare off-diagonal finite intervals with standard persistence;
- compare essential intervals;
- test plateau-heavy examples explicitly;
- test that local ids map back to valid `Simplex_tree` handles;
- keep the public API independent of the MorseFrames prototype namespace.

The experimental motivation is that, on plateau-heavy inputs, the Morse frame
can be much smaller than the original filtration. The first GUDHI conversation
should focus on API shape, correctness checks, and maintainability before making
strong performance claims.

## Current Limitations

The current GUDHI-facing prototype computes ordinary persistence over `Z2`.
Prime fields are implemented in the broader MorseFrames prototype, but they are
not yet promoted to the GUDHI-facing API because the first patch should stay
small.

The implementation currently targets simplicial complexes through
`Gudhi::Simplex_tree<>`. Cubical complexes are feasible in principle, but would
need a separate complex-view adapter and coefficient-incidence handling, so they
should not be part of the first contribution.

## Suggested First Discussion

The most useful first maintainer discussion is probably:

- Does `Gudhi::morse_persistence` look like an acceptable namespace and module
  boundary?
- Is a header-only module acceptable for this first version?
- Should the result expose local simplex ids plus a handle map, or should it
  expose handles directly in intervals?
- Should all four strategies be public initially, or should the first API expose
  only `F_MAX` plus one simpler baseline?
- What level of documentation is needed before a pull request?
