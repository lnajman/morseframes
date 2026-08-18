# Morse Sequence Strategies

MorseFrames separates the construction of a Morse sequence from the later
reference/coreference and persistence computations. This is intentional: the
same public pipeline can be used with several sequence constructors, and new
constructors should be easy to add.

All strategies exposed here produce a simplex-wise Morse sequence compatible
with the input filtration. Regular pairs are restricted to a single filtration
level. Most methods work directly on complexes with plateaus;
`process-lower-stars` deliberately has the stricter classical lower-star input
contract described below.

A flooding sequence is an `F`-sequence whose order is globally nondecreasing
with respect to the filtration: once a simplex of value `lambda` has appeared,
no later simplex has a smaller filtration value. The saturated, same-level,
plateau-greedy, and named `flooding-*` strategies are flooding constructions in
this sense. The `f-max` and `f-min` strategies are valid `F`-sequence
constructors, but they are global seed-and-expand strategies and are not
required to be flooding.

## Canonical Strategy Names

The Python API accepts these canonical names:

```python
"saturated"
"f-max"
"process-lower-stars"
"process-lower-stars-parallel"
"f-min"
"same-level-reduction"
"plateau-greedy"
"flooding-max"
"flooding-min"
"flooding-reduction-kernel"
"flooding-reduction-kernel-parallel"
"flooding-minmax"
"flooding-maxmin"
```

The default is `"saturated"`.

Some legacy aliases are still accepted for compatibility, but documentation and
new examples should use the canonical names above.

## Common Vocabulary

The implementation uses the following operational vocabulary.

- A simplex is **fillable** in an increasing construction when all its boundary
  faces have already been inserted.
- A simplex is **pairable** in an increasing construction when it has exactly
  one missing boundary face, and that missing face lies in the same filtration
  level.
- A **coperforation** inserts a fillable simplex as critical.
- A decreasing construction uses the dual language: remaining cofaces,
  reductions, and perforations.

The words "min" and "max" refer to the direction and priority convention of the
sequence constructor, not to the coefficient field or to a choice of homology
versus cohomology.

## Saturated

Canonical name:

```python
"saturated"
```

The saturated strategy scans filtration levels in increasing order. Inside each
level, it repeatedly performs all currently available same-level regular pairs.
When no such pair is available, it marks one fillable simplex critical and then
continues. This is the simplest flooding construction exposed by the package.

Informally:

```text
pair until stuck,
mark one fillable simplex critical,
continue.
```

This is a good default because it is simple, plateau-aware, and usually gives a
small critical complex. Its exact result depends on deterministic tie-breaking
inside a filtration plateau.

## F-Max

Canonical name:

```python
"f-max"
```

`f-max` implements a forward seed-and-expand construction. It scans the simplex
filtration in increasing order. When the current admissible pair queue is empty,
the next not-yet-inserted simplex, often a `0`-simplex, becomes a critical seed.
The builder then gives priority to all same-level coreduction-like pairs
unlocked by the current state: an upper simplex `tau` is paired with its unique
remaining boundary face `sigma` when both lie at the same level. Once no such
pair remains, the next seed is chosen.

Informally:

```text
prefer same-level coreduction-like pairs,
otherwise add the next available critical seed.
```

This is the implementation corresponding to the `Max(S,F)` style used in our
experiments. It is a valid `F`-sequence constructor, but it is not necessarily a
flooding construction.

## ProcessLowerStars

Canonical name:

```python
"process-lower-stars"
```

This is the simplicial ProcessLowerStars construction. It orders the vertices
by increasing filtration value and partitions every simplex into the lower star
of its unique maximum vertex. Each lower star is then processed independently
with the forward one-missing-face expansion used by `f-max`. Candidate cells
are ordered by the Robins key: the ranks of their vertices in decreasing order,
compared lexicographically.

The current implementation intentionally enforces the classical generic-input
case:

- every cell is a nonempty simplex;
- vertex filtration values are injective; and
- every simplex has the same filtration level as its maximum vertex.

An attachment face of a lower star belongs to an earlier lower star and is
treated as already inserted. Lower stars are processed in increasing maximum-
vertex order, so the emitted result is a globally filtration-monotone Morse
sequence. Invalid inputs raise an error instead of being silently refined or
tie-broken.

The strategy is available from the C++ and Python APIs but is not part of the
default automatic strategy portfolio, because that portfolio accepts arbitrary
filtered complexes.

The parallel variant has canonical name:

```python
"process-lower-stars-parallel"
```

It computes disjoint lower-star kernels concurrently with the package's bounded
executor. Lower stars are assigned greedily from largest to smallest, using
their simplex cardinality as a deterministic work estimate, so workers receive
more even loads than with equal vertex-count chunks. Profile metrics report the
number and maximum size of lower stars plus the minimum and maximum scheduled
task loads. Local events are then replayed in increasing maximum-vertex order.
Consequently it produces exactly the same deterministic sequence as
`process-lower-stars`. In Python, `max_workers` limits the complete worker
budget, including the calling thread. The pure-Python fallback accepts the same
name and preserves the result but currently executes the kernels sequentially.

## F-Min

Canonical name:

```python
"f-min"
```

`f-min` is the decreasing dual of `f-max`. It scans from the high end of the
filtration, often from top-dimensional simplexes. When no same-level
reduction-like pair is available, the next not-yet-removed simplex becomes a
critical seed. The builder then removes same-level reduction-like pairs when a
lower simplex has a unique remaining coface. The decreasing events are finally
reversed to produce the increasing Morse sequence used by the rest of the
pipeline.

Informally:

```text
work from high to low,
prefer same-level reduction-like pairs,
reverse the events into an increasing sequence.
```

This is the implementation corresponding to the `Min(S,F)` style used in our
experiments. It is a valid `F`-sequence constructor, but it is not necessarily a
flooding construction.

## Same-Level Reduction

Canonical name:

```python
"same-level-reduction"
```

This strategy processes each filtration plateau independently. Within one level
it collapses same-level free-face pairs until no same-level free face remains.
The remaining active simplexes in that level are then marked critical. Finally,
the collapse pairs are emitted in reverse order so that the resulting Morse
sequence is an increasing sequence. Because it exhausts one filtration value
before moving to the next, it is a flooding construction.

Informally:

```text
collapse until no free face,
mark all leftovers critical,
emit collapsed pairs in reverse order.
```

This strategy is useful as a plateau-local reduction baseline. It is also the
strategy behind the older `"coreduction"` alias, but that alias is deliberately
not used in new documentation because it is ambiguous.

## Plateau-Greedy

Canonical name:

```python
"plateau-greedy"
```

`plateau-greedy` follows the same flooding recurrence as `"saturated"`, but it
uses a more strategic choice when it is stuck and must mark a fillable simplex
critical. Among the fillable candidates, it scores each candidate by how many
currently almost-pairable same-level cofaces it would unlock.

Informally:

```text
pair until stuck,
mark one strategically chosen fillable simplex critical,
continue.
```

This is experimental. It is meant to test whether local plateau information can
reduce the number of critical simplexes or the later reducer work.

## Flooding Reduction Kernel

Canonical name:

```python
"flooding-reduction-kernel"
```

This strategy is the sequential reference implementation of the
reduction-kernel algorithm. It processes one filtration section at a time. In
each kernel round it computes the current facets, protects the core of every
facet cell, greedily reduces each cell to a deterministic local kernel, and
merges the independent local reductions. Kernel rounds repeat until stable;
only then is one current facet perforated. Decreasing events are reversed to
produce the increasing flooding sequence used by the rest of MorseFrames.

The implementation uses Proposition 1's attachment characterization: a face of
a facet cell belongs to its protected core exactly when it is contained in a
different current facet. When several local kernels are possible, candidates
are selected in the level bucket's dimension/lexicographic order.

This version is intentionally sequential and serves as the deterministic
reference for every parallel execution policy.

The mutable active-set computation lives in `ReductionKernelWorkspace`, behind
the read-only `ComplexView` interface. Frame profiles report the number of
processed levels, kernel rounds, facet kernels, reductions, and perforations,
as well as time spent finding facets, computing facet incidence, constructing
cells, computing local reductions, aggregating facet results, and merging a
round.

Round merging visits only the simplices named by accepted reduction events.
Their conflict markers are reset as they are removed, avoiding bucket-wide
marker clearing and removal scans while preserving deterministic event order.
Sequential facet discovery compacts an ordered active-simplex list in place.
That list also bounds facet-incidence reset and uncached low-dimensional cell
construction, so later rounds do not revisit already removed simplices.
The ordinary metrics-free sequential path consumes facet results immediately,
keeping only their compact event stream for conflict checking and merging.
Diagnostic and parallel executions retain explicit facet results for phase
accounting and concurrent result collection.
The metrics-free result type omits all diagnostic fields at compile time;
parallel metrics-free collection therefore also moves smaller result objects.

For repeated sequential gradients on an owning `FilteredComplex`, callers may
invoke `complex_.prepare_reduction_kernel_cache()` once. ReductionKernel then
reads immutable, precomputed same-level closure ranges and coboundary adjacency
instead of reconstructing or filtering them for every sequence. Cache
construction and memory are reported by that explicit call and are therefore
excluded from subsequent gradient timings. The generic `ComplexView` contract
is unchanged, and multiworker execution continues to use worker-local topology.

The companion strategy `"flooding-reduction-kernel-parallel"` uses the same
workspace and local-kernel routine. Facet cells in a round are evaluated in
bounded batches against one immutable active-set snapshot. A reusable task pool
is shared by level and facet work; waiting tasks cooperatively execute queued
work, allowing nested parallelism without deadlock or repeated thread creation.
The current facets are discovered in parallel in level-bucket chunks. Facet
incidence for each active face is then computed in parallel once per round and
saturated at two; incidence greater than one identifies the protected core
directly, without rescanning every other facet inside every local kernel. Facet
results are combined by an order-preserving binary tree.

The configured worker count is a strict global budget. Independent level
streams are assembled in increasing level order, and every tree node appends
its right event stream after its left stream, so the parallel and sequential
strategies produce the same deterministic Morse sequence while implementing
the parallel phases of Algorithms 1 and 2. The pure-Python fallback preserves
the same semantics sequentially.

The native worker budget defaults to the available hardware concurrency and can
be fixed for reproducible runs:

```python
sequence = compute_morse_sequence(
    complex_,
    algorithm="flooding-reduction-kernel-parallel",
    max_workers=4,
)
```

## Flooding Variants

Canonical names:

```python
"flooding-max"
"flooding-min"
"flooding-minmax"
"flooding-maxmin"
```

The named flooding strategies process one filtration level at a time while
maintaining both active boundary and active coboundary counts. They remove
simplexes from the active plateau using four local operations:

- **coreduction:** remove a pair from the maximal/increasing side;
- **coperforation:** mark a maximal-side critical simplex;
- **reduction:** remove a pair from the minimal/decreasing side;
- **perforation:** mark a minimal-side critical simplex.

The variants differ only in priority:

```text
flooding-max:    coreduction, then coperforation
flooding-min:    reduction, then perforation
flooding-minmax: reduction, coreduction, perforation, coperforation
flooding-maxmin: coreduction, reduction, coperforation, perforation
```

The `minmax` and `maxmin` variants are intentionally different. The former
starts from the minimal/decreasing side, while the latter starts from the
maximal/increasing side.

Flooding is experimental in MorseFrames. It is useful for comparing local
plateau behavior, but it is not currently the default strategy.

## Relation To Coreduction Terminology

The word "coreduction" is used in several nearby literatures, and it can be a
source of confusion.

In a Mrozek-Batko-style coreduction primitive, the computation is naturally
viewed as a deletion process: a coface can be removed together with a unique
remaining boundary face. This is closest in spirit to decreasing or active-set
operations such as `f-min` and the reduction/coreduction sides of the flooding
strategies.

The MorseFrames strategy named `"same-level-reduction"` is not meant to claim
identity with every published coreduction algorithm. It performs plateau-local
free-face collapses and then reverses the collapse order into an increasing
Morse sequence. Older code accepted `"coreduction"` as an alias for this
strategy, but new code should prefer `"same-level-reduction"`.

## Relation To CAM

MorseFrames is not an implementation of CAM. CAM is a useful comparison point
because it also exploits coreduction-style simplifications before or during
persistence computations. In MorseFrames, the main object is instead a Morse
sequence and the associated reference/coreference maps. The persistence reducer
then works on the critical frame rather than on the full filtration matrix.

This distinction matters for benchmarking: comparisons with GUDHI or CAM-like
methods should state which parts of the pipeline are included, for example
sequence construction, reference/coreference construction, reducer time, and
total end-to-end time.

## Relation To The GUDHI Adapter

The experimental GUDHI adapter does not define new mathematics. It exposes the
same MorseFrames sequence/reference/persistence pipeline on a GUDHI
`Simplex_tree` view. Its purpose is to test whether the interface can be made
compatible with GUDHI-style data structures and, later, whether a smaller
upstream contribution is realistic.

For now, the native Python and C++ MorseFrames APIs are the reference
implementation. The GUDHI-facing API should be treated as experimental.
