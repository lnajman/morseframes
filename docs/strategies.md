# Morse Sequence Strategies

MorseFrames separates the construction of a Morse sequence from the later
reference/coreference and persistence computations. This is intentional: the
same public pipeline can be used with several sequence constructors, and new
constructors should be easy to add.

All strategies exposed here produce a simplex-wise Morse sequence compatible
with the input filtration. Regular pairs are restricted to a single filtration
level, so the methods work directly on complexes with plateaus. No lower-star
refinement is required by the API.

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
"f-min"
"same-level-reduction"
"plateau-greedy"
"flooding-max"
"flooding-min"
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
