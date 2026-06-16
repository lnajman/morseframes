# Python API Guide

This page summarizes the current public Python interface. The package name is
`morseframes`, usually imported as:

```python
import morseframes as mf
```

The Python package uses the native C++ backend when `_morse_core` is available
and falls back to a pure-Python implementation otherwise.

```python
print(mf.__version__)
print(mf.cpp_backend_available())
```

To force the pure-Python fallback:

```sh
MORSEFRAMES_DISABLE_CPP_BACKEND=1 python3 script.py
```

## Complexes

The main input type is `FilteredComplex`, a finite filtered abstract simplicial
complex with contiguous simplex ids.

```python
complex_ = mf.FilteredComplex.from_simplices([
    ([0], 0.0),
    ([1], 0.0),
    ([0, 1], 1.0),
])
```

Useful constructors:

```python
mf.FilteredComplex.from_simplices(items)
mf.FilteredComplex.from_facets(facets, filtration=...)
mf.FilteredComplex.from_lower_star(facets, vertex_values, dimension_offset=0.0)
mf.FilteredComplex.from_graph(num_vertices, edges, vertex_filtration=..., edge_filtration=...)
mf.FilteredComplex.from_rips_distance_matrix(distances, threshold=..., max_dimension=...)
mf.FilteredComplex.from_gudhi_simplex_tree(simplex_tree)
```

Useful queries:

```python
complex_.size
complex_.num_levels
complex_.dimension(simplex_id)
complex_.filtration(simplex_id)
complex_.vertices(simplex_id)
complex_.simplex_id([0, 1])
[0, 1] in complex_
complex_.boundary(simplex_id)
complex_.coboundary(simplex_id)
complex_.boundary_simplices([0, 1, 2])
complex_.coboundary_simplices([0, 1])
```

For incremental construction, use `SimplexTreeBuilder`:

```python
builder = mf.SimplexTreeBuilder()
builder.insert([0, 1, 2], 2.0, include_faces=True)
complex_ = builder.finalize()
```

## Sequence Strategies

Most high-level functions accept an `algorithm` keyword. The currently exposed
strategies are:

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

## Morse Sequences

Compute a Morse sequence:

```python
sequence = mf.compute_morse_sequence(complex_, algorithm="f-max")
```

Important fields and helpers:

```python
sequence.steps
sequence.critical_simplices
sequence.paired_with
sequence.critical_index(simplex_id)
sequence.steps_as_simplices(complex_)
sequence.critical_simplices_as_simplices(complex_)
```

Each step is either `"critical"` or `"regular_pair"`.

## Reference And Coreference Maps

Reference maps and coreference maps are annotation tables indexed by simplex id.
An annotation is a tuple of critical labels. Use the simplex helpers to convert
labels back to simplex vertex tuples.

```python
references = mf.compute_reference_map(complex_, sequence)
coreferences = mf.compute_coreference_map(complex_, sequence)

reference_simplices = mf.reference_map_as_simplices(complex_, sequence, references)
coreference_simplices = mf.coreference_map_as_simplices(complex_, sequence, coreferences)
```

The sequence and reference map can be computed in one pass:

```python
frame = mf.compute_morse_sequence_and_reference_map(complex_, algorithm="f-max")
sequence = frame.sequence
references = frame.references
```

The dual coreference frame is also available:

```python
dual_frame = mf.compute_morse_sequence_and_coreference_map(
    complex_,
    algorithm="same-level-reduction",
)
coreferences = dual_frame.coreferences
```

## Morse Complexes

The Morse complex is represented on critical simplex ids.

```python
morse_complex = mf.compute_morse_complex(complex_, sequence, references)
```

Useful fields and helpers:

```python
morse_complex.critical_simplices
morse_complex.boundary(critical_id)
morse_complex.boundary_as_simplices(complex_, critical_id)
```

The reference and coreference complexes can also be computed explicitly:

```python
reference_complex = mf.compute_reference_complex(complex_, sequence, references)
coreference_complex = mf.compute_coreference_complex(complex_, sequence, coreferences)
```

## Persistence

The main entry point is:

```python
diagram = mf.compute_morse_persistence(
    complex_,
    sequence=sequence,
    references=references,
)
```

If no sequence or reference map is supplied, they are computed internally:

```python
diagram = mf.compute_morse_persistence(complex_, algorithm="f-min")
```

The result is a `PersistenceDiagram`:

```python
diagram.finite_pairs
diagram.essential
diagram.finite_barcode()
diagram.finite_barcode(include_zero=True)
diagram.essential_barcode()
diagram.intervals_as_simplices(complex_)
```

Ordinary full-complex persistence is available as a validation baseline:

```python
standard = mf.compute_standard_persistence(complex_)
assert diagram.finite_barcode() == standard.finite_barcode()
assert diagram.essential_barcode() == standard.essential_barcode()
```

## Prime Fields

Persistence over a prime field `F_p` is available by passing `modulus=p`.
Composite moduli are rejected.

```python
diagram_f3 = mf.compute_morse_persistence(
    complex_,
    algorithm="f-max",
    modulus=3,
)

standard_f3 = mf.compute_standard_persistence(complex_, modulus=3)
```

Explicit signed mod-p reference and coreference maps are also exposed:

```python
references_f3 = mf.compute_reference_map_modp(complex_, sequence, modulus=3)
diagram_f3 = mf.compute_morse_persistence_modp(
    complex_,
    sequence=sequence,
    references=references_f3,
    modulus=3,
)
```

## Profiling And Selection

The profiling API builds the Morse sequence and compact reference-reduction
input without running the full pivot reduction.

```python
profile = mf.profile_morse_reference_frame(complex_, algorithm="f-max")
print(profile.critical_ratio)
print(profile.estimated_reducer_work)
```

To compare strategies:

```python
profiles = mf.profile_morse_sequence_algorithms(
    complex_,
    algorithms=("saturated", "f-max", "f-min", "same-level-reduction"),
)

best = mf.select_morse_sequence_profile(
    complex_,
    selection_metric="profile_total_work",
)
```

For a measured, timing-based comparison:

```python
benchmarks = mf.benchmark_sequence_algorithms(
    complex_,
    algorithms=("saturated", "f-max", "f-min"),
    repeats=3,
)
```

The adaptive entry point can choose between Morse persistence and ordinary
persistence:

```python
result = mf.compute_persistence_adaptive(
    complex_,
    sequence_algorithm="auto",
)
print(result.method)
print(result.finite_barcode())
```

## Plotting Helpers

Plotly helpers live in `morseframes.plotting`.

```python
from morseframes import plotting

field = plotting.build_noisy_sine_square(size=12, noise=0.05, seed=0)
fig = plotting.plot_morse_surface_with_persistence(field, algorithm="f-max")
fig.show()
```

Plotly is optional and can be installed with:

```sh
python3 -m pip install -e ".[plotting]"
```

