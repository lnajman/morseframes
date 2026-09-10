# C++ complex-view API

The C++ Morse kernels are written against a lightweight `ComplexView` API rather
than against one storage class. This lets the same sequence, reference-map, and
Morse-persistence code run on the compact prototype complex or directly on a
GUDHI-style `Simplex_tree` view.

## Required operations

A complex view is a const-readable object whose simplex ids are contiguous in
`[0, size())`. It must expose:

```cpp
std::size_t size() const;
std::uint16_t dimension(SimplexId simplex) const;
LevelId level(SimplexId simplex) const;
double filtration(SimplexId simplex) const;

const std::vector<VertexId>& vertices(SimplexId simplex) const;
const std::vector<SimplexId>& boundary(SimplexId simplex) const;
const std::vector<SimplexId>& coboundary(SimplexId simplex) const;

const std::vector<SimplexId>& filtration_order() const;
const std::vector<SimplexId>& simplices_of_level(LevelId level) const;
const std::vector<double>& level_values() const;
std::size_t num_levels() const;
```

The helper trait `morseframes::is_complex_view_v<T>` checks this interface at
compile time. The sequence and reference-persistence templates also use
`static_assert` guards so missing methods fail with an explicit message.

## Semantic conventions

The `vertices(simplex)` vector is sorted and canonical. Boundaries and
coboundaries are expressed with the same local simplex ids. Filtration is
monotone on faces, and `level(simplex)` indexes the sorted unique values returned
by `level_values()`.

Tie-breaking inside sequence strategies is representation-independent: when two
candidates have the same algorithmic score, they are ordered by filtration level,
dimension, then vertex tuple. This matters for comparing a direct GUDHI view with
the compact imported complex.

## Public headers

Use `gudhi/Morse_persistence.h` for the GUDHI-shaped public API. It exposes
`Gudhi::morse_persistence::compute_morse_persistence`,
`compute_morse_sequence_and_reference_map`, `Simplex_tree_view`, and
`Morse_sequence_strategy`.

The lower-level MorseFrames headers remain available under `morseframes/...`:
`morseframes/morse_reference_api.hpp` contains the generic `ComplexView` entry
points, and `morseframes/simplex_tree_morse.hpp` contains the direct
Simplex-tree adapter used by the public wrapper.

## Owning-complex construction diagnostics

`FilteredSimplicialComplex::finalize()` preserves lexicographic simplex IDs,
filtration ordering, boundary order and coboundary order. The finalized simplex
array itself is the sorted lookup index: a compact array records the start of
each first-vertex range, and lookup binary-searches the range's vertex key and
then its simplex records. No second tree or duplicate vertex-key storage is
needed. The range keys are sparse (not an array indexed by vertex ID), and
range offsets remain valid when the complex is copied or moved. Const lookup
does not lazily mutate the index.

Finalization also reserves simplex/boundary storage, reuses a temporary face
buffer and reads each pending filtration value directly. These shared changes
apply to every strategy using this owning complex and to arbitrary monotone
filtrations; direct external complex views are unaffected. Pending insertion
still uses the ordered map, and newly inserted simplices become findable only
after finalization, as before. Boundary lookup uses the same compact index and
still checks missing faces and filtration monotonicity.

For a separate diagnostic run:

```cpp
morseframes::ComplexConstructionMetrics metrics;
complex.finalize_with_metrics(metrics);
```

The fields (seconds) are `reset_seconds`, `index_and_simplices_seconds`,
`levels_seconds`, `boundaries_seconds`, `coboundaries_seconds` and
`orders_and_buckets_seconds`. Boundary construction includes face-closure and
filtration-monotonicity checks. Metrics are reset for each call. Input insertion
is not included: measure it separately. Ordinary `finalize()` compiles the same
implementation without internal clocks. Do not mix diagnostic samples with
uninstrumented performance measurements.

`tools/benchmark_complex_construction.py` compares identical resident input
arrays against two header snapshots, checks exact reference complexes and
gradients, and measures construction, gradient execution and fresh-process peak
memory separately. Snapshots without the bulk header use the legacy per-face adapter: its
per-insertion diagnostic clocks perturb enumeration, whose residual includes
clock overhead. Current headers use the bulk path below with separate named
phase diagnostics. Neither diagnostic mode contributes performance samples.

## Bulk construction from a vertex function

For a complex supplied as cells and a function on vertices:

```cpp
#include <morseframes/lower_star_complex.hpp>

morseframes::FilteredSimplicialComplex complex;
std::vector<double> values{0, 1, 2, 3};
std::vector<std::vector<morseframes::VertexId>> cells{{2, 0, 1}, {3}};
morseframes::add_lower_star_cells(complex, values, cells);
complex.finalize();
```

The shared C++ helper inserts every nonempty face of the supplied cells, using
the maximum vertex value as its filtration. IDs index `values`; isolated
vertices must appear as singleton cells. Mixed-dimensional, duplicated and
unordered cells are supported. No dimension-three restriction is imposed:
faces of cardinality one through four use compact fixed-size temporary records,
and larger faces use variable-length records. One dimension is enumerated,
sorted and deduplicated at a time, then its unique faces are merged into the
ordered pending map. Temporary memory scales with the number of generated
faces in the largest dimensional batch, not only with the unique faces.

The helper does not finalize the complex or run any gradient algorithm. All
gradient strategies can consume the resulting complex; RK itself need not use
a max-vertex filtration. This is an opt-in C++ path, adopted by the resident
construction/gradient benchmarks, not an automatic Python-constructor change.

Existing entries retain their first filtration value when the new value agrees
within the existing `1e-12` duplicate tolerance; conflicting values throw.
Max-vertex ties preserve the first input occurrence, including the sign of zero,
as in the legacy per-cell enumeration. Empty cells, repeated vertices in a cell,
missing vertex values and NaN values are rejected before insertion; infinite
values are allowed. An empty cell collection adds nothing. Insertion/allocation
errors can leave partial additions, just like a sequence of `add_simplex`
calls. Finalization still performs closure and monotonicity checks. Input
arrays are not modified.

For separate diagnostics, pass a `LowerStarConstructionMetrics` object to
`add_lower_star_cells_with_metrics`. It is reset on entry and reports
`validation_seconds` (including canonical cell copying), `enumeration_seconds`,
`sort_and_dedup_seconds`, `insertion_seconds` (including temporary cleanup),
`generated_faces` and `unique_faces_submitted`. The last counter includes faces
already in the pending map. The ordinary helper has no internal clock reads.

## Lightweight ReductionKernel initialization

For RK-only sequence construction, include
`morseframes/reduction_kernel_sequence.hpp`:

```cpp
#include <morseframes/reduction_kernel_sequence.hpp>

morseframes::ReductionKernelSequenceBuilder builder(complex_view);
auto sequence = builder.build_flooding_reduction_kernel_parallel(8);
```

This builder exposes the sequential and parallel RK methods, their step-callback
variants, and the execution-options variant. It reuses the existing RK
implementation; it does not introduce another gradient algorithm. The view must
remain valid and unchanged while the builder is used. Optional sequence metrics
and the detailed/coarse profiling flag have the same meaning as on
`FSequenceBuilder`.

Initialization still validates that `filtration_order()` is a permutation of
`[0, size())`, rejecting incorrect lengths, invalid IDs, and duplicates. It uses
a temporary byte per simplex for that check, then releases it. It does not build
the general-purpose rank, level, and dimension caches: RK reads the view's level
buckets and metadata directly. On a 64-bit build, this avoids retaining 14 bytes
per simplex of cache payload (excluding vector/allocator overhead). This is not
a claim about the peak memory of the entire algorithm.

The shared RK reference-frame and compact persistence-input entry points, and
the native Python RK sequence, profiling, and persistence paths use this builder.
`FSequenceBuilder` retains its existing eager initialization and all strategy
methods, including its backward-compatible RK methods. Callers that construct
`FSequenceBuilder` directly must opt into the new RK-only type to obtain the
initialization saving. Other strategies, including F-Max, are unchanged.

No lazy mutable cache is introduced. Repeated const builds use independent
workspace state; concurrent calls require a thread-safe immutable view and no
shared writable metrics object, as with the existing kernels.

## GUDHI-style Simplex_tree entry point

Include `gudhi/Morse_persistence.h` to use the direct Simplex-tree path:

```cpp
#include <gudhi/Simplex_tree.h>
#include <gudhi/Morse_persistence.h>

namespace mp = Gudhi::morse_persistence;

Gudhi::Simplex_tree<> st;
// Insert simplices, assign filtrations, then initialize the filtration cache.
st.initialize_filtration();

auto result = mp::compute_morse_persistence(
    st,
    mp::Morse_sequence_strategy::F_MAX);

for (const auto& pair : result.off_diagonal_intervals()) {
  auto birth_handle = result.simplex_tree_handle(pair.birth);
  auto death_handle = result.simplex_tree_handle(pair.death);
  (void)birth_handle;
  (void)death_handle;
}
```

The result owns the `Simplex_tree_view` so local simplex ids in the returned
sequence and diagram can still be mapped back to Simplex-tree handles. The
Simplex tree itself must outlive the result because the view stores handles into
that tree.

The file `examples/gudhi_simplex_tree_morse.cpp` is a minimal complete example
for this integration path. The file
`examples/example_morse_persistence_from_simplex_tree.cpp` is written in a
candidate upstream GUDHI example style.
