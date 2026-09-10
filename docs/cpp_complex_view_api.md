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
