# Cubical Complexes

Cubical support is experimental and is developed on the
`feature/cubical-complexes` branch. The first backend is
`morseframes::CubicalGrid2DComplex`, a rectangular two-dimensional cubical grid
constructed from vertex filtration values.

```cpp
#include "morseframes/cubical_complex.hpp"
#include "morseframes/morse_sequence.hpp"
#include "morseframes/reference_persistence.hpp"
#include "morseframes/standard_persistence.hpp"

std::vector<double> values = {
    0.0, 1.0, 0.0,
    1.0, 2.0, 1.0,
    0.0, 1.0, 0.0,
};

morseframes::CubicalGrid2DComplex grid(3, 3, values);
auto sequence = morseframes::FSequenceBuilder(grid).build_saturated();
auto diagram = morseframes::compute_morse_reference_persistence(grid, sequence);
```

The constructor receives the number of grid vertices in the `x` and `y`
directions, followed by row-major vertex values. A cell filtration value is the
maximum of the values on its vertices, so the resulting filtration is monotone on
faces and preserves equal-value plateaus.

The view exposes the same structural operations as the simplicial complex view:
`size`, `dimension`, `level`, `filtration`, `vertices`, `boundary`,
`coboundary`, `filtration_order`, `simplices_of_level`, `level_values`, and
`num_levels`. The `vertices` method returns a canonical grid-vertex key for the
cell; for example, a square returns its four corner vertex ids.

For prime fields, the grid also exposes the optional coefficient hook:

```cpp
std::uint32_t boundary_coefficient(SimplexId cell,
                                   std::size_t boundary_index,
                                   std::uint32_t modulus) const;
```

The signs follow the product orientation. Horizontal edges are oriented from
left to right, vertical edges from bottom to top, and a square has boundary
right vertical edge, left vertical edge, top horizontal edge, bottom horizontal
edge with signs `+`, `-`, `-`, `+`.

This first backend is intentionally explicit: boundaries and coboundaries are
precomputed and each local incidence list has fixed inline capacity. The natural
next optimization is an implicit grid view for large images or volumes, but the
explicit view is the right correctness target before that step.
