# Cubical Complexes

Cubical support is experimental and is developed on the
`feature/cubical-complexes` branch. The first backend is
`morseframes::CubicalGrid2DComplex`, a rectangular two-dimensional cubical grid
constructed from vertex filtration values.

## Python

The Python interface exposes the same backend through
`morseframes.CubicalGrid2DComplex` when the native extension is available.

Run the complete tutorial from the repository root with:

```sh
PYTHONPATH=python python3 python/examples/cubical_grid_tutorial.py \
  --width 8 \
  --height 8 \
  --algorithm f-max \
  --modulus 3
```

When the optional `gudhi` Python package is importable, the tutorial also builds
the same lower-star grid with `gudhi.CubicalComplex`, checks that the barcodes
agree, and prints a GUDHI cubical timing. Pass `--skip-gudhi` to keep the run
limited to MorseFrames.

```python
import morseframes as mf

grid = mf.CubicalGrid2DComplex.from_vertex_values(
    3,
    3,
    [
        0.0, 1.0, 0.0,
        1.0, 2.0, 1.0,
        0.0, 1.0, 0.0,
    ],
)

sequence = mf.compute_morse_sequence(grid, algorithm="f-max")
diagram = mf.compute_morse_persistence(grid, sequence, modulus=3)
standard = mf.compute_standard_persistence(grid, modulus=3)
gudhi = mf.gudhi_cubical_barcode(grid, modulus=3)

assert diagram.finite_barcode() == standard.finite_barcode()
assert gudhi[0] == standard.finite_barcode()
```

The object has the same inspection helpers as the simplicial Python complex:
`size`, `level_values`, `filtration_order`, `simplex_records()`, `vertices`,
`boundary`, `coboundary`, `filtration`, and `simplices_of_level`. It also adds
grid-coordinate helpers: `vertex(x, y)`, `horizontal_edge(x, y)`,
`vertical_edge(x, y)`, and `square(x, y)`. The `cell_type(cell)` method returns
`"vertex"`, `"horizontal_edge"`, `"vertical_edge"`, or `"square"`.

## C++

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
