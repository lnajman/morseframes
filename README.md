# MorseFrames

MorseFrames is an experimental C++ and Python library for working with
Morse sequences, reference and coreference maps, Morse complexes, and
Morse-based persistent homology.

The current codebase contains:

- a header-only C++ core for filtered simplicial complexes;
- several Morse sequence strategies, including saturated, F-Min, F-Max,
  same-level reduction, plateau-greedy, and flooding variants;
- reference and coreference map construction;
- Morse-reference and Morse-coreference persistence;
- ordinary persistence for validation;
- prime-field coefficients `F_p` in the Python interface;
- a Python API with a pure-Python fallback and optional nanobind backend;
- an experimental GUDHI-facing adapter for `Simplex_tree`;
- tests, examples, and reproducible benchmark drivers.

The GUDHI adapter is included to make future upstream integration easier. It
should be considered experimental.

## Build the C++ Tests

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

When the C++ GUDHI and Boost headers are available, CMake also builds optional
GUDHI adapter tests, examples, and benchmarks. They can be supplied explicitly:

```sh
cmake -S . -B build \
  -DMORSEFRAMES_GUDHI_INCLUDE_DIR=/path/to/gudhi/include \
  -DMORSEFRAMES_BOOST_INCLUDE_DIR=/path/to/boost/include
```

To skip the optional GUDHI tools:

```sh
cmake -S . -B build -DMORSEFRAMES_BUILD_GUDHI_TOOLS=OFF
```

## Python Quick Start

The Python package is installable in editable mode. This uses
`scikit-build-core` and `nanobind` to build the optional native backend:

```sh
python3 -m pip install -e .
python3 -c "import morseframes as mf; print(mf.__version__, mf.cpp_backend_available())"
```

For development, install the test dependency as well:

```sh
python3 -m pip install -e ".[dev]"
```

Example:

```python
import morseframes as mf

complex_ = mf.FilteredComplex.from_simplices([
    ([0], 0.0),
    ([1], 0.0),
    ([0, 1], 1.0),
])

sequence = mf.compute_morse_sequence(complex_, algorithm="f-max")
references = mf.compute_reference_map(complex_, sequence)
diagram = mf.compute_morse_persistence(complex_, sequence, references)

print(diagram.finite_barcode())
print(diagram.essential_barcode())
```

Prime-field coefficients are available by passing a prime modulus:

```python
diagram = mf.compute_morse_persistence(complex_, algorithm="f-min", modulus=3)
```

## Python Tests

```sh
MORSEFRAMES_DISABLE_CPP_BACKEND=1 python3 -m unittest discover -s python/tests -p "test_*.py"
python3 python/examples/quickstart.py
```

The pure-Python fallback can still be exercised directly from the source tree
without building the native backend, or by disabling the backend explicitly:

```sh
PYTHONPATH=python python3 -c "import morseframes as mf; print(mf.cpp_backend_available())"
MORSEFRAMES_DISABLE_CPP_BACKEND=1 python3 -m unittest discover -s python/tests -p "test_*.py"
```

## Documentation

- `docs/python_api.md` summarizes the public Python API.
- `docs/strategies.md` explains the Morse sequence strategies and terminology.
- `docs/python_prime_field_tutorial.md` explains persistence over `F_p`.
- `docs/cpp_complex_view_api.md` describes the C++ complex-view interface.
- `docs/gudhi_contribution_design_note.md` summarizes the experimental GUDHI adapter.
- `docs/benchmark_reproduction.md` explains how to regenerate benchmark tables.
- Generated benchmark table fragments live in `docs/*_table.tex`; the scripts
  that regenerate them live in `tools/`.

## Status

This is research code. The public API is useful for experimentation, but names
and interfaces may still change while the paper and GUDHI integration mature.
