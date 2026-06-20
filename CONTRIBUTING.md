# Contributing

MorseFrames is research software. Contributions are welcome, especially small
tests, documentation fixes, reproducibility improvements, and focused bug
reports.

## Development Setup

From the repository root, build and test the C++ core:

```sh
cmake -S . -B build -DMORSEFRAMES_BUILD_GUDHI_TOOLS=OFF
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Install the Python package in editable mode:

```sh
python3 -m pip install -e ".[dev]"
python3 -c "import morseframes as mf; print(mf.__version__, mf.cpp_backend_available())"
```

Run the pure-Python fallback tests:

```sh
MORSEFRAMES_DISABLE_CPP_BACKEND=1 python3 -m unittest discover -s python/tests -p "test_*.py"
```

Run the quickstart example:

```sh
python3 python/examples/quickstart.py
```

## Documentation

Install the documentation dependencies and build the Sphinx site:

```sh
python3 -m pip install -r docs/requirements.txt
MORSEFRAMES_DISABLE_CPP_BACKEND=1 python3 -m sphinx -W --keep-going -b html docs docs/_build/html
```

The public documentation is hosted at:

<https://morseframes.readthedocs.io/>

## Optional GUDHI Tools

The GUDHI-facing examples, tests, and benchmarks are enabled when GUDHI and
Boost headers are available:

```sh
cmake -S . -B build \
  -DMORSEFRAMES_GUDHI_INCLUDE_DIR=/path/to/gudhi/include \
  -DMORSEFRAMES_BOOST_INCLUDE_DIR=/path/to/boost/include
```

## Benchmarks

Benchmark and table-generation scripts live in `tools/`. The public
reproduction instructions are in:

```text
docs/benchmark_reproduction.md
```

Use benchmark results as comparative signals rather than absolute machine
independent timings. When changing algorithms or reducers, prefer reporting
ratios, dataset metadata, and repeated-run summaries.

## Pull Request Checklist

- Keep changes focused.
- Add or update tests when behavior changes.
- Run the C++ tests, Python fallback tests, and documentation build when
  touching the corresponding surface.
- Keep generated scratch files out of commits.
- Do not commit private manuscript drafts or private discussion notes to this
  public repository.
