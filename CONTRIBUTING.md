# Contributing

MorseFrames is currently research software. Contributions are welcome, but the
API and repository layout may still evolve.

Before opening a substantial pull request, please open an issue describing the
change. This is especially useful for changes to Morse sequence strategies,
coefficient fields, GUDHI integration, or benchmark methodology.

## Local Checks

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
PYTHONPATH=python python3 -m pytest python/tests
```

Optional GUDHI adapter checks require C++ GUDHI and Boost headers.

