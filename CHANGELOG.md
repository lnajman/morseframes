# Changelog

All notable changes to MorseFrames are recorded here.

## 0.1.0-alpha.1 - 2026-06-20

Initial public alpha release.

### Added

- Header-only C++ core for filtered simplicial complexes.
- Python interface with pure-Python fallback and optional nanobind backend.
- Morse sequence strategies: saturated, F-Min, F-Max, same-level reduction,
  plateau-greedy, and flooding variants.
- Reference and coreference map construction.
- Morse-reference and Morse-coreference persistence over `Z_2` and prime
  fields `F_p`.
- Simplex-tree-like builder for constructing complexes.
- Experimental GUDHI-facing adapter for `Simplex_tree`.
- Plotly demonstration for inspecting Morse critical data on a triangulated
  square.
- Reproducible benchmark scripts and generated public benchmark table
  fragments.
- Read the Docs documentation, including overview, quickstart, minimal
  examples, and Morse-data inspection guide.
- GitHub Actions CI for C++, Python, and documentation builds.

### Notes

- This is research software. Public APIs are useful for experiments, but may
  still evolve while the paper and GUDHI integration mature.
