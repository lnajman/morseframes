# Morse persistence discussion package

This note is the current hand-off map for discussing the prototype and deciding
what to do next.  It separates the research/backend result from the current
GUDHI-integration result, because those are different claims.

## Main artifact

Start with the viewer-friendly report:

```text
docs/experiments_morse_persistence_viewer.pdf
```

The main source is:

```text
docs/experiments_morse_persistence.tex
```

The most important section is now:

```text
Backend and integration comparisons with GUDHI
```

It states the two claims explicitly:

- Backend comparison: the compact C++ Morse backend, currently driven through
  the Python-facing benchmark pipeline, is faster than GUDHI on the selected
  Roadmap averages.
- Integration comparison: the current native GUDHI-facing path starts from
  `Gudhi::Simplex_tree`.  After inline local storage for small vertex,
  boundary, coboundary, and annotation lists, dense low-dimensional face lookup
  on compact vertex handles, rank-based Morse candidate queues,
  fixed-dimension native-view ordering, and the byte-valued hot state used
  during sequence/reference construction, plus the specialized `Z_2` reducer
  cleanup, it is faster than GUDHI in the repeated
  large `f-min` validation.  The quick multi-strategy integration table remains
  mixed, but the grid `f-max` and `f-min` rows now cross parity.

The compact overhead diagnosis is:

```text
Read-only view construction: 39.0%
Morse frame construction:   33.6%
Reducer kernel:             27.2%
All pre-reducer work:       72.6%
```

This is the present explanation for why the backend result and the integration
result still differ by family and strategy: the repeated large `f-min`
validation now crosses parity on all listed cases, but quick strategy rows still
include several below-parity measurements because most of the native path is
spent before the reducer starts.

## Reproducible checks

The helper script:

```text
tools/reproduce_discussion_package.py
```

runs the smoke checks used for discussion:

```sh
python3 tools/reproduce_discussion_package.py
```

By default it:

- configures the C++ build directory if needed;
- builds the C++ targets;
- runs the C++ tests;
- runs the prime-field Python tutorial over `F_3`.

To also rerun the native GUDHI quick benchmark without overwriting the canonical
report CSVs:

```sh
python3 tools/reproduce_discussion_package.py --with-native-gudhi-benchmark
```

This writes reproduction artifacts under `docs/reproduction_native_gudhi_*`.
Those files are for sanity checking.  The canonical report tables remain:

```text
docs/native_gudhi_view_quick.csv
docs/native_gudhi_stage_profile_quick.csv
docs/native_gudhi_large_fmin_repeats.csv
```

Short native benchmark runs are timing-sensitive; use the canonical CSVs above
for the report claims, and use reproduction runs mainly to check that the local
pipeline still executes.

## Python interface demonstration

Prime-field persistence is demonstrated in:

```text
docs/python_prime_field_tutorial.md
python/examples/prime_field_tutorial.py
```

Typical commands:

```sh
python3 python/examples/prime_field_tutorial.py --modulus 3
python3 python/examples/prime_field_tutorial.py --modulus 5 --algorithm f-max
```

The tutorial checks that Morse persistence, Morse coreference persistence, and
ordinary full-complex persistence agree over the chosen prime field.

## GUDHI integration notes

The GUDHI-facing design notes are:

```text
docs/gudhi_contribution_design_note.md
docs/gudhi_upstream_patch_map.md
```

The current recommendation for a first upstream discussion is deliberately
small:

- `Z2` first for the GUDHI-facing API;
- direct plateau handling, no lower-star refinement requirement;
- stable strategy set: `same-level`, `f-max`, `f-min`, and `plateau-greedy`;
- keep flooding and prime-field API as prototype features until the first
  integration question is settled.

## What not to claim yet

Do not claim that the native GUDHI plugin is universally faster end-to-end.  The
repeated large `f-min` validation is now faster on all listed cases
(mean `GUDHI/Direct` between `1.08` and `1.60`), with individual flag/clique
samples still close to parity.  The quick multi-strategy table still ranges
from `0.73` to `1.33`, so broader repeated validation is needed before making a
strong integration claim.

Do claim that the compact Morse backend is often faster than GUDHI on the
selected Roadmap benchmark averages, and that the current GUDHI-facing path
identifies a concrete engineering bottleneck before the reducer.

## Recommended next discussion

The next scientific/engineering question is:

```text
Should we keep optimizing the GUDHI-facing view construction and broadening the
native integration benchmark, or should the next research step be new Morse
sequence algorithms and broader benchmark families?
```

My recommendation is now to keep both tracks open.  The GUDHI-facing path has
become competitive enough to justify another focused pass on native integration
costs, but for the paper direction the more informative next work is still a
clearer strategy comparison over a larger benchmark family.
