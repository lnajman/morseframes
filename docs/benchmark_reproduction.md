# Benchmark Reproduction

This page explains how to regenerate the public benchmark artifacts in this
repository. It is meant for software reproducibility: manuscript text and
discussion notes live outside the public repository, in the private manuscript
workspace until a public preprint or published version exists.

The current resident-array comparison reports loading and native construction
separately and compares the remaining algorithm work. Complete resident-to-gradient
totals are retained, starting from a common in-memory mesh and vertex function.
The latest [PLS key-storage study](pls_key_arena_benchmark.md) gives commands
for fine profiling, the controlled storage optimization and its fresh TTK rerun.
The [direct optimized-PLS/TTK study](pls_ttk_gradient_benchmark.md) adds PLS to
the resident-array harness and measures all four methods together, with balanced
execution order, a confirmation run and unchanged timing boundaries.
The [simplicial gradient study](simplicial_gradient_benchmark.md) records the
latest PLS workspace optimization and two controlled gradient-only comparisons
of RK, PLS and F-Max through dimension seven, including critical counts and
separate loading/construction times. TTK was not rerun for that study.
The [boundary and ordering reuse study](native_order_reuse_benchmark.md) records
the latest shared-constructor change, repeated construction measurements and
the refreshed RK/F-Max/TTK comparison with timing caveats.
The [compact finalization study](native_finalization_benchmark.md) records the
preceding shared-index change and its construction/memory comparison.
The preceding [bulk face-construction study](native_bulk_construction_benchmark.md) records
the opt-in shared C++ construction path, its controlled comparison,
memory costs and gradient-timing caveats. The earlier
[native construction study](native_construction_benchmark.md) records the
shared-constructor optimization on revision `d7cddc7`, its separate timing
and memory checks, and the refreshed 16-configuration RK/F-Max/TTK comparison.
The `5598963` direct RK/F-Max/TTK tables, including lightweight
RK initialization, remain available below as a pre-construction-optimization
snapshot. The lightweight RK initialization section records the latest
kernel-specific initialization update.
The closure-based incidence, controlled packed-coface A/B, direct TTK comparison,
and earlier phase-profile sections retain historical snapshots; those timings
must not be treated as fresh measurements of subsequent changes.

Run commands from the repository root.

```sh
cd morseframes
```

The examples below write raw CSV, Markdown summaries, and diagnostic prose to
`../work/`. The public repository tracks the scripts and selected LaTeX table
fragments, but not the manuscript prose built from them.

Some rendering commands below intentionally write tracked files under `docs/`.
Use those commands when you want to refresh the public table fragments. If you
only want to test the workflow on a local machine, redirect the table outputs to
`../work/` or restore the tracked table fragments afterward.

## Output Policy

Tracked public artifacts:

- `docs/*_table.tex`: LaTeX table fragments used to report benchmark results.
- `tools/*.py`: benchmark, validation, and table-rendering scripts.
- `benchmarks/benchmark_gudhi_view.cpp`: native GUDHI-view benchmark.
- `benchmarks/benchmark_reduction_kernel_ab.cpp`: native worker for controlled
  comparisons between ReductionKernel revisions.

Local or private artifacts:

- `../work/*.csv`, `../work/*.md`, `../work/*.json`: raw benchmark outputs and
  summaries.
- `docs/*_prose.tex`: generated prose fragments. These are ignored by Git and
  should be copied into the private notes repository only when needed.
- report PDFs and manuscript drafts: private-note material, not public package
  documentation.

## Quick Validation

These checks are the fastest way to confirm that the source tree is usable.

```sh
MORSEFRAMES_DISABLE_CPP_BACKEND=1 \
  python3 -m unittest discover -s python/tests -p "test_*.py"

PYTHONPATH=python python3 python/examples/quickstart.py
PYTHONPATH=python python3 python/examples/prime_field_tutorial.py --modulus 3
```

To include the native C++ backend, install the package in editable mode:

```sh
python3 -m pip install -e ".[dev]"
python3 -c "import morseframes as mf; print(mf.__version__, mf.cpp_backend_available())"
```

The C++ smoke tests are:

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

## Synthetic Morse vs Standard Benchmarks

The main synthetic runner is `tools/benchmark_persistence.py`. It can run one
strategy, or the default strategy portfolio with `--sequence-algorithm
portfolio`.

Small smoke run:

```sh
mkdir -p ../work
PYTHONPATH=python python3 tools/benchmark_persistence.py \
  --preset smoke \
  --sequence-algorithm portfolio \
  --format summary \
  --output ../work/benchmark_smoke_summary.txt
```

Regenerate the public synthetic scale table:

```sh
mkdir -p ../work
PYTHONPATH=python python3 -c "import morseframes as mf; print(mf.cpp_backend_available())"

PYTHONPATH=python python3 tools/benchmark_persistence.py \
  --families lower-star plateau rips \
  --sizes 48 \
  --seeds 0 1 2 \
  --repeats 3 \
  --sequence-algorithm portfolio \
  --validation-mode core \
  --format csv \
  --output ../work/synthetic_scale_size48_portfolio.csv

PYTHONPATH=python python3 tools/render_synthetic_scale_table.py \
  --input ../work/synthetic_scale_size48_portfolio.csv \
  --table-output docs/synthetic_scale_table.tex \
  --prose-output ../work/synthetic_scale_prose.tex
```

The table reports `Std/Morse`, so values above `1` mean the Morse pipeline is
faster than ordinary full-complex persistence for that row.

CSV and JSON rows also retain gradient quality, rather than reporting timing
alone. The fields include the total critical count,
`critical_simplices_by_dimension`, `critical_ratio`, `num_regular_pairs`, and
`sequence_seconds_per_eliminated_simplex`. Timing is split into sequence and
downstream persistence phases, with `total_seconds` recording their complete
end-to-end pipeline. This makes comparisons between strategies meaningful even
when a faster constructor produces a larger Morse complex.

The tracked synthetic table is a native-backed core-mode benchmark. Before
replacing it, make sure the backend check above prints `True`; otherwise the CSV
will contain `cpp_backend=False` rows and the timing will describe the
pure-Python fallback instead of the optimized C++ backend.

## ProcessLowerStars Scaling

The focused ProcessLowerStars runner constructs two controlled simplicial
families. The balanced family gives every anchor vertex the same triangle-fan
workload. The skewed family keeps the same number of simplices but concentrates
most of that work in one lower star. All vertex values are injective and every
simplex uses the exact max-vertex lower-star extension.

```sh
PYTHONPATH=python python3 tools/benchmark_process_lower_stars.py \
  --anchors 16 \
  --balanced-fans 8 32 128 \
  --light-fan 2 \
  --workers 1 2 4 8 \
  --repeats 7 \
  --warmups 2 \
  --format csv \
  --output ../work/process_lower_stars_scaling.csv
```

For each scale, the heavy fan is selected automatically so both families have
the same total fan size and therefore the same simplex count. For example,
`16 * 8 = 98 + 15 * 2`. This isolates load imbalance from input size. Each
parallel row is checked against the sequential algorithm
for the exact Morse-step sequence and against ordinary persistence for the
barcode. The output records estimated task loads, critical counts by dimension,
construction and downstream persistence times, speedup, and parallel
efficiency. A `cpp_backend=False` row is a correctness run of the sequential
fallback, not a parallel-performance measurement.

The current controlled run on an Apple M1 Max (10 CPU cores), using five
measured repetitions after one warm-up, is shown below. At 12,304 simplices the
balanced case reaches 2.26x sequence speedup and 1.68x end-to-end speedup. The
matched skewed case reaches only 1.43x and 1.29x, respectively, while its
estimated eight-worker load ratio rises to 6.79. The two cases have identical
critical counts `(2048, 2032, 0)`, so this gap is attributable to scheduling
imbalance rather than a different Morse complex.

![ProcessLowerStars scaling](process_lower_stars_scaling.svg)

The exact paper-ready values are generated in
`docs/process_lower_stars_scaling_table.tex`. These measurements are a local
scheduler study, not yet the comparison with Robins' implementation; that
external benchmark remains a separate stage.

## Gradient Comparison with Construction Reported Separately

This is the current primary comparison. Native construction is reported separately,
not hidden or optimized by the scheduling changes. The common input is still
resident vertex values and maximal-cell connectivity, and every repetition uses
fresh native objects. Current results use the lightweight RK revision `5598963`;
the earlier `591b658` and `d466174` measurements are retained below as historical data.

The partition is explicit:

| Algorithm | Construction, reported separately | Algorithm time, used for comparison |
| --- | --- | --- |
| F-Max | MorseFrames simplex enumeration, representation/filtration and finalization | Fresh builder setup and F-Max gradient construction |
| RK | MorseFrames simplex enumeration, representation/filtration and finalization | Fresh builder setup, workspace/pool creation, level processing and replay |
| TTK | Native object and cell representation setup, connectivity preconditioning | Vertex ordering and `buildGradient`, including lower-star construction and matching |

MorseFrames assigns the max-vertex extension while constructing its filtered
complex; this remains in its combined construction phase. TTK's function-dependent
vertex ordering is **not** removed from algorithm time. Its lower stars are built
inside `buildGradient` and remain counted there. Structural connectivity preparation
is included in TTK construction, as it is in MorseFrames finalization. This is a
comparison conditional on the listed constructed native representations, not a
claim that the algorithms have identical input structures or output formats.

Each non-profiled performance repetition now records outer phase boundaries.
For each individual sample, **resident-to-gradient total = construction + algorithm**.
TTK's ordering occurs before connectivity preparation in its existing call order,
so its algorithm time sums the ordering and gradient intervals from the same run.
All phases are retained; none are estimated by subtracting separately aggregated
medians. The benchmark uses only four (MorseFrames) or six (TTK) outer clock reads,
without internal gradient profiling. Fine diagnostic samples remain separate and
never enter the comparison. Per-phase performance values and a narrower gradient-only
value remain available alongside the construction-separated primary metric.

File loading/parsing/validation is recorded once per native process and shown
separately as shared input loading, not charged independently to each algorithm.
It includes the adapter's coordinate-array initialization and is not a cold-cache
disk benchmark. Input generation remains outside all recorded intervals. As before,
reference checks, post-readiness destruction and persistence are excluded.

The new native/study schemas are `resident-gradient-v2` and
`resident-gradient-study-v2`. Readers retain historical v1 total/diagnostic support,
but reject v1 evidence for construction-separated comparisons. The default v2
table compares algorithm time; `--comparison total` explicitly selects full totals.
All components, ratios and speedups come from raw non-profiled performance samples.

```sh
LC_ALL=C tools/build_ttk_gradient_benchmark.sh \
  ../work/ttk-benchmark morseframes_resident_gradient_benchmark

OMP_WAIT_POLICY=PASSIVE LC_ALL=C python3 tools/benchmark_resident_gradients.py \
  --benchmark ../work/ttk-benchmark/build-f4ffd1a1049d0ccf6e8f3eb4f7c096a6cc251ba0/morseframes_resident_gradient_benchmark \
  --terrain-sizes 16 64 --volume-sizes 8 16 --seeds 0 2 --workers 1 2 4 8 \
  --repeats 12 --diagnostics 6 --warmups 2 \
  --input-dir ../rk-ab-inputs --output ../rk-resident-5598963-main.json

python3 tools/render_resident_gradients.py \
  --input ../rk-resident-5598963-main.json \
  --table-output docs/resident_gradient_algorithm_table.tex \
  --construction-output docs/resident_gradient_construction_table.tex \
  --phases-output docs/resident_gradient_unprofiled_phases_table.tex
```

The total and construction medians in separate tables need not add to the algorithm
median: the additive identity is checked before aggregation, at the individual
repetition level. F-Max is sequential at every displayed worker setting; it is
remeasured alongside RK and TTK. Raw input/binary/source hashes, revision, wait policy,
orders, exact reference checks and critical counts remain recorded as described below.

### Previous construction-separated results (`5598963`)

The main study reruns the original eight inputs at 1/2/4/8 workers: 32
configurations with 12 non-profiled performance repetitions, six separate
diagnostic repetitions and two warmups per mode. The confirmation reruns terrain
`n=64` and volume `n=16`, both seeds, at 1/8 workers with 18 performance
repetitions, six diagnostics and two warmups. A separate extension measures
volume `n=32`, both seeds, at 1/8 workers with the confirmation settings.
It is not mixed into the original eight-input comparison.

All three studies started from clean revision
`5598963f1a7a975e4ee14552e194729ffc05d64c` on the Apple M1 Max, using the native ARM
release build, pinned TTK revision `f4ffd1a1049d0ccf6e8f3eb4f7c096a6cc251ba0`, and
`OMP_WAIT_POLICY=PASSIVE`. The TTK checkout was clean. Build flags are Apple Clang
15, `-O3 -DNDEBUG -arch arm64`, C++17 and OpenMP for TTK's classic backend. The
driver uses the lightweight RK builder; F-Max remains sequential even in rows
labelled with multiple workers.

Provenance shared by the three raw files:

- Header SHA-256: `03636fd51de667e246fe3c8f9629dad2200b8a38171fa6a3d1d66d6388d3aa7b`.
- Driver SHA-256: `a1ca9b1f76203e376a5b8b8d612ddd9cba0d7aee007d0f12e65aa76ecff0086b`.
- Executable SHA-256: `2c543a11098a788859d53bb945c0b132dacfc0d1faa4477316a1f2a6fbe2c990`.

Selected main-run **algorithm times**, excluding native construction (milliseconds,
medians over the two seed-specific sample medians):

| Input | Workers | Sequential F-Max | RK | TTK |
| --- | ---: | ---: | ---: | ---: |
| 2D terrain, `n=64` | 1 | 1.341 | 1.008 | 1.044 |
| 2D terrain, `n=64` | 8 | 1.361 | 0.670 | 0.626 |
| 3D volume, `n=16` | 1 | 8.839 | 3.774 | 10.626 |
| 3D volume, `n=16` | 8 | 9.397 | 1.799 | 3.212 |

RK has lower paired-median algorithm time than both TTK and F-Max in all 16 main
volume configurations. On the terrains, RK is ahead of TTK in six of 16
configurations, and ahead of F-Max in ten of 16. Small terrains still make
parallel setup costs visible. These are observed rankings on these synthetic
inputs, not claims of universal superiority or statistical significance.

```sh
OMP_WAIT_POLICY=PASSIVE LC_ALL=C python3 tools/benchmark_resident_gradients.py \
  --benchmark ../work/ttk-benchmark/build-f4ffd1a1049d0ccf6e8f3eb4f7c096a6cc251ba0/morseframes_resident_gradient_benchmark \
  --terrain-sizes 64 --volume-sizes 16 --seeds 0 2 --workers 1 8 \
  --repeats 18 --diagnostics 6 --warmups 2 \
  --input-dir ../rk-ab-inputs --output ../rk-resident-5598963-confirmation.json

OMP_WAIT_POLICY=PASSIVE LC_ALL=C python3 tools/benchmark_resident_gradients.py \
  --benchmark ../work/ttk-benchmark/build-f4ffd1a1049d0ccf6e8f3eb4f7c096a6cc251ba0/morseframes_resident_gradient_benchmark \
  --terrain-sizes --volume-sizes 32 --seeds 0 2 --workers 1 8 \
  --repeats 18 --diagnostics 6 --warmups 2 \
  --input-dir ../rk-ab-inputs --output ../rk-resident-5598963-volume32.json
```

Confirmed eight-worker **algorithm times** on the original larger inputs:

| Input | Sequential F-Max (ms) | RK (ms) | TTK (ms) | Paired RK/TTK |
| --- | ---: | ---: | ---: | ---: |
| 2D terrain, `n=64` | 1.323 | 0.703 | 0.649 | 1.044 |
| 3D volume, `n=16` | 9.252 | 1.660 | 2.728 | 0.600 |

The confirmation has RK ahead of both competitors in all four volume
configurations. On terrain, RK is ahead of F-Max in all four configurations and
ahead of TTK in one of four. The eight-worker terrain aggregate retains a small
TTK lead; the eight-worker volume aggregate favors RK. Main/confirmation paired
RK/TTK ratios are respectively 1.080/1.044 on terrain and 0.555/0.600 on volume.
Ratios are formed per repetition before aggregation and need not equal quotients
of the displayed time medians.

Selected main-run **construction and shared loading times**, reported separately:

| Input | Workers | F-Max construction (ms) | RK construction (ms) | TTK construction (ms) | Shared loading (ms) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2D terrain, `n=64` | 8 | 21.799 | 21.860 | 0.523 | 3.582 |
| 3D volume, `n=16` | 8 | 132.161 | 130.463 | 4.676 | 10.757 |

F-Max and RK share the same construction implementation; differences between
their construction columns are measurement variation. Shared loading is one
read/parse/validation per native invocation, not a per-algorithm charge or a
cold-cache disk benchmark. Every raw performance sample satisfies construction
plus algorithm equals its full total. Separately aggregated medians need not add.

The three primary tracked fragments are regenerated from
`rk-resident-5598963-main.json`: `resident_gradient_algorithm_table.tex`,
`resident_gradient_construction_table.tex` and
`resident_gradient_unprofiled_phases_table.tex`. The phase fragment reports all
outer phases for terrain `n=64` and volume `n=16` at 1/8 workers. Raw JSON also
retains fine F-Max diagnostics, RK workspace/pool, level-processing and replay
timings, narrower gradient-only timings and full resident-to-gradient totals.
TTK lower-star construction and matching remain combined inside its gradient
phase; vertex ordering stays in algorithm time. No persistence is timed.

Other desktop applications remained active; our builds and tests did not overlap
the measurements. All samples, including outliers, are retained. Share these
results with the same-machine/synthetic-input caveat. Differences from older
absolute timings are not a controlled estimate of the initialization improvement;
the separate RK-only A/B study below provides that comparison.

#### Larger-volume extension

The separate `rk-resident-5598963-volume32.json` study contains 792,051 simplices
and 32,768 vertices per input. Its algorithm times are:

| Workers | Sequential F-Max (ms) | RK (ms) | TTK (ms) | Paired RK/TTK |
| --- | ---: | ---: | ---: | ---: |
| 1 | 122.039 | 44.182 | 92.643 | 0.464 |
| 8 | 126.266 | 14.148 | 21.359 | 0.658 |

RK has lower paired-median algorithm time than both competitors in all four
seed/worker configurations. This is an additional size, not a repeated-session
confirmation of the `n=32` ranking. Full timings remain conditional on the
specified native representations and output formats.

The extension's construction and loading costs are still reported separately:

| Workers | F-Max construction (ms) | RK construction (ms) | TTK construction (ms) | Shared loading (ms) |
| --- | ---: | ---: | ---: | ---: |
| 1 | 1141.253 | 1134.132 | 56.648 | 92.993 |
| 8 | 1188.479 | 1208.015 | 35.944 | 103.807 |

TTK has lower paired-median full resident-to-gradient time than both MorseFrames
methods in **all 44 configurations**, including this extension. The construction
gap is therefore important for end-to-end use even though it is excluded from
the requested primary algorithm comparison.

To regenerate the extension's complete phase and comparison tables:

```sh
python3 tools/render_resident_gradients.py \
  --input ../rk-resident-5598963-volume32.json \
  --table-output docs/resident_gradient_volume32_algorithm_table.tex \
  --construction-output docs/resident_gradient_volume32_construction_table.tex \
  --phases-output docs/resident_gradient_volume32_unprofiled_phases_table.tex
```

For example, the eight-worker RK builder takes 0.634 ms and its gradient phase
13.398 ms; F-Max's corresponding phases are 6.252 and 120.347 ms. TTK vertex
ordering takes 2.127 ms and combined lower-star construction/matching 19.107 ms.
These are independently aggregated non-profiled phase medians, so their sums
need not equal the aggregate algorithm medians above. The raw sample identities
are exact. Detailed diagnostic clocks remain separate.

#### Critical counts and validation

All 44 configurations pass exact within-algorithm reference checks, including
TTK's cell pairings and MorseFrames' complete sequence fields. Each method agrees
with its own one-worker reference, not necessarily with another method's gradient.
Critical counts agree across all three algorithms and worker counts on every
input, and the original eight inputs retain their previous counts:

| Input | Seed | Critical vertices | Critical edges | Critical triangles | Critical tetrahedra |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2D terrain, `n=16` | 0 | 6 | 6 | 1 | -- |
| 2D terrain, `n=16` | 2 | 8 | 9 | 2 | -- |
| 2D terrain, `n=64` | 0 | 126 | 224 | 99 | -- |
| 2D terrain, `n=64` | 2 | 176 | 334 | 159 | -- |
| 3D volume, `n=8` | 0 | 12 | 21 | 13 | 3 |
| 3D volume, `n=8` | 2 | 16 | 26 | 13 | 2 |
| 3D volume, `n=16` | 0 | 21 | 38 | 29 | 11 |
| 3D volume, `n=16` | 2 | 26 | 54 | 36 | 7 |
| 3D volume, `n=32` | 0 | 128 | 455 | 379 | 51 |
| 3D volume, `n=32` | 2 | 152 | 649 | 570 | 72 |

The three completed raw studies were audited for input/source/header/binary
hashes, input and worker coverage, balanced algorithm orders, repetition counts,
nonnegative phase partitions, critical counts and Euler characteristic. Stored
medians, paired ratios and worker speedups were recomputed from raw samples;
the table renderer also recomputes its values instead of trusting stored summaries.
The local audit companion is `../rk-resident-5598963-validation.py`; the public
reproduction/validation routines remain in `tools/benchmark_resident_gradients.py`
and `tools/render_resident_gradients.py`. The six table fragments were checked
against fresh renders from their respective raw files.

Validation assessment: **share with caveats**. These are preliminary synthetic
grid results on one desktop, with two seeds per size and other applications
active. Outliers are not trimmed. The original larger cases have a separate
repeat, while `n=32` has one study. No universal ranking, higher-dimensional
comparison or full-persistence speedup is established. The native benchmark
integration and all 147 Python tests pass; the fallback suite runs 147 tests with
nine native-only skips. Strict documentation and generated-table checks pass.

### Construction-separated snapshot (`591b658`)

The refreshed main study covers eight inputs at 1/2/4/8 workers: 32 configurations,
12 performance repetitions, six separate diagnostic repetitions, and two warmups
per mode. The confirmation repeats terrain `64` and volume `16`, both seeds, at
1/8 workers with 18 performance repetitions and otherwise identical settings.
Both runs started from a clean `591b65870b4a314164778b89f990b19b2d2fe1d4` checkout,
on the Apple M1 Max with the native ARM release build, pinned TTK revision
`f4ffd1a1049d0ccf6e8f3eb4f7c096a6cc251ba0`, and `OMP_WAIT_POLICY=PASSIVE`.

The checked CMake build uses Apple Clang 15, `-O3 -DNDEBUG -arch arm64`, C++17,
and OpenMP for TTK. The measured executable SHA-256 is
`cba26795d6347f2ecfffe450c5acd44407d20a2577079e9c935efe7cf88a8e26`;
the header digest is
`babf46699673de2aee8dee6c9091893f52c257f62b5a9efc995a09c2a47fd927`.
Both raw JSON files retain matching binary, driver, input and header hashes.
The following main-run **algorithm times** exclude native construction and are
milliseconds, aggregated as medians over the two seed-specific medians:

| Input | Workers | F-Max | RK | TTK |
| --- | ---: | ---: | ---: | ---: |
| 2D terrain, `n=64` | 1 | 1.451 | 1.134 | 1.034 |
| 2D terrain, `n=64` | 8 | 1.465 | 0.894 | 0.693 |
| 3D volume, `n=16` | 1 | 10.407 | 4.663 | 11.264 |
| 3D volume, `n=16` | 8 | 10.894 | 2.105 | 3.472 |

RK has lower paired-median algorithm time than TTK and F-Max in all 16 volume
configurations. On terrains it is faster than TTK in one of 16 configurations,
and faster than F-Max in nine of 16. Small terrains can lose to sequential
F-Max because parallel setup costs are included. This is a workload-dependent
comparison, not a universal ranking.

```sh
OMP_WAIT_POLICY=PASSIVE LC_ALL=C python3 tools/benchmark_resident_gradients.py \
  --benchmark ../work/ttk-benchmark/build-f4ffd1a1049d0ccf6e8f3eb4f7c096a6cc251ba0/morseframes_resident_gradient_benchmark \
  --terrain-sizes 64 --volume-sizes 16 --seeds 0 2 --workers 1 8 \
  --repeats 18 --diagnostics 6 --warmups 2 \
  --input-dir ../rk-ab-inputs --output ../rk-resident-591b658-confirmation.json
```

Confirmed eight-worker algorithm times are:

| Input | Sequential F-Max (ms) | RK (ms) | TTK (ms) | Paired RK/TTK |
| --- | ---: | ---: | ---: | ---: |
| 2D terrain, `n=64` | 1.515 | 0.763 | 0.708 | 1.117 |
| 3D volume, `n=16` | 10.756 | 2.136 | 3.182 | 0.687 |

The confirmation has RK ahead of both competitors in all four volume
configurations, and TTK ahead of RK in all four terrain configurations. RK is
ahead of F-Max in all eight. On the eight-worker volume, paired RK/TTK ratios
are 0.733 in the main run and 0.687 in confirmation; the corresponding terrain
ratios are 1.220 and 1.117. Ratios are computed per repetition before aggregation
and need not equal quotients of the displayed time medians.

Other desktop applications remained active, while our local builds and tests
did not overlap either study. All samples are retained, including timing outliers:
for example, the main seed-0 eight-worker volume has TTK algorithm samples from
2.206 to 37.872 ms. The repeated direction supports a preliminary comparison,
not a quiet-machine performance guarantee or a precisely established universal
speed ratio. These sessions are **not** an A/B comparison against `d466174`;
differences from historical absolute timings must not be attributed solely to
the RK optimization. The controlled plateau A/B study remains separate below.

Selected main-run **construction and shared loading times** are reported separately:

| Input | Workers | F-Max construction (ms) | RK construction (ms) | TTK construction (ms) | Shared loading (ms) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2D terrain, `n=64` | 8 | 24.416 | 23.645 | 0.555 | 3.716 |
| 3D volume, `n=16` | 8 | 133.701 | 141.281 | 4.966 | 11.257 |

F-Max and RK use the same construction implementation; differences between their
columns reflect separate measurements, not separate construction algorithms.
TTK remains substantially faster when this native construction is included.
Shared loading is one read/parse per invocation, summarized across seeds, not a
per-algorithm cost or repeated cold-cache experiment. Construction plus algorithm
equals the resident-to-gradient total **for every raw sample**, but separately
aggregated medians and rounded table values need not add.

All 40 configurations pass exact within-algorithm reference checks. Each
algorithm's critical-count vector is invariant across worker counts; the three
algorithms also agree on all tested inputs. Counts by dimension are:

| Input | Seed | Critical vertices | Critical edges | Critical triangles | Critical tetrahedra |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2D terrain, `n=16` | 0 | 6 | 6 | 1 | -- |
| 2D terrain, `n=16` | 2 | 8 | 9 | 2 | -- |
| 2D terrain, `n=64` | 0 | 126 | 224 | 99 | -- |
| 2D terrain, `n=64` | 2 | 176 | 334 | 159 | -- |
| 3D volume, `n=8` | 0 | 12 | 21 | 13 | 3 |
| 3D volume, `n=8` | 2 | 16 | 26 | 13 | 2 |
| 3D volume, `n=16` | 0 | 21 | 38 | 29 | 11 |
| 3D volume, `n=16` | 2 | 26 | 54 | 36 | 7 |

Matching counts do not imply identical gradients across different algorithms.
Input coverage, hashes, per-sample phase sums, stored medians, paired ratios and
critical counts were checked against the raw evidence. At that revision, tracked fragments
were `resident_gradient_algorithm_table.tex`, `resident_gradient_construction_table.tex`
and `resident_gradient_unprofiled_phases_table.tex`. They were regenerated from the
main study (the historical versions remain in Git); the latter reported every
outer phase for the larger inputs at 1/8
workers. Detailed RK/F-Max timings and full totals remain in both raw JSON files.
TTK lower-star construction and matching remain combined inside its gradient
phase. No persistence work is timed.

### Construction-separated results (historical `d466174`)

These measurements are the `d466174` snapshot, preceding the facet-discovery
scheduling update below. Their timing definition remains the primary comparison;
the refreshed `5598963` results above supersede them as the current measurements.

The main run uses the same Apple M1 Max, native ARM release build, pinned TTK,
eight inputs, 1/2/4/8 workers, 12 performance repetitions, six diagnostic
repetitions and two warmups per mode under `OMP_WAIT_POLICY=PASSIVE`. Every
one of the 32 configurations passes exact within-algorithm reference checks;
critical counts agree across algorithms. The following **algorithm times**
exclude the separately reported construction costs. Values are milliseconds,
medians over the two seed-specific medians:

| Input | Workers | F-Max | RK | TTK |
| --- | ---: | ---: | ---: | ---: |
| 2D terrain, `n=64` | 1 | 1.611 | 1.242 | 1.183 |
| 2D terrain, `n=64` | 8 | 1.579 | 0.961 | 0.878 |
| 3D volume, `n=16` | 1 | 12.580 | 6.016 | 13.379 |
| 3D volume, `n=16` | 8 | 12.609 | 2.572 | 3.983 |

RK has a lower paired-median algorithm time than TTK in all 16 volume
configurations and five of 16 terrain configurations. It is also faster than
F-Max in all 16 volume configurations and ten of 16 terrain configurations;
small parallel terrains can lose to sequential F-Max. This is not a universal
RK ranking or an optimization gain: the algorithm code is unchanged, and the
reported interval now excludes native construction.

Selected **native construction times**, from the same performance repetitions,
are reported separately:

| Input | Workers | F-Max construction | RK construction | TTK construction | Shared loading |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2D terrain, `n=64` | 8 | 28.930 | 28.366 | 0.624 | 3.929 |
| 3D volume, `n=16` | 8 | 151.317 | 150.301 | 5.356 | 11.672 |

Construction values are medians over seed-specific repetition medians.
Shared loading values are medians over one read/parse per seed/worker process;
they do not describe repeated or cold-cache loading. F-Max and RK use the same
construction implementation, measured independently in their balanced run order.
Differences between those construction columns are timing variation, not different
construction algorithms. Full performance totals and all critical-count vectors
remain in `../rk-construction-split-main.json`.

A fresh confirmation in `../rk-construction-split-confirmation.json` repeats
terrain `64` and volume `16`, both seeds, and workers `1 8`, using 18 performance
repetitions, six diagnostic repetitions and two warmups. All eight configurations
again pass exact reference and critical-count checks. RK is faster than TTK in all
four volume configurations; TTK is faster in all four terrain configurations.
At eight workers, the confirmed algorithm medians are F-Max/RK/TTK
`1.778/1.127/0.871` ms for terrain and `13.333/2.981/3.931` ms for volume.
The paired RK/TTK ratio on the larger eight-worker volume is `0.663` in the main
run and `0.734` in confirmation. These preliminary results support the observed
direction, not a universal speed ratio; detailed internal profiles and historical
end-to-end totals are not mixed into this comparison.

## Resident-Array Totals (Historical v1 Comparison)

The measurements and commands in this historical section correspond to the
driver at `b637fb4`. The current driver emits v2 evidence; use the new section
above and fresh output paths for new measurements. Existing v1 JSON and tables
remain unchanged and readable.

This earlier comparison prioritized **one complete gradient from a common in-memory
complex and function**. The starting representation consists of vertex values,
vertex coordinates, and maximal-cell vertex arrays. No vertex ranks,
lower-star lists, full simplex enumeration, incidence arrays, library-specific
triangulation, or builder is prepared for free. The current adapter covers
pure 1D--3D meshes; the generated comparison cases are 2D terrains and 3D
tetrahedral volumes. The function on simplices is the max-vertex extension.

Every repetition starts from these same resident arrays and constructs fresh
native objects. Timing stops when the algorithm's native gradient is ready.
It includes native representation conversion, filtration extension, ordering,
connectivity preparation, builder setup, internal workspace/task-pool creation,
and gradient construction. Input generation, file I/O, reference validation,
destruction after gradient readiness, and persistence are excluded. Temporary
cleanup performed internally before a build method returns remains included.
Warmups prime the process/runtime, not a reusable prepared complex or gradient.

This timing boundary differs from both earlier benchmarks:

| Benchmark | Starting point | Native preparation |
| --- | --- | --- |
| Resident-array comparison | Common mesh arrays and vertex values | Included for every algorithm |
| RK implementation A/B | Finalized MorseFrames filtered complex | Complex construction excluded; fresh builder included |
| Earlier prepared-kernel TTK comparison | Native topology, TTK offsets/connectivity, MorseFrames builder | Excluded from headline kernel times |

The new native executable does not change the historical executable or its
output schema. The new Python runner and table renderer reject legacy timing
schemas, incomplete studies, missing phases, nonfinite values, failed exact
reference checks, and double-counted nested phase times.

### Global and phase timings

The main performance runs have only outer clock reads and no local diagnostic
instrumentation. F-Max, RK, and TTK execute in all six permutations, balanced
over repetition counts divisible by six. Every raw total is retained; the
headline statistics are medians and IQRs, with paired per-repetition ratios.
The worker setting applies to RK and TTK; F-Max remains sequential and is
remeasured alongside them. Native preparation may itself remain sequential.

Separate diagnostic repetitions report these outer phases, which partition
each individual diagnostic total:

- F-Max and RK: native complex construction and filtration preparation;
  builder setup; gradient construction.
- TTK: native object initialization; vertex ordering using TTK's own
  `preconditionOrderArray`; representation setup; connectivity preconditioning;
  gradient construction with the cache bypassed. The same worker budget is
  supplied to TTK's ordering, triangulation, and gradient routines.

Within RK's gradient phase, coarse diagnostics further report workspace/pool
setup, level processing, replay, and an explicit residual. F-Max's existing
finer instrumentation reports workspace initialization, candidate seeding,
candidate selection, emission/updates, callbacks, and a residual. These are
**children of the gradient phase**, not additional total-time components.
The unmodified pinned TTK backend constructs lower stars inside `buildGradient`;
lower-star construction and matching are included but reported together.
Separating those two internal TTK operations would require additional
instrumentation and is not claimed here.

Diagnostic samples never enter the headline speed comparisons. F-Max's fine
timers can noticeably perturb its gradient phase. Phase shares use each
diagnostic sample's own total, not an uninstrumented median; medians over
samples or seeds and rounded phase values need not sum to the displayed
median total. No cumulative worker duration is presented as elapsed time.

### Reproduction

The build helper's second argument selects the new target; omitting it still
builds the historical prepared-kernel executable. Rebuild before measuring.
The TTK revision is embedded in the executable and checked by the runner.
Input, binary, driver, and MorseFrames header hashes plus source state, worker
orders, raw measurements, and the OpenMP wait policy are saved with the study.
Source hashes describe the source state at run time, not a substitute for
rebuilding the executable after edits.

```sh
LC_ALL=C tools/build_ttk_gradient_benchmark.sh \
  ../work/ttk-benchmark morseframes_resident_gradient_benchmark

OMP_WAIT_POLICY=PASSIVE LC_ALL=C python3 tools/benchmark_resident_gradients.py \
  --benchmark ../work/ttk-benchmark/build-f4ffd1a1049d0ccf6e8f3eb4f7c096a6cc251ba0/morseframes_resident_gradient_benchmark \
  --terrain-sizes 16 64 --volume-sizes 8 16 --seeds 0 2 --workers 1 2 4 8 \
  --repeats 12 --diagnostics 6 --warmups 2 \
  --input-dir ../rk-ab-inputs --output ../rk-resident-gradients-main.json

python3 tools/render_resident_gradients.py \
  --input ../rk-resident-gradients-main.json \
  --table-output docs/resident_gradient_comparison_table.tex \
  --phases-output docs/resident_gradient_phases_table.tex
```

The main run covers eight inputs and 32 worker configurations on the Apple M1
Max, using the release native ARM build and `OMP_WAIT_POLICY=PASSIVE`. Every
timed and diagnostic result, including warmups, matches its own algorithm's
one-worker reference: full Morse sequence fields for F-Max/RK and per-cell
critical/upward/downward pairing data for TTK. MorseFrames reference sequences
are also validated independently of timing. Critical counts by dimension are
reported separately for every algorithm; unlike same-algorithm reference
differences, differences between algorithms are allowed and explicitly flagged.
All counts agree in this run.

Selected **uninstrumented total times**, milliseconds, are medians over the
two seed-specific medians:

| Input | Workers | F-Max | RK | TTK |
| --- | ---: | ---: | ---: | ---: |
| 2D terrain, `n=64` | 1 | 24.61 | 24.28 | 1.54 |
| 2D terrain, `n=64` | 8 | 25.23 | 24.72 | 1.36 |
| 3D volume, `n=16` | 1 | 151.71 | 138.13 | 18.26 |
| 3D volume, `n=16` | 8 | 141.94 | 136.32 | 8.02 |

TTK has the lower paired-median total in all 32 configurations. In separate
diagnostics, representation/filtration preparation accounts for about 96%
of RK's total on the largest eight-worker terrain and 98% on the largest
volume. For the latter, RK's diagnostic gradient phase is about 1.67 ms and
TTK's is 3.05 ms, but the complete measured workflows have the opposite
ranking. These diagnostic phase times are not an uninstrumented kernel
comparison. The result identifies the current MorseFrames input adapter and
native preparation as the dominant cost, not an intrinsic limitation of RK.
The adapter currently enumerates each maximal cell's faces and inserts them
through `FilteredSimplicialComplex::add_simplex`, then finalizes the complex;
bulk/compact input construction is a separate optimization opportunity.

A fresh confirmation repeats the largest terrain (`64`) and volume (`16`),
both seeds, and workers `1 8`, with `--repeats 18 --diagnostics 6 --warmups 2`.
It is saved as `../rk-resident-gradients-confirmation.json`. All eight
configurations retain exact reference agreement and matching critical counts,
and TTK again has the lower paired-median total in every configuration.
The confirmed eight-worker RK/TTK total medians are 24.03/1.35 ms for the
terrain and 143.21/8.90 ms for the volume (medians over the two seed medians).
These independent sessions support the direction of the observed ranking,
not a universal speed ratio or an inference about other input representations.

Raw evidence stays outside the public repository. The selected total and
outer-phase tables are tracked; nested per-algorithm diagnostics and all
critical-count vectors remain in the JSON. These results are local to this
implementation, machine, workload family, and wait policy.

The optional native integration tests exercise graph, tied triangle, and
shared tetrahedral inputs, fresh runs at one/four workers, phase accounting,
and invalid-input rejection. Enable them with:

```sh
OMP_WAIT_POLICY=PASSIVE \
MORSEFRAMES_RESIDENT_BENCHMARK=../work/ttk-benchmark/build-f4ffd1a1049d0ccf6e8f3eb4f7c096a6cc251ba0/morseframes_resident_gradient_benchmark \
  python3 -m pytest -q python/tests/test_resident_gradient_benchmark.py
```

## TTK ProcessLowerStars Reference (Prepared Kernels)

The external reference benchmark uses TTK's classic `DiscreteGradient`
backend, which implements the Robins ProcessLowerStars algorithm for explicit
1D--3D triangulations. TTK is pinned to revision
`f4ffd1a1049d0ccf6e8f3eb4f7c096a6cc251ba0`. It is built separately from the
package, without VTK, ParaView, or the standalone applications:

```sh
TTK_BENCHMARK=$(tools/build_ttk_gradient_benchmark.sh)

PYTHONPATH=python python3 tools/benchmark_ttk_process_lower_stars.py \
  --ttk-benchmark "$TTK_BENCHMARK" \
  --terrain-sizes 16 32 64 \
  --volume-sizes 4 8 12 16 \
  --seeds 0 1 2 \
  --workers 1 2 4 8 \
  --repeats 5 \
  --warmups 1 \
  --output ../work/ttk_process_lower_stars.csv

python3 tools/render_ttk_process_lower_stars_table.py \
  --input ../work/ttk_process_lower_stars.csv \
  --output docs/ttk_process_lower_stars_table.tex
```

The driver exports the same injective terrain and tetrahedral-volume cases used
by the unified gradient comparison. For every case and worker count it checks
TTK's critical-simplex count, including the count in each dimension, against
both F-Max and MorseFrames ProcessLowerStars. A mismatch aborts the run.

`ttk_process_lower_stars_seconds` measures only gradient construction after
TTK's explicit triangulation has been preconditioned. Every timed call uses
TTK's cache-bypass path, so it recomputes the gradient rather than fetching a
cached result. `ttk_setup_seconds` and `ttk_precondition_seconds` are reported
separately. This keeps the algorithm comparison focused on gradient
construction while preserving the otherwise hidden topology-preparation cost.

On the Apple M1 Max reference run, all 84 configurations have matching
critical counts by dimension. Sequentially, MorseFrames ProcessLowerStars takes
median times of 7.35 and 5.18 times F-Max in 2D and 3D, while TTK takes 0.51 and
1.18 times F-Max. Thus the current MorseFrames implementation takes 13.46 times
TTK's kernel time in 2D and 4.24 times in 3D. At eight workers, the median
MorseFrames/TTK gaps are 9.42 in 2D and 5.39 in 3D. These results show that the
earlier ProcessLowerStars slowdown is an implementation cost rather than an
inherent property of the Robins algorithm.

The aggregate ratios are recorded in
`docs/ttk_process_lower_stars_table.tex`. Ratios below one denote a faster time
than the denominator. ``TTK setup+kernel'' is a diagnostic estimate formed by
adding one measured explicit-triangulation setup and preconditioning pass to
the prepared gradient-kernel time. The other columns compare prepared gradient
kernels only.

### Controlled ReductionKernel A/B comparison

Use this workflow to evaluate an implementation change. It compiles the same
native driver against two snapshots of the headers using identical compiler
flags. Each process keeps the same topology resident; only one process runs a
gradient at a time while the other waits for a command. It alternates baseline
and candidate order in balanced blocks and saves every timing sample, paired
block ratios, interquartile ranges and a descriptive bootstrap interval.
Run it without other tests, builds or benchmarks in progress.

```sh
python3 tools/benchmark_reduction_kernel_ab.py \
  --baseline d32f71a --candidate WORKTREE \
  --sizes 16 24 32 --seeds 0 1 2 --workers 1 2 4 8 \
  --blocks 12 --repeats 3 --warmups 2 \
  --input-dir ../rk-ab-inputs --output ../rk-packed-cofaces-ab.json
```

By default, inputs use the existing injective tetrahedral-volume generator.
`--family terrain` selects triangulated terrains, and `--filtration plateau`
sets every vertex value to zero while retaining the generated topology. Both
the family and filtration mode are saved with each case. The timed
scope includes a fresh builder, workspace, task-pool setup/teardown and event
replay, with no diagnostic instrumentation or persistence computation.
Topology construction, sequence validation, protocol I/O and destruction of
the returned sequence are outside the timer. The two native sequential
sequences are compared exactly through exported step records, and every timed
sequence is checked field by field against its build's sequential reference.
All checks must pass before a case is recorded. Header, input, driver and
binary hashes, compiler target and flags are retained in the JSON. Choose a
new output path for each run; existing evidence is never overwritten.

The controlled native ARM run on this MacBook Pro (Apple clang 15, C++17,
`-O3 -DNDEBUG -pthread`) covered 36 configurations, including 32-cubed volumes
with 792,051 simplices. The table gives median reductions in total gradient
time across nine size/seed cases per worker count. Positive values mean faster.

| Workers | Packed facets: `b01eaa0` to `ebddcc4` | Packed core: `ebddcc4` to `d32f71a` | Packed cofaces: relative to `d32f71a` |
| ---: | ---: | ---: | ---: |
| 1 | 8.2% | 10.4% | 42.9% |
| 2 | 5.8% | 10.1% | 35.7% |
| 4 | 7.7% | 10.7% | 27.6% |
| 8 | 6.1% | 9.0% | 23.2% |

Packed facets were faster in 34/36 configurations; packed cores were faster in
36/36. All exact sequence comparisons passed. These paired results replace the
earlier small, separate-run estimates of 1.3% and 5.6% for packed facets. Raw
evidence is saved locally as `../rk-packed-facets-ab-1.json` and
`../rk-packed-core-ab-1.json`. To reproduce the first comparison, use
`--baseline b01eaa0 --candidate ebddcc4`; for the second, use
`--baseline ebddcc4 --candidate d32f71a`. Bootstrap intervals describe the
observed session and do not establish uncertainty across machines or sessions.
These measurements do not update the separate TTK comparison below.

A fresh confirmation session used sizes 16 and 32, seeds 0 and 2, workers 1
and 8, and 16 balanced blocks of three repetitions. All eight configurations
again favored the packed core; median reductions were 10.2% sequentially and
8.0% with eight workers. The 32-cubed, seed-zero sequential case was nearly
tied (1.1% reduction with an interval spanning parity), so the benefit is not
uniform. The raw confirmation is `../rk-packed-core-ab-confirmation.json`.
Use the corresponding `--sizes`, `--seeds`, `--workers` and `--blocks` options
with a new output filename to repeat it.

Packed cofaces additionally remove adjacency-list traversal from packed facet
discovery and the local unique-coface test. Their mask construction and active
mask updates are included in the measured total. All 36 configurations favored
this change, and each case's paired-block bootstrap interval was below parity.
All exact sequence checks passed. Raw samples are saved as
`../rk-packed-cofaces-ab-1.json`.
A fresh run on sizes 16 and 32, seeds 0 and 2, workers 1 and 8, with 16 blocks
of three repetitions again favored packed cofaces in all eight configurations.
Median reductions were 42.7% sequentially and 25.6% with eight workers, and all
paired-block intervals were below parity. Its exact sequence checks also
passed; raw evidence is `../rk-packed-cofaces-ab-confirmation.json`.

### Direct TTK versus parallel ReductionKernel (Prepared Kernels)

The earlier prepared-kernel comparison runs F-Max, TTK ProcessLowerStars, and the
parallel ReductionKernel in the same native process. The execution order is
balanced over all six permutations of the algorithms. Each CSV row retains
all timing samples and their interquartile ranges. Reported times are medians;
ratios between algorithms are medians of their paired per-repetition ratios.
This replaces the former independently selected best-time estimates.
Both TTK and MorseFrames topology construction are outside the gradient timing,
and TTK's gradient cache is bypassed on every call. Unlike the total-construction
A/B benchmark above, this comparison reuses the MorseFrames builder outside
the timer. It measures the prepared gradient kernels, including RK's task pool
and event replay. The returned sequence is assigned to a reusable output.

The main run uses `OMP_WAIT_POLICY=PASSIVE` so idle OpenMP workers do not spin
while another algorithm is being measured. The selected policy is recorded
in every CSV row. Timing and especially small-input scaling depend on this
runtime setting; comparisons must use the same policy and timing definition.

```sh
OMP_WAIT_POLICY=PASSIVE PYTHONPATH=python python3 tools/benchmark_ttk_reduction_kernel.py \
  --ttk-benchmark "$TTK_BENCHMARK" \
  --terrain-sizes 16 32 64 \
  --volume-sizes 4 8 12 16 24 32 \
  --seeds 0 1 2 \
  --workers 1 2 4 8 \
  --repeats 18 \
  --warmups 3 \
  --output ../ttk_reduction_kernel_packed_cofaces.csv

python3 tools/render_ttk_reduction_kernel_table.py \
  --input ../ttk_reduction_kernel_packed_cofaces.csv \
  --output docs/ttk_reduction_kernel_table.tex
```

All 108 configurations have identical critical-simplex counts by dimension.
On 2D terrains, TTK is faster in all 36 configurations: the median
ReductionKernel/TTK ratios are 1.26, 1.61, 2.03, and 1.95 at one, two, four,
and eight workers. RK's median ratios to sequential F-Max are 0.79, 0.78, 0.83
and 1.07. Additional workers do not benefit these small RK terrain workloads.

On 3D volumes, the median ReductionKernel/TTK ratios are 0.43, 0.53, 0.62,
and 0.67 at one, two, four, and eight workers. RK wins in 65 of 72
configurations, including all 18 sequential cases. The median eight-worker
ratios to F-Max are 0.20 for RK and 0.36 for TTK. RK reaches a median 2.22-fold
eight-worker speedup from its own one-worker time, while TTK reaches 3.45-fold.
TTK scales more strongly but starts from a slower one-worker time on these
volumes. Restricting to the 24-cubed and 32-cubed inputs, eight-worker RK has
median ratios of 0.12 to F-Max and 0.67 to TTK.

These are workload- and runtime-specific prepared-kernel comparisons, not a
claim that RK is uniformly fastest. The aggregate table is stored in
`docs/ttk_reduction_kernel_table.tex`. This refresh changes the statistic,
OpenMP policy and input-size range from the older `ebddcc4` table; the controlled
A/B experiment above is the evidence for the isolated implementation gain.

A runtime-policy check reran terrains of sizes 32/64 and volumes of sizes
16/32, seeds 0/2, with 1/8 workers and `OMP_WAIT_POLICY` unset. All 16 critical
count comparisons passed. On its four eight-worker volume cases, RK/TTK was
0.52 (versus 0.65 on the matching passive-policy cases), and RK won every case
under both policies. The eight-worker terrain ranking reversed: RK/TTK was
0.34 with the default policy, versus 1.98 with passive waiting. Thus the 3D
advantage survived this check, while the 2D parallel comparison is sensitive
to runtime behavior and should not be presented as an unconditional algorithm
ranking. These were separate sessions, not an interleaved policy experiment.
The control samples are in `../ttk_reduction_kernel_cofaces_default_policy.csv`.
To reproduce them, use `env -u OMP_WAIT_POLICY` with the same command and
`--terrain-sizes 32 64 --volume-sizes 16 32 --seeds 0 2 --workers 1 8`, writing
to a separate output file.

## Unified Gradient-Only Strategy Comparison

This is the central internal benchmark for discrete-gradient construction. It
compares F-Max, sequential and parallel ProcessLowerStars, and sequential and
parallel ReductionKernel on both triangulated 2D terrains and tetrahedral 3D
volumes. The timed path calls only `profile_morse_sequence`: it does not build a
reference map, a reduction plan, or a persistence diagram. Parallel sequences
are constructed once outside the timing loop and must agree exactly with their
sequential counterpart. Critical simplices are recorded by dimension and
compared with F-Max. ReductionKernel uses its default uncached path here, which
is the honest comparison for one gradient; repeated-gradient cache reuse is
measured separately below.

```sh
PYTHONPATH=python python3 tools/benchmark_gradient_strategies.py \
  --terrain-sizes 16 32 64 \
  --volume-sizes 4 8 12 16 \
  --seeds 0 1 2 \
  --workers 1 2 4 8 \
  --repeats 11 \
  --warmups 1 \
  --format csv \
  --output ../work/gradient_strategy_benchmark.csv

MPLCONFIGDIR=../work/matplotlib-cache \
  python3 tools/render_gradient_strategy_benchmark.py \
  --input ../work/gradient_strategy_benchmark.csv \
  --figure-output docs/gradient_strategy_comparison.svg \
  --table-output docs/gradient_strategy_comparison_table.tex
```

The run contains 21 complexes and 231 measured rows. All parallel gradients
match their sequential counterpart exactly, and all five approaches have zero
critical-count difference from F-Max in every case. In 2D, sequential
ProcessLowerStars and ReductionKernel take median times of 6.19 and 1.45 times
the F-Max time. At eight workers these ratios fall to 3.47 and 0.99, so the
parallel ReductionKernel is the faster of the two and is effectively tied with
F-Max. Four workers give the best aggregate 2D ReductionKernel ratio, 0.90;
additional scheduling overhead outweighs useful work at eight workers. In 3D,
the corresponding sequential ratios are 4.71 and 1.41; at eight workers they
fall to 1.69 and 0.49. ReductionKernel is therefore the faster parallel method
in both dimensions and is about twice as fast as F-Max over the aggregate 3D
corpus. At grid sizes 12 and 16 it reaches 0.39 and 0.33 times the F-Max time,
respectively. ReductionKernel scales from one to eight workers by 1.40-fold in
2D and 2.99-fold in 3D, versus 1.66-fold and 2.80-fold for ProcessLowerStars.
The lower 2D ReductionKernel scaling reflects saturation at four workers, not
a regression in absolute performance relative to F-Max.

The profiler times an uninstrumented construction and collects detailed phase
counters in a separate diagnostic run. This prevents high-frequency timing
calls inside ReductionKernel facet tasks from biasing the comparison. The
optimized kernel caches compact same-level face closures for triangular and
tetrahedral sections, reuses per-worker level scratch across filtration levels
and facet-discovery and result buffers across rounds, and stores the small
per-facet cells, removal masks, and events inline. Facet events are consumed
directly in deterministic order instead of rebuilding intermediate event
vectors. Per-level events are written into disjoint slices of one preallocated
arena, eliminating one allocation per filtration level without changing
reverse-per-level replay. Each event stores only its lower and upper simplex;
an invalid upper simplex denotes a perforation. The 16-byte events are
placement-constructed only when emitted, so unused arena capacity is neither
initialized nor read. Higher-dimensional cells transparently fall back to
dynamic storage. Relative to the preceding cached-closure implementation, the
median ReductionKernel/F-Max ratio falls from 5.51 to 3.04 sequentially in
2D and from 3.74 to 2.06 sequentially in 3D; the eight-worker 3D ratio falls
from 0.94 to 0.68.
Dynamic level claiming then lowers the eight-worker ratio from 1.89 to 1.33 in
2D and from 0.68 to 0.63 in 3D by eliminating level sorting and static
simplex-count partitions.
Reusing the level scratch owned by each worker subsequently lowers the
sequential ratios from 3.08 to 2.06 in 2D and from 2.14 to 1.81 in 3D. The
eight-worker ratios also fall from 1.33 to 1.31 and from 0.63 to 0.62,
respectively.
The shared event arena then lowers those eight-worker ratios from 1.31 to 1.13
in 2D and from 0.62 to 0.56 in 3D; the sequential ratios fall from 2.06 to 1.86
and from 1.81 to 1.74.
In a seven-repeat paired comparison against the initial 24-byte,
value-initialized arena, the compact uninitialized representation improves
median eight-worker time by 6.5 percent in 2D and 5.9 percent in 3D. Its
one-worker improvements are 1.3 and 2.8 percent, respectively.
The ordinary construction path now instantiates a compile-time metrics-free
kernel: clocks, diagnostic counters, per-facet diagnostic arrays, and the
per-level metrics vector remain available to the separate profiling run but
are absent from timed gradient construction. In a seven-repeat paired run over
all 105 ReductionKernel configurations, this lowers median time by 7.7 percent
and wins 89 comparisons. Median improvements are 7.7 percent sequentially and
15.2 percent at eight workers in 2D, and 5.9 and 6.4 percent, respectively, in
3D. The few regressions are concentrated in the smallest sub-millisecond cases.
Kernel-round merging is event-driven: only simplices named by accepted facet
reductions are visited and cleared, avoiding two full level-bucket passes per
round. Against the preceding topology-cache implementation, this reduces
median cached gradient time by 7.1 percent in 2D and 5.7 percent in 3D.
Sequential facet discovery also maintains an ordered compact list of active
simplices. Removed entries are discarded while discovering the next facets,
and the same list bounds incidence reset and low-dimensional cell construction.
In the initial seven-repeat comparison against event-driven merging alone,
cached median time falls by another 7.8 percent in 2D and 11.9 percent in 3D.
The metrics-free sequential path consumes each facet result immediately and
retains only its compact reduction events for the coordinator merge. Removing
the array of large intermediate facet-result objects lowers cached median time
by a further 4.9 percent in 2D and 7.7 percent in 3D.
Metrics-free facet results are separately specialized to contain only the
inline event buffer; diagnostic counters exist only in the instrumented type.
Across two confirmation runs this lowers default uncached median time by 3.2
percent in 2D and 4.7 percent in 3D, while cached time remains neutral in 2D
and improves by 1.9 percent in 3D.

### Reusable ReductionKernel Topology Cache

The owning `FilteredComplex` can explicitly precompute immutable same-level
closure ranges and coboundary adjacency for repeated sequential ReductionKernel
gradients. The focused benchmark constructs separate cached and uncached
copies, alternates their measurement order, verifies identical sequences, and
reports cache build time and memory separately:

```sh
PYTHONPATH=python python3 tools/benchmark_reduction_kernel_cache.py \
  --terrain-sizes 16 32 64 \
  --volume-sizes 4 8 12 16 \
  --seeds 0 1 2 \
  --repeats 11 \
  --warmups 2 \
  --output ../work/reduction_kernel_cache.csv
```

Across all 21 cases, every cached gradient exactly matches its uncached
counterpart and every cached run is faster. Median sequential speedup is
1.43-fold on terrains and 1.62-fold on tetrahedral volumes. Median cache build
costs are 0.23 ms and 1.45 ms, with median allocated footprints of 0.27 MiB and
1.65 MiB, respectively. The build cost is recovered after median counts of
1.57 terrain gradients and 1.34 volume gradients. The largest `n=16` volume
cache occupies 5.80 MiB. Multiworker ReductionKernel deliberately retains
worker-local topology construction because shared-cache access did not improve
its wall time.

![Gradient-only strategy comparison](gradient_strategy_comparison.svg)

The aggregate values are generated in
`docs/gradient_strategy_comparison_table.tex`; the raw CSV retains the exact
critical counts by dimension for every complex and worker count.

## Simplicial Strategy Comparison

The first internal comparison uses connected triangulated terrains rather than
independent synthetic lower stars. A smooth random field is sampled on each
grid, the vertices are strictly ranked, and every higher-dimensional simplex
receives the value of its maximum vertex. This gives an injective lower-star
filtration while retaining nontrivial terrain topology.

```sh
PYTHONPATH=python python3 tools/benchmark_simplicial_strategies.py \
  --sizes 16 32 64 \
  --seeds 0 1 2 \
  --parallel-workers 8 \
  --repeats 5 \
  --warmups 1 \
  --format csv \
  --output ../work/simplicial_strategy_benchmark.csv

MPLCONFIGDIR=../work/matplotlib-cache \
  python3 tools/render_simplicial_strategy_benchmark.py \
  --input ../work/simplicial_strategy_benchmark.csv \
  --figure-output docs/simplicial_strategy_comparison.svg \
  --table-output docs/simplicial_strategy_comparison_table.tex \
  --kernel-table-output docs/reduction_kernel_metrics_table.tex
```

Every measured Morse pipeline is checked against ordinary persistence for the
same barcode. The output records critical counts by dimension, regular-pair
counts, exact sequence agreement with ProcessLowerStars, and both construction
and end-to-end timings. Dedicated reduction-kernel fields record its levels,
rounds, facet kernels, reductions, perforations, parallel batches, concurrency,
and instrumented phase work. The phase durations are cumulative across tasks;
for a parallel row they are not wall-clock durations. Each reported case uses
the fastest of five measured runs after one warm-up.

The initial Apple M1 Max run covers nine cases (three sizes by three seeds).
ProcessLowerStars, reduction kernels, both eight-worker versions, F-Max, F-Min,
and Saturated produce the same critical count in every case. Same-level
reduction produces a median of 3.95 times as many critical simplices. Relative
to F-Max, median end-to-end time is 2.35 times as large for sequential
ProcessLowerStars and 2.39 times as large for sequential reduction kernels.
Eight-worker ProcessLowerStars improves to 1.29 times the F-Max time, while
the optimized eight-worker reduction kernel reaches 1.09 times the F-Max time.
The structural counters are unchanged: the runtime improvement comes from
coarse scheduling rather than a different Morse complex.

![Simplicial strategy comparison](simplicial_strategy_comparison.svg)

The grid-size aggregates are generated in
`docs/simplicial_strategy_comparison_table.tex`, and the reduction-kernel
operation and scheduling counters in `docs/reduction_kernel_metrics_table.tex`.
This is a MorseFrames-internal comparison; it does not replace the planned
external benchmark against Robins' ProcessLowerStars implementation.

## Tetrahedral Strategy Comparison

The three-dimensional companion uses a conforming Freudenthal triangulation:
each grid cube is split into the six tetrahedra defined by the permutations of
the coordinate axes. The same smooth-field ranking makes all vertex values
distinct, and the lower-star extension assigns every edge, triangle, and
tetrahedron the value of its maximum vertex. Sizes 4, 8, and 12 contain 883,
9,843, and 36,851 simplices, respectively, so they span approximately the same
range as the two-dimensional corpus.

```sh
PYTHONPATH=python python3 tools/benchmark_tetrahedral_strategies.py \
  --sizes 4 8 12 16 \
  --seeds 0 1 2 \
  --parallel-workers 8 \
  --repeats 5 \
  --warmups 1 \
  --format csv \
  --output ../work/tetrahedral_strategy_benchmark.csv

MPLCONFIGDIR=../work/matplotlib-cache \
  python3 tools/render_simplicial_strategy_benchmark.py \
  --input ../work/tetrahedral_strategy_benchmark.csv \
  --figure-output docs/tetrahedral_strategy_comparison.svg \
  --table-output docs/tetrahedral_strategy_comparison_table.tex \
  --kernel-table-output docs/tetrahedral_reduction_kernel_metrics_table.tex \
  --title "MorseFrames strategies on injective tetrahedral volumes"
```

The runner preserves the two-dimensional timing boundaries and correctness
checks. In particular, every strategy must reproduce the ordinary-persistence
barcode, while sequential and parallel versions of each new construction must
produce the same dimension-wise critical counts.

On the Apple M1 Max, the nine-case run finds identical critical counts for
ProcessLowerStars, reduction kernels, both eight-worker versions, F-Max, F-Min,
and Saturated. Same-level reduction creates a median of 7.57 times as many
critical simplices. Relative to F-Max, median end-to-end time is 2.27 times as
large for sequential ProcessLowerStars and 3.31 times as large for sequential
reduction kernels. The eight-worker versions reduce those ratios to 1.16 and
1.19, respectively. Parallel ProcessLowerStars is therefore slightly faster on
this rerun, although both remain close to the highly optimized sequential F-Max
path in three dimensions.

![Tetrahedral strategy comparison](tetrahedral_strategy_comparison.svg)

The grid-size aggregates are generated in
`docs/tetrahedral_strategy_comparison_table.tex`, with reduction-kernel counters
in `docs/tetrahedral_reduction_kernel_metrics_table.tex`.

## Tetrahedral Worker Scaling

The focused 3D scaling study times only discrete-gradient construction for the
parallel ProcessLowerStars and reduction-kernel implementations at one, two,
four, and eight workers. Each algorithm uses its own one-worker parallel
execution as the speedup baseline; every measured sequence is checked for exact
agreement with the corresponding sequential implementation. No reference map
or persistence computation is performed.

```sh
PYTHONPATH=python python3 tools/benchmark_tetrahedral_worker_scaling.py \
  --sizes 4 8 12 \
  --seeds 0 1 2 \
  --workers 1 2 4 8 \
  --repeats 5 \
  --warmups 1 \
  --format csv \
  --output ../work/tetrahedral_worker_scaling.csv

MPLCONFIGDIR=../work/matplotlib-cache \
  python3 tools/render_tetrahedral_worker_scaling.py \
  --input ../work/tetrahedral_worker_scaling.csv \
  --figure-output docs/tetrahedral_worker_scaling.svg \
  --table-output docs/tetrahedral_worker_scaling_table.tex
```

Across the nine cases, eight-worker ProcessLowerStars reaches a median
gradient-construction speedup of 2.61, versus 3.09 for the reduction kernel.
The corresponding median construction times are 1.39 ms and 0.45 ms, so the
optimized ReductionKernel is now faster in absolute time. Median eight-worker
efficiencies are 0.33 and 0.39, respectively.

![Tetrahedral worker scaling](tetrahedral_worker_scaling.svg)

The grid-size aggregates are generated in
`docs/tetrahedral_worker_scaling_table.tex`.

## Lightweight RK Initialization

`ReductionKernelSequenceBuilder` reuses the existing RK kernel without allocating
the general-purpose builder's rank, level, and dimension caches. Initialization
still validates the complete filtration-order permutation with a temporary byte
per simplex. Generic `FSequenceBuilder` initialization is unchanged. The shared
reference-frame and compact persistence-input paths and native Python RK entry
points select the lightweight builder. See the
[C++ API contract](cpp_complex_view_api.md#lightweight-reductionkernel-initialization).

This is an implementation change on the current development branch, not a change
to the gradient algorithm or a migration of `main`. It applies to RK's gradient
stage within full persistence computations; total persistence speedup must be
measured separately and is not inferred from gradient timings.

The A/B driver selects the lightweight builder when its header is available and
the legacy general-purpose builder otherwise. The same driver is compiled against
both snapshots. Both versions time fresh builder initialization, workspace, worker
pool, matching, replay and internal teardown. Neither moves preparation outside
the timer. Native complex construction, file loading, correctness checks and
returned-gradient destruction remain excluded. The resident-array three-way driver
also selects the new RK builder, with loading and native construction still
reported separately. Its current tables are refreshed on `5598963` above; the
controlled A/B study below measures the initialization change separately.

The controlled baseline is `afeaa9bf8d1034864f88661e4d48223258ceafae`.
Main runs use eight alternating blocks of three unprofiled measurements, two
warmups, workers 1/2/4/8, and the following families. Reuse the input directory but
choose fresh output paths when reproducing these commands:

```sh
python3 tools/benchmark_reduction_kernel_ab.py \
  --baseline afeaa9b --family terrain --filtration lower-star \
  --sizes 16 64 --seeds 0 2 --workers 1 2 4 8 --blocks 8 --repeats 3 --warmups 2 \
  --input-dir ../rk-ab-inputs --output ../rk-lean-builder-terrain-ab.json
python3 tools/benchmark_reduction_kernel_ab.py \
  --baseline afeaa9b --family volume --filtration lower-star \
  --sizes 16 32 --seeds 0 2 --workers 1 2 4 8 --blocks 8 --repeats 3 --warmups 2 \
  --input-dir ../rk-ab-inputs --output ../rk-lean-builder-volume-ab.json
python3 tools/benchmark_reduction_kernel_ab.py \
  --baseline afeaa9b --family terrain --filtration plateau \
  --sizes 16 32 --seeds 0 --workers 1 2 4 8 --blocks 8 --repeats 3 --warmups 2 \
  --input-dir ../rk-ab-inputs --output ../rk-lean-builder-terrain-plateau-ab.json
python3 tools/benchmark_reduction_kernel_ab.py \
  --baseline afeaa9b --family volume --filtration plateau \
  --sizes 4 8 12 --seeds 0 --workers 1 2 4 8 --blocks 8 --repeats 3 --warmups 2 \
  --input-dir ../rk-ab-inputs --output ../rk-lean-builder-volume-plateau-ab.json
```

The scripts retain raw timings, alternating orders, source/header/binary/input
hashes, exact sequence checks and critical counts. Bootstrap intervals describe
paired blocks within one desktop session, not independent-machine uncertainty.
Diagnostic phase timings are separate from performance estimates. No builds or
test runs overlap the measurements; unrelated desktop activity is not suppressed.

### Controlled results and validation

The main study contains 52 configurations. A second study repeats all 52 with
12 blocks of five measurements; use the same four commands above with
`--blocks 12 --repeats 5` and replace each output suffix `-ab.json` with
`-confirmation.json`. Both studies use two warmups. Every configuration passes
exact old/new gradient checks; input hashes, critical counts, repetition counts,
all raw-sample medians and paired summaries have also been checked independently.

Both studies measured the same candidate header snapshot, SHA-256
`03636fd51de667e246fe3c8f9629dad2200b8a38171fa6a3d1d66d6388d3aa7b`,
against baseline header digest
`babf46699673de2aee8dee6c9091893f52c257f62b5a9efc995a09c2a47fd927`.
The common driver digest is
`e86b237a81d057e2bac530c9ce65aee78728989c679272673acd503278bf7fd5`.
The candidate was an uncommitted worktree based on `afeaa9b`; an immutable rerun
should select the revision with the matching header digest via `--candidate`.
The native ARM builds use Apple Clang 15 and `-std=c++17 -O3 -DNDEBUG -pthread`.
Each raw file records its actual binary digest and timing scope.

For the largest ordinary volume (`n=32`, 792,051 simplices), every one of the
eight seed/worker configurations has a below-one paired bootstrap interval in
both studies. Confirmation times are:

| Workers | Previous RK (ms) | Lightweight RK (ms) |
| --- | ---: | ---: |
| 1 | 61.865 | 55.045 |
| 2 | 38.937 | 32.691 |
| 4 | 29.524 | 21.383 |
| 8 | 25.264 | 17.619 |

Entries are medians over the two seed-specific sample medians, with native
construction excluded. At eight workers, the main-study paired ratios are
0.696 and 0.701; confirmation ratios are 0.628 and 0.734 (approximately 27--37%
less gradient time). Confirmation intervals are [0.580, 0.677] and [0.701, 0.822].
At one worker, confirmation paired reductions are about 12--13%. Paired ratios
are computed per alternating block and need not equal quotients of table entries.

Ordinary volume `n=16` also improves at eight workers in both studies, but its
estimated gain varies between sessions. Terrain and constant-value plateau
results do not establish a broad, repeatable improvement. No configuration has
an above-one interval in either study; that is not proof of equivalence or absence
of regressions on other inputs. Small inputs remain sensitive to scheduling and
desktop noise. All samples, including outliers, are retained.

### Separate initialization profile

```sh
LC_ALL=C python3 tools/benchmark_reduction_kernel_phases.py \
  --input-dir ../rk-ab-inputs --output ../rk-phases-lean-builder.json
```

This covers the same 13 inputs at 1/2/4/8 workers, with 24 unprofiled, eight coarse,
and three detailed samples per configuration. Every raw summary and speedup was
recomputed. Input/sequence hashes, critical counts, and 16 deterministic work
counters match `rk-phases-facet-work.json` in all 52 configurations. In particular,
kernel rounds, facet work and scheduling counts have not changed.

On volume `n=32` at eight workers, coarse builder initialization is now
0.632/0.788 ms for seeds 0/2, about 3.7/4.0% of coarse total time. The earlier
profile recorded 7.743/7.275 ms and 31.5/30.8%. These profiles explain the removed
work but are separate sessions and do not supply the performance speedup claims;
those come from the unprofiled A/B samples above. The new phase binary digest is
`f68fb800518de42ffe4e8db91a43eb629f960c0e3f669d4395da54f8801d59a1`.

Normal and address/undefined-behavior-sanitized C++ tests pass, including malformed
permutations, generic views, repeated/concurrent const calls, callbacks, cached
and uncached kernels, exact gradients and persistence barcodes. Native Python:
147 tests pass, including the rebuilt resident TTK integration. The fallback
suite runs 147 tests with nine native-only skips; both backend corpus checks pass
all nine cases. Strict documentation and generated-summary checks pass. Full
persistence correctness is verified, but full persistence runtime is not measured
by this gradient-only study.

## Work-Aware Facet Execution

Local facet execution now uses a work estimate to select the task count, rather
than always dispatching up to one task per worker whenever two facets remain.
For each current facet, the estimate counts:

- sparse/cached closure entries, including inactive entries the kernel still scans;
- live closure bits in the packed representation;
- the active-list length for the graph fallback, whose local kernel scans that list.

The task budget is `min(workers, facet_count, floor(estimated_work / 1024))`.
Budgets below two use sequential execution. Estimation stops once the maximum
budget is reached and uses saturating arithmetic. It does not reconstruct local
cells, inspect every sparse closure entry, or run the kernel twice. The 1,024-unit
threshold is an implementation heuristic, not a hardware-independent optimum or
a prediction of local reduction cost. Unequal closures and repeated reduction
searches can still produce unequal task costs.

The same decision feeds metrics-free, coarse-profiled and detailed execution.
Detailed facet-execution wall time includes estimation. Sequential policy and
independent-level tasks bypass estimation; level parallelism and the shared worker
budget are unchanged. The existing dynamic chunk claiming, canonical result slots,
ordered merge and exception draining are retained. No gradient choice changes.

Tests cover work totals 2,047/2,048/2,049, partial worker budgets, contraction to
cheap rounds, cached/uncached topology, graph scans, packed masks, uneven facets,
and 9-dimensional simplices with overflowing inline storage. They compare exact
sequence fields and work counters with sequential execution; failure-injection
tests still require all submitted tasks to drain before an exception escapes.

Local validation passes the C++ suite both normally and with AddressSanitizer /
UndefinedBehaviorSanitizer, all 146 native Python tests (including the rebuilt
optional resident TTK benchmark integration), and the fallback suite (146 run,
nine native-only tests skipped). The nine-case correctness corpus passes with
both backends at 1/2/4 workers. Strict Sphinx and generated-summary checks pass.

### Controlled comparison

Baseline: `012d804`. Candidate: `591b658`, the work-aware implementation identified in the
raw outputs by header SHA-256
`babf46699673de2aee8dee6c9091893f52c257f62b5a9efc995a09c2a47fd927`.
The native ARM Clang 15 flags are `-std=c++17 -O3 -DNDEBUG -pthread` on the
Apple M1 Max. Each timing includes fresh builder, workspace, pool, gradient,
replay and internal teardown from an already finalized native complex.
Loading/construction, validation, protocol I/O and returned-sequence destruction
are excluded equally. There is no persistence computation or new TTK/F-Max run.

```sh
LC_ALL=C python3 tools/benchmark_reduction_kernel_ab.py \
  --baseline 012d804 --candidate 591b658 \
  --family terrain --filtration plateau --sizes 16 32 --seeds 0 \
  --workers 1 2 4 8 --blocks 8 --repeats 3 --warmups 2 \
  --input-dir ../rk-ab-inputs \
  --output ../rk-facet-work-terrain-plateau-ab.json
```

Run the other families with the same options and the following substitutions:

| Family | Filtration | Sizes | Seeds | Output basename |
| --- | --- | --- | --- | --- |
| volume | plateau | 4, 8, 12 | 0 | `rk-facet-work-volume-plateau-ab.json` |
| terrain | lower-star | 16, 64 | 0, 2 | `rk-facet-work-terrain-lower-star-ab.json` |
| volume | lower-star | 16, 32 | 0, 2 | `rk-facet-work-volume-lower-star-ab.json` |

The confirmation repeats **every** configuration with `--blocks 16 --repeats 5`
and output suffix `-confirmation.json` in place of `-ab.json`. Both studies
alternate revision order within each configuration and retain every raw sample.
No benchmark timings overlap our test or build processes. Other desktop
applications were active; the paired-block bootstrap describes these sessions,
not independent-machine uncertainty or a quiet-machine performance guarantee.

Confirmed eight-worker plateau results are:

| Input | Before (ms) | After (ms) | Paired after/before ratio | 95% paired-block interval |
| --- | ---: | ---: | ---: | --- |
| 2D terrain, `n=16` | 0.844 | 0.445 | 0.541 | [0.504, 0.644] |
| 2D terrain, `n=32` | 2.807 | 2.236 | 0.805 | [0.753, 0.825] |
| 3D volume, `n=4` | 0.504 | 0.303 | 0.597 | [0.552, 0.705] |
| 3D volume, `n=8` | 2.191 | 1.934 | 0.898 | [0.885, 0.947] |
| 3D volume, `n=12` | 8.456 | 7.966 | 0.941 | [0.906, 0.975] |

Times are medians of raw samples; ratios are medians of paired block ratios,
not quotients of those displayed medians. These are within-RK improvements over
the facet-discovery revision, not RK speedups relative to TTK or F-Max.

All 52 main and 52 confirmation configurations preserve every gradient-sequence
field and critical count. Raw timings, paired summaries, input hashes and the
candidate header digest were independently checked. All five eight-worker plateau
paired medians improve in both runs; all five confirmation intervals lie below
one. Smaller worker budgets have mixed or inconclusive effects: the larger
terrain's two-worker ratio changes from 1.180 in the main run to 1.001 in the
confirmation, and the larger volume's from 1.029 to 0.987.

No ordinary lower-star control has an above-one bootstrap interval in both runs.
The conspicuous main-run sequential `n=32`, seed-2 volume slowdown (ratio 1.286)
does not repeat (0.988). Eight-worker median control ratios are 0.975 / 0.978 for
terrains and 1.011 / 1.030 for volumes (main / confirmation). This supports the
plateau improvement with a variability caveat, not proven equivalence on ordinary
inputs, universal speedup, or a claim that all facet time was scheduling overhead.
The 1,024-unit threshold was not tuned between these two studies.

### Separate phase and scheduling check

```sh
LC_ALL=C python3 tools/benchmark_reduction_kernel_phases.py \
  --input-dir ../rk-ab-inputs --output ../rk-phases-facet-work.json

python3 tools/render_reduction_kernel_phases.py \
  --input ../rk-phases-facet-work.json \
  --table-output docs/reduction_kernel_facet_work_phases_table.tex
```

This records 13 inputs at 1/2/4/8 workers, with 24 uninstrumented, eight coarse
and three detailed samples per configuration. All 52 exact-sequence checks pass.
Against `rk-phases-discovery.json`, sequence hashes, kernel/round counts,
incidence/candidate/coface visits, overflow counts and discovery task counts are
unchanged. Submitted local facet tasks at eight workers decrease as follows:

| Plateau input | Before tasks | After tasks | Before parallel rounds | After parallel rounds |
| --- | ---: | ---: | ---: | ---: |
| 2D terrain, `n=16` | 121 | 9 | 16 | 4 |
| 2D terrain, `n=32` | 249 | 113 | 32 | 20 |
| 3D volume, `n=4` | 47 | 2 | 6 | 1 |
| 3D volume, `n=8` | 95 | 46 | 12 | 7 |
| 3D volume, `n=12` | 143 | 94 | 18 | 13 |

The smaller terrain uses at most three facet tasks per round, and the smallest
volume uses two despite an eight-worker budget. Local execution remains about
67% of detailed level wall time on the larger terrain, and about 31% on the
larger volume. These separate diagnostic shares are not causal decompositions
of the uninstrumented A/B improvement; cumulative core/local times nest inside
facet execution and must not be added to it.

Parallel RK still does not consistently beat sequential RK on plateaus. The
unprofiled phase samples give 2.089 ms versus 1.377 ms on the larger terrain,
and 8.226 ms versus 8.160 ms on the larger volume (eight workers versus one).
Ordinary largest-volume inputs retain roughly 2.62-fold eight-worker speedup
in this profile; builder initialization is about 31% of coarse total time.
This remains a distinct optimization target, not part of the scheduling change.

## Facet-Discovery Scheduling Update (Historical `3c18fa7`)

Sparse facet discovery now dispatches work according to the **remaining active
simplex count**, not the original level-bucket size. It submits
`min(workers, floor(active_count / 4096))` tasks only if that number is at least
two and intra-level parallelism is allowed. Each task receives a contiguous,
balanced range of at least 4,096 active candidates. Smaller rounds use the
existing ordered sequential scan. This threshold is an implementation heuristic
for the cheap immediate-coface test, not a claim of optimal scheduling on every
machine or complex.

Before parallel discovery, the coordinator stably compacts the active list.
Workers write disjoint flag/counter slots for that list; unused buffer capacity
is neither cleared nor scanned. The coordinator collects facets in the original
canonical order after joining every task. Submission or accessor exceptions drain
all submitted work before captured buffers can leave scope. Cached same-level
coboundaries use the same policy. Packed-mask discovery, local facet-task batching,
independent-level scheduling, and the global worker budget are unchanged.

Native tests exercise 8,191/8,192/8,193 active simplices, the three-task threshold,
an eight-task case, cached/uncached topology, and contraction to one/two active
vertices. They require exact gradient fields and equal incidence/coface-visit
counts against sequential RK. Failure-injection tests check that discovery tasks
finish before an exception propagates and the executor remains usable. Python
tests also check the dispatch count on both sides of the threshold. Graph-only
local-cell construction is not changed by this scheduling update.

The controlled A/B baseline is `d466174`. Its gradient headers are the same as
the preceding `ec408dc` facet-task batching revision. The measured scheduling
candidate is `3c18fa7`. Timing starts from an
already finalized native complex and includes fresh builder/workspace/pool
creation, gradient computation, replay and internal teardown. Input loading and
native complex construction remain outside this algorithm comparison. This RK
implementation A/B study is not a fresh comparison against TTK or F-Max.

```sh
LC_ALL=C python3 tools/benchmark_reduction_kernel_ab.py \
  --baseline d466174 --candidate 3c18fa7 \
  --family volume --filtration plateau --sizes 4 8 12 --seeds 0 \
  --workers 1 2 4 8 --blocks 8 --repeats 3 --warmups 2 \
  --input-dir ../rk-ab-inputs \
  --output ../rk-discovery-volume-plateau-ab.json
```

Complementary main runs use the same options with these substitutions:

| Family | Filtration | Sizes | Seeds | Output basename |
| --- | --- | --- | --- | --- |
| terrain | plateau | 16, 32 | 0 | `rk-discovery-terrain-plateau-ab.json` |
| volume | lower-star | 16, 32 | 0, 2 | `rk-discovery-volume-lower-star-ab.json` |
| terrain | lower-star | 16, 64 | 0, 2 | `rk-discovery-terrain-lower-star-ab.json` |

The main run covers 52 configurations on the Apple M1 Max with native ARM Clang
15 and `-std=c++17 -O3 -DNDEBUG -pthread`. The candidate header digest is
`52d2a43d2330a2eaaf7f422e6d23589cb880ba7a127c9a36c07f60adb4fd57e8`.
Every sequence field and critical count agrees across revisions and worker
counts. All 15 multiworker plateau configurations have a below-one paired
median ratio; 12 have a paired-bootstrap interval wholly below one. Selected
eight-worker results are:

| Plateau input | Before (ms) | After (ms) | Paired after/before ratio |
| --- | ---: | ---: | ---: |
| 2D terrain, `n=32` | 4.554 | 2.935 | 0.683 |
| 3D volume, `n=8` | 2.910 | 2.389 | 0.788 |
| 3D volume, `n=12` | 11.326 | 9.777 | 0.895 |

The latter two main-run intervals include one. A fresh confirmation uses twelve
blocks of three repetitions: plateau volumes `8 12`, terrain `32`, seed `0`,
workers `1 8`, and all 32 ordinary lower-star controls above. Output basenames
replace `-ab.json` with `-confirmation.json`. Confirmed eight-worker results are:

| Plateau input | Before (ms) | After (ms) | Paired after/before ratio | 95% paired-block interval |
| --- | ---: | ---: | ---: | --- |
| 2D terrain, `n=32` | 4.172 | 3.093 | 0.807 | [0.651, 0.881] |
| 3D volume, `n=8` | 3.196 | 2.489 | 0.793 | [0.646, 0.921] |
| 3D volume, `n=12` | 9.866 | 9.192 | 0.923 | [0.807, 0.962] |

These are roughly 19%, 21% and 8% reductions in paired median algorithm time,
respectively, not universal speedup factors. Times are medians over raw samples;
ratios are medians over paired block ratios and need not equal ratios of the
displayed medians. All 38 confirmation configurations preserve exact gradients.
Construction/loading are excluded from both versions equally and no diagnostic
samples enter these A/B estimates.

The main ordinary-volume controls initially suggested a 10--20% slowdown at
`n=32`. A dedicated twelve-block repeat at one/eight workers, both seeds,
(`rk-discovery-v1-volume-lower-star-confirmation.json`) did not reproduce that
shift: sequential ratios are 0.983/0.996 and eight-worker ratios 1.038/1.009.
The broader confirmation also does not reproduce it. No ordinary configuration
has an above-one paired-bootstrap interval in both main and broad confirmation
sessions. Eight-worker median control ratios are 1.073 / 1.011 for volumes and
1.003 / 0.956 for terrains (main / broad confirmation). This supports retaining
the change with a variability caveat, not a proof of performance equivalence,
a general lower-star improvement, or an established compiler-inlining diagnosis.
All 94 A/B configurations, including the dedicated repeat, pass exact checks.

A separate candidate phase profile covers the same 13 inputs at 1/2/4/8 workers,
with 24 uninstrumented, eight coarse and three detailed samples per configuration.
All 52 configurations preserve the exact sequence. Compared with the historical
`ec408dc` profile, local facet-task counts, facet-kernel counts, incidence visits,
and discovery coface/mask work are identical. Discovery task counts at eight
workers change as follows:

| Plateau input | Previous discovery tasks | New discovery tasks |
| --- | ---: | ---: |
| 2D terrain, `n=16` | 144 | 0 |
| 2D terrain, `n=32` | 272 | 0 |
| 3D volume, `n=4` | 64 | 0 |
| 3D volume, `n=8` | 112 | 4 |
| 3D volume, `n=12` | 160 | 38 |

On the larger terrain, discovery is about 12% of detailed level wall time and
facet execution about 70%; on the larger volume these shares are about 28%
and 36%. Diagnostic timings are separate snapshots, not the basis of the A/B
speedup estimate. In particular, the larger volume's detailed discovery time
does not fall relative to the earlier snapshot despite fewer tasks; profiling
overhead and session variability must not be mistaken for a causal phase result.

Parallel RK is **not consistently faster than sequential RK on plateaus**. The
confirmed larger volume takes 9.192 ms with eight workers versus 8.782 ms
sequentially, while the separate phase run gives 8.061 versus 8.531 ms. The
larger terrain remains slower in parallel. Ordinary `n=32` volumes retain about
2.46-fold eight-worker speedup in the phase run, with builder initialization
about 30% of coarse total time. Further local-facet granularity and builder work
remain distinct possible targets; they are not changed here.

```sh
LC_ALL=C python3 tools/benchmark_reduction_kernel_phases.py \
  --input-dir ../rk-ab-inputs --output ../rk-phases-discovery.json

python3 tools/render_reduction_kernel_phases.py \
  --input ../rk-phases-discovery.json \
  --table-output docs/reduction_kernel_discovery_phases_table.tex
```

### F-Max diagnostic timer correction

Final integration checks exposed an existing overlap between F-Max's emission
and callback timers: the callback interval began before emission's final clock
read. A scheduling pause between those reads could be charged to both phases.
A separate profiling-only correction closes emission before opening callback
timing, for both critical insertions and regular pairs. It changes no gradient
choices or callbacks and leaves any timer bookkeeping gap in the residual.
Repeated tiny native diagnostics retain the strict no-double-counting check.

The RK A/B evidence above belongs to `3c18fa7`. Building its unprofiled RK driver
before and after this F-Max correction gives **identical executable bytes**
(SHA-256 `82cbc63155869348b18aff704560b00ae9c34e567799572d06bd355f79664910`),
so the correction does not alter that measured RK executable. Historical
non-profiled comparisons do not depend on the overlapping counters. Older
F-Max fine diagnostic components, however, should not be treated as a strictly
disjoint partition even when their sum happens to fit inside the outer total.

## Facet-Task Batching Update (Historical `ec408dc`)

Following the `f466631` incidence update, intra-level RK now submits at most one
long-lived facet task per configured worker in each round. Tasks claim chunks
of consecutive facets and write disjoint preallocated result slots; merging
still consumes results in canonical order. Chunk size is
`max(1, facet_count / task_count / 4)` with integer division. This removes the
individual facet futures and the barrier after each worker-sized wave, without
changing local kernel choices, sequential execution, independent-level
scheduling, or the global worker budget. Submitted tasks are drained before
any local-kernel or submission exception propagates.

The diagnostic `parallel_batches` counter now counts parallel facet rounds,
not the old worker-sized waves. The new `facet_parallel_tasks` counts submitted
facet-worker tasks; it is at most workers times parallel rounds. These counters
are exposed through the C++ and native Python profiling APIs. Coarse profiling
and uninstrumented execution do not collect detailed facet counters.

The controlled comparison uses baseline `f466631`, eight alternating blocks
of three uninstrumented repetitions and two warmups, on the Apple M1 Max with
native ARM Clang 15 and `-std=c++17 -O3 -DNDEBUG -pthread`. Timing includes a
fresh builder, workspace, worker pool, replay, and internal teardown; topology,
validation, protocol I/O, and returned-sequence destruction are excluded.
There is no persistence computation. The candidate header digest recorded by
the harness is `ee0fc1448aa1296752f589af1c9c1a4838696d742fcf9847a3d5fda577d74676`.

```sh
LC_ALL=C python3 tools/benchmark_reduction_kernel_ab.py \
  --baseline f466631 --candidate WORKTREE \
  --family volume --filtration plateau --sizes 4 8 12 --seeds 0 \
  --workers 1 2 4 8 --blocks 8 --repeats 3 --warmups 2 \
  --input-dir ../rk-ab-inputs \
  --output ../rk-facet-batch-volume-plateau-ab.json
```

Complementary runs use the same options with these substitutions:

| Family | Filtration | Sizes | Seeds | Output basename |
| --- | --- | --- | --- | --- |
| terrain | plateau | 16, 32 | 0 | `rk-facet-batch-terrain-plateau-ab.json` |
| volume | lower-star | 16, 32 | 0, 2 | `rk-facet-batch-volume-lower-star-ab.json` |
| terrain | lower-star | 16, 64 | 0, 2 | `rk-facet-batch-terrain-lower-star-ab.json` |

All 52 configurations preserve every sequence field and critical count across
revisions and worker counts. Selected eight-worker plateau results are:

| Input | Before (ms) | After (ms) | Paired after/before ratio |
| --- | ---: | ---: | ---: |
| 2D terrain, `n=32` | 66.49 | 4.49 | 0.0749 |
| 3D volume, `n=8` | 26.61 | 2.72 | 0.0957 |
| 3D volume, `n=12` | 122.61 | 10.16 | 0.0829 |

Times are medians of all raw samples; ratios are medians of paired block
ratios, so they need not equal ratios of the displayed times. All 15
multiworker plateau configurations improve, with paired-bootstrap intervals
below one. For the largest volume, two/four/eight-worker paired ratios are
0.200, 0.097, and 0.083. The sequential candidate median is 8.53 ms: batching
greatly reduces the parallel penalty but **does not establish an advantage
over sequential RK on these plateaus**.

A fresh confirmation uses twelve blocks of three repetitions: plateau volumes
`8 12`, terrain `32`, seed `0`, workers `1 8`, and all ordinary lower-star
controls above. Output basenames replace `-ab.json` with `-confirmation.json`.
All 38 additional configurations preserve exact gradients. Confirmed
eight-worker paired ratios are 0.0979 and 0.0610 for the two volumes, and
0.0599 for the terrain; each is below one in all twelve blocks. The confirmed
largest volume takes 8.43 ms at eight workers versus 7.43 ms sequentially.

Ordinary lower-star controls show mixed timing shifts. Median eight-worker
paired ratios across inputs are 1.003 initially / 0.974 in confirmation for
volumes and 1.018 / 1.013 for terrains. No ordinary configuration has an
above-one paired-bootstrap interval in both sessions. Two appear only in
confirmation (volume `n=32`, seed 0, four workers; terrain `n=16`, seed 0,
two workers). This is a nonregression check, not proof of performance
equivalence or a general lower-star speedup. Session-level bootstrap intervals
do not cover machine or between-session variability.

A separate candidate phase run covers 13 inputs at 1/2/4/8 workers, retaining
24 uninstrumented, eight coarse, and three detailed samples per configuration,
plus warmups. All 52 configurations preserve the exact gradient. The largest
volume still evaluates 40,806 facet kernels and visits 580,792 incidence
entries, but at eight workers now submits only 143 facet tasks in 18 parallel
rounds. The historical implementation used 5,108 worker-sized waves. The
`n=32` terrain similarly drops from 2,744 waves to 249 tasks in 32 rounds,
with its 21,854 facet kernels unchanged.

For the largest volume, eight-worker detailed facet execution including
dispatch/wait is about 2.62 ms and 29.6% of level wall time, versus about
139 ms / 95.4% in the earlier diagnostic snapshot. Facet discovery accounts
for 28.1%, closure construction 23.7%, and incidence 13.1%. For the terrain,
facet execution is 49.9% and facet discovery 35.4%. These are median per-run
elapsed shares on single-level inputs; cumulative local worker times must not
be added to them. The fresh A/B tests above, not differences between separate
diagnostic sessions, support the performance claim.

The phase run's uninstrumented largest-volume medians are 7.63, 8.74, 8.00,
and 8.67 ms at 1/2/4/8 workers. More workers still do not guarantee a faster
plateau computation. Ordinary `n=32` lower stars retain about 2.81-fold
eight-worker speedup in this separate phase run, with builder initialization
about 31% of total coarse time. Phase-specific task granularity, especially
facet discovery on modest active buckets, is the next plateau target; generic
builder initialization remains separate ordinary-lower-star work. The discovery
scheduling update above addresses the former target; these batching measurements
remain the historical `ec408dc` snapshot.

```sh
LC_ALL=C python3 tools/benchmark_reduction_kernel_phases.py \
  --input-dir ../rk-ab-inputs --output ../rk-phases-facet-batch.json

python3 tools/render_reduction_kernel_phases.py \
  --input ../rk-phases-facet-batch.json \
  --table-output docs/reduction_kernel_facet_batch_phase_table.tex
```

The earlier phase tables remain historical. All raw timing samples, input and
header hashes, intervals, counters, and exact-check results remain in the
local JSON files under the repository output policy.

Unit tests check exact gradients under uninstrumented, coarse, and detailed
execution at 1/2/4/8 workers, with uneven chunk tails, mixed facet sizes,
dimensions through five, inline/event overflow, and shared sparse faces.
Structural bounds ensure at most one facet task per worker per parallel round.
A throwing complex-view test verifies that all submitted facet tasks finish
before the exception returns, with the executor still alive and reusable.
Existing packed/sparse boundary tests and independent-level checks remain.

## Closure-Based Incidence Update (Historical: `f466631`)

This update implements the incidence optimization identified by the `8a2bd06`
profile below. Sparse levels now accumulate saturated incidence by traversing
cached same-level facet closures for every execution policy. This removes the
parallel all-pairs simplex/facet scan. Levels containing only vertices and
edges enumerate their closures directly. Incidence remains a coordinator step
before facet tasks begin; independent levels and facet reductions retain their
parallel execution. Per-facet scheduling and generic builder initialization
are deliberately unchanged in this isolated optimization.

The controlled A/B comparison uses `8a2bd06` as the baseline, the same native
driver against both header snapshots, alternating execution order, eight
blocks of three repetitions, and two warmups. It includes fresh builder,
workspace, worker-pool setup/teardown, and replay, without instrumentation or
persistence; topology, validation, protocol I/O, and returned-sequence
destruction remain outside timing.

```sh
LC_ALL=C python3 tools/benchmark_reduction_kernel_ab.py \
  --baseline 8a2bd06 --candidate f466631 \
  --family volume --filtration plateau --sizes 4 8 12 --seeds 0 \
  --workers 1 2 4 8 --blocks 8 --repeats 3 --warmups 2 \
  --input-dir ../rk-ab-inputs \
  --output ../rk-closure-incidence-volume-plateau-ab.json
```

The complementary runs use the same options with these substitutions:

| Family | Filtration | Sizes | Seeds | Output basename |
| --- | --- | --- | --- | --- |
| terrain | plateau | 16, 32 | 0 | `rk-closure-incidence-terrain-plateau-ab.json` |
| volume | lower-star | 16, 32 | 0, 2 | `rk-closure-incidence-volume-lower-star-ab.json` |
| terrain | lower-star | 16, 64 | 0, 2 | `rk-closure-incidence-terrain-lower-star-ab.json` |

The four runs total 52 worker configurations on the Apple M1 Max using native
ARM Clang 15, C++17, and `-O3 -DNDEBUG -pthread`. All exact sequences and
critical counts agree across revisions and worker counts. Selected larger
plateau timings at eight workers, in milliseconds, are:

| Plateau | Before | After |
| --- | ---: | ---: |
| 2D terrain, `n=32` | 153.63 | 61.70 |
| 3D volume, `n=8` | 86.97 | 24.60 |
| 3D volume, `n=12` | 1,462.70 | 130.56 |

The columns are medians of raw uninstrumented timings. For the largest volume,
median paired candidate/baseline ratios are 0.0155, 0.0490, and 0.0895 at two,
four, and eight workers, respectively. Its sequential ratio is 0.987, with
the paired-bootstrap interval spanning one. The eight-worker update is about
11 times faster than the old parallel path, but it is **still much slower
than sequential RK on this plateau** (candidate median 7.77 ms). Removing the
incidence problem does not remove fine-grained facet-task overhead.

A fresh confirmation uses twelve blocks of three repetitions. For plateaus it
selects volumes `8 12`, terrain `32`, seed `0`, and workers `1 8`; for ordinary
lower stars it repeats all configurations in the table above. Output names
replace `-ab.json` with `-confirmation.json`. These 38 additional configurations
also preserve every exact sequence. The confirmed eight-worker paired ratios
are 0.2725 and 0.0858 for the `n=8` and `n=12` volumes and 0.3732 for the
`n=32` terrain. Thus the large-plateau improvements reproduce independently.

The 32 ordinary lower-star controls in each session have mixed timing shifts.
At eight workers, median paired ratios for volumes are 0.970 initially and
0.984 in confirmation; for terrains they are 0.931 and 1.021. No ordinary
case has a paired-bootstrap interval wholly above one in both sessions.
This is a regression check, not a claim of a general lower-star speedup or
formal performance equivalence; absolute timings and small differences remain
sensitive to system state.

A separate `f466631` phase profile confirms that the all-pairs work is
gone. Every worker count now records identical sparse incidence-entry visit
counts. For the `n=12` volume this is 580,792 closure-entry visits, the same as
the sequential implementation, replacing the old parallel path's 663,665,205
containment tests. Eight-worker diagnostic incidence time is about 1.33 ms.
Facet execution (including dispatch/wait) takes about 139 ms, or 95.4% of the
level-processing interval. These are separate diagnostic measurements, not
the uninstrumented timings in the A/B table. This motivated the facet-task
batching update above; builder initialization remains a separate target for
ordinary lower stars. Check out `f466631` before reproducing this historical
phase snapshot:

```sh
LC_ALL=C python3 tools/benchmark_reduction_kernel_phases.py \
  --input-dir ../rk-ab-inputs --output ../rk-phases-closure-incidence.json

python3 tools/render_reduction_kernel_phases.py \
  --input ../rk-phases-closure-incidence.json \
  --table-output docs/reduction_kernel_incidence_phase_table.tex
```

This historical phase run contains 13 inputs at 1/2/4/8 workers, with separate
uninstrumented and diagnostic samples and exact sequence checks throughout.
The preceding `docs/reduction_kernel_phase_table.tex` remains the historical
pre-incidence table. Unit tests additionally bound incidence-entry visits by
facet-closure size for sparse graph, triangle, and dimension-four inputs,
including multiple levels, shared faces with incidence greater than two, and
isolated vertices. Existing packed/sparse boundary tests remain in place.

Raw paired-block samples, IQRs, bootstrap intervals, input/header hashes, and
exact-check results remain in the JSON files. A bootstrap interval describes
one session; it does not capture independent-session or machine variability.

## ReductionKernel Phase Profile Before Closure-Based Incidence

The September 9, 2026 profile of the packed-coface algorithm (`e1aaeaa`, with
profiling-only additions) identifies two different optimization targets:
large single plateaus suffer from the parallel incidence algorithm, while
ordinary lower stars increasingly expose serial builder initialization.
These measurements concern gradient construction only.

```sh
python3 tools/benchmark_reduction_kernel_phases.py \
  --terrain-sizes 16 64 --volume-sizes 16 32 \
  --plateau-terrain-sizes 16 32 --plateau-volume-sizes 4 8 12 \
  --seeds 0 2 --workers 1 2 4 8 \
  --blocks 8 --repeats 3 --coarse-blocks 8 --detailed-blocks 3 --warmups 2 \
  --input-dir ../rk-ab-inputs --output ../rk-phases-cofaces-1.json

python3 tools/render_reduction_kernel_phases.py \
  --input ../rk-phases-cofaces-1.json \
  --table-output docs/reduction_kernel_phase_table.tex
```

The run covers 13 complexes and 52 worker configurations on an Apple M1 Max
(eight performance cores, two efficiency cores), using native ARM Clang 15,
C++17, `-O3 -DNDEBUG -pthread`. Ordinary lower-star inputs use both seeds;
plateaus use the first seed's topology with all vertex values set to zero.
Each configuration retains 24 uninstrumented samples, eight coarse samples,
and three detailed samples, plus separate warmups. Worker-count order rotates
and reverses between blocks. All input generation and compilation finish before
timing. Every timed and diagnostic sequence agrees field-for-field with a
validated sequential sequence on the same input. No persistence is computed.

The timer includes a fresh builder, workspace, worker pool, gradient, ordered
replay, and internal teardown. Topology construction, protocol I/O, validation,
and destruction of the returned sequence are excluded. Unlike the direct TTK
table, builder construction is included here, so times from the two tables
must not be substituted for each other.

Coarse profiling disables local diagnostic timers and counters, retaining
outer wall-clock phases and elapsed duration/load counters for long-lived
level tasks. Detailed runs separately record mask/closure construction,
facet discovery, protected-face incidence, local reductions, facet execution
including dispatch/wait, aggregation, and merging. The JSON preserves every
sample, quartiles, input/header/driver hashes, and exact-check status. The field
`overhead_ratio` is only the ratio of diagnostic to uninstrumented medians from
separate run groups, **not** a causal estimate of instrumentation overhead.
Reported phase shares use each coarse run's own measured total. Residual time
includes workspace/builder destruction and worker-pool shutdown. Detailed
subphase durations accumulate across levels/workers, and core/local durations
nest inside facet execution; they cannot be summed with outer wall times.

At the largest ordinary 3D size (`n=32`, 792,051 simplices), median time across
the two input medians is 60.44 ms sequentially and 22.87 ms at eight workers;
the median per-input speedup is 2.64. At eight workers, coarse phase shares
are approximately 30.7% builder initialization, 5.2% setup, 52.1% level
processing, 6.9% replay, and 3.8% residual. Shares are separately aggregated
medians and need not sum exactly to 100%. Elapsed level-task overlap is
7.43 and 7.66 for the two seeds at eight workers. This is **not CPU utilization**
and does not separate memory contention from OS scheduling, but gives no
reason to prioritize another level-load-balancing change. Approximately 99.8%
of simplices use packed-eligible levels, with maximum bucket size 75.
Small workloads remain sensitive to worker-pool overhead: the `n=16` ordinary
terrains take roughly 0.06--0.07 ms sequentially and 0.27--0.29 ms at eight
workers. Increasing the worker budget is therefore not always beneficial.

The constant-field plateaus exercise the sparse path and spatial parallelism.
For the `n=32` terrain (5,891 simplices), time rises from 1.38 ms sequentially
to 155.69 ms at eight workers. For the `n=12` volume (36,851 simplices), it rises
from 8.56 ms to 1,501.23 ms. Thus the earlier lower-star results must not be
generalized to large plateaus. These large slowdowns are not gradient errors:
the sequences are identical at every worker count.

Code inspection and diagnostic counters identify the mechanism. In
`ReductionKernelWorkspace::compute_facet_incidence`, the sequential sparse
path walks cached facet closures, but the intra-level parallel path tests
active simplices against current facets. On the `n=12` volume this changes
580,792 closure-entry visits into 663,665,205 simplex/facet containment tests.
These count different operations, not interchangeable units. Protected-face
incidence alone accounts for about 1,295 ms of the detailed eight-worker
run. On the smaller `n=8` volume, the corresponding counts are 104,650 and
35,725,927; incidence takes about 61 ms and facet execution about 29 ms.
The latter includes thousands of individually submitted facet tasks; this is
a second cost, distinct from the all-pairs scan.

A fresh six-complex, 24-configuration confirmation uses the same command with
`--terrain-sizes 64 --volume-sizes 32 --plateau-terrain-sizes 32
--plateau-volume-sizes 8` and output
`../rk-phases-cofaces-confirmation.json`. All exact checks pass again. The
ordinary `n=32` volume speedup is 2.62 and the eight-worker builder share is
29.7%. The plateau terrain takes 1.34 ms sequentially versus 156.62 ms at eight
workers; the `n=8` plateau volume takes 1.50 ms versus 91.43 ms. The large
`n=12` plateau is not included in this independent confirmation. Raw timing
variability remains available; these are one-machine observations, not
cross-machine confidence bounds or proof of the gain from an unimplemented
optimization.

An additional uninstrumented A/B control against `e1aaeaa` checks the profiling
additions themselves. The first eight-configuration run and a fresh
four-configuration repeat do not show a repeatable performance shift: for
example, the large seed-0 sequential candidate/baseline ratio changes from
1.138 to 0.917 between runs. All exact comparisons pass. These controls are
stored in `../rk-phases-uninstrumented-control.json` and
`../rk-phases-uninstrumented-control-2.json`. They are not evidence for a
speed improvement or a formal equivalence margin; they reinforce the need
for independent sessions alongside within-session timing distributions.

The profile recommended **removing the parallel all-pairs incidence
scan**, using cached closures for all execution policies while retaining
saturated counts (zero, one, or multiple incident facets). A simple linear
closure traversal is an important baseline; any parallel version should use
race-free local accumulation and deterministic merging. Then address
per-facet task batching if it remains material. For ordinary lower stars,
inspect whether the generic builder's rank/level/dimension tables can be
avoided or deferred for RK without weakening input validation. Ordered replay
is not the first target in this profile. None of these algorithmic changes is
included in the profiling-only update. The closure-based incidence update
above subsequently implements the first recommendation; each change needs a
controlled A/B test and exact-gradient checks on both packed levels and large
sparse plateaus.

## Gradient-only Tetrahedral Phase Profile (Historical)

This benchmark times only construction of the discrete gradient. It invokes
`FSequenceBuilder` with a no-op callback: no reference map, reduction plan, or
persistence computation occurs in the measured path. The native profiler
separates ProcessLowerStars into builder initialization, global lower-star
setup, local-star processing, and ordered event replay. It also records
critical-simplex counts by dimension and checks every parallel gradient against
the corresponding sequential result.

```sh
PYTHONPATH=python python3 tools/benchmark_tetrahedral_phase_profile.py \
  --sizes 4 8 12 16 \
  --seeds 0 1 2 \
  --workers 1 8 \
  --repeats 5 \
  --format csv \
  --output ../work/tetrahedral_phase_profile.csv

python3 tools/render_tetrahedral_phase_profile.py \
  --input ../work/tetrahedral_phase_profile.csv \
  --table-output docs/tetrahedral_phase_profile_table.tex
```

Across the twelve cases, eight-worker ProcessLowerStars reaches a 2.23-fold
gradient-construction speedup; ReductionKernel reaches 1.90-fold. Their median
eight-worker construction times are 4.53 ms and 1.09 ms, respectively. At
eight workers, ProcessLowerStars spends 28.9 percent in global setup, 42.9
percent in parallel local-star processing, 3.4 percent in ordered replay, and
4.5 percent in builder initialization. ReductionKernel spends 82.4 percent of
its diagnostic wall time processing levels; setup and replay account for 13.2
and 3.1 percent. The CSV also records ReductionKernel chunk counts, per-worker
level and simplex loads, task time, and effective task parallelism, alongside
kernel rounds, facet kernels, topology scans, membership tests, and inline
buffer overflows. All measured triangular and
tetrahedral kernels report zero inline-buffer overflows. For the `n=16`
volumes, the medians include roughly
458,000 local candidate visits but only 70,000 membership tests, confirming
that preserving the short-circuit local scan is preferable to eagerly building
all local incidences. These figures deliberately make no claim about
persistence performance.

The phase and concurrency aggregates are generated in
`docs/tetrahedral_phase_profile_table.tex`.

## Reduction-Kernel Scaling

The optimized parallel scheduler launches one long-lived task per worker. Each
task dynamically claims a small chunk of filtration levels from an atomic
counter, amortizing counter contention while retaining enough chunks for load
balancing. This avoids both thousands of tiny executor tasks and the former
sorting and static simplex-count partition. Each task retains its facet flags,
closure tables, incidence counters, and result buffers between claimed levels.
It also writes into a preassigned slice of the shared event arena. It disables
nested facet tasks while multiple levels are running concurrently. A single
large plateau still uses the facet-parallel reduction-kernel path.

The diagnostic profile records claimed chunks, processed levels and simplices,
and task durations per worker. Two- and four-worker runs are generally well
balanced, while eight-worker runs vary substantially with operating-system
scheduling. Experiments with a start rendezvous and reserved initial chunks
produced contradictory repeat measurements and were therefore not retained.
The remaining eight-worker gap cannot yet be attributed solely to either load
imbalance or level-kernel cost.

```sh
PYTHONPATH=python python3 tools/benchmark_reduction_kernel_scaling.py \
  --sizes 16 32 64 \
  --seeds 0 1 2 \
  --workers 1 2 4 8 \
  --repeats 5 \
  --warmups 1 \
  --format csv \
  --output ../work/reduction_kernel_scaling.csv

MPLCONFIGDIR=../work/matplotlib-cache \
  python3 tools/render_reduction_kernel_scaling.py \
  --input ../work/reduction_kernel_scaling.csv \
  --figure-output docs/reduction_kernel_scaling.svg \
  --table-output docs/reduction_kernel_scaling_table.tex
```

This older Python-level pipeline benchmark includes sequence-object
materialization, reference maps, and persistence, so it is retained only as an
overhead diagnostic rather than evidence for gradient speed. Across the nine
terrains, median native-sequence speedups are 1.05x, 1.11x, and 1.07x at two,
four, and eight workers; median end-to-end speedup is 1.02x at eight workers.
Every worker count is checked for the exact sequential reduction-kernel
sequence and the standard barcode. The gradient-only results above are the
relevant measurements for the present study.

![Reduction-kernel scaling](reduction_kernel_scaling.svg)

The grid-size aggregates are generated in
`docs/reduction_kernel_scaling_table.tex`.

## Roadmap and External Data

The benchmark runner also has Roadmap and CAM-style families:

```text
cam-s4-rips
roadmap-rips
```

Roadmap datasets are cached under `../work/roadmap-data` by default. Missing
Roadmap files are not downloaded unless requested explicitly:

```sh
PYTHONPATH=python python3 tools/benchmark_persistence.py \
  --preset roadmap \
  --sequence-algorithm portfolio \
  --download-roadmap-data \
  --format csv \
  --output ../work/roadmap_portfolio.csv
```

Use this only when network access is acceptable.

## Native GUDHI-View Benchmark

The native GUDHI benchmark compares three in-process paths on the same
`Gudhi::Simplex_tree` input:

- `Direct`: MorseFrames through a read-only `Simplex_tree` view.
- `Import`: copy into the compact owning MorseFrames complex first.
- `GUDHI`: GUDHI persistent cohomology on the original `Simplex_tree`.

This benchmark is optional because it needs GUDHI and Boost headers. Configure
them explicitly when CMake cannot find them:

```sh
cmake -S . -B build-gudhi \
  -DMORSEFRAMES_GUDHI_INCLUDE_DIR=/path/to/gudhi/include \
  -DMORSEFRAMES_BOOST_INCLUDE_DIR=/path/to/boost/include

cmake --build build-gudhi --target morseframes_benchmark_gudhi_view
```

Quick run:

```sh
mkdir -p ../work
./build-gudhi/morseframes_benchmark_gudhi_view \
  --quick \
  --repeats 3 \
  > ../work/native_gudhi_view_quick.csv

PYTHONPATH=python python3 tools/render_native_gudhi_view_table.py \
  --input ../work/native_gudhi_view_quick.csv \
  --output docs/native_gudhi_view_quick_table.tex \
  --summary

PYTHONPATH=python python3 tools/render_native_gudhi_stage_profile.py \
  --input ../work/native_gudhi_view_quick.csv \
  --table-output docs/native_gudhi_stage_profile_quick_table.tex \
  --prose-output ../work/native_gudhi_stage_profile_quick_prose.tex \
  --summary
```

Default-size repeat run:

```sh
./build-gudhi/morseframes_benchmark_gudhi_view \
  --repeats 30 \
  > ../work/native_gudhi_view_default_r30.csv

PYTHONPATH=python python3 tools/render_native_gudhi_view_table.py \
  --input ../work/native_gudhi_view_default_r30.csv \
  --output docs/native_gudhi_view_default_r30_table.tex \
  --caption-title "Native \\texttt{Gudhi::Simplex\\_tree} default benchmark." \
  --label tab:native-gudhi-view-default-r30 \
  --summary
```

Larger lean run:

```sh
./build-gudhi/morseframes_benchmark_gudhi_view \
  --large \
  --lean \
  --repeats 30 \
  > ../work/native_gudhi_large_lean_r30.csv

PYTHONPATH=python python3 tools/render_native_gudhi_view_table.py \
  --input ../work/native_gudhi_large_lean_r30.csv \
  --output docs/native_gudhi_large_lean_r30_table.tex \
  --caption-title "Native \\texttt{Gudhi::Simplex\\_tree} larger lean benchmark." \
  --label tab:native-gudhi-large-lean-r30 \
  --summary
```

In these tables, `GUDHI/Direct < 1` means GUDHI is faster end-to-end, while
`GUDHI/Reducer > 1` means the Morse reducer kernel alone is faster than GUDHI
persistence after the Morse input has already been built.

## Prime-Field Overhead

Prime-field coefficient experiments are generated by
`tools/benchmark_prime_field_overhead.py`.

Quick local run:

```sh
mkdir -p ../work
PYTHONPATH=python python3 tools/benchmark_prime_field_overhead.py \
  --families lower-star plateau rips \
  --sizes 8 12 16 \
  --seeds 0 1 \
  --algorithms saturated f-max same-level-reduction \
  --primes 3 5 \
  --repeats 5 \
  --output-csv ../work/prime_field_overhead_quick.csv \
  --output-md ../work/prime_field_overhead_quick.md
```

Composite moduli are intentionally rejected by the barcode API; these reducers
work over fields `F_p`.

## Profile-Selection Validation

The profile-selection scripts compare cheap strategy-selection metrics against
measured portfolio timings. These runs are more expensive than the smoke tests.

Preview the commands without executing them:

```sh
PYTHONPATH=python python3 tools/run_fair_profile_validation.py \
  --validation-preset report \
  --dry-run
```

Regenerate the public validation table from fresh timings:

```sh
mkdir -p ../work
PYTHONPATH=python python3 tools/run_fair_profile_validation.py \
  --validation-preset report \
  --output-dir ../work \
  --table-output docs/profile_metric_fair_validation_table.tex \
  --prose-output ../work/profile_metric_fair_validation_prose.tex \
  --manifest-output ../work/fair_profile_validation_manifest.md
```

If CSVs already exist in `../work`, summaries can be regenerated without
rerunning timings:

```sh
PYTHONPATH=python python3 tools/run_fair_profile_validation.py \
  --validation-preset report \
  --output-dir ../work \
  --summaries-only \
  --table-output docs/profile_metric_fair_validation_table.tex \
  --prose-output ../work/profile_metric_fair_validation_prose.tex \
  --manifest-output ../work/fair_profile_validation_manifest.md
```

Selector decision and feature diagnostic tables are rendered from the validation
CSVs:

```sh
PYTHONPATH=python python3 tools/summarize_selector_decisions.py \
  --table-output ../work/profile_selector_decision_summary.txt \
  --csv-output ../work/profile_selector_decision_summary.csv \
  --latex-output docs/profile_selector_decision_summary_table.tex \
  --prose-output ../work/profile_selector_decision_summary_prose.tex

PYTHONPATH=python python3 tools/analyze_selector_features.py \
  --table-output ../work/selector_feature_diagnostic.txt \
  --csv-output ../work/selector_feature_diagnostic.csv \
  --json-output ../work/selector_feature_diagnostic.json \
  --latex-output docs/selector_feature_diagnostic_table.tex \
  --prose-output ../work/selector_feature_diagnostic_prose.tex
```

## Benchmark Summary Page

The visible compact tables in `docs/benchmark_summary.md` are generated from the
tracked LaTeX table fragments:

```sh
python3 tools/render_benchmark_summary.py
```

CI checks that this generated block is up to date:

```sh
python3 tools/render_benchmark_summary.py --check
```

## Before Committing Regenerated Results

Before committing regenerated table fragments, run:

```sh
git diff -- docs tools benchmarks
git diff --check
python3 tools/render_benchmark_summary.py --check
MORSEFRAMES_DISABLE_CPP_BACKEND=1 \
  python3 -m unittest discover -s python/tests -p "test_*.py"
```

Commit only public artifacts that are meant to be reproducible from this
repository. Keep manuscript text, discussion packages, generated prose, and PDFs
in the private notes repository.
