# Reusing simplex order during finalization

Revision `ba57c79` speeds up shared native construction by reusing information
already encoded in lexicographic simplex IDs. Across two controlled studies,
**39/40 paired construction medians improve**. All 16 intervals for volume n=16/32
(two seeds, two worker settings, two studies) are below one. Smaller-case
results are less conclusive. Main largest-volume construction falls from
**662.947 to 588.546 ms**, with peak memory essentially unchanged.
This is a construction improvement, not a claim
of faster gradient matching or full persistence.

This follows the [compact lookup study](native_finalization_benchmark.md).
Loading and construction remain separately reported and excluded from the
primary gradient comparison, as agreed. Earlier snapshots are preserved.

## Implementation and correctness

The change is in `include/morseframes/filtered_complex.hpp`:

- Boundary construction tracks the current first-vertex range while walking
  lexicographic IDs. Deleting any vertex except the first keeps the face in that
  range, so its prefix-index search can be skipped. Deleting the first vertex
  still uses the general lookup.
- Filtration sorting still compares level, then dimension, then lexicographic
  simplex order. Its final comparison now uses IDs instead of vertex vectors;
  IDs are already positions in the lexicographically sorted simplex array.

No stored index, mutable cache or additional heap allocation is introduced.
Construction remains sequential. Every strategy using this owning complex can
benefit, including persistence callers; direct external complex views are
unchanged. The implementation supports arbitrary monotone filtrations and is
not restricted to vertex lower stars or dimensions two/three.

Simplex IDs, filtrations, levels, boundary deletion order, coboundary order,
filtration order and level buckets are unchanged. Missing faces and
non-monotone filtrations are still rejected. Independent map-based tests check
these contracts through dimension seven, with sparse vertex IDs up to the
maximum, mixed dimensions, plateaus, non-lower-star filtrations, shuffled
insertion, repeated finalization and every facet-deletion position. Existing
copy/move, concurrent const lookup and pending-insertion tests remain in place.

Full cross-version dumps compare complexes and sequential F-Max/RK sequences
byte for byte. Every measured run also checks fingerprints, including parallel
RK, and critical counts are unchanged. Local ordinary C++ and ASan/UBSan tests,
151 native Python tests (including both benchmark adapters), the fallback suite
(151 run, 10 skipped), and both nine-case RK corpora at workers 1/2/4 pass.

## Controlled construction comparison

Baseline: `4b9f1495e9053eccf8849f77fa0d2cd4b2bc50f2`.
Candidate: `ba57c79e8b8b71275d1d8a51596b1cf0446e56f5`.
Both use the same shared bulk helper, driver, compiler flags and input arrays.

Each study covers terrain n=16/64 and volume n=8/16/32, seeds 0/2 and worker
settings 1/8: ten inputs and 20 configurations. The repeat reverses both input
and worker order. Each configuration uses eight balanced version-order blocks,
two repetitions per version/block and two warmups. F-Max/RK execution order
also alternates. Each input/version has three separate diagnostics and three
fresh-process memory samples. Worker counts apply to the subsequent RK run,
not to the sequential constructor.

Both studies start clean on the native ARM Apple M1 Max, Apple Clang 15,
`-std=c++17 -O3 -DNDEBUG -pthread`. Builds and tests do not overlap measurements.
Other desktop applications remain active, and all outliers are retained.

The paired-block bootstrap interval is below one in **11/20** main and **14/20**
repeat configurations: **25/40 overall**. The other 15 intervals include one;
none is wholly above one. Only the main terrain n=16/seed 0/one-worker paired
median exceeds one (1.002). These are session-level intervals, not independent
machine uncertainty or universal speedup guarantees.

Main construction times (ms), selecting the eight-worker setting and taking
medians over the two seed-specific medians:

| Input | Previous | Order reuse | Paired new/previous |
| --- | ---: | ---: | ---: |
| 2D terrain, n=16 | 0.586 | 0.556 | 0.941 |
| 2D terrain, n=64 | 12.647 | 11.605 | 0.918 |
| 3D volume, n=8 | 6.235 | 5.701 | 0.879 |
| 3D volume, n=16 | 66.149 | 56.629 | 0.831 |
| 3D volume, n=32 | 662.947 | 588.546 | 0.895 |

Largest-volume repeat ratios range from 0.860 to 0.906; all four intervals
are below one. Main inconclusive intervals concern terrain n=16/64 and volume
n=8; repeat inconclusive intervals concern all four terrain n=16 configurations
and terrain n=64/seed 2 at both worker settings. Full intervals remain in the
raw files. Do not infer cumulative speedups by comparing absolute times from
different historical studies.

## Phase attribution and memory

At the main n=32/eight-worker setting, uninstrumented finalization falls from
**463.715 to 379.981 ms**. Enumeration/insertion, whose code is unchanged, is
**197.025 / 200.091 ms**. These separately aggregated medians need not sum to
the construction medians.

Separate n=32 diagnostic phase times (ms; medians within seed, then across seeds):

| Phase | Previous | Order reuse |
| --- | ---: | ---: |
| Index and simplex records | 32.384 | 31.000 |
| Filtration levels | 58.346 | 57.789 |
| Boundary lookup and checks | 165.324 | 102.077 |
| Coboundaries | 65.267 | 65.487 |
| Ordering and level buckets | 124.879 | 103.186 |

The reductions are concentrated in the two targeted phases. This is a combined
change, not a separate ablation of each optimization. Diagnostic timings are
not performance samples; changes in untouched phases should not be attributed
to new code. Reset and bulk subphases remain in the raw evidence.

Main-study fresh-process peak RSS (MiB; medians over seed medians):

| Input | Previous | Order reuse |
| --- | ---: | ---: |
| 2D terrain, n=16 | 1.648 | 1.617 |
| 2D terrain, n=64 | 7.477 | 7.477 |
| 3D volume, n=8 | 4.000 | 4.031 |
| 3D volume, n=16 | 25.469 | 25.531 |
| 3D volume, n=32 | 195.961 | 195.859 |

There is no material memory saving to claim. These peaks are measured before
gradients and validation and include the runtime, resident inputs, allocator
and construction temporaries. They are not live heap sizes.

## Gradient timing caveats

Gradient code and exact outputs are unchanged, but timing equivalence is not
established:

- Main volume n=32/seed 0/one-worker F-Max: paired ratio **1.097**, interval
  **[1.029, 1.143]**, medians **132.071 → 144.456 ms**. The repeat ratio is
  **1.017**, interval **[0.884, 1.081]**.
- Main volume n=32/seed 0/eight-worker RK: ratio **1.095**, interval
  **[1.014, 1.195]**, medians **13.293 → 14.156 ms**. The repeat ratio is
  **1.068**, interval **[0.988, 1.216]**. Its repeat median still points toward
  a slowdown; an interval including one is not evidence of no effect.
- Repeat terrain n=16/seed 0/eight-worker F-Max: ratio **1.125**, interval
  **[1.044, 1.193]**, medians **0.081 → 0.083 ms**. It was not an above-one
  interval in the main study. The paired ratio is not the quotient of those
  medians, especially after rounding.

No identical gradient configuration has an above-one interval in both studies.
Nevertheless, neither these runs nor the fresh comparison below prove the
absence of a gradient regression. Desktop activity and session variation are
limitations; no particular scheduling, allocation or code-layout cause has
been isolated. Retain the earlier study's timing caveats as historical evidence
too. The supported result is improved construction, not uniformly faster
gradient execution.

## Refreshed gradient comparison

The clean-start resident-array study at `ba57c79` covers terrain n=64, volume
n=8/16/32, seeds 0/2 and workers 1/8: eight inputs and 16 configurations.
It uses 12 performance repetitions, six separate diagnostics and two warmups
per mode, with balanced six-way algorithm order. TTK stays at clean pinned
revision `f4ffd1a1049d0ccf6e8f3eb4f7c096a6cc251ba0`, classic backend, native
ARM release, OpenMP and `OMP_WAIT_POLICY=PASSIVE`.

Eight-worker algorithm times (ms; medians over seed-specific medians):

| Input | Sequential F-Max | RK | TTK |
| --- | ---: | ---: | ---: |
| 2D terrain, n=64 | 1.476 | 0.757 | 0.721 |
| 3D volume, n=8 | 0.807 | 0.460 | 0.605 |
| 3D volume, n=16 | 11.379 | 2.012 | 3.217 |
| 3D volume, n=32 | 134.542 | 15.053 | 23.796 |

RK leads both competitors by paired-median algorithm ratios in all 12 volume
configurations, including one-worker cases. TTK leads the eight-worker terrain
aggregate. Critical counts agree across all three methods in all 16
configurations, and each method passes its exact reference checks. Dimensionwise
count vectors are retained per input; for n=32 they are **(128, 455, 379, 51)**
at seed 0 and **(152, 649, 570, 72)** at seed 2, ordered by dimension 0–3.
This is a fresh comparison, not a causal comparison against earlier absolute
gradient timings and not a replacement for the controlled warnings above.

Corresponding native construction costs:

| Input | F-Max construction (ms) | RK construction (ms) | TTK construction (ms) |
| --- | ---: | ---: | ---: |
| 2D terrain, n=64 | 11.293 | 11.233 | 0.550 |
| 3D volume, n=8 | 5.512 | 5.558 | 0.746 |
| 3D volume, n=16 | 58.863 | 58.823 | 4.850 |
| 3D volume, n=32 | 594.322 | 590.581 | 39.257 |

F-Max and RK use identical construction code; differences between their
construction measurements reflect variation. TTK retains the lower full
resident-to-gradient time in all 16 configurations. On n=32/eight workers,
full RK/TTK medians are **605.834 / 62.852 ms**.

Loading is reported separately as one shared read/parse/validation per native
invocation, not as repeated per-algorithm work or cold-cache I/O. Native
construction is excluded from the algorithm comparison. Fresh MorseFrames
builders and TTK vertex ordering/lower-star processing remain included in
algorithm time. Complete one/eight-worker algorithm, construction/loading and
outer-phase tables are `resident_gradient_ba57c79_algorithm_table.tex`,
`resident_gradient_ba57c79_construction_table.tex` and
`resident_gradient_ba57c79_phases_table.tex`.

## Reproduction and provenance

Raw files are preserved outside the repository beside the earlier studies:

- `../construction-order-main.json`: UTC 2026-09-10 12:10:30 to 12:14:38;
  SHA-256 `45d5dff20ed57acf79c6b4ea8c72cf7495f1c3318b3c83d79d8735219e0806a0`.
- `../construction-order-confirmation.json`: UTC 2026-09-10 12:14:51 to 12:19:05;
  SHA-256 `a2692aaa1a8c582c0aedbbce4bae57a3c82f388e5f8e406a2ea1dcaf318646e4`.
- `../resident-construction-order-ba57c79.json`: UTC 2026-09-10 12:19:25 to 12:22:42;
  SHA-256 `02e36d3ff0bbcec8b9bc1547b069474533a12de4c2f11ab628a4ff593f1271cb`.
- `../construction-order-pilot.json`: exploratory dirty-checkout pilot, not
  pooled with either controlled study.

The files retain all raw samples, input/revision/binary hashes, diagnostics,
critical counts and memory measurements. Build/measurement flags are unchanged.

- Baseline headers: `7b1309f114613cdf1025c862c82871855520f9098cd3b62ddf4538aa39721c5e`.
- Candidate headers: `f161815c05f2822aaacd821fbbc8fc2f0ad9753d2872f9b7694d69aca3bc30cf`.
- Construction driver: `a5359ed40f71fa68e7d7cc76eb3cdfc594437d476418e1ac40ad02623cf1eebb`.
- Construction runner: `4e0887fd2816a46249a2dc2eaef8d0b383a1e8ce621ac8e1727efacbb28a71e1`.
- Resident driver: `393c3af7f2df2d973031cc38ef7c5c44dcffecd11faa2434f8e03b97d9ab0bc8`.
- Resident executable: `091e314e16206697f4f5a2cd9face66128e437a7701d2fc1c5d45757b4b474b2`.

Build the pinned native resident driver and generate inputs as described in
[benchmark reproduction](benchmark_reproduction.md). From the repository root:

```sh
LC_ALL=C python3 tools/benchmark_complex_construction.py \
  --baseline 4b9f149 --candidate ba57c79 \
  --inputs ../rk-ab-inputs/terrain-n16-seed0.txt ../rk-ab-inputs/terrain-n16-seed2.txt \
    ../rk-ab-inputs/terrain-n64-seed0.txt ../rk-ab-inputs/terrain-n64-seed2.txt \
    ../rk-ab-inputs/volume-n8-seed0.txt ../rk-ab-inputs/volume-n8-seed2.txt \
    ../rk-ab-inputs/volume-n16-seed0.txt ../rk-ab-inputs/volume-n16-seed2.txt \
    ../rk-ab-inputs/volume-n32-seed0.txt ../rk-ab-inputs/volume-n32-seed2.txt \
  --workers 1 8 --blocks 8 --repeats 2 --warmups 2 --diagnostics 3 \
  --output ../construction-order-main.json

LC_ALL=C python3 tools/benchmark_complex_construction.py \
  --baseline 4b9f149 --candidate ba57c79 \
  --inputs ../rk-ab-inputs/volume-n32-seed2.txt ../rk-ab-inputs/volume-n32-seed0.txt \
    ../rk-ab-inputs/volume-n16-seed2.txt ../rk-ab-inputs/volume-n16-seed0.txt \
    ../rk-ab-inputs/volume-n8-seed2.txt ../rk-ab-inputs/volume-n8-seed0.txt \
    ../rk-ab-inputs/terrain-n64-seed2.txt ../rk-ab-inputs/terrain-n64-seed0.txt \
    ../rk-ab-inputs/terrain-n16-seed2.txt ../rk-ab-inputs/terrain-n16-seed0.txt \
  --workers 8 1 --blocks 8 --repeats 2 --warmups 2 --diagnostics 3 \
  --output ../construction-order-confirmation.json

OMP_WAIT_POLICY=PASSIVE LC_ALL=C python3 tools/benchmark_resident_gradients.py \
  --benchmark ../work/ttk-benchmark/build-f4ffd1a1049d0ccf6e8f3eb4f7c096a6cc251ba0/morseframes_resident_gradient_benchmark \
  --terrain-sizes 64 --volume-sizes 8 16 32 --seeds 0 2 --workers 1 8 \
  --repeats 12 --diagnostics 6 --warmups 2 --input-dir ../rk-ab-inputs \
  --output ../resident-construction-order-ba57c79.json

python3 tools/render_resident_gradients.py \
  --input ../resident-construction-order-ba57c79.json \
  --table-output docs/resident_gradient_ba57c79_algorithm_table.tex \
  --construction-output docs/resident_gradient_ba57c79_construction_table.tex \
  --phases-output docs/resident_gradient_ba57c79_phases_table.tex
```

Choose unused output names for reruns; the runner refuses to overwrite evidence.

## Validation scope

Assessment: share with timing caveats. The benchmark sources, case coverage,
input/header hashes, exact-reference checks, timing partitions, diagnostic
counters and memory bounds are validated. Tables are checked against recomputed
summaries. Separate phase medians need not sum to total medians, and paired
ratios are not quotients of displayed time medians. Performance evidence here
is limited to the tested 2D/3D grids; correctness coverage is broader.
Local audit scripts `../validate_construction_studies.py` and
`../validate_order_report.py` recompute summaries and check all five rendered
tables and the three exact TeX fragments. The public runners and renderer
retain their own schema, timing-partition and reference checks.
