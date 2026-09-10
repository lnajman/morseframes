# Native construction: profiling and first optimization

This is the `d7cddc7` snapshot. The subsequent
[bulk face-construction study](native_bulk_construction_benchmark.md) records
the next optimization and its gradient-timing caveats.

The shared owning-complex constructor takes about **20% less time on the tested grids**
after revision `d7cddc7`. This is a construction improvement, not a new RK
algorithm. Gradient outputs are preserved, but gradient timings are not uniformly
unchanged: a small sequential-RK regression is reported below.

## What changed

The initial profile identified substantial costs in creating the finalized
simplex index and looking up boundary faces. An exploratory hash-table version
helped smaller cases but regressed on the largest volume; it was discarded.
The retained implementation keeps both ordered maps and all ordering conventions:

- Insert into the finalized index with an end hint, since the pending simplices
  are already lexicographically sorted.
- Read each filtration value directly from the pending entry instead of looking
  it up again.
- Reserve the final simplex array and each boundary vector.
- Reuse one temporary boundary-face vertex buffer across the traversal.

The changes are in `include/morseframes/filtered_complex.hpp`, shared by F-Max,
RK and other users of `FilteredSimplicialComplex`. They do not change the
benchmark's face-enumeration adapter, the gradient kernels, their spatial/level
scheduling, or direct external complex views. Construction remains sequential.
Full-persistence performance was not measured.

## Controlled construction comparison

Baseline `4eeb905597f7b1d49f5aff864b3474a97f2264c8` introduces only construction
diagnostics and the comparison harness, before the optimization. Candidate
`d7cddc7ad16ab0e41051ca5c4b991f45bc5a365e` is the optimized constructor.
Both studies started from a clean candidate checkout on the Apple M1 Max, using
native ARM Apple Clang 15, `-std=c++17 -O3 -DNDEBUG -pthread`.

The main study covers terrain sizes 16/64 and volume sizes 8/16/32, seeds 0/2,
and RK worker counts 1/8: ten inputs and 20 configurations. Each configuration
has eight alternating version-order blocks, two performance repetitions per
version/block and two warmups. F-Max/RK execution order also alternates.
There are three separate diagnostic runs and three fresh-process memory samples
per version/input. Worker counts apply to RK execution after construction, not
to the sequential constructor.

The repeat covers volume sizes 8/32, both seeds and both worker counts, with
the same settings: eight further configurations. All **28 construction paired-block
bootstrap intervals lie below one**. These intervals describe the sampled
sessions; they do not establish independent-machine uncertainty or universal
speedups. No builds or tests overlapped performance measurements; other desktop
applications remained active and all outliers were retained.

Main-study construction times (milliseconds), selecting the eight-worker
configurations and taking medians over the two seed-specific medians:

| Input | Before | After | Paired after/before |
| --- | ---: | ---: | ---: |
| 2D terrain, n=16 | 1.132 | 0.900 | 0.799 |
| 2D terrain, n=64 | 21.167 | 17.087 | 0.804 |
| 3D volume, n=8 | 11.545 | 9.058 | 0.789 |
| 3D volume, n=16 | 121.516 | 95.755 | 0.799 |
| 3D volume, n=32 | 1121.474 | 926.838 | 0.818 |

Ratios use paired block medians before aggregation; they need not equal quotients
of the displayed time medians. The largest-volume repeat's ratios range from
0.781 to 0.824 across its four seed/worker configurations.

## Refreshed RK / F-Max / TTK comparison

A separate clean-start study at `d7cddc7` measures terrain n=64 and volume
n=8/16/32, seeds 0/2, workers 1/8: eight inputs and 16 configurations. It uses
12 performance repetitions, six separate diagnostics and two warmups per mode,
with balanced six-way algorithm order. Loading and native construction remain
separate; TTK vertex ordering and lower-star processing remain in algorithm time.
F-Max is always sequential. These are new same-machine observations, not a
controlled estimate obtained by subtracting older tables.

Eight-worker algorithm times, excluding construction (milliseconds; medians of
the two seed-specific medians):

| Input | Sequential F-Max | RK | TTK |
| --- | ---: | ---: | ---: |
| 2D terrain, n=64 | 1.305 | 0.620 | 0.597 |
| 3D volume, n=8 | 0.736 | 0.392 | 0.507 |
| 3D volume, n=16 | 8.497 | 1.434 | 2.553 |
| 3D volume, n=32 | 114.164 | 12.768 | 19.090 |

RK has lower paired-median algorithm time than both competitors in all 12
volume configurations, including the one-worker cases. The eight-worker terrain
aggregate still gives TTK a small lead. Critical counts agree among all three
methods in all 16 configurations, and match the earlier input-specific results.
Each method also passes its own exact reference checks.

The corresponding eight-worker native construction costs remain substantial:

| Input | F-Max construction (ms) | RK construction (ms) | TTK construction (ms) |
| --- | ---: | ---: | ---: |
| 2D terrain, n=64 | 17.058 | 17.050 | 0.495 |
| 3D volume, n=8 | 9.255 | 9.298 | 0.646 |
| 3D volume, n=16 | 97.478 | 95.717 | 4.188 |
| 3D volume, n=32 | 914.886 | 914.737 | 34.622 |

F-Max and RK use the same construction implementation; their column differences
are measurement variation. TTK has lower full resident-to-gradient time than
both MorseFrames methods in all 16 configurations. For n=32 at eight workers,
the full RK/TTK medians are 929.507 / 53.486 ms. Separately aggregated phase
medians need not add to these total medians.

The complete one/eight-worker algorithm, construction/loading and outer-phase
tables are `resident_gradient_d7cddc7_algorithm_table.tex`,
`resident_gradient_d7cddc7_construction_table.tex` and
`resident_gradient_d7cddc7_phases_table.tex`. Raw evidence is
`../resident-construction-d7cddc7.json` (UTC 2026-09-10 10:47:12 to 10:51:46).
TTK uses the clean pinned revision
`f4ffd1a1049d0ccf6e8f3eb4f7c096a6cc251ba0`, classic backend, OpenMP and
`OMP_WAIT_POLICY=PASSIVE`; the native ARM release build is unchanged apart from
the MorseFrames header update. Its resident-driver SHA-256 is
`a1ca9b1f76203e376a5b8b8d612ddd9cba0d7aee007d0f12e65aa76ecff0086b`, executable
SHA-256 is `ac1acbd95ce5b09a094d060b7cdabc0112ef79f16a4ebbb0caea6631cee0ddf1`,
and the headers hash is the candidate hash recorded below.

```sh
OMP_WAIT_POLICY=PASSIVE LC_ALL=C python3 tools/benchmark_resident_gradients.py \
  --benchmark ../work/ttk-benchmark/build-f4ffd1a1049d0ccf6e8f3eb4f7c096a6cc251ba0/morseframes_resident_gradient_benchmark \
  --terrain-sizes 64 --volume-sizes 8 16 32 --seeds 0 2 --workers 1 8 \
  --repeats 12 --diagnostics 6 --warmups 2 --input-dir ../rk-ab-inputs \
  --output ../resident-construction-d7cddc7.json

python3 tools/render_resident_gradients.py \
  --input ../resident-construction-d7cddc7.json \
  --table-output docs/resident_gradient_d7cddc7_algorithm_table.tex \
  --construction-output docs/resident_gradient_d7cddc7_construction_table.tex \
  --phases-output docs/resident_gradient_d7cddc7_phases_table.tex
```

## Phase attribution and memory

Largest-volume diagnostic phase medians (milliseconds, first within each seed,
then across seeds):

| Phase | Before | After |
| --- | ---: | ---: |
| Adapter: face enumeration and insertion, inclusive | 456.998 | 454.551 |
| Reset | <0.001 | <0.001 |
| Finalized index and simplex records | 197.614 | 107.276 |
| Filtration levels | 55.187 | 53.651 |
| Boundaries and filtration validation | 287.741 | 168.509 |
| Coboundaries | 59.048 | 56.980 |
| Filtration ordering and level buckets | 122.511 | 120.477 |

The adapter diagnostics additionally split canonicalization/deduplication
(295.491 / 293.378 ms before/after) from an enumeration residual
(161.611 / 161.188 ms). That residual includes per-insertion clock overhead;
it is **not an uninstrumented estimate of pure enumeration time**. These two
subphases are nested inside the adapter total, not additional work. Separately
aggregated medians need not sum. Performance claims use only the uninstrumented
paired runs, not this phase table.

Fresh-process peak RSS after construction, before any gradient/reference or
validation allocation (MiB; medians across three samples per seed, then seeds):

| Input | Before | After |
| --- | ---: | ---: |
| 2D terrain, n=16 | 1.703 | 1.602 |
| 2D terrain, n=64 | 11.234 | 9.000 |
| 3D volume, n=8 | 5.594 | 4.320 |
| 3D volume, n=16 | 39.844 | 30.516 |
| 3D volume, n=32 | 256.109 | 252.961 |

RSS includes the native runtime, resident input, allocator state and constructor
temporaries. It is not live heap size or an isolated allocation count. The
largest-case saving is only about 3 MiB (1.2%); this is primarily a time
optimization, not a solution to large-input memory growth.

## Gradient checks and timing caveats

For each version/input, the harness writes and compares complete reference
complexes (vertices, filtration, levels, IDs, boundaries, coboundaries and
orders), plus complete F-Max and RK sequences. Every subsequent performance
sample checks fingerprints of the complex and each algorithm's sequence
against its reference, including RK at eight workers. Diagnostic fingerprints
and structural counters also agree. Each method is checked against its own
reference; equal critical counts do not imply identical cross-method gradients.

Critical counts agree between versions and between F-Max and RK on these inputs.
However, the constructor changes allocation behavior, so unchanged combinatorial
outputs do not imply unchanged execution times:

- Sequential RK on volume n=8, seed 2 takes 0.371 -> 0.416 ms in the main study
  and 0.394 -> 0.430 ms in the repeat. Paired ratios are 1.147 and 1.080, with
  above-one intervals in both studies. The roughly 0.04 ms regression is much
  smaller than the approximately 2.4 ms construction saving for this input.
- The main study also flags sequential RK on n=8, seed 0; that interval does not
  remain above one in the repeat.
- The repeat flags F-Max on n=32, seed 2 in the one-worker configuration
  (106.275 -> 109.355 ms; paired ratio 1.036). Its main-study interval does not
  exclude one. No general claim of gradient-time equivalence is made.

Allocation/cache effects are a possible explanation for these timing changes,
not an established cause. The next construction target is the unchanged
enumeration/insertion stage, especially repeated shared-face submissions.

## Reproduction and provenance

Run from the repository root. Inputs are the same resident-array files described
in [benchmark reproduction](benchmark_reproduction.md); input generation and
loading are outside these construction timers. Both header snapshots must
support the construction-metrics API introduced in `4eeb905`.

```sh
LC_ALL=C python3 tools/benchmark_complex_construction.py \
  --baseline 4eeb905 --candidate d7cddc7 \
  --inputs ../rk-ab-inputs/terrain-n16-seed0.txt ../rk-ab-inputs/terrain-n16-seed2.txt \
    ../rk-ab-inputs/terrain-n64-seed0.txt ../rk-ab-inputs/terrain-n64-seed2.txt \
    ../rk-ab-inputs/volume-n8-seed0.txt ../rk-ab-inputs/volume-n8-seed2.txt \
    ../rk-ab-inputs/volume-n16-seed0.txt ../rk-ab-inputs/volume-n16-seed2.txt \
    ../rk-ab-inputs/volume-n32-seed0.txt ../rk-ab-inputs/volume-n32-seed2.txt \
  --workers 1 8 --blocks 8 --repeats 2 --warmups 2 --diagnostics 3 \
  --output ../construction-ordered-main.json
```

Repeat with only volume n=8/32, both seeds, and a new output path
`../construction-ordered-confirmation.json`. The runner refuses to overwrite
previous evidence. The exploratory raw files are `../construction-hash-pilot.json`
and `../construction-ordered-pilot.json`; neither supplies the headline results.

- Main UTC window: 2026-09-10 10:36:37 to 10:42:01.
- Repeat UTC window: 2026-09-10 10:42:07 to 10:46:59.
- Baseline headers SHA-256: `8908e3cd77d2778bfc344e2bb256f7c5ee18531e4f16fb89a307f57a4f1f7c18`.
- Candidate headers SHA-256: `5ee7b89e47a48fc76370743819b4907e44557d00ea8798fd63dc75038862df4d`.
- Driver SHA-256: `d05a68ad1ec947ed14a5667be16a47d24fea6810d1e0ba6f4ae2ef3703a2a94d`.
- Runner SHA-256: `df14899d94f5bc396a808fa4f400c5a6842d705a5e1bc163df7d58050a2cecc9`.

Raw JSON retains compiler flags, executable/input hashes, sample order,
individual measurements, counters, critical counts, memory samples and
diagnostics. The local read-only audit companion is
`../validate_construction_studies.py`; it recomputes summaries and checks input
hashes, matrix coverage, phase partitions and reference metadata. Raw files
and this local companion are outside the public repository; the reusable
benchmark runner, tests and selected tables are tracked.

Validation: normal C++ and ASan/UBSan suites passed; 151 native Python tests
passed, the fallback suite ran 151 tests with 10 skips, and both nine-case
native/fallback corpora passed. New C++ tests cover shuffled insertion, large
vertex IDs, dimensions 1–7, duplicate tolerance, copies/moves, repeated
finalization, cache invalidation and concurrent const lookups. The standalone
driver also passed its native integration test. Strict Sphinx and benchmark
renderer checks passed. The three generated table fragments were checked
byte-for-byte against fresh rendering, and the rendered report's five tables
were checked for expected row/column counts and current values.

Assessment: **share with caveats**—two seeds, regular synthetic grids, one
desktop machine, limited memory instrumentation and the gradient-timing
trade-offs above. Persistence was checked for correctness, not performance.
