# Compact lookup during finalization

Revision `5812849` removes the duplicate finalized simplex tree index. In the
main study, per-size construction summaries take about **10–22% less time**;
largest-volume peak RSS falls from **253.195 to 195.891 MiB**. The memory
reduction is more consistent than the largest-volume timing gain. Several
gradient-timing warnings remain, so this is not a claim of uniformly faster
gradient execution.

This follows the [bulk construction study](native_bulk_construction_benchmark.md).
Loading and native construction remain separate from gradient comparison.
The primary deliverable is a shared-library implementation with reproducible
benchmark evidence; no full-persistence speedup is claimed.

## What changed

Finalized simplex records already appear in lexicographic vertex order, and
their array positions are their IDs. The new lookup reuses that sorted array
instead of building another tree containing copied vertex keys and IDs. A small
sorted array records the start of each first-vertex range. Lookup first finds
the range, then searches its simplex records; the equal first vertex can be
skipped during comparisons inside that range.

The range keys are sparse: a large vertex ID does not require a large dense
array. Offsets, rather than pointers, preserve copy/move behavior. Const lookup
does not modify a cache. Missing and empty queries, duplicate-vertex rejection,
face-closure checks and filtration-monotonicity checks retain their behavior.
Insertion still uses the pending ordered map; new pending simplices become
findable only after finalization, as before.

The change is in `include/morseframes/filtered_complex.hpp`. It applies to
every caller using this owning complex, not just the bulk helper or RK. It
supports arbitrary monotone filtrations and is not restricted to dimensions
two/three. Direct external complex views and their indexes are unchanged.
Construction remains sequential. Gradient code, simplex IDs, vertex values,
boundary/coboundary order, filtration order and level buckets are unchanged.

## Controlled construction comparison

Baseline: `8b88ffad974d2cf7ec878d03e729c39fcdb0abb7` (bulk construction with the
duplicate finalized tree). Candidate: `5812849c0e08edb7e35cb7dfa707746df11aa96a`
(compact lookup). Both versions use the same bulk insertion helper, driver,
compiler flags and existing resident input arrays.

The main study covers terrain n=16/64 and volume n=8/16/32, seeds 0/2,
workers 1/8: ten inputs and 20 configurations. The repeat covers terrain n=64
and volume n=8/16/32, both seeds and worker counts, in reversed input/worker
order: eight inputs and 16 configurations. Each configuration has eight
alternating version-order blocks, two repetitions per version/block and two
warmups. F-Max/RK execution order also alternates. Each input/version has three
separate diagnostics and three fresh-process memory samples.

Both studies start clean on the native ARM Apple M1 Max, Apple Clang 15,
`-std=c++17 -O3 -DNDEBUG -pthread`. No builds or tests overlap measurements;
other desktop applications remain active and all outliers are retained. Worker
counts apply to RK after construction, not to the sequential constructor.

All **36 paired construction medians are below one**. The paired-block bootstrap
interval excludes one in 18/20 main configurations and 12/16 repeat configurations:
**30/36 overall**. These intervals describe the sampled sessions, not independent
machine uncertainty or universal speedups. The six inconclusive intervals are:

- Main: volume n=16/seed 0/eight workers and n=32/seed 2/eight workers.
- Repeat: volume n=32/seed 0 at both worker counts, terrain n=64/seed 0/eight
  workers, and volume n=8/seed 2/one worker.

Main construction times (ms), selecting the eight-worker setting and taking
medians over the two seed-specific medians:

| Input | Previous | Compact lookup | Paired new/previous |
| --- | ---: | ---: | ---: |
| 2D terrain, n=16 | 0.757 | 0.579 | 0.781 |
| 2D terrain, n=64 | 14.009 | 11.162 | 0.799 |
| 3D volume, n=8 | 6.887 | 5.534 | 0.806 |
| 3D volume, n=16 | 74.159 | 60.624 | 0.822 |
| 3D volume, n=32 | 716.155 | 660.705 | 0.894 |

Ratios aggregate paired block ratios, not quotients of the displayed time
medians. Largest-volume repeat ratios range from 0.831 to 0.967; two of its
four intervals include one. Do not estimate a combined improvement by
subtracting absolute times from earlier studies.

## Phase attribution and memory

On n=32 at the main eight-worker setting, uninstrumented finalization falls
from **528.743 to 464.590 ms**. Enumeration/insertion, whose code is unchanged,
is **187.281 / 197.219 ms**. Separately aggregated phase medians need not sum
to total medians.

Separate n=32 finalization diagnostics (ms; medians within seed, then across seeds):

| Phase | Previous | Compact lookup |
| --- | ---: | ---: |
| Index and simplex records | 100.703 | 34.256 |
| Filtration levels | 57.184 | 59.685 |
| Boundary lookup and checks | 173.458 | 167.019 |
| Coboundaries | 57.944 | 58.830 |
| Ordering and level buckets | 119.997 | 119.049 |

The clearest phase saving is removal of duplicate index construction; the
aggregate boundary improvement is smaller. These are separate diagnostic
samples, not performance measurements. Reset time and all bulk subphases remain
in the raw files. Boundary lookup and ordering still deserve attention, but
this experiment does not isolate a further optimization for either.

Main-study fresh-process peak RSS (MiB; medians over seed medians):

| Input | Previous | Compact lookup |
| --- | ---: | ---: |
| 2D terrain, n=16 | 1.734 | 1.625 |
| 2D terrain, n=64 | 9.305 | 7.414 |
| 3D volume, n=8 | 4.672 | 3.969 |
| 3D volume, n=16 | 32.578 | 25.500 |
| 3D volume, n=32 | 253.195 | 195.891 |

The largest-volume reduction is about 57.3 MiB (22.6%). Per-input median peak
RSS is lower for every input in both studies. These fresh-process measurements
occur before gradients and validation; they include the runtime, resident
input, allocator and construction temporaries. They are not live heap sizes.

## Correctness and timing caveats

Full cross-version reference dumps compare the complex and sequential F-Max/RK
sequences byte for byte. Every measured run checks its complex and gradient
fingerprints, including parallel RK. Critical counts are unchanged. Additional
C++ tests exercise sparse IDs through the maximum `VertexId`, present/missing
subsets, empty and duplicate queries, absent first-vertex ranges, missing faces,
insertion before/after finalization, copies/moves, repeated finalization and
concurrent const lookup. The existing mixed-dimensional tests extend through
dimension seven; performance evidence here remains limited to the tested grids.

Gradient timings are not uniformly improved:

- Main: volume n=32/seed 2/eight workers gives above-one F-Max and RK intervals.
  Paired ratios are 1.224 and 1.158; time medians are 118.884 → 137.157 ms and
  13.279 → 14.238 ms, respectively.
- Repeat of that configuration: F-Max's ratio is 1.079 with an interval including
  one; RK's ratio is 0.830 with its interval below one. The earlier RK slowdown
  is not reproduced. This does not establish F-Max timing equivalence.
- Other repeat warnings: one-worker F-Max on volume n=32/seed 0 has ratio 1.097
  (118.526 → 142.535 ms), and eight-worker F-Max on volume n=16/seed 2 has ratio
  1.085 (10.100 → 10.687 ms). Neither had an above-one interval in the main study.

No identical gradient configuration has an above-one interval in both studies,
but this is not proof of no regression. In particular, F-Max's large-volume
timings warrant caution. Allocation/cache state, scheduling and code layout are
possible influences, not isolated causes. Retain this as a measured construction
and memory improvement with timing caveats, not a blanket gradient speedup.

Local checks pass: ordinary C++ tests, ASan/UBSan, 151 native Python tests
(including both native benchmark adapters), fallback suite (151 run, 10 skipped),
and the nine-case RK corpus on both backends at workers 1/2/4.

## Refreshed gradient comparison

A separate clean-start resident-array study at `5812849` covers terrain n=64,
volume n=8/16/32, seeds 0/2 and workers 1/8: eight inputs and 16 configurations.
It uses 12 performance repetitions, six separate diagnostics and two warmups
per mode, with balanced six-way algorithm order. TTK remains at clean pinned
revision `f4ffd1a1049d0ccf6e8f3eb4f7c096a6cc251ba0`, classic backend, native
ARM release, OpenMP and `OMP_WAIT_POLICY=PASSIVE`.

Eight-worker algorithm-only times (ms; medians over seed-specific medians):

| Input | Sequential F-Max | RK | TTK |
| --- | ---: | ---: | ---: |
| 2D terrain, n=64 | 1.307 | 0.712 | 0.632 |
| 3D volume, n=8 | 0.740 | 0.445 | 0.660 |
| 3D volume, n=16 | 8.953 | 1.702 | 2.764 |
| 3D volume, n=32 | 119.798 | 12.338 | 19.228 |

RK leads both competitors by paired-median algorithm ratios in all 12 volume
configurations, including one-worker cases. TTK leads the eight-worker terrain
aggregate. Critical counts agree across all three methods in all 16 configurations,
and each method passes its exact reference checks. This new comparison does not
erase the controlled A/B warnings or establish a causal gradient improvement.

Corresponding native construction costs:

| Input | F-Max construction (ms) | RK construction (ms) | TTK construction (ms) |
| --- | ---: | ---: | ---: |
| 2D terrain, n=64 | 11.331 | 11.254 | 0.506 |
| 3D volume, n=8 | 5.696 | 5.790 | 0.740 |
| 3D volume, n=16 | 59.764 | 60.943 | 4.188 |
| 3D volume, n=32 | 639.118 | 636.944 | 35.755 |

F-Max/RK construction differences reflect variation between measurements of
the same implementation. TTK has lower full resident-to-gradient time in all
16 configurations; n=32/eight-worker full RK/TTK medians are **649.972 / 54.744 ms**.
Separately aggregated phase medians need not sum to those totals.

Loading is reported separately as one shared read/parse/validation per invocation,
not as a cost repeated for each method or a cold-cache disk benchmark. Native
construction is excluded from algorithm comparison, while fresh MorseFrames
builders and TTK vertex ordering/lower-star processing remain included. Complete
one/eight-worker algorithm, construction/loading and outer-phase tables are
`resident_gradient_5812849_algorithm_table.tex`,
`resident_gradient_5812849_construction_table.tex` and
`resident_gradient_5812849_phases_table.tex`. Earlier snapshots are preserved.

## Reproduction and provenance

Raw evidence remains outside the repository alongside the historical studies:

- `../construction-prefix-main.json`: UTC 2026-09-10 11:36:23 to 11:40:43;
  SHA-256 `e6bb6df504f81cef38f3bf47de0de9c40bce6e3d2947900381cb5bb4029e5d95`.
- `../construction-prefix-confirmation.json`: UTC 2026-09-10 11:41:54 to 11:46:24;
  SHA-256 `34611c5377827424082024201334756c9295f06e78a5bb0016b19e4c9a33ecc2`.
- `../construction-prefix-pilot.json`: exploratory dirty-checkout pilot,
  retained but not pooled with the controlled studies.
- `../resident-construction-prefix-5812849.json`: UTC 2026-09-10 11:47:05 to 11:50:33;
  SHA-256 `d1a77a0c5fd77ac25a5e298abe1e2ea26c4da0e449dca4dd0da7cf37eb344213`.

The controlled files record revisions, input/header/binary hashes, all samples,
diagnostics and memory measurements. Summaries and timing partitions were
recomputed; input hashes, complete case/worker coverage, reference agreement,
critical counts, diagnostic counters and memory bounds were checked.

- Baseline headers: `355c07e166b57f3f32d2bce3fdef14f45359f6015fcd589d3c392d77d75ca5f7`.
- Candidate headers: `7b1309f114613cdf1025c862c82871855520f9098cd3b62ddf4538aa39721c5e`.
- Construction driver: `a5359ed40f71fa68e7d7cc76eb3cdfc594437d476418e1ac40ad02623cf1eebb`.
- Construction runner: `4e0887fd2816a46249a2dc2eaef8d0b383a1e8ce621ac8e1727efacbb28a71e1`.
- Resident driver: `393c3af7f2df2d973031cc38ef7c5c44dcffecd11faa2434f8e03b97d9ab0bc8`.
- Resident executable: `902dfc1f17e9f004adc89a8f3ae62a5dd9d151b699c0ed05510acb58cd616b8b`.

```sh
LC_ALL=C python3 tools/benchmark_complex_construction.py \
  --baseline 8b88ffa --candidate 5812849 \
  --inputs ../rk-ab-inputs/terrain-n16-seed0.txt ../rk-ab-inputs/terrain-n16-seed2.txt \
    ../rk-ab-inputs/terrain-n64-seed0.txt ../rk-ab-inputs/terrain-n64-seed2.txt \
    ../rk-ab-inputs/volume-n8-seed0.txt ../rk-ab-inputs/volume-n8-seed2.txt \
    ../rk-ab-inputs/volume-n16-seed0.txt ../rk-ab-inputs/volume-n16-seed2.txt \
    ../rk-ab-inputs/volume-n32-seed0.txt ../rk-ab-inputs/volume-n32-seed2.txt \
  --workers 1 8 --blocks 8 --repeats 2 --warmups 2 --diagnostics 3 \
  --output ../construction-prefix-main.json

LC_ALL=C python3 tools/benchmark_complex_construction.py \
  --baseline 8b88ffa --candidate 5812849 \
  --inputs ../rk-ab-inputs/volume-n32-seed2.txt ../rk-ab-inputs/volume-n32-seed0.txt \
    ../rk-ab-inputs/volume-n16-seed2.txt ../rk-ab-inputs/volume-n16-seed0.txt \
    ../rk-ab-inputs/terrain-n64-seed2.txt ../rk-ab-inputs/terrain-n64-seed0.txt \
    ../rk-ab-inputs/volume-n8-seed2.txt ../rk-ab-inputs/volume-n8-seed0.txt \
  --workers 8 1 --blocks 8 --repeats 2 --warmups 2 --diagnostics 3 \
  --output ../construction-prefix-confirmation.json

OMP_WAIT_POLICY=PASSIVE LC_ALL=C python3 tools/benchmark_resident_gradients.py \
  --benchmark ../work/ttk-benchmark/build-f4ffd1a1049d0ccf6e8f3eb4f7c096a6cc251ba0/morseframes_resident_gradient_benchmark \
  --terrain-sizes 64 --volume-sizes 8 16 32 --seeds 0 2 --workers 1 8 \
  --repeats 12 --diagnostics 6 --warmups 2 --input-dir ../rk-ab-inputs \
  --output ../resident-construction-prefix-5812849.json

python3 tools/render_resident_gradients.py \
  --input ../resident-construction-prefix-5812849.json \
  --table-output docs/resident_gradient_5812849_algorithm_table.tex \
  --construction-output docs/resident_gradient_5812849_construction_table.tex \
  --phases-output docs/resident_gradient_5812849_phases_table.tex
```

Use unused output names for reruns; the runner rejects overwriting evidence.
