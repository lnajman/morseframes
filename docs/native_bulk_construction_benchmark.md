# Bulk face construction

This is the `10f641f` snapshot. The subsequent
[compact finalization study](native_finalization_benchmark.md) records the next
shared-index optimization and its timing caveats.

Revision `10f641f` reduces total native construction time by roughly **15–25%**
on the tested grids, relative to the preceding optimized constructor. It adds
an opt-in shared C++ construction path, not a new gradient kernel. A small
sequential-terrain gradient regression remains and is reported below.

This follows the [first native-construction optimization](native_construction_benchmark.md).
The timing boundary is unchanged: common vertex/function/cell arrays are already
resident, native construction is reported separately, and gradient comparison
includes fresh builders (and TTK ordering/lower-star processing).

## Implementation and scope

`add_lower_star_cells` in `include/morseframes/lower_star_complex.hpp` constructs
the closure of supplied cells with max-vertex filtration. It canonicalizes the
input cells once, enumerates one dimension at a time into contiguous temporary
records, sorts and deduplicates those records, and merges unique faces into the
ordered pending map using insertion hints. On volume n=32, there are still
2,681,190 generated faces, but only 792,051 unique faces are submitted to the map.

The finalized representation, lexicographic IDs, incidence order and gradient
algorithms are unchanged. Signed-zero ties preserve the legacy input-order
convention. Mixed-dimensional, duplicated and unordered cells are supported;
isolated vertices use singleton cells. Compact records cover dimensions 0–3,
with variable-length records above dimension three. Compatibility tests cover
dimensions through seven, but higher-dimensional performance was not measured.

Both current resident benchmark drivers adopt the helper. Other C++ callers can
opt in; existing `add_simplex` callers, Python construction and direct external
complex views do not automatically use it. RK itself remains applicable beyond
max-vertex filtrations. Construction remains sequential. Full-persistence
performance was not measured. See the [API contract](cpp_complex_view_api.md)
for validation, duplicate tolerance, diagnostics and partial-insertion behavior.

## Controlled comparison

Baseline: `9320f003b1cb48bbf7baebfa7c5f6371acaf054f` (the preceding constructor).
Candidate: `10f641f931d7fbaa15f85022f1bbf7cdb7ab7ab7` (bulk construction).
The same current driver is compiled against each header snapshot; it selects
the legacy per-face adapter when the new header is absent.

The main study covers terrain n=16/64 and volume n=8/16/32, seeds 0/2, with
RK workers 1/8: ten inputs and 20 configurations. Each configuration has eight
alternating version-order blocks, two performance repetitions per version/block,
and two warmups. F-Max/RK order also alternates. Three separate diagnostics and
three fresh-process memory samples are collected per version/input.

The repeat covers terrain n=64 and volume n=8/32, both seeds and both worker
counts, with the same settings and reversed input/worker order: 12 further
configurations. All **32 construction paired-block bootstrap intervals lie
below one**. Largest-volume repeat ratios range from 0.725 to 0.773.

Both studies started from clean candidate checkouts on the Apple M1 Max, with
native ARM Apple Clang 15 and `-std=c++17 -O3 -DNDEBUG -pthread`. No builds or
tests overlapped measurements. Other desktop applications remained active;
all outliers were retained. Intervals describe these same-machine sessions,
not independent-machine uncertainty or universal speedups. Worker counts apply
to RK after construction; they do not parallelize the constructor.

Main construction times (ms), using the eight-worker configurations and taking
medians across the two seed-specific medians:

| Input | Previous | Bulk | Paired bulk/previous |
| --- | ---: | ---: | ---: |
| 2D terrain, n=16 | 0.896 | 0.756 | 0.830 |
| 2D terrain, n=64 | 17.336 | 14.918 | 0.853 |
| 3D volume, n=8 | 9.690 | 7.163 | 0.751 |
| 3D volume, n=16 | 95.361 | 71.302 | 0.743 |
| 3D volume, n=32 | 968.074 | 733.088 | 0.763 |

Ratios aggregate paired block ratios, not quotients of the displayed time
medians. Historical absolute times should not be subtracted from these new
observations to estimate an additional speedup.

## Phases and memory

On n=32 in the main eight-worker comparison, uninstrumented enumeration and
insertion fall from **415.793 to 188.340 ms**. Finalization is approximately
unchanged: **546.379 to 548.742 ms**. It is now the dominant construction cost.
Separately aggregated phase medians need not sum to the total median.

Separate bulk diagnostics on n=32 (ms; medians within seed, then across seeds):

| Bulk phase | Time |
| --- | ---: |
| Validation and canonical cell copying | 6.380 |
| Enumeration | 20.363 |
| Sorting and deduplication | 109.756 |
| Unique insertion and temporary cleanup | 64.242 |

These are attribution diagnostics, not performance samples. The old adapter's
per-insertion clocks perturb its runtime substantially, so its instrumented
adapter total is not directly comparable to the lightly instrumented bulk
phases. All finalization subphases remain in the raw evidence.

Main-study fresh-process peak RSS (MiB; medians over the seed medians):

| Input | Previous | Bulk |
| --- | ---: | ---: |
| 2D terrain, n=16 | 1.617 | 1.688 |
| 2D terrain, n=64 | 8.953 | 9.250 |
| 3D volume, n=8 | 4.297 | 4.656 |
| 3D volume, n=16 | 30.539 | 32.531 |
| 3D volume, n=32 | 252.953 | 253.297 |

The largest per-input median increase is 2 MiB. RSS is sampled before gradients
and validation, and includes the executable/runtime, resident input, allocator
and construction temporaries; it is not live heap size. The small measured
peak increase does not imply that batching needs no extra temporary memory.
Its face buffer scales with generated faces in the largest dimensional batch.

## Correctness and gradient-timing caveats

Cross-version reference dumps compare the full complex and sequential F-Max/RK
sequences byte for byte. Each measured run checks its complex and gradient
fingerprints against those references, including parallel RK. Critical-simplex
counts are unchanged on every input. Unit tests additionally exercise mixed
dimensions, shuffled/repeated cells, existing entries and duplicate tolerance,
signed zeros, infinities, rejected inputs, copies/moves, cache invalidation,
repeated finalization and combinatorial count overflow.

Construction speedups do **not** establish unchanged gradient timing:

- Terrain n=64, seed 2, one-worker RK has an above-one interval in both studies.
  Its paired ratio is 1.122 in the main study and 1.052 in the repeat. The
  repeat's medians are 0.946 → 0.994 ms, while construction falls from
  16.903 → 13.739 ms. In the main study the unpaired medians are nearly equal
  (0.976 → 0.977 ms); the paired-ratio statistic answers a different question.
- The repeat also gives above-one F-Max intervals on one-worker terrain n=64:
  seed 0 has paired ratio 1.036 (1.313 → 1.347 ms), and seed 2 has ratio 1.066
  (1.239 → 1.290 ms). These intervals did not exclude one in the main study.
- A short exploratory pilot showed a much larger terrain F-Max slowdown, which
  was not reproduced at that magnitude in the controlled studies. Its raw file
  is retained but is not pooled with them.

The gradient code and logical representation are unchanged; allocation/cache
state or compiled-code layout could affect timings, but their contribution has
not been isolated. Retain the path as a construction optimization with these
caveats, not as a blanket improvement to gradient execution.

Local validation passes: ordinary C++ tests, ASan/UBSan, 151 native Python tests
(including both native benchmark adapters), fallback tests (151 run, 10 skipped),
and the nine-case RK corpus with native and fallback backends at workers 1/2/4.

## Refreshed RK / F-Max / TTK comparison

A separate clean-start study at `10f641f` covers terrain n=64 and volume
n=8/16/32, seeds 0/2, workers 1/8: eight inputs and 16 configurations. It uses
12 performance repetitions, six separate diagnostics and two warmups per mode,
with balanced six-way algorithm order. Both MorseFrames methods use the same
new construction helper. TTK uses the classic backend at clean pinned revision
`f4ffd1a1049d0ccf6e8f3eb4f7c096a6cc251ba0`, native ARM release, OpenMP and
`OMP_WAIT_POLICY=PASSIVE`.

Eight-worker algorithm-only times (ms; medians over seed-specific medians):

| Input | Sequential F-Max | RK | TTK |
| --- | ---: | ---: | ---: |
| 2D terrain, n=64 | 1.334 | 0.676 | 0.650 |
| 3D volume, n=8 | 0.762 | 0.434 | 0.578 |
| 3D volume, n=16 | 9.375 | 1.679 | 2.869 |
| 3D volume, n=32 | 123.329 | 13.793 | 20.051 |

RK leads both competitors by paired-median algorithm ratios in all 12 volume
configurations, including the one-worker cases. TTK retains a small lead on
the eight-worker terrain aggregate. Critical counts agree among all three
methods in all 16 configurations, with exact within-method reference checks.
These are current same-machine observations, not a controlled gradient speedup
estimated by subtracting a previous study.

The corresponding native construction costs remain much larger in MorseFrames:

| Input | F-Max construction (ms) | RK construction (ms) | TTK construction (ms) |
| --- | ---: | ---: | ---: |
| 2D terrain, n=64 | 14.249 | 14.128 | 0.531 |
| 3D volume, n=8 | 7.177 | 7.277 | 0.679 |
| 3D volume, n=16 | 73.493 | 72.737 | 4.350 |
| 3D volume, n=32 | 755.084 | 749.837 | 36.398 |

F-Max/RK construction differences are variation between measurements of the
same implementation. TTK wins full resident-to-gradient time against both
MorseFrames methods in all 16 configurations. At n=32 and eight workers, full
RK/TTK medians are **765.554 / 55.996 ms**. Separately aggregated phase medians
need not sum to those totals.

Loading remains a separate shared cost, not one repeated for each algorithm.
It is file read/parse/validation under the current cache state, not a cold-cache
disk benchmark. Construction excludes TTK vertex ordering and lower-star
processing: those remain in algorithm time. The complete one/eight-worker
algorithm, construction/loading and outer-phase tables are
`resident_gradient_10f641f_algorithm_table.tex`,
`resident_gradient_10f641f_construction_table.tex` and
`resident_gradient_10f641f_phases_table.tex`. Historical fragments are retained.

## Reproduction and provenance

Raw evidence is kept outside the repository, preserving earlier studies:

- `../construction-bulk-main.json`: UTC 2026-09-10 11:10:39 to 11:15:32;
  SHA-256 `2b63fba99c5dd43de982c49a9cce8a0fd3854a684dda2eeba29feb9ea8522faa`.
- `../construction-bulk-confirmation.json`: UTC 2026-09-10 11:15:59 to 11:20:19;
  SHA-256 `8f782914526105b10b38ea5984dbb1d609e57091dbe9646650203615facf2be0`.
- `../construction-bulk-pilot.json`: exploratory, dirty-checkout pilot only.
- `../resident-construction-bulk-10f641f.json`: UTC 2026-09-10 11:20:52 to 11:24:48;
  SHA-256 `72749e6c613b629c1b5c775440a97b349bdfe517dcd2e9e23a986a2bbf2a89a4`.

Both controlled files record input and binary hashes, complete samples, phase
diagnostics, memory samples and critical counts. Their summaries and timing
partitions were recomputed; input hashes, complete case/worker coverage,
diagnostic counts, reference agreement and memory bounds were checked.

- Previous headers SHA-256: `5ee7b89e47a48fc76370743819b4907e44557d00ea8798fd63dc75038862df4d`.
- Bulk headers SHA-256: `355c07e166b57f3f32d2bce3fdef14f45359f6015fcd589d3c392d77d75ca5f7`.
- Construction driver SHA-256: `a5359ed40f71fa68e7d7cc76eb3cdfc594437d476418e1ac40ad02623cf1eebb`.
- Construction runner SHA-256: `4e0887fd2816a46249a2dc2eaef8d0b383a1e8ce621ac8e1727efacbb28a71e1`.
- Resident driver SHA-256: `393c3af7f2df2d973031cc38ef7c5c44dcffecd11faa2434f8e03b97d9ab0bc8`.
- Resident executable SHA-256: `f3e2277d93309d395d1063b67d95541c7699a03834bef5bb7f33a5453afd1c02`.

```sh
LC_ALL=C python3 tools/benchmark_complex_construction.py \
  --baseline 9320f00 --candidate 10f641f \
  --inputs ../rk-ab-inputs/terrain-n16-seed0.txt ../rk-ab-inputs/terrain-n16-seed2.txt \
    ../rk-ab-inputs/terrain-n64-seed0.txt ../rk-ab-inputs/terrain-n64-seed2.txt \
    ../rk-ab-inputs/volume-n8-seed0.txt ../rk-ab-inputs/volume-n8-seed2.txt \
    ../rk-ab-inputs/volume-n16-seed0.txt ../rk-ab-inputs/volume-n16-seed2.txt \
    ../rk-ab-inputs/volume-n32-seed0.txt ../rk-ab-inputs/volume-n32-seed2.txt \
  --workers 1 8 --blocks 8 --repeats 2 --warmups 2 --diagnostics 3 \
  --output ../construction-bulk-main.json

LC_ALL=C python3 tools/benchmark_complex_construction.py \
  --baseline 9320f00 --candidate 10f641f \
  --inputs ../rk-ab-inputs/volume-n32-seed2.txt ../rk-ab-inputs/volume-n32-seed0.txt \
    ../rk-ab-inputs/terrain-n64-seed2.txt ../rk-ab-inputs/terrain-n64-seed0.txt \
    ../rk-ab-inputs/volume-n8-seed2.txt ../rk-ab-inputs/volume-n8-seed0.txt \
  --workers 8 1 --blocks 8 --repeats 2 --warmups 2 --diagnostics 3 \
  --output ../construction-bulk-confirmation.json

OMP_WAIT_POLICY=PASSIVE LC_ALL=C python3 tools/benchmark_resident_gradients.py \
  --benchmark ../work/ttk-benchmark/build-f4ffd1a1049d0ccf6e8f3eb4f7c096a6cc251ba0/morseframes_resident_gradient_benchmark \
  --terrain-sizes 64 --volume-sizes 8 16 32 --seeds 0 2 --workers 1 8 \
  --repeats 12 --diagnostics 6 --warmups 2 --input-dir ../rk-ab-inputs \
  --output ../resident-construction-bulk-10f641f.json

python3 tools/render_resident_gradients.py \
  --input ../resident-construction-bulk-10f641f.json \
  --table-output docs/resident_gradient_10f641f_algorithm_table.tex \
  --construction-output docs/resident_gradient_10f641f_construction_table.tex \
  --phases-output docs/resident_gradient_10f641f_phases_table.tex
```

Use unused output names for reruns; the runner rejects overwriting evidence.
