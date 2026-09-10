# RK level timelines: one large level limits scaling

## Result and recommendation

**The difficult cases are dominated by one long level, not principally by late
scheduling.** At eight workers, the longest level occupies about 91% of the
level-processing wall time in 7D seed 3 and 89% in 6D seed 4. The same levels
dominate every coarse trace in both measurement sessions.

The next algorithm experiment should therefore target **parallel work within
large levels**, especially closure preparation, while retaining level parallelism
for small levels and the existing bounded executor. Reordering levels alone is
unlikely to remove this bottleneck. This is a recommendation, not an implemented
optimization or a promised speedup.

There is also an unresolved regression-control warning: volume-32 has slower
unprofiled RK medians in both old/new sessions. The optional trace adds no vector
allocation or per-level timestamp calls to ordinary entry points, but these
measurements do **not** establish zero performance impact from the refactor.
Before adopting a new performance baseline, isolate this warning with an
ordinary-only driver/binary comparison and repeat the volume control.

**Follow-up:** the [ordinary-binary isolation study](rk_binary_isolation_benchmark.md)
has now completed. The larger slowdown is not consistently reproduced, and
same-executable process controls show comparable variability. A small sequential
cost remains possible; the original observations above are preserved, not
reclassified as proof of zero overhead.

## Implementation and timing scope

Candidate `8edb1ec7ee459c6036a4817f7e2e52907b20b5cc` adds a separate optional
native level-profile entry point, JSON commands, a frozen runner and audits.
The RK scheduler, kernels, level order and deterministic replay are unchanged.
The ordinary builder instantiates the non-tracing implementation; the traced
version records task IDs, start offsets, durations, level sizes, event counts,
and optionally existing detailed phases/counters. Task IDs identify persistent
tasks/scratch slots, not physical CPU cores or OS threads.

Start from the same finalized native complex in memory. Loading and construction
are separate. Fresh builder, pools, all algorithm-specific preparation, local
work, replay and in-method cleanup remain in gradient time. No persistence is
computed. TTK is not rerun; the previous
[direct TTK comparison](pls_key_arena_benchmark.md) is unchanged.

The [protocol](rk_level_profile_protocol.md) specifies four synthetic inputs:
7D/side-2 seeds 3 and 1, and 6D/side-2 seeds 4 and 1. Within each dimension the
geometry is identical and only the injective vertex ordering changes. These
single-hypercube triangulations are diagnostic cases, not a representative
higher-dimensional corpus. Each input runs at 1, 4 and 8 workers, with two
warmups per mode and ten balanced blocks of five profiling modes. Confirmation
reverses input, worker and block order. No local build/test work overlaps timing.

## Coarse traces identify the bottleneck

Times below are medians in milliseconds. Start is relative to the level-phase
origin, not the whole gradient call. Percentages are medians of within-trace
duration/wall fractions; do not recompute them from separate displayed medians.
W denotes worker count; share is the longest level's duration divided by level wall time.

| Case | Study | W | Start | Duration | Wall | Share |
| --- | --- | --- | --- | --- | --- | --- |
| 7D seed 3 | Main | 4 | 5.292 | 31.743 | 37.317 | 85.74% |
| 7D seed 3 | Confirmation | 4 | 5.099 | 33.514 | 39.120 | 86.39% |
| 7D seed 3 | Main | 8 | 3.138 | 32.298 | 35.528 | 91.19% |
| 7D seed 3 | Confirmation | 8 | 3.021 | 32.124 | 35.060 | 91.54% |
| 6D seed 4 | Main | 4 | 0.226 | 1.283 | 1.522 | 84.86% |
| 6D seed 4 | Confirmation | 4 | 0.230 | 1.254 | 1.503 | 84.43% |
| 6D seed 4 | Main | 8 | 0.151 | 1.341 | 1.507 | 88.80% |
| 6D seed 4 | Confirmation | 8 | 0.142 | 1.267 | 1.429 | 89.05% |

Level 116 contains 69,072 of 189,171 simplices in 7D seed 3; level 53 contains
5,274 of 18,731 in 6D seed 4. Both are owned by **vertex 0**, but neither is
filtration level 0. The runner explicitly maps vertex-indexed star sizes through
the input scalar ordering.

Holding the measured level durations fixed, even ideal ordering cannot reduce
level wall time below the longest duration. This leaves only about 8.5–8.8%
of the eight-worker level-phase time removable by scheduling in the 7D case,
and 10.9–11.2% in the 6D case. Whole-gradient gains would be smaller with other
phases fixed. These are diagnostic bounds under fixed durations, not bounds on
real scheduler experiments: contention, cache behavior and scratch reuse may
change the durations themselves.

The matched controls differ: the longest-level share at eight workers is
51.77%/54.04% for 7D seed 1 and 61.92%/59.57% for 6D seed 1 (main/confirmation).
The longest level can change between repetitions. Thus the finding supports
selective treatment of expensive levels, not parallelizing every small level.

## Closure preparation is the largest detailed phase

Detailed instrumentation perturbs execution and cannot replace unprofiled
benchmark times. It nevertheless identifies consistent phase structure in the
dominant levels. The entries below are median per-level phase fractions.

| Case | Study | Closure | Facet execution | Discovery | Incidence |
| --- | --- | --- | --- | --- | --- |
| 7D seed 3 | Main | 52.18% | 31.95% | 10.58% | 4.18% |
| 7D seed 3 | Confirmation | 52.49% | 32.26% | 9.78% | 4.10% |
| 6D seed 4 | Main | 45.79% | 35.73% | 10.13% | 6.04% |
| 6D seed 4 | Confirmation | 45.80% | 35.27% | 10.20% | 6.19% |

In 7D seed 3, the detailed dominant level takes 39.931/40.218 ms. Its closure
phase takes 20.777/21.056 ms, including traversal 8.084/7.834 ms,
materialization 4.684/4.573 ms, boundary indexing 3.684/3.960 ms and sorting
2.487/2.458 ms. These are separate medians, not an additive exact partition.
The 6D dominant level takes 1.735/1.875 ms, with closure 0.779/0.829 ms.

Core and local reduction are children of facet execution, not additional
phases. In the 7D dominant level their median times are 5.545/5.621 ms and
4.884/4.871 ms. Parallelizing only the final local reductions would therefore
leave most of the measured level cost serial. A useful next design must examine
closure preparation as well as independent facet kernels, with disjoint task
storage and deterministic merging. No such change is made in this pass.

## Instrumentation cost is reported separately

Across all twelve input/worker configurations, coarse trace/plain paired-median
ratios range from 0.964–1.026 in main and 0.877–1.045 in confirmation. Ratios
below one are not evidence that timestamps make RK intrinsically faster.
Detailed trace/plain ratios are 1.139–1.313 and 1.159–1.441 respectively.

For the difficult eight-worker 7D case, coarse trace/plain is
0.995 [0.962, 1.048] / 0.983 [0.948, 1.056], whereas detailed trace/plain is
1.241 [1.196, 1.280] / 1.223 [1.172, 1.315]. Brackets are 95% paired-block
bootstrap intervals. The incremental coarse trace/existing-coarse ratio is
1.019 [1.005, 1.061] in main, but 0.968 [0.793, 1.051] in confirmation.
For 6D seed 1, incremental detailed trace/existing-detailed is
1.156 [1.005, 1.179] in confirmation. Do not claim negligible or zero overhead
universally. All six mode comparisons at all worker counts are in the
[numerical appendix](rk_level_profile_tables.md).

## Unprofiled regression controls remain qualified

Baseline `876354e777fc491c63d4274e78e6af62045c91cd` and the candidate were tested
with the same current driver on seven inputs at 1/8 workers, six balanced blocks
of two repetitions per algorithm, two warmups, and a reversed confirmation.
Complete ordered old/new reference dumps match, as do the audited work counters.

The 7D seed-3 eight-worker RK ratio is 1.049 [1.026, 1.210] in main but
0.972 [0.941, 0.996] in confirmation: that slowdown does not repeat.
Volume-32 seed 0 remains a warning:

| Workers | Main RK candidate/baseline | Confirmation RK candidate/baseline |
| --- | --- | --- |
| 1 | 1.085 [0.973, 1.110] | 1.044 [1.012, 1.114] |
| 8 | 1.022 [0.986, 1.038] | 1.075 [0.9996, 1.1325] |

No RK configuration has a slowdown interval wholly above one in both sessions,
but that is not an equivalence test. The unchanged PLS control on volume-32 at
one worker is also slower in both sessions: 1.094 [1.036, 1.130] and
1.064 [1.023, 1.103]. This cautions against attributing all differences to RK
itself, without establishing their cause. Binary layout, allocation history and
run conditions have not been isolated. Preserve the preceding frozen benchmark
as the headline comparison until this warning is resolved. The appendix retains
all algorithms and all controls, including negative and inconclusive results.

## Correctness, reproduction and provenance

All three algorithms have identical critical-count vectors on the four focused
inputs; they are unchanged across worker counts and sessions:

| Case | Critical counts by dimension, starting at dimension 0 |
| --- | --- |
| 7D seed 3 | `[5, 17, 23, 10, 0, 0, 0, 0]` |
| 7D seed 1 | `[4, 4, 1, 0, 0, 0, 0, 0]` |
| 6D seed 4 | `[1, 3, 4, 1, 0, 0, 0]` |
| 6D seed 1 | `[4, 4, 1, 0, 0, 0, 0]` |

The detailed dominant levels contain 34,536 reductions, zero perforations,
8 rounds and 29,631 facet kernels (7D), and 2,637 reductions, zero perforations,
7 rounds and 1,916 facet kernels (6D). All per-level counters agree across
repetitions and worker counts. Traces cover every level/simplex/event, task
intervals do not overlap, and phase nesting/global sums pass the raw audits.
Exact ordered gradients are compared within each algorithm, not across algorithms.

Local validation: ARM64 C++ reference suite and ARM64 AddressSanitizer/UBSan
suite pass (sanitizer run 8.23 s); rebuilt native Python suite has 173 passes and
4 optional benchmark skips; fallback suite runs 177 tests with 15 skips.
The first new C++ fixture incorrectly tried to finalize an empty complex and
was corrected before freezing. The prior x86-64 test process stalled while
reporting that failure; final C++ validation uses explicit ARM64 builds.

Machine: Apple M1 Max, 64 GiB, native ARM64, Apple Clang 15. Diagnostic main ran
2026-09-10 17:11:34–17:11:54 UTC; confirmation 17:12:07–17:12:27.
Regression-control sessions ran 17:13:00–17:14:05 and 17:14:15–17:15:20 UTC
(including their setup/build time). Raw JSON and retained diagnostic binaries
and ordered dumps are local artifacts in the parent work directory, not published
with this documentation. The appendix records the four full SHA-256 hashes.

Reaudit and check the generated appendix from the repository root:

```bash
LC_ALL=C python tools/profile_reduction_kernel_levels.py --audit \
  ../rk-level-profile-main.json ../rk-level-profile-confirmation.json
LC_ALL=C python tools/validate_simplicial_gradient_ab.py \
  ../rk-level-profile-controls-main.json ../rk-level-profile-controls-confirmation.json \
  --reversed-confirmation
LC_ALL=C python tools/render_rk_level_profile.py --check docs/rk_level_profile_tables.md
```
