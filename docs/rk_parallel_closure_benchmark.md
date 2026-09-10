# Parallel RK closure preparation: useful on one large level, still opt-in

## Decision

**Implemented and correctness-tested, but disabled by default.** Parallel closure
preparation helps the difficult 7D seed-3 case, but the initial work gate is too
broad: the easier 7D case and difficult 6D case regress when the switch is enabled.
The two precommitted benefit targets were not both met, so the broad gate is not
adopted as a default optimization. Existing level/facet parallelism remains active;
only the new closure-parallel feature is opt-in.

At eight workers, the median of four independent process-pair estimates gives
**13.3% less gradient time** on 7D seed 3 versus the same binary with the switch
off, and **12.9% less** versus the previous implementation. The 7D seed-1 control
instead takes 6.8% more time with the switch on. The 6D seed-4 median is 10.8%
slower, with substantial variation between process pairs.

The next experiment should test a more selective work/level-size gate, preserving
these unfavourable controls. The largest level contains 69,072 simplices in the
successful case, versus 10,308 in 7D seed 1 and 5,274 in 6D seed 4. This is a
useful distinction to test, not a validated threshold or a proven scheduling cause.
Allocation/copying and dispatch/wait costs remain possible contributors. No new
gate was tuned after inspecting these results.

## Eight-worker comparison

Ratios below one mean less complete-gradient time. Each entry is the median of
four paired-block estimates, followed by their **observed range**, not a
confidence interval. Both candidate modes run in the same process and executable.
The old/new comparison uses separate processes. Every measured call is unprofiled.

| Input | On / off | On / previous implementation |
| --- | --- | --- |
| 7D seed 3 | 0.8667 [0.8376, 0.8692] | 0.8713 [0.8533, 0.8866] |
| 7D seed 1 | 1.0684 [1.0445, 1.1216] | 1.0124 [1.0025, 1.0921] |
| 6D seed 4 | 1.1081 [0.9848, 1.2434] | 1.0837 [1.0042, 1.1838] |
| 6D seed 1 | 0.9699 [0.9517, 1.0386] | 1.0247 [0.9678, 1.0793] |
| Volume-32 seed 0 | 0.9941 [0.9936, 1.0012] | 1.0041 [0.9575, 1.0374] |
| Volume-16 seed 0 | 0.9617 [0.9529, 0.9903] | 0.9863 [0.9833, 1.0118] |

All four 7D seed-3 on/off block-bootstrap intervals lie below one, as do all four
on/previous intervals. In 7D seed 1, three of four on/off intervals lie above one.
For 6D seed 4, three point estimates are slower and one is faster; two intervals
lie above one and two span it. Do not describe its penalty as four identical
outcomes. Individual 95% within-process block intervals, all three comparison
pairs, F-Max/PLS controls, and builder/kernel/total phases are retained in the
[numerical appendix](rk_parallel_closure_tables.md).

Volume levels have at most 75 simplices and use the unchanged packed path;
the new sparse closure dispatch does not run there. Accordingly, an apparent
volume on/off improvement must not be attributed to parallel closure work.
At one worker the switch also does not dispatch tasks. Differences in these
controls illustrate residual measurement/context variation, not new algorithmic
speedups. Six inputs on one machine do not establish a universal improvement.

## Sequential version controls remain qualified

Turning the feature off is a default-policy decision, **not proof that introducing
its code has zero performance cost**. In sequential 7D seed 3, candidate-off /
previous has median 1.0617 and pair range [1.0101, 1.1212]. Its two main-session
intervals are above one; the two confirmation intervals span one. In sequential
volume-32 the corresponding median is 1.0570, range [0.8819, 1.1076], and main
and confirmation differ strongly. These warnings remain visible rather than
being classified as harmless noise or a demonstrated universal regression.

The same-executable runtime toggle is the cleanest evidence for this feature's
effect. It still retains cross-call allocator/cache history and interleaves all
three methods in balanced orders. Old/new process effects, code generation and
context variation are not individually isolated. Do not replace existing
paper-facing timing tables with this experiment; TTK is not rerun or re-ranked.

## Implementation and use

Sparse closures are still prepared lazily, only for newly exposed facets.
Independent tasks read the immutable boundary index, maintain private visited
flags/traversal vectors/entry arenas, and write disjoint facet-range slots.
The coordinator concatenates contiguous facet chunks in their original order
and fixes the range offsets before incidence/local reduction reads them.
Inactive faces remain in each immutable closure, just as before.

Tasks share the existing bounded executor; there is no extra pool or oversubscription.
Scratch is indexed by task, not OS thread, because cooperative waits may execute
another level task on the same thread. Submitted tasks are drained on exceptions
before captured scratch/event buffers can disappear, in both ordinary and traced
level builds. The existing level order, boundary-index construction, facet
reduction and replay algorithms are unchanged.

The prototype gate submits at most one task per available executor worker and
requires an estimated 8,192 closure entries per task, with at least two tasks.
The estimate is a saturated simplex-face upper bound; it can overestimate actual
same-level work. This heuristic is retained for reproducibility, not endorsed
as the final adaptive scheduler. The low-level boolean formerly named
`allow_intra_level_parallelism` is now `allow_parallel_facet_operations`;
independent closure tasks are controlled separately by execution options.

The new C++ opt-in is:

```cpp
morseframes::ReductionKernelExecutionOptions options;
options.policy = morseframes::ReductionKernelExecutionPolicy::Parallel;
options.max_workers = 8;
options.parallel_closure_preparation = true; // Experimental; default is false.
auto sequence = builder.build_flooding_reduction_kernel_with_execution_options(
    options, [](const auto&, const auto&) {});
```

No new Python strategy is introduced. The native benchmark explicitly selects
on or off through `run` and `run_closure_serial`; therefore its experimental
`run` mode must not be confused with the library's final default.

Detailed profiling reports parallel closure wall time, merge time, task/batch
counts and cumulative worker traversal/sort/materialization times. Worker times
may overlap and must not be added to elapsed closure time. Existing serial child
phases preserve their previous meaning. Validators accept historical profiles
without these fields and reject incomplete new global/per-level coverage.

## Memory cost

Supplemental measurements reuse the exact frozen ordinary binaries in 36 fresh
processes, three observations per build/input/worker combination. The entries
below are median **whole-process peak RSS after gradient construction**, in MiB.
They include the complex and any construction high-water mark, not just live
gradient allocations. Returned gradient/builders are still alive at the sample.

| Input | Workers | Previous peak MiB | Experimental-on peak MiB |
| --- | ---: | ---: | ---: |
| 7D seed 3 | 1 | 84.83 | 85.05 |
| 7D seed 3 | 8 | 90.95 | 104.77 |
| 6D seed 4 | 1 | 8.83 | 9.17 |
| 6D seed 4 | 8 | 10.16 | 10.42 |
| Volume-32 seed 0 | 1 | 214.48 | 214.45 |
| Volume-32 seed 0 | 8 | 214.62 | 214.66 |

The successful 7D case thus also has a visible memory tradeoff: approximately
13.8 MiB more process peak at eight workers. These observations are not a
dedicated allocation profile or cross-machine memory guarantee. The runner and
its exact script contents/hash are retained with the supplemental raw data.

## Timing scope, verification and frozen evidence

The [protocol](rk_parallel_closure_protocol.md) was committed before measurement.
Main and reversed sessions each have 24 fresh process pairs: six selected inputs,
1/8 workers, two independent pairs per configuration. Each pair uses two warmups
per mode and six balanced blocks of six repetitions, covering all mode and
algorithm permutations. There are 7,776 timed gradient calls per session.
Confirmation reverses input/worker order and pair schedules. No builds/tests
overlap the timed sessions. The machine is the same Apple M1 Max, native ARM64,
macOS 14.6, Apple Clang 15 with `-std=c++17 -O3 -DNDEBUG -pthread`.

Loading and native construction are reported separately. Gradient time includes
fresh builder, algorithm-specific preparation, pools, local work, replay and
in-method cleanup; returned-object destruction and reference validation are
outside the timer. No persistence computation or timed instrumentation.

Complete ordered sequential-reference dumps agree across all processes/builds
and both sessions, with identical recorded critical counts. Every measured mode
and worker-count run checks its ordered fingerprint against that reference.
Native tests additionally compare shared
high-dimensional/multilevel cases against an independently built eager-cache
sequence, with 1/2/4/8 workers, both switch states, scratch reuse and failures.
AddressSanitizer passed with leak detection disabled because the macOS runtime
does not support it. All 1,296 phase summaries are independently recomputed from
raw repetitions, not just regenerated with the runner's statistics helper.

- Baseline headers: `b2f2ac7` (identical to `8edb1ec`).
- Frozen prototype and runner: `e25c86016447645d44b41069828d77adce28d0a3`.
- Main SHA-256: `a9534795e9a48e20c483e41a21ae72eec3737a080c6ad2881e05cb885cc4a9d7`.
- Confirmation SHA-256: `45b75e3ca7526cc98098e4c4df5fc91c9922d68c7d090b3208efeb44ce0c7800`.
- Supplemental memory SHA-256: `7a650d41bfc6b97ae079e92302eadf4b49af40f384da51dc67497e918c49d741`.

The default was changed to false **after** these frozen measurements; no work
gate was changed. Raw JSON, binaries, symbol listings and reference dumps remain
in the parent `work` directory locally, not yet in a public artifact archive.
An extra unfrozen exploratory profile is retained but excluded from the evidence:
its worktree header metadata did not describe its prebuilt binary's header snapshot.
The headline studies and supplemental memory data use verified frozen binaries.

The data-validation skill led to independent recomputation, same-executable
controls, explicit process-level replication and retention of negative outcomes.
The implementation is ready for scoped experiments; the broad default scheduler
needs revision before adoption.

## Reproduction

From a clean checkout of `e25c860`, run the following, then repeat with `--reverse`
and a fresh confirmation output filename:

```sh
LC_ALL=C OMP_WAIT_POLICY=PASSIVE python3 tools/benchmark_rk_parallel_closure.py \
  --baseline b2f2ac7 --candidate e25c860 \
  --inputs ../higher-dimensional-inputs/grid-d7-n2-seed3.txt \
           ../higher-dimensional-inputs/grid-d7-n2-seed1.txt \
           ../higher-dimensional-inputs/grid-d6-n2-seed4.txt \
           ../higher-dimensional-inputs/grid-d6-n2-seed1.txt \
           ../rk-ab-inputs/volume-n32-seed0.txt ../rk-ab-inputs/volume-n16-seed0.txt \
  --output ../rk-parallel-closure-main.json
```

The current checkout can audit and regenerate the original results using recorded
git snapshots and retained artifacts:

```sh
python3 tools/benchmark_rk_parallel_closure.py --audit \
  ../rk-parallel-closure-main.json ../rk-parallel-closure-confirmation.json
python3 tools/render_rk_parallel_closure.py --inputs \
  ../rk-parallel-closure-main.json ../rk-parallel-closure-confirmation.json \
  --output docs/rk_parallel_closure_tables.md --check
python3 tools/check_rk_parallel_closure_memory.py \
  --study ../rk-parallel-closure-main.json --output ../fresh-closure-memory.json
```

Output paths must be fresh. Absolute paths in raw provenance refer to the
originating workspace; do not silently rewrite frozen evidence when relocating it.
