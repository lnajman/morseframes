# Selective parallel RK closures: benefit retained, still opt-in

## Decision

The selective gate retains a useful benefit on the difficult 7D case while
avoiding closure dispatch on the previously penalized smaller levels. At eight
workers, 7D seed 3 uses **14.9% less complete-gradient time** versus the same
executable with closure parallelism off, and **13.3% less** versus the previous
implementation with its opt-in feature off. All four process-pair estimates and
their conditional block-bootstrap intervals favour the selective-on path for
both comparisons.

The experiment uses one fixed, data-informed size threshold. The feature remains
**opt-in**, regardless of this selected-case follow-up's outcome; it is not a
default-policy change or a universal scheduling result.

The previous [broad-gate experiment](rk_parallel_closure_benchmark.md) is preserved
as a separate study. This candidate requires at least **32,768 simplices in an
individual level**, as well as the original minimum estimated work per task.
The threshold was committed before measurements and was not retuned afterwards.

## Eight-worker results and acceptance

Each entry is the median of four independent process-pair estimates, followed by
their **observed range**, not a confidence interval. Ratios below one mean less
time. On/off uses one resident candidate process; baseline comparisons use
separate old/new processes. Every primary measured call is unprofiled.

| Input | Selective on / off | Selective on / baseline |
| --- | --- | --- |
| 7D seed 3 | 0.8511 [0.8461, 0.8689] | 0.8672 [0.8418, 0.8795] |
| 7D seed 1 | 0.9894 [0.9535, 1.0254] | 1.0084 [0.9660, 1.0162] |
| 6D seed 4 | 0.9633 [0.9489, 0.9830] | 0.9887 [0.9626, 0.9977] |
| 6D seed 1 | 1.0143 [0.9759, 1.0874] | 1.0253 [0.9351, 1.1660] |
| Volume-32 seed 0 | 1.0118 [0.9861, 1.0258] | 1.0126 [0.9719, 1.0355] |
| Volume-16 seed 0 | 0.9922 [0.9082, 1.0705] | 0.9983 [0.8971, 1.1206] |
| 7D seed 6 | 1.0024 [0.9661, 1.0643] | 1.0248 [0.9554, 1.0773] |
| 7D seed 7 | 1.0005 [0.9666, 1.0284] | 0.9862 [0.9391, 1.0457] |

The primary target passes the precommitted median-on/off < 0.95 criterion with a
benefit in both sessions. No tested on/off control has a four-pair median penalty
above 5%, at either worker count. In particular, the previous 7D seed-1 and 6D
seed-4 median on/off penalties do not recur in this selective-gate study.
This does **not** establish equivalence: individual control estimates still vary
substantially. For example, volume-16 at eight workers ranges from 0.9082 to
1.0705. All inputs except 7D seed 3 dispatch zero new closure tasks, so their
apparent gains must not be attributed to closure parallelism.

The large case's four candidate-off median gradient times are 37.51, 37.52,
38.45 and 38.57 ms; selective-on gives 32.68, 31.76, 32.51 and 32.63 ms, in
main-pair-0/1 then confirmation-pair-0/1 order. The paired ratios above are
computed from individual blocks, not by dividing these displayed medians.
The [numerical appendix](rk_selective_closure_tables.md) retains every method,
worker count, process pair, phase and conditional interval.

## Sequential and version warnings remain

Passing the runtime-toggle controls is not proof that the new code has zero
cost. In sequential 7D seed 1, candidate-off/baseline has median 1.0515, range
[0.9951, 1.0988]; two of four conditional intervals are above one. In sequential
7D seed 7 it is 1.0534 [0.9042, 1.1235], with three intervals above one and one
below one. These version/process effects are not resolved by the gate.

Sequential 6D seed 4 also has a selective-on/baseline median of 1.0622
[1.0305, 1.0963], although only one conditional interval excludes one. Its
on/off median is 1.0120 [0.9727, 1.0397], and no closure tasks run. These warnings
remain visible rather than being classified as harmless noise or evidence of a
universal algorithmic regression. The detailed off/baseline comparisons and
unchanged F-Max/PLS controls remain in the appendix.

## Dispatch checks

The independent diagnostic pass uses the same frozen candidate, but a full
profiling binary rather than the ordinary-only timing binaries. All ten detailed
level traces per configuration agree on the task counts below. These counts are
cumulative over closure batches, not simultaneous worker counts.

| Input | Largest level | Closure tasks at 8 workers |
| --- | ---: | ---: |
| 7D seed 3 | 69,072 | 34 |
| 7D seed 1 | 10,308 | 0 |
| 6D seed 4 | 5,274 | 0 |
| 6D seed 1 | 1,576 | 0 |
| Volume-32 seed 0 | 75 | 0 |
| Volume-16 seed 0 | 75 | 0 |
| 7D seed 6 | 17,156 | 0 |
| 7D seed 7 | 15,236 | 0 |

Only the 69,072-simplex level dispatches. At one worker all configurations have
zero closure tasks. Every traced dispatched level meets the threshold; several
smaller simultaneous levels do not qualify by summing their sizes.

Seeds 6 and 7 were selected and generated before inspecting their profiles or
timings. Both turn out to be below the threshold. They add controls, **not an
independent large-level benefit test**. The successful-case selection therefore
remains a material limitation for generalization.

## Implementation and correctness

```cpp
morseframes::ReductionKernelExecutionOptions options;
options.policy = morseframes::ReductionKernelExecutionPolicy::Parallel;
options.max_workers = 8;
options.parallel_closure_preparation = true; // Still opt-in; default false.
options.parallel_closure_min_level_size = 32768; // Inclusive experimental gate.
```

Setting `parallel_closure_min_level_size` to zero restores the original broad
gate. This does not remove the existing 8,192-estimated-entries-per-task rule or
the requirement for at least two tasks. The size check precedes the work scan
and private task-buffer allocation. It does not change closure contents,
traversal/merge ordering, the level scheduler, boundary indexing or replay.
Packed, cached and graph-specialized paths remain unchanged.

The benchmark's per-level profiler now explicitly enables the same opt-in policy
as `run` and `rk_plain`. Previously it could inherit the library's disabled
default while other experimental benchmark modes enabled the feature. Ordinary-only
timing binaries exclude this profiling code. This correction does not change the
library default or rewrite historical evidence.

Native tests compare complete ordered sequences with an independently prepared
eager-cache oracle. They exercise zero/default/equal-to-size/size-plus-one
thresholds, 1/2/4/8 workers, both switch states, coarse/detailed builds, shared
faces, concurrent levels, scratch reuse and unchanged work counters. Existing
exception-draining tests remain active. AddressSanitizer passed with leak
detection disabled because the macOS runtime does not support it. The Python
suite passed with **187 passed, 4 skipped**, including native on/off and
global/per-level gate checks.

## Measurement contract and limitations

The [protocol](rk_selective_closure_protocol.md) fixes eight synthetic inputs,
1/8 workers, two fresh process pairs per configuration per session, and a second
session with reversed input/worker/pair orders. Every pair has two warmups per
mode and six blocks of six repetitions, covering all six mode and algorithm
permutations. There are four independent pair estimates per configuration and
10,368 timed gradient calls per session. No builds, tests or diagnostic profiling
overlap the primary timing sessions.

The baseline is `12ee94a` with its opt-in closure feature **disabled**, using
`--baseline-closure-serial`. Candidate-on and candidate-off share an executable
and resident process; the former enables the selective gate. This is not a direct
new selective-on versus old broad-on comparison. Candidate-off/baseline retains
version/process-context warnings rather than assuming zero overhead.

Report the median and observed range of four process-pair estimates; their range
is not a confidence interval. Individual 95% paired-block bootstrap intervals
(2,000 resamples, seed 0) are conditional on one process pair. Neither a small
change nor an interval spanning one demonstrates equivalence. Cross-call
allocator/cache history and active-desktop interference remain possible.

The native ARM64 machine is the same Apple M1 Max (8 performance plus 2
efficiency cores, 64 GiB), macOS 14.6, Apple Clang 15. It was not a dedicated idle
benchmark host: desktop activity was observed before timing. Recorded one-minute
host load averages changed from 38.8 to 9.8 in the main session and 53.7 to 9.8
in confirmation. These include the benchmark itself and other activity; they
are not an isolated background-CPU measure or a correction for interference.

Loading/native construction are separate per-process observations. Gradient time
includes a fresh builder, algorithm-specific preparation, worker pools, local
work, replay and in-method cleanup. Returned-object destruction and reference
checking are outside timing. The numerical appendix retains absolute gradient
times, separate builder/kernel/total ratios, F-Max/PLS controls and critical counts.
Detailed phase traces are diagnostic evidence, not interchangeable with
uninstrumented timing. No persistence or TTK rerun; paper-facing rankings are not
replaced by this scheduling experiment.

Peak RSS is not remeasured here. The previous broad-on study's memory increase
on the large 7D case remains a relevant caveat, not a new selective-gate estimate.
Broader large-level inputs and intermediate worker counts remain necessary
before considering default enablement. The data-validation skill informs these
controls, independent calculations and limitations within the repository workflow.

## Frozen evidence

- Baseline: `12ee94a` (closure parallelism disabled).
- Candidate, protocol and runner: `1cf3066dff82daaea7a0855d6bdb41d4aa4c7388`.
- Main timing SHA-256: `687eaa1d106b910194ddf8d45a9b9a668be2259e0bcf0d2fb743e3b4f9c3fc8b`.
- Confirmation SHA-256: `0cc0ba274cc342fcce65032cd550cf8e4375c1894733ad9c37aebcdc23bdf267`.
- Diagnostic profile SHA-256: `de68ca303f32bed00265583b9e63d8f36f1a93ba11364a964cedb024413f1de2`.

Both timing artifacts passed source/header/input/binary and schedule audits;
all **1,728 phase summaries** were independently recomputed from raw repetitions.
All ordered sequential reference dumps agree across builds, processes and
sessions, and every measured call checks its ordered fingerprint against that
reference. Critical counts and Euler/step accounting agree throughout.

Raw JSON, header snapshots, ordinary binaries, symbol listings and exact reference
dumps remain in the parent `work` directory locally, not in a public archive.
The diagnostic artifact is `rk-selective-closure-profile.json`; its detailed
volume traces make it approximately 3 GiB. Its internal source/header/binary
provenance was audited against the frozen commit. All task-gate claims above were
also checked against its complete detailed level traces.

## Reproduction

Use a clean checkout of `1cf3066`. The six original input paths match the prior
experiment. Generate the added seeds with the frozen helper before timing:

```python
from pathlib import Path
import sys
sys.path.insert(0, "tools")
from benchmark_simplicial_gradients import grid_input
folder = Path("../rk-selective-closure-inputs")
folder.mkdir(exist_ok=True)
for seed in (6, 7):
    path = folder / f"grid-d7-n2-seed{seed}.txt"
    contents = grid_input(7, 2, seed)
    if path.exists():
        assert path.read_text() == contents
    else:
        path.write_text(contents)
```

Run the following with a fresh output path, then repeat with `--reverse` and a
fresh confirmation path. The baseline switch is essential to this comparison:

```sh
LC_ALL=C OMP_WAIT_POLICY=PASSIVE python3 tools/benchmark_rk_parallel_closure.py \
  --baseline 12ee94a --candidate 1cf3066 --baseline-closure-serial \
  --inputs ../higher-dimensional-inputs/grid-d7-n2-seed3.txt \
           ../higher-dimensional-inputs/grid-d7-n2-seed1.txt \
           ../higher-dimensional-inputs/grid-d6-n2-seed4.txt \
           ../higher-dimensional-inputs/grid-d6-n2-seed1.txt \
           ../rk-ab-inputs/volume-n32-seed0.txt ../rk-ab-inputs/volume-n16-seed0.txt \
           ../rk-selective-closure-inputs/grid-d7-n2-seed6.txt \
           ../rk-selective-closure-inputs/grid-d7-n2-seed7.txt \
  --output ../rk-selective-closure-main.json
```

Audit the retained runs and verify the generated appendix:

```sh
python3 tools/benchmark_rk_parallel_closure.py --audit \
  ../rk-selective-closure-main.json ../rk-selective-closure-confirmation.json
python3 tools/render_rk_parallel_closure.py --inputs \
  ../rk-selective-closure-main.json ../rk-selective-closure-confirmation.json \
  --report rk_selective_closure_benchmark.md \
  --protocol rk_selective_closure_protocol.md \
  --output docs/rk_selective_closure_tables.md --check
```

The optional diagnostic uses `tools/profile_reduction_kernel_levels.py` with
the same eight `--inputs`, `--workers 1 8`, and a fresh `--output`. It compiles
its own verified full binary from the clean frozen checkout; never substitute
untracked prebuilt binaries for the recorded header snapshot. Do not run it
alongside performance measurements.
