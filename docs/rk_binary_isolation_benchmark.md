# RK binary isolation: the earlier slowdown is not consistently reproduced

## Outcome

**The larger volume-32 warning does not repeat consistently with fresh processes
and ordinary-only binaries.** The reversed session gives an ordinary RK new/old
difference of +1.9% sequentially and approximately zero at eight workers. Both
intervals span no change. The sequential point estimate remains positive in both
sessions, so a small cost is **not excluded** and this is not an equivalence claim.

A post-hoc control runs the **exact same retained executable** in independent
processes. Its sequential RK paired estimates range from -4.8% to +9.2% across
three fresh pairs. Thus non-version variation can reach the size of the earlier
warning on this machine. This does not identify a CPU, allocator, clock or binary
layout mechanism, nor does it prove that all old/new differences are noise.

No RK algorithm, scheduler, library header or gradient semantics changed in this
pass. There is no demonstrated repeatable ordinary-only regression that justifies
an algorithm fix here. We can resume a **separate, measured experiment in
within-level closure preparation**, retaining this volume control and adding
independent process replication to future performance guards. Keep existing
paper-facing benchmark results frozen: this is a diagnostic study, not a new
TTK/PLS/RK ranking or a claim of zero tracing-refactor overhead.

## Volume-32 results

These are candidate/baseline paired-block ratios for complete gradient time;
greater than one means slower. Brackets are 95% block-bootstrap intervals.
Ordinary executables exclude profiling code at compile time. Full executables
contain diagnostic commands, but **all measured calls are unprofiled**.

| Executable | Workers | Main | Reversed confirmation |
| --- | ---: | --- | --- |
| Ordinary-only | 1 | 1.0524 [0.9625, 1.1300] | 1.0188 [0.9825, 1.0410] |
| Ordinary-only | 8 | 0.9713 [0.9451, 1.0514] | 1.0000 [0.9619, 1.0136] |
| Full | 1 | 0.9994 [0.9771, 1.0333] | 1.0077 [0.9833, 1.0285] |
| Full | 8 | 0.9819 [0.8680, 1.0281] | 1.0342 [0.9583, 1.0994] |

For all six input/worker configurations, ordinary RK new/old intervals include
one in both sessions. These within-process intervals do not establish equivalence
or incorporate between-process variation. In particular, the ordinary sequential
volume-32 estimates are +5.2% and +1.9%, not two estimates of exactly zero.

The [full numerical appendix](rk_binary_isolation_tables.md) retains all methods,
all four comparisons, builder/kernel/total phases, critical counts, preparation
observations, ordered-reference hashes, executable hashes and timing windows.

## Controls prevent a premature causal claim

F-Max and PLS algorithm sources do not change between these header revisions.
Nevertheless, some executable/process comparisons change their timings too:

- In main, 7D seed 3 at one worker has full-binary new/old ratios of 1.1113 for
  RK, 1.1290 for F-Max and 1.1212 for PLS; each interval is above one. None of
  these slowdowns repeats with an interval wholly above one in confirmation.
- On volume-16 at one worker, ordinary new/old F-Max and PLS ratios are below
  one with intervals below one in both sessions. Do not call these RK-induced
  algorithmic improvements: those two algorithms were not modified.
- Within unchanged **baseline headers**, volume-16 sequential full/ordinary RK
  is 0.9487 [0.9227, 0.9891] and 0.9095 [0.8849, 0.9416]. Executable mode
  and process conditions can matter even without the tracing refactor. These
  measurements do not separate their individual contributions.

No unfavourable control is discarded. The appendix also retains the main-only
volume-32 baseline full/ordinary penalty and the larger volume-16 F-Max changes.
The executable-mode contrast bundles available code and compile-time
specialization; it is not an isolated instruction-cache experiment.

## Same-executable A/A control

This control was added **after** the two frozen four-build sessions, prompted by
changes in unmodified methods. It reuses the main study's candidate ordinary
binary, without recompilation. Labels A and B are separate processes, not code
versions. Each row below is a different fresh process pair; do not pool these
within-pair bootstrap intervals into an across-process interval.

| Workers | Pair | RK B/A [95% block interval] |
| --- | ---: | --- |
| 1 | 0 | 0.9516 [0.8923, 1.0270] |
| 1 | 1 | 0.9866 [0.9544, 1.0179] |
| 1 | 2 | 1.0923 [1.0000, 1.1314] |
| 8 | 0 | 1.0132 [0.9715, 1.0441] |
| 8 | 1 | 1.0006 [0.9679, 1.0213] |
| 8 | 2 | 1.0007 [0.9783, 1.0104] |

The sequential pair-2 lower endpoint rounds to 1.0000; do not interpret that
rounded value as an exact boundary test. F-Max and PLS also show A/A differences,
including intervals excluding one in some pairs. These cannot be code-version
effects because both processes execute the same binary. Three process pairs per
worker count are a diagnostic check, not an estimate of the machine's full
performance distribution. One A/A construction observation is 1.553 seconds
versus about 0.53 seconds in several neighbouring processes; it is retained,
outside gradient timing, without attributing a cause or excluding the run.

## Scope, correctness and provenance

Measurements ran on 2026-09-10 on the Apple M1 Max, native ARM64, macOS 14.6,
Apple Clang 15, `-std=c++17 -O3 -DNDEBUG -pthread`. There are ten logical CPUs;
the tested worker counts are 1 and 8. F-Max remains sequential at either setting.
Inputs are the selected synthetic volume-32/seed-0, volume-16/seed-0 and
7D/side-2/seed-3 complexes, not a representative corpus or cross-machine study.

The [precommitted protocol](rk_binary_isolation_protocol.md) crosses old/new
headers with ordinary/full executables. For every input/worker configuration,
four fresh processes construct references, run two unprofiled warmups, then eight
balanced blocks of six repetitions, covering all algorithm order permutations.
Confirmation reverses input, worker, launch and block order. There are 3,456
timed gradient calls per four-build session. A/A adds 864 calls, four balanced
blocks per independent pair. No builds or tests overlap timing.

Loading and native construction are separately recorded. The comparison starts
from the finalized complex in memory and includes fresh builder, all
algorithm-specific preparation, pools, local work, replay and in-method cleanup.
Returned-object destruction and reference validation are outside the timer.
There is no persistence computation and no timed profiling. The preparation
observations are not a statistically controlled construction benchmark.

All builds agree on complete ordered reference dumps, critical counts, complex
identity and Euler accounting. All six A/A references match the original
volume-32 reference. Native ordinary/full tests also check 1/8-worker execution
and rejection of profiling commands by the ordinary executable. Raw symbol lists
have no traced RK instantiation in ordinary builds; the candidate full binary
has 32 matching trace symbols. Symbol counts are diagnostic evidence, not a
proof that the executable-mode change is a single-variable intervention.

Frozen header revisions:

- Baseline: `876354e777fc491c63d4274e78e6af62045c91cd`.
- Candidate: `8edb1ec7ee459c6036a4817f7e2e52907b20b5cc`.
- Four-build runner: `9414d4674fe4a0f1cec86233f3d6dc92518695ec`.
- A/A runner: `18a402f80b24a68e26041ec30fff7c67aef0779a`.

Raw JSON and `.artifacts/` directories are retained in the parent `work`
directory locally; they are **not committed or publicly archived**. Audits use
recorded git revisions, input/header hashes, retained binaries/symbols and dump
hashes. The independently implemented renderer recomputes all **486 phase
summaries**, including medians, quartiles, block ratios and bootstrap intervals,
from raw repetitions without importing the runners' statistics functions.

The metric-diagnostics and validation skills informed the separation of code
version, executable capability and process controls, plus independent numerical
QA. They do not supply the experimental evidence or establish a causal mechanism.

## Reproduce and audit

From a clean checkout of runner commit `9414d46`, with the three saved input
arrays available, run the main study and repeat with `--reverse` and a fresh
output filename:

```sh
LC_ALL=C OMP_WAIT_POLICY=PASSIVE python3 tools/benchmark_rk_binary_isolation.py \
  --baseline 876354e777fc491c63d4274e78e6af62045c91cd \
  --candidate 8edb1ec7ee459c6036a4817f7e2e52907b20b5cc \
  --inputs ../rk-ab-inputs/volume-n32-seed0.txt \
           ../rk-ab-inputs/volume-n16-seed0.txt \
           ../higher-dimensional-inputs/grid-d7-n2-seed3.txt \
  --workers 1 8 --blocks 8 --repeats 6 --warmups 2 \
  --output ../rk-binary-isolation-main.json
```

From the clean A/A runner commit, point to the completed main study and its
retained executable:

```sh
LC_ALL=C OMP_WAIT_POLICY=PASSIVE python3 tools/benchmark_rk_same_binary.py \
  --binary-study ../rk-binary-isolation-main.json \
  --input ../rk-ab-inputs/volume-n32-seed0.txt \
  --output ../rk-binary-isolation-aa.json
```

The current tools can audit the original frozen artifacts without requiring the
current checkout's driver to equal the historical source snapshot:

```sh
python3 tools/benchmark_rk_binary_isolation.py --audit \
  ../rk-binary-isolation-main.json ../rk-binary-isolation-confirmation.json
python3 tools/benchmark_rk_same_binary.py --audit ../rk-binary-isolation-aa.json
python3 tools/render_rk_binary_isolation.py --inputs \
  ../rk-binary-isolation-main.json ../rk-binary-isolation-confirmation.json \
  ../rk-binary-isolation-aa.json --output docs/rk_binary_isolation_tables.md --check
```

Reruns require fresh output and artifact paths. Absolute paths in raw provenance
refer to the originating local workspace; relocating artifacts requires an
explicit relocation procedure rather than editing the frozen evidence silently.
