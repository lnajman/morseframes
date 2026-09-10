# RK ordinary-binary regression isolation protocol

## Question and fixed scope

The [per-level diagnostic study](rk_level_profile_benchmark.md) left a warning:
volume-32 seed 0 had slower unprofiled RK medians after the optional trace
refactor. The original measurements used diagnostic-capable binaries and mixed
unprofiled and profiled calls within each resident process. This pass asks whether
the warning repeats in **ordinary-only executables**, and whether it depends on
executable capability or process history. It does not change library headers,
RK kernels, scheduling, native construction, or gradient semantics.

Baseline headers: `876354e777fc491c63d4274e78e6af62045c91cd`.
Candidate headers: `8edb1ec7ee459c6036a4817f7e2e52907b20b5cc`.
The current branch's headers are identical to the latter candidate; this pass
changes only benchmark/test/documentation tooling.

## Four-build comparison

Compile the same native driver with each header set, once normally and once with
`MORSEFRAMES_BENCHMARK_ORDINARY_ONLY`. This yields old/new ordinary-only and
old/new diagnostic-capable (full) binaries. Every measured call is unprofiled,
even in full binaries. No profiling commands run between configurations either.

The ordinary-only macro removes profiling command handlers and diagnostic
helpers, and fixes metrics to null. It retains the same `Run` layout, ownership,
returned objects and timing boundaries. Symbol evidence must show no traced RK
instantiation or `run_levels` function in ordinary-only binaries. This contrasts
a bundle of executable differences (available code and compile-time
specialization), not a single CPU-cache mechanism. Do not attribute a difference
to instruction layout, allocation history or scheduling without further evidence.

For each input/worker configuration, launch all four processes afresh. Each builds
the same native complex and validated sequential references, then performs two
unprofiled warmups. Keep all four resident while timing one at a time. This removes
profile-call and previous-worker-count history from the new sessions; reference
validation still establishes a common, nonempty allocation history.

Eight blocks use four cyclic executable orders and their reversals, balancing
position and pair order. Each executable runs six repetitions per block, covering
all six F-Max/PLS/RK order permutations. Reversed confirmation changes input,
worker, process-launch and within-block executable order. Run both 1 and 8 workers;
F-Max remains sequential. Primary input is volume-32 seed 0; controls are
volume-16 seed 0 and 7D/side-2 seed 3. These are selected synthetic controls,
not a representative population of complexes.

Loading/common construction remain separate. Gradient timing includes fresh
builder, algorithm-specific setup, pools, local work, replay and natural in-method
cleanup. Returned builder/gradient teardown and reference validation are outside
timing. No persistence, precomputed algorithm cache or timed profiling. Complete
all builds before measurements; do not run tests or CPU-heavy analysis concurrently.

## Comparisons and acceptance criteria

Report paired-block candidate/baseline ratios separately for ordinary and full
binaries, plus full/ordinary within each header version. Keep full raw repetitions,
builder/kernel/total phase times, metadata, critical counts and bootstrap intervals.
Do not form a paired comparison from a ratio of two independently reported medians.
Uncertainty comes from paired blocks in one session, not across machines.

Recompute key ratios independently from raw runs. Check exact ordered reference
dumps across all four builds and sessions, all configurations/repetitions,
phase accounting, process launch order, method permutations, input/header/source
hashes, and absence of traced instantiations in ordinary-only binaries. Preserve
binaries, symbol listings and dumps so the compilation comparison can be audited.

An interval spanning one is inconclusive, not proof of equivalence. A single
non-repeating slowdown does not establish a regression. A repeatable ordinary-only
slowdown keeps the warning open and warrants checking generated code before further
optimization. If only full binaries show the warning, executable effects become
a stronger lead; this still does not prove a specific mechanism. If neither
repeats, say the prior warning was not reproduced under the fresh-process protocol,
not that it never occurred. Report every control, including slower outcomes.

`tools/benchmark_rk_binary_isolation.py` implements the frozen run and audit.
This protocol is committed before performance measurement. The metric-diagnostics
and validation skills informed the separate comparisons, independent recomputation,
and distinction between observations and causal explanations; the repository's
benchmark artifacts remain the controlling evidence and requested deliverable.
