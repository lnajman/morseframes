# RK per-level diagnostic protocol

This pass diagnoses the repeatable large-level imbalance in the
[frozen worker-scaling study](rk_worker_scaling_benchmark.md). It does not change
the scheduler, level order, local kernels, or replay order.

## Scope and instrumentation

The optional C++ `build_flooding_reduction_kernel_with_level_profile` entry point
returns the ordinary gradient and fills a `ReductionKernelLevelProfile` supplied
by the caller. Ordinary entry points instantiate the non-tracing implementation:
no trace vector or per-level clock calls. The new diagnostics are not exposed in
Python's gradient API; the standalone native benchmark emits JSON.

Each row identifies a filtration level, persistent task/scratch slot, simplex
count, event count, start offset and duration. Task IDs are **not physical thread
or CPU-core identities**. Timestamp origin is the start of the level-processing
phase. A task's intervals must not overlap. Copying a completed row is outside its
duration but inside the phase wall time. Trace allocation is inside setup time.
The caller must not read or reuse the profile while the build is running. The
completion flag remains false on failure, including replay callback exceptions;
submitted level tasks are drained before an exception is propagated.

With the builder's ordinary detailed-metrics setting, each row additionally
copies the existing phase times and selected work counters. Without that setting,
these fields stay zero. Closure, facet discovery, incidence/essential work, facet
execution, aggregation and merge are disjoint phase measurements. Core and local
reduction are children of facet execution, not additional elapsed phases. Closure
initial preparation, traversal, sorting and materialization are children of
closure; packed preparation and boundary indexing are nested inside initial
preparation. If facets run in parallel within a single level, child task times
can overlap. The focused injective multi-level corpus uses level parallelism,
with nested facet parallelism disabled by the existing implementation.

## Frozen experiment

Primary cases: 7D side-2 seed 3 and 6D side-2 seed 4. Matched controls: seed 1 of
each geometry. These are synthetic triangulated hypercubes with shuffled
injective vertex values, not representative real-world datasets. Do not equate
vertex IDs with filtration-level ranks: map using the input values.

For 1, 4 and 8 workers, run five modes: unprofiled RK, existing coarse/detailed
RK, and coarse/detailed level traces. Two warmups per mode precede ten measured
blocks. Five cyclic orders and their reversals balance each mode's position and
each pair's order. Repeat with inputs, worker counts and block sequence reversed.
Use `tools/profile_reduction_kernel_levels.py`; it requires a clean frozen commit,
builds the native driver before measuring, and retains raw rows, binary, reference
dumps, source/input hashes, metadata, critical counts and paired-overhead summaries.

Loading and common native construction remain separate. Fresh builder, pool
startup, all algorithm-specific preparation, local work and replay stay within
gradient time. Diagnostic timings are **not headline algorithm comparisons**.
Report coarse/detailed overhead relative to unprofiled RK, and incremental trace
overhead relative to matching existing profiles. Bootstrap intervals describe
the small local paired-block experiment, not population guarantees.

The longest measured level is a lower bound on the level-phase duration for an
unchanged indivisible-level implementation. Compare its start offset and duration
with phase wall time to distinguish late scheduling from a long indivisible
task. This is a diagnostic bound, not a prediction of a new scheduler's speed:
changing scheduling can change contention, cache behavior and scratch reuse.

## Validation and regression guard

C++ checks compare exact ordered gradients across cached/uncached complexes,
flat/tied/injective weights, multiple dimensions and worker counts. Native JSON
checks verify level coverage, vertex-rank mapping, event/simplex accounting,
nonoverlapping task timelines, phase nesting and per-level/global sums.
The runner compares each timed gradient fingerprint to a validated reference;
complete ordered reference dumps are retained and compared across sessions.

Separately use `tools/benchmark_simplicial_gradients.py` to compare unprofiled
baseline `876354e777fc491c63d4274e78e6af62045c91cd` with the frozen candidate.
Use both difficult cases, both matched controls, 5D seed 0, 4D seed 2 and volume-32
seed 0; 1/8 workers, six blocks of two repetitions, two warmups, three profiles,
and a reversed confirmation. Audit exact old/new reference dumps and unchanged
work counters. Retain slower and inconclusive controls as well as improvements.
No tests, compilation, or CPU-heavy analysis may overlap performance measurements.
