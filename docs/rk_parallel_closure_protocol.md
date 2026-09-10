# Parallel RK closure preparation: frozen experiment protocol

## Change and correctness contract

Prepare independent, previously unseen sparse facet closures in private task
buffers; merge in canonical facet order before incidence/local kernels read them.
Retain the exact same closures, including inactive faces, and exact event order.
Use the existing bounded executor without new pools or threads. Task-indexed
scratch is necessary because cooperative waits may execute another level task on
the same thread. Both submission and execution failures must drain live captures.

Packed levels, precomputed-cache use, sequential construction, boundary-index
construction, facet discovery/execution and replay algorithms are unchanged.
The level scheduler can now run independent closure tasks alongside other levels;
it still disables nested facet operations. The low-level boolean is renamed to
`allow_parallel_facet_operations` to make that distinction explicit. The new
execution option `parallel_closure_preparation` can disable the experiment.

The work gate uses an upper bound on nonempty simplex faces, saturated at 8,192
per unprepared facet. Submit at most `min(workers, facets, estimated_work/8192)`
tasks; fewer than two stays serial. This is a heuristic, not a measured-time
prediction. Contiguous facet chunks and ordered arena concatenation preserve
determinism. Do not tune this gate after viewing confirmatory results without
recording a separate candidate/study.

Native tests cover 1/2/4/8 workers, enabled/disabled closure preparation, coarse
and detailed execution, shared high-dimensional faces, several simultaneous large
levels, scratch reuse on later graph levels, and exact agreement with an
independently prepared eager-cache sequence. Exception tests cover ordinary and
traced level tasks. Sanitizer checks precede performance measurement.

Parallel closure elapsed time includes dispatch, waits and merge. Its worker
traversal/sort/materialization fields are **cumulative**, not additional elapsed
phases. Existing serial closure child fields keep their previous meaning.
Profiling output and validators retain backward compatibility with old artifacts.

## Fixed comparison

Baseline: the previous completed branch state, `b2f2ac7` (headers identical to
`8edb1ec`). Freeze the candidate and runner in a clean commit before timing.
Compile the same current driver against both snapshots with ordinary-only mode,
`clang++ -std=c++17 -O3 -DNDEBUG -pthread`. Retain binaries, header hashes and
raw symbol evidence confirming the absence of traced RK instantiations.

Three modes: baseline; candidate with closure parallelism disabled; candidate
with it enabled. Both candidate modes use the **same resident process and same
executable** through `run_closure_serial` and `run`. This directly compares the
runtime switch while controlling executable/process identity, but does retain
cross-call allocator/cache history. All modes compute F-Max/PLS/RK, so unchanged
algorithms remain controls. No profiling is interleaved with timed calls.

Inputs: 7D/side-2 seeds 3 and 1, 6D/side-2 seeds 4 and 1, volume-32 seed 0,
volume-16 seed 0. These selected synthetic cases are not a representative corpus.
Workers: 1 and 8. For every configuration, use two independent fresh process
pairs, each with two warmups per mode and six blocks of six repetitions. Cover
all six mode permutations and all six algorithm permutations. Alternate launch
and mode order between process pairs. Confirmation reverses input/worker order
and pair schedules. This yields four independent process pairs per configuration
across the two sessions. Complete all compilation/tests before timing.

Compare complete gradient time from a finalized complex in memory: fresh builder,
all algorithm-specific preparation, worker pools, local work, replay and in-method
cleanup are included. Loading/native construction are recorded separately;
returned-object destruction and reference checking are outside timing. No
persistence computation. TTK and existing paper-facing rankings are not rerun.

## Acceptance and evidence

Retain every run, exact reference dump, critical count, phase and preparation
observation. Audit source/header/input/binary provenance, configurations, order,
method permutations, timing partitions and exact gradient agreement. Independently
recompute all consequential ratios and intervals from raw runs.

Report candidate-on/baseline, candidate-off/baseline and candidate-on/off paired
block ratios, plus each fresh process pair separately. Block-bootstrap intervals
are conditional on a process pair. Across four independent pairs, report the
descriptive median and range of pair estimates; do not pool blocks and label that
an across-process confidence interval. Neither an interval spanning one nor a
small observed change demonstrates equivalence.

Primary benefit targets are the difficult 7D seed-3 and 6D seed-4 cases at eight
workers. Retain default enablement only with a repeatable useful benefit (target
median on/off ratio below 0.95) and without a consistent material regression on
volume or easy-case controls. A consistent penalty above 5% needs investigation
or narrower/default-disabled gating. Do not discard unfavourable controls; if the
evidence is inconclusive, report that qualification and avoid a universal claim.
Any gate/default follow-up is a separate, recorded decision, not a replacement
for the frozen results.

The data-validation skill informs the independent recomputation and distinction
between within-process and between-process evidence. The deliverable remains the
requested repository implementation and its reproducible benchmark notes.
