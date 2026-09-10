# Frozen gradient scaling: measurement protocol

Protocol fixed before measurements on 10 September 2026. This pass changes no
algorithm, native data structure, scheduling policy, or timing boundary.
Headers are frozen at `232d65069125b2fd3959f9eee16da32f8ef0a33a`.

The question is whether RK's parallel advantage over our simplicial PLS is
robust to vertex-order seeds, and whether existing worker diagnostics support
prioritizing another scheduling change. F-Max remains a sequential control.
TTK, persistence, memory benchmarking, and scheduler implementation are out of scope.

## Inputs and fixed sample counts

- Freudenthal grids 5D/side 3, 6D/side 2, 7D/side 2, seeds 0–5: 18 inputs.
- Controls: 4D/side 4, volume16 and volume32, each seeds 0 and 2: six inputs.
- Main and reversed confirmation: all 24 inputs, worker counts 1/2/4/8.
- Eight blocks per input; every block includes every worker count, each with
  six repetitions of all three methods. Each method order occurs once per
  count/block. Two initial warmups per count.
- Worker positions and relative order of every pair are balanced within each
  study. Confirmation reverses input order and each block's worker order.
- Three separate PLS profiles, three coarse RK profiles, and three detailed
  RK profiles per input/count, after that input's headline measurements.
- A separate old/new boundary-index control check repeats the six controls
  with baseline `a0665429b6245e7803cc5dbac6cb37b63148a75f` and current headers,
  workers 1/8, six blocks of two repetitions, two warmups, three profiles,
  plus reversed confirmation. This is needed to distinguish scaling/noise
  from an optimization regression; a current-only sweep cannot do so.

Use fresh output files for all runs, retain every observation, and do not
choose extra repetitions based on favorable results. A separate small smoke
test validates the runner but is not performance evidence. No builds/tests
or concurrent benchmark tasks during measured intervals.

## Timing and interpretation

Common finalized native complex resident. Count fresh builders, every
algorithm-specific preparation step, worker creation/startup, local work,
replay and natural in-method cleanup. Returned builders and gradients remain
alive at the stop. Loading and construction are separate startup observations.

Worker ratios are medians of block ratios `t(w)/t(1)`; speedup is their
reciprocal. Method ratios pair methods within each block/count. Bootstrap
95% intervals resample eight blocks and describe within-session variation,
not uncertainty over machines. Report seed-level evidence and cross-seed
ranges; no pooled replication claim or multiple-testing-adjusted inference.
Use both sessions to assess repeatability. F-Max does not acquire workers;
variation across worker configurations measures context/noise, not F-Max scaling.

Exact ordered reference-dump hashes must match across studies. Native runs
check their ordered-gradient fingerprints against sequential references.
Check Euler identities, critical counts by dimension, raw coverage/order,
phase accounting, immutable RK work counts, source/header/input hashes, and
recompute summaries. Methods may legitimately have different critical counts.

Existing RK task durations cover each persistent worker task from its start
until it finishes claiming levels. They are elapsed lifetimes, not CPU time,
individual-level timings, or hardware utilization. Their sum overlaps in wall
time. `sum(task lifetimes)/(workers * max(task lifetime))` is only a lifetime
balance indicator. Simplex and level extrema have no paired task identity and
are not necessarily proportional to actual work. Large indivisible levels,
start delays, heterogeneous cores, and contention remain possible explanations.
Do not claim a cause this instrumentation cannot distinguish.

## Delivery and decision

Keep the established repository/Sphinx benchmark documentation as the primary
artifact. The technical report structure is summary, scoped evidence, protocol,
robustness/caveats, next step and open questions. Exact numeric tables are
preferred for seed/count/interval lookup; do not imply smooth interpolation
between worker counts. Validate rendered tables against raw data.

Recommend a scheduling experiment only if the diagnostics support it;
otherwise prioritize a measured preparation/local-work bottleneck. No new
algorithm optimization is authorized by this diagnostic pass itself.
