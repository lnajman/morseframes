# PLS priority-key storage: profiling and optimization

**Contiguous priority-key storage makes PLS faster without changing its ordered
gradient.** All 56 paired configuration medians improve across two A/B studies;
53 of their 56 within-session intervals are wholly below one. On the largest
eight-worker volume, paired median time reductions are **27.6% and 30.4%** in
the main and confirmation studies. TTK remains faster on every tested 2D/3D
configuration, while optimized parallel PLS now beats RK on both tested 7D inputs.

This is one bounded implementation change to the existing simplicial PLS/F-Max
variant. It does not change RK, native complex construction, persistence, the
matching algorithm, priority order, scheduling policy or ordered event replay.
The [previous direct PLS/TTK study](pls_ttk_gradient_benchmark.md) remains a
historical snapshot; its timings are not used as an optimization baseline here.

## Diagnosis and change

Revision `c64a403b5588304e9d196f253b21fbcdefdf74b9` adds fine diagnostic timers.
They split setup into output allocation, vertex ordering, executor startup,
storage initialization, owner/key construction and star partitioning. Local
wall time is split into scheduling and execution (including dispatch/wait on
the parallel path). Natural reverse-order destruction is measured through
checkpoints for events/indices, priority keys, membership storage, executor
shutdown and vertex maps. Cleanup remains inside the gradient call; profiling
does not explicitly free buffers earlier than ordinary execution.

The exploratory one/eight-worker profiles identify per-simplex priority-key
allocation/destruction as an avoidable cost, not thread shutdown or replay.
Revision `21a6c48500ff0d9a74e3c78f76e57e04ce2d2f4b` replaces the vector of
separately allocated keys with one rank array and an offset array. It preserves
each complete variable-length key, its descending rank order, lexicographic
comparison and simplex-ID tie break. Workers fill disjoint ranges; completed
keys are read-only during local matching. There is no fixed dimension cap.

Both A/B versions contain the same instrumentation. Performance calls have
internal profiling disabled; separate diagnostics identify the cost change.
Main-study diagnostic medians for volume n=32, eight workers, in milliseconds:

| Phase | Before | After |
| --- | ---: | ---: |
| Vertex preparation | 15.563 | 14.611 |
| Storage initialization | 1.852 | 4.474 |
| Owner/key construction | 8.674 | 4.233 |
| Star partitioning | 14.631 | 13.896 |
| Scheduling | 3.049 | 2.703 |
| Local execution and waiting | 28.288 | 27.035 |
| Priority-key cleanup | 27.517 | 0.573 |
| Total cleanup, including priority keys | 38.168 | 9.624 |

The confirmation key-cleanup medians are 27.437 / 0.486 ms and total cleanup
37.765 / 9.953 ms. Storage initialization becomes more expensive because it
also computes offsets and initializes contiguous storage; that cost is counted.
The total-cleanup row contains the key-cleanup row and must not be added to it.
These are medians of each seed's diagnostic medians, not an additive partition
of independently measured headline medians. Per-sample residuals remain explicit.
No new peak-memory measurements are made; fewer allocations is not a measured
peak-RSS claim.

## Controlled before/after results

The A/B timer starts with the same finalized native complex resident and includes
a fresh builder plus the entire gradient call, including preparation and cleanup.
Loading/construction are reported separately. The returned gradient and builder
remain alive at the stop boundary. Every method's own sequential reference is
validated before timing; full ordered old/new complex and gradient dumps agree.

Eight-worker main-study PLS algorithm times, milliseconds:

| Input | Before | After | Paired after/before |
| --- | ---: | ---: | ---: |
| 2D terrain, n=64 | 4.451 | 3.646 | 0.868 |
| 3D volume, n=16 | 12.323 | 8.832 | 0.713 |
| 3D volume, n=32 | 119.475 | 86.122 | 0.724 |
| 4D grid, n=4 | 2.022 | 1.525 | 0.756 |
| 5D grid, n=3 | 6.220 | 4.408 | 0.722 |
| 6D grid, n=2 | 2.488 | 1.778 | 0.747 |
| 7D grid, n=2 | 27.296 | 20.103 | 0.733 |

Times are medians of seed-specific medians. Ratios aggregate paired block-median
ratios across seeds, not quotients of the displayed marginal time medians.
Per-configuration main reductions range from 9.7–38.0% sequentially and
8.6–29.6% with eight workers; confirmation ranges are 5.8–27.7% and 15.1–34.9%.
Three optimization intervals overlap one: main volume32/seed0/one worker,
main terrain64/seed0/eight workers, and confirmation 7D/seed0/one worker.
Their positive median changes remain inconclusive individually.

On the 7D inputs, the main aggregate eight-worker times are **20.103 ms PLS**,
25.734 ms RK and 28.832 ms sequential F-Max. Confirmation times are
21.347 / 27.248 / 29.090 ms. Each 7D seed has a PLS/RK and PLS/F-Max interval
favoring PLS in both studies. This strengthens PLS as a higher-dimensional
baseline; it does not establish a universal ranking. Sequential F-Max still
has lower paired-median times than sequential PLS on every tested input.

## Fresh direct comparison with TTK

The resident-array harness measures F-Max, RK, PLS and TTK together, using
fresh native objects per call. Its primary comparison excludes native
construction but includes MF builder setup and all gradient preparation;
TTK vertex ordering, lower-star construction and matching are included.
The complete [timing contract](pls_ttk_gradient_benchmark.md#timing-contract)
is unchanged. In particular, MF returns an ordered Morse sequence and TTK its
native pairings; converting TTK output to a sequence is not timed.

Main-study algorithm times, milliseconds; F-Max is sequential in every row:

| Input | Workers | F-Max | RK | PLS | TTK | Paired PLS/TTK |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2D terrain, n=64 | 1 | 1.309 | 0.990 | 3.757 | 1.038 | 3.565 |
| 2D terrain, n=64 | 8 | 1.288 | 0.670 | 3.266 | 0.627 | 5.158 |
| 3D volume, n=16 | 1 | 9.606 | 4.026 | 18.163 | 10.936 | 1.686 |
| 3D volume, n=16 | 8 | 9.097 | 1.547 | 8.238 | 2.763 | 3.054 |
| 3D volume, n=32 | 1 | 124.463 | 46.318 | 202.361 | 98.581 | 2.064 |
| 3D volume, n=32 | 8 | 125.473 | 12.418 | 91.324 | 19.502 | 4.632 |

TTK has lower paired medians in all 24 configurations across both studies;
all 24 PLS/TTK intervals are above one. TTK is faster in 287 of 288 individual
paired repetitions; the remaining sample is retained. Confirmation aggregate
PLS/TTK ratios at one/eight workers are 3.706/5.088 (terrain), 1.666/3.102
(volume16), and 2.045/4.509 (volume32). RK has lower paired-median times than
all alternatives on every volume configuration in both direct studies.

Use the A/B studies to estimate the optimization effect, and the direct studies
to compare with TTK. They have different native-object lifetimes, allocator
histories and executables; do not mix their times into cross-harness ratios.
The resident construction/outer-phase tables are the generated
`resident_gradient_pls_arena_21a6c48_*_table.tex` fragments. For example,
eight-worker volume32 native construction is 558.587 ms for PLS and 36.522 ms
for TTK, kept outside the gradient comparison. Shared loading is 95.058 ms,
one read/parse per native invocation, not per algorithm or a cold-cache I/O test.

## Correctness, design and limitations

The independent rescanning oracle now covers dimensions 0–9, shuffled injective
functions, sparse vertex IDs, non-pure complexes, repeated builders, callbacks
and workers 1/2/4/8. Tests retain invalid-filtration rejection. Ordinary core
C++ and ASan/UBSan tests pass. The rebuilt native Python suite passes all
163 tests; the fallback suite runs 163 with 13 skips. Native benchmark tests
validate fine phase accounting and the TTK adapter in 1D–3D. Six historical
v2/v3 TeX tables reproduce byte for byte with the current renderer.

Every A/B run checks its own method's sequential fingerprint; complete old/new
reference dumps agree for all 14 inputs. Critical counts and Euler values are
checked independently. All four methods' critical counts agree on the six
TTK inputs. Higher-dimensional differences are accepted: on 4D seed0, F-Max
has counts [19, 46, 43, 16, 1], versus [19, 46, 42, 15, 1] for PLS and RK.
These counts and all ordered gradients are unchanged by the optimization.

All four studies start clean at `21a6c48` on the native ARM Apple M1 Max
(10 cores), Apple Clang 15, `-O3 -DNDEBUG`, with `OMP_WAIT_POLICY=PASSIVE`.
A/B snapshots use C++17 and pthreads. Direct TTK studies retain the pinned
unmodified classic backend `f4ffd1a1049d0ccf6e8f3eb4f7c096a6cc251ba0`, explicit
triangulation, OpenMP and its existing `TTK_ENABLE_KAMIKAZE` configuration.
No builds or tests overlap benchmark timing; desktop applications remain active.

Each A/B study uses 14 inputs, seeds 0/2, workers 1/8, six balanced old/new
blocks with two repetitions per block, two warmups and three separate profiles.
Every three-method execution permutation appears twice. The confirmation
reverses input and worker order. Each direct TTK study uses six inputs,
workers 1/8, twelve performance repetitions, four diagnostics and two warmups
per mode. Its four-round Williams design balances positions and directed
adjacent pairs; confirmation reverses input/workers and shifts the cycle by two.

Intervals resample paired blocks (A/B) or paired repetitions (TTK) within a
session. Six A/B blocks are a small sample; intervals do not measure
independent-machine uncertainty. No timing samples are removed. Unchanged
F-Max/RK controls vary across versions (ratios approximately 0.89–1.13);
their changes are not credited to this PLS optimization. The grids are small
synthetic workloads, not a representative application corpus or a controlled
dimension-only scaling experiment. These are gradient results, not persistence
benchmarks.

## Reproduction and evidence

Use the deterministic 2D/3D inputs and pinned native TTK build from the
[previous resident-array study](pls_ttk_gradient_benchmark.md#reproduction-and-artifacts).
From the repository root:

```sh
OMP_WAIT_POLICY=PASSIVE python3 tools/benchmark_simplicial_gradients.py \
  --baseline c64a403 --candidate 21a6c48 \
  --inputs ../rk-ab-inputs/terrain-n64-seed0.txt ../rk-ab-inputs/terrain-n64-seed2.txt \
    ../rk-ab-inputs/volume-n16-seed0.txt ../rk-ab-inputs/volume-n16-seed2.txt \
    ../rk-ab-inputs/volume-n32-seed0.txt ../rk-ab-inputs/volume-n32-seed2.txt \
  --grids 4:4 5:3 6:2 7:2 --seeds 0 2 --input-dir ../higher-dimensional-inputs \
  --workers 1 8 --blocks 6 --repeats 2 --warmups 2 --profiles 3 \
  --memory-repeats 0 --output ../pls-arena-ab-main.json
```

Repeat with `--workers 8 1 --reverse-inputs` and a fresh confirmation output.
For TTK, rebuild `morseframes_resident_gradient_benchmark` at `21a6c48` and run
the previous resident-array commands with output paths
`../pls-arena-ttk-main.json` and `../pls-arena-ttk-confirmation.json`; retain
`--include-pls`, twelve repetitions, four diagnostics, two warmups, and the
confirmation's reversed input/workers and `--order-offset 2`.

Optional fine profiles are exported as `pls_profile_seconds`. Their children
overlap coarse parents and are validated separately by `tools/pls_phase_profile.py`.
Historical records without this optional field remain valid. The standalone
`tools/profile_process_lower_stars.py` can profile a frozen native simplicial
worker via `--binary`, `--inputs`, `--workers` and a fresh `--output`.

Raw study artifacts are retained locally; repository tools and selected tables
are tracked. SHA-256 hashes:

- A/B main: `2194660e3791c9c8a76a3cfc1b8e08e0dfb296be5a766f1bcc3d2b0e1b211449`.
- A/B confirmation: `7ab603f5366b1efd6dc248354b5594d6f4ff9e6eec0dc83da8e662545b3422c3`.
- TTK main: `9c29b87cedec6d95b3a731bfeb7b0850e31094d5a1275f019a23578779908da1`.
- TTK confirmation: `d7b8abbc0f81172f74218663bb7e1b81678d3ae85457d7d5a47e78f031a0eeae`.
- Candidate headers: `ae60887d5e62946c02e086f2999ec8d6cd52916fd265b0ea800019db34246494`.
- Resident native binary: `cc21b0266f25fcbf5b5bc6ad8b263f45c516adaeb2e32f298a7b580a1ad92df3`.

The post-run audit recomputes summaries, intervals, order balance and phase
partitions from the raw samples, checks hashes and exact identities, and
verifies the reversed confirmation schedules. Earlier studies are preserved.
