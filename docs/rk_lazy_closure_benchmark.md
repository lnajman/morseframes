# RK: preparing closures only for exposed facets

**On-demand closure preparation reduces eight-worker RK time by about 13–20%
on the two tested 7D inputs, in two clean A/B studies.** Ordered gradients and
critical counts are unchanged. RK is now statistically indistinguishable from
PLS on 7D seed 0, but remains slower on seed 2. Sequential F-Max still leads
clearly in 7D. The 5D parallel control does not establish an improvement.

This is one bounded RK implementation change. PLS, F-Max, native construction,
the filtration, facet scheduling, reduction decisions, and event replay are
unchanged. The [PLS key-storage study](pls_key_arena_benchmark.md) remains the
latest direct TTK comparison; **TTK was not rerun for this RK change**.

## Diagnosis: closure preparation is a major 7D cost

The initial diagnostic pass covers 5D and 7D, seeds 0 and 2, at one and eight
workers. In sequential 7D, closure preparation takes about 55–61 ms, while the
level-processing wall time is about 100–112 ms. Global setup and replay each
take less than 0.4 ms. Local reduction is the next substantial component.
Worker task durations are uneven, but this pass changes neither their schedule
nor the number of tasks. Summed worker durations are wall-clock measurements,
not pure CPU time or a direct measure of memory bandwidth.

Two preliminary closure-union alternatives were evaluated and retained as
unsuccessful or superseded exploratory evidence. Sorting occupied bitset words
instead of simplex indices did not improve RK. Retaining compressed words for
repeated face-closure unions helped sequential 7D, but parallel gains were
uneven. Neither variant is in the final implementation.

The retained change is simpler: **do not prepare a sparse closure until its
simplex actually becomes a facet**. The old implementation prepared every
simplex's same-level closure at the start of its level. The new implementation:

1. Initializes empty ranges for large, uncached levels.
2. After facet discovery, traverses the same-level boundary of each newly
   exposed facet, visiting each face once and sorting by canonical bucket rank.
3. Retains that immutable closure for subsequent rounds of the level.

Inactive faces remain in the stored closure and are filtered at the existing
consumption points. Thus cell order, incidence, local pair selection, and replay
are preserved. All closure writes finish before local facet tasks read them;
concurrent levels have separate scratch. Small packed levels (at most 128
simplices), graph-only handling, and precomputed closure caches retain their
existing paths. There is no new dimension limit or dense closure matrix.

Separate detailed diagnostic medians from the main study, milliseconds:

| Input | Workers | Closure before | Closure after |
| --- | ---: | ---: | ---: |
| 7D, seed 0 | 1 | 67.564 | 46.580 |
| 7D, seed 2 | 1 | 58.492 | 31.834 |
| 7D, seed 0 | 8 | 79.118 | 46.826 |
| 7D, seed 2 | 8 | 69.090 | 40.402 |

These diagnostics support a reduction in closure-preparation work, not an
additive explanation of independently measured headline medians. Parallel
closure timers sum across levels/tasks and may exceed elapsed time. Detailed
core and local-reduction timers are nested inside facet execution and must not
be added to that parent. Detailed runs also materialize diagnostic results;
they are not performance runs. The candidate's closure timer includes the new
per-round preparation checkpoints, even on paths that return immediately.
No peak-memory study or bandwidth measurement is claimed.

## Timing contract and controlled comparison

Both versions use the same benchmark driver and instrumentation. Baseline
`ff32e4427c700523a6fb4adc9fed61222a62d180` has the diagnostic commands but the
old RK algorithm. Candidate `8957918333d1c8c056206fe5a83b07374ec21c3c` contains
the on-demand change and its tests. The worktree was clean for both studies.

Measurements on 10 September 2026 use this Apple M1 Max, native ARM64,
Apple Clang 15, C++17, `-O3 -DNDEBUG -pthread`. Each study has 14 inputs:
2D terrain n=64, 3D volumes n=16/32, and 4D/5D/6D/7D grids of side 4/3/2/2,
each with seeds 0 and 2. Each input is tested at one and eight workers.

Timing starts with the common finalized native complex resident. It includes
the fresh builder and complete gradient call: method-specific preparation,
thread startup, local computation, replay, and natural in-method cleanup.
The returned gradient and builder remain alive at the stop boundary. There is
no persistence computation. Loading and native construction are separate.

Each configuration has two warmups, six paired blocks of two repetitions, and
three separate PLS profiles plus three coarse and three detailed RK profiles
per version. Old/new order is balanced; all six method orders occur equally.
The confirmation reverses input and worker-count order. No compilation or
tests run during timing. Exact old/new complex and complete ordered gradient
dumps are compared before measurements; timed/profiled runs match each method's
own validated sequential reference.

Ratios below are medians of paired block ratios, not quotients of displayed
marginal time medians. Intervals are 95% paired-block bootstrap intervals within
one session, not independent-machine confidence intervals.

| 7D seed | Workers | Main after/before [interval] | Confirmation after/before [interval] |
| --- | ---: | --- | --- |
| 0 | 1 | 0.830 [0.819, 0.873] | 0.822 [0.765, 0.919] |
| 2 | 1 | 0.760 [0.721, 0.771] | 0.756 [0.737, 0.763] |
| 0 | 8 | 0.851 [0.830, 0.880] | 0.873 [0.858, 0.876] |
| 2 | 8 | 0.796 [0.752, 0.835] | 0.809 [0.758, 0.886] |

All eight 7D intervals favor the candidate. Sequential reductions are about
17–24%; eight-worker reductions are about 13–20%. These findings reproduce on
the frozen candidate, independently of the exploratory variant selection.

### Controls: useful gains are not universal

In 5D, sequential paired reductions are 5.2–8.5%, with all four intervals below
one. At eight workers the main ratios are 0.964/1.004 and the confirmation
ratios 0.987/1.009; all four intervals include one. There is no established
parallel gain in this control.

In 6D, eight-worker paired reductions range from 9.8–16.0%. Three of four
intervals favor the candidate; main seed 0 includes one. Sequential 6D improves
in median, but the main seed-2 interval is very wide, [0.481, 1.720]. That
outlier-sensitive configuration is retained, not discarded or relabeled a win.

The 2D/3D/4D controls have mixed small changes. In the main run, volume16/seed0
and 4D/seed2 are slower sequentially with intervals above one. Neither slowdown
has an interval wholly above one in the confirmation, but their positive
median ratios remain visible in the raw results. Unchanged F-Max/PLS controls
also fluctuate. This supports neither a universal speedup nor proof that every
small-input regression has been excluded. There is no configuration with an
RK slowdown interval wholly above one in both studies.

## Where RK now stands against PLS and F-Max

Eight-worker configurations, gradient time in milliseconds. F-Max remains
sequential in these same three-method comparisons.

| Study | 7D seed | RK before | RK after | PLS | F-Max |
| --- | ---: | ---: | ---: | ---: | ---: |
| Main | 0 | 26.502 | 22.103 | 23.018 | 31.123 |
| Main | 2 | 25.756 | 20.330 | 17.592 | 28.468 |
| Confirmation | 0 | 25.741 | 22.429 | 22.238 | 27.305 |
| Confirmation | 2 | 25.517 | 21.035 | 17.442 | 26.298 |

For seed 0, RK/PLS paired ratios are 1.007 [0.902, 1.049] and
1.028 [0.959, 1.090]: no resolved winner. For seed 2 they are
1.141 [1.082, 1.244] and 1.144 [1.081, 1.305]: PLS retains an advantage.
Eight-worker RK beats sequential F-Max on both 7D seeds in both studies;
one-worker RK remains substantially slower than F-Max.

These are small synthetic high-dimensional grids, not a representative corpus
or evidence that RK is intrinsically the fastest algorithm. A sensible next
investigation is the remaining local reduction work and level imbalance on
7D seed 2, retaining the present timing boundary and low-dimensional controls.

## Construction and critical counts remain separate

Main-study candidate loading/construction times, milliseconds, shown as the
median of the two seed observations per family. These are single startup
observations per input/version, not a dedicated construction benchmark.

| Input | Loading | Native construction |
| --- | ---: | ---: |
| 2D terrain, n=64 | 3.476 | 10.451 |
| 3D volume, n=16 | 13.166 | 50.764 |
| 3D volume, n=32 | 91.561 | 530.062 |
| 4D grid, n=4 | 1.225 | 9.838 |
| 5D grid, n=3 | 2.584 | 45.639 |
| 6D grid, n=2 | 0.651 | 18.731 |
| 7D grid, n=2 | 4.059 | 296.706 |

All method-specific preparation, including the newly deferred closure work,
remains inside the reported gradient time. Nothing is moved into native
construction or an untimed warmup cache.

Critical counts by dimension are unchanged in every old/new comparison, as
are the complete ordered gradients. The three methods agree on 13 of 14
inputs. The existing 4D/seed0 difference remains legitimate: F-Max produces
`[19, 46, 43, 16, 1]`, while PLS and RK produce `[19, 46, 42, 15, 1]`.
The 7D counts are `[5, 10, 7, 1, 0, 0, 0, 0]` for seed 0 and
`[3, 2, 0, 0, 0, 0, 0, 0]` for seed 2, for all three methods. Raw identities
retain every count vector; all satisfy the complex's Euler characteristic.

## Verification and reproduction

The C++ suite passes normally and under address/undefined-behavior sanitizers.
New tests compare against the unchanged eager cache in dimensions 4, 7, and 9,
with plateau, tied, and injective vertex filtrations, sparse vertex IDs,
shared facets, graph-only levels, and one/two/four/eight workers. Existing
packed/sparse boundary, callback-failure, and concurrency tests also pass.
The rebuilt native Python suite has 160 passes and four optional resident
benchmark skips (no TTK adapter supplied); fallback runs 164 tests with 13
skips. No TTK result is inferred from these tests.

Run from the repository root, using fresh output paths:

```sh
OMP_WAIT_POLICY=PASSIVE LC_ALL=C python3 tools/benchmark_simplicial_gradients.py \
  --baseline ff32e44 --candidate 8957918 \
  --inputs ../rk-ab-inputs/terrain-n64-seed0.txt ../rk-ab-inputs/terrain-n64-seed2.txt \
    ../rk-ab-inputs/volume-n16-seed0.txt ../rk-ab-inputs/volume-n16-seed2.txt \
    ../rk-ab-inputs/volume-n32-seed0.txt ../rk-ab-inputs/volume-n32-seed2.txt \
  --grids 4:4 5:3 6:2 7:2 --seeds 0 2 --input-dir ../higher-dimensional-inputs \
  --workers 1 8 --blocks 6 --repeats 2 --warmups 2 --profiles 3 --rk-profiles 3 \
  --memory-repeats 0 --output ../rk-lazy-closure-ab-main.json
```

Repeat with `--workers 8 1 --reverse-inputs` and the fresh confirmation output.
Audit both frozen studies without new timings:

```sh
LC_ALL=C python3 tools/validate_simplicial_gradient_ab.py \
  ../rk-lazy-closure-ab-main.json ../rk-lazy-closure-ab-confirmation.json \
  --reversed-confirmation
```

The audit recomputes summaries, intervals, phase partitions, order balance,
critical counts, unchanged detailed RK work counters, and source/input hashes.
It verifies the reversed confirmation schedule and both frozen header sets.
The standalone `tools/profile_reduction_kernel_simplicial.py` also accepts
`--binary`, `--inputs`, `--workers`, and a fresh `--output` for diagnostic-only
passes; its coarse residual includes uninstrumented cleanup and bookkeeping.

Raw artifacts remain local, with repository tools and selected results tracked.
SHA-256 identities:

- Main: `a641876d1bdc2f9de97b9a87e3b7baf23bbcf4a2c066519d24c3605688117cd2`.
- Confirmation: `4aebb95ef6ee3ceb103e9ae03943fa49476e7eae97c9827bc620ce7bbe2e1cbd`.
- Candidate headers: `d6ec19547461a2dcb88869ad5e0144cadcafbb96557a9559ab318b1b8ea7848f`.
- Initial `rk-7d-profile-baseline.json`: `480ffb050efbcc6b7a02a80ed52eb1f03289864866f35a4181e48f5f9109b080`.
- Rejected `rk-word-union-exploratory.json`: `e7baaeae8ff1fa0eea7c6b969905b12dc22dcfe4232ce72274305f19aae0bae6`.
- Superseded `rk-compressed-union-exploratory.json`: `294ba787976e768f71e6d1bb12aeea0c7143044246ad3808372f5e2e70039b6b`.
- Retained-candidate exploration `rk-lazy-closure-exploratory.json`: `5b4593dff6f04f947e4f1c82b5cb5c990b0a3252ef78648a942ab4b9f60c8ae7`.
