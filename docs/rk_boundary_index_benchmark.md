# RK: level-local boundary indexing

**Sequential RK is 9–15% faster on both tested 7D inputs across two frozen
A/B studies.** Eight-worker paired medians improve by 6–15%, but the
confirmation intervals include one, so the parallel improvement remains
unresolved. In the same-run method comparisons, RK beats PLS on both 7D
inputs at one and eight workers in both studies.

This pass adds closure diagnostics and one optimization: prepare same-level
boundary indices once per large uncached level, then reuse them across facet
closures. PLS, F-Max, scheduling, construction, canonical ordering, local
pair selection, and replay are unchanged. The preceding
[eligible-candidate study](rk_eligible_candidates_benchmark.md) remains
historical evidence. TTK is not rerun; the [PLS key-storage study](pls_key_arena_benchmark.md)
remains the latest direct TTK comparison on its documented earlier RK revision.

## Diagnosis: traversal, not sorting

The new detailed profile separates initial setup, packed preparation,
sparse traversal/deduplication, sorting, and storing the cells.
**Packed preparation and boundary-index preparation are children of initial
setup**, not additional phases to add to it. Sparse traversal, sorting, and
storage are separate children of the overall closure timer. Unattributed
loop/clock overhead remains in that parent timer.

The initial frozen diagnostic revision is `f469358`. It uses three repetitions
after two warmups on six inputs at one/eight workers. Selected sequential
medians show why traversal was chosen:

| Initial profile | Initial setup (ms) | Traversal (ms) | Sorting (ms) | Storage (ms) |
| --- | --- | --- | --- | --- |
| 5D, seed 0 | 0.277 | 2.565 | 0.652 | 0.580 |
| 5D, seed 2 | 0.286 | 2.612 | 0.663 | 0.613 |
| 7D, seed 0 | 0.381 | 26.751 | 3.345 | 4.520 |
| 7D, seed 2 | 0.518 | 23.874 | 3.027 | 4.151 |

The baseline visits 7,711,028 and 7,197,925 complex-boundary records while
constructing the two 7D inputs' sparse cells. Frontier/output buffers grow only
3/5 and 3/6 times respectively in the sequential diagnostic calls; sorting and
storage are smaller measured components than traversal. These are instrumented
observations, not unprofiled performance claims. Timing and counting each facet
can perturb the relative costs.

## One bounded implementation change

For each uncached level whose bucket exceeds 128 simplices, build an offset
array and a contiguous array of **same-level boundary ranks**. This filters
out other levels and translates simplex IDs to canonical bucket indices once.
Subsequent facet traversals use those ranks directly.

Boundary order is preserved within each simplex. The same visited flags remove
duplicates, and the existing sort restores canonical cell order. Inactive
faces remain in the immutable boundary table and closures; no face is removed
merely because an earlier reduction made it inactive. Face-first validation is
performed while preparing the table. The existing externally cached path is
unchanged and remains an independent reference.

Packed levels and graph-only levels do not build this table. Scratch buffers
belong to one level worker, are reused across levels, and are rebuilt before
a new sparse level reads them. The coordinator still finishes preparation
before local facet tasks start. There is no new dimension bound or scheduling
change.

This trades extra temporary storage for fewer repeated indirect accesses:
`O(n + m)` indices for a level with `n` simplices and `m` same-level boundary
incidences. Each worker retains buffer capacity for reuse during the call.
No peak-memory benchmark was run, so this is a structural bound, not a measured
memory improvement.

For 7D seed 0, preparation reads 991,415 original boundary records and stores
803,969 local boundary entries across levels; subsequent closure traversals
read 5,810,365 local entries. Seed 2 uses 989,688 original reads, 802,552 stored
entries, and 5,440,380 subsequent local reads. Do not compare only the smaller
preparation count against the entire old traversal: the local reads still
exist. Sparse cell counts and contents, duplicate hits, sorting input sizes,
and all local reduction search/scan counts remain unchanged.

## Frozen total-gradient comparison

Baseline `a0665429b6245e7803cc5dbac6cb37b63148a75f` contains the diagnostics
but no boundary index. Candidate `cb04b26253b578ee1be292623b6e38c82022a783`
adds the index. The optional diagnostic groups let the current driver compile
against historical headers; this was checked against `d4c7212`.

Both studies ran on 10 September 2026 on the Apple M1 Max, native ARM64,
Apple Clang 15, C++17, `-O3 -DNDEBUG -pthread`. Main:
16:02:50–16:04:32 UTC. Confirmation: 16:04:43–16:06:24 UTC.
Each used clean sources and the same 14 inputs: terrain n=64, volumes n=16/32,
and 4D/5D/6D/7D grids of side 4/3/2/2, seeds 0 and 2, at one/eight workers.

Each configuration uses two warmups and six paired blocks of two repetitions.
Version order is balanced, all six method orders occur equally, and confirmation
reverses input and worker order. Each version/configuration also has three
separate PLS profiles and three coarse plus three detailed RK profiles.
No compilation or tests run during timings. All raw outliers are retained on
this interactive workstation.

Start with the same finalized native complex resident in memory. Include a
fresh builder, all algorithm-specific preparation (including the new index),
worker startup, local work, replay, and natural in-method cleanup. Returned
gradients and builders remain alive at the stop boundary. Loading and native
construction stay separate; no persistence or untimed algorithm-specific cache
is introduced.

Ratios below are medians of paired-block ratios, not ratios of marginal
medians. Their 95% bootstrap intervals describe within-session variability,
not independent-machine uncertainty.

| 7D seed | Workers | Main after/before [95% interval] | Confirmation after/before [95% interval] |
| --- | --- | --- | --- |
| 0 | 1 | 0.848 [0.834, 0.868] | 0.876 [0.867, 0.888] |
| 2 | 1 | 0.913 [0.889, 0.926] | 0.857 [0.841, 0.897] |
| 0 | 8 | 0.853 [0.784, 1.042] | 0.914 [0.848, 1.005] |
| 2 | 8 | 0.865 [0.846, 0.920] | 0.941 [0.832, 1.050] |

All four sequential 7D intervals favor the candidate: 8.7–15.2% reductions.
The parallel medians favor it, but both seed-0 intervals and the seed-2
confirmation interval include one. Do not claim a repeatable parallel
speedup from these results alone.

### Controls and deviations

No 5D or 6D configuration has an interval wholly below one in both studies.
Some individual runs favor the candidate; others are inconclusive.
The only interval wholly above one is confirmation volume16/seed0 with one
worker: 1.031 [1.002, 1.104]. Main gives 0.971 [0.848, 1.017].
That slowdown is retained, not dismissed because it does not repeat.

Confirmation 4D/seed2 at eight workers has a large paired ratio:
1.733 [0.990, 2.285], versus 0.946 [0.718, 1.018] in main.
The unchanged PLS control also rises in that confirmation configuration
(paired median 1.292). This is compatible with session interference, but does
not establish the cause. Confirmation 5D/seed0 at eight workers is also noisy:
0.999 [0.777, 2.515].

The previous volume32/seed2 sequential risk remains worth watching: ratios
1.025 [0.987, 1.062] and 1.061 [0.949, 1.101]. Both medians are slower, even
though both intervals include one. These controls prohibit a blanket
no-regression claim. The strongest evidence is the sequential 7D result.

## Does preparation pay for itself?

Separate sequential diagnostic medians, including the new table's preparation:

| Study | 7D seed | Total closure before → after (ms) | Traversal before → after (ms) | Index preparation after (ms) |
| --- | --- | --- | --- | --- |
| Main | 0 | 39.636 → 29.577 | 28.862 → 9.740 | 9.148 |
| Main | 2 | 31.936 → 28.074 | 22.408 → 8.997 | 8.780 |
| Confirmation | 0 | 36.779 → 29.303 | 26.199 → 9.553 | 9.531 |
| Confirmation | 2 | 33.399 → 25.974 | 23.474 → 8.612 | 7.940 |

The new index costs roughly 8–10 ms in these calls, but total closure
preparation still falls on all four sequential 7D comparisons.
The nested index timer must not be added again to total closure time.

Parallel diagnostics sum work across tasks rather than measuring global
elapsed time. In particular, main 7D/seed0 cumulative closure time rises from
60.703 to 73.756 ms even while the unprofiled total median falls. This
discrepancy is retained; cumulative task times, instrumented overhead and
scheduling variability are not an additive prediction of wall-clock gains.

## Current RK / PLS / F-Max comparison

Eight-worker configurations; F-Max remains sequential:

| Study | 7D seed | RK before (ms) | RK after (ms) | PLS (ms) | F-Max (ms) |
| --- | --- | --- | --- | --- | --- |
| Main | 0 | 19.999 | 16.939 | 21.292 | 27.426 |
| Main | 2 | 17.048 | 14.838 | 17.808 | 27.772 |
| Confirmation | 0 | 18.797 | 17.025 | 21.499 | 27.516 |
| Confirmation | 2 | 16.711 | 15.335 | 17.289 | 25.423 |

Paired RK/PLS ratios for seed 0 are 0.788 [0.681, 0.865] and
0.777 [0.697, 0.870]; for seed 2, 0.843 [0.734, 0.891] and
0.898 [0.846, 0.981]. All four favor RK. The separate sequential RK/PLS
comparisons also favor RK on both seeds in both studies.
Eight-worker RK beats sequential F-Max here, but sequential F-Max is still
roughly 2.1–2.2 times faster than sequential RK.

These synthetic cases do not establish that RK is generally the fastest
algorithm. Do not multiply ratios from different optimization studies to
claim a measured cumulative gain.

## Remaining phases and construction

Coarse RK diagnostic medians for 7D seed 2, one worker:

| RK phase, 7D seed 2, one worker | Main before (ms) | Main after (ms) | Confirmation before (ms) | Confirmation after (ms) |
| --- | --- | --- | --- | --- |
| Builder | 0.184 | 0.193 | 0.152 | 0.222 |
| Workspace/pool setup | 0.170 | 0.154 | 0.157 | 0.167 |
| Level processing | 61.158 | 59.339 | 68.949 | 52.831 |
| Replay | 0.335 | 0.327 | 0.419 | 0.339 |
| Unattributed remainder | 0.146 | 0.161 | 0.156 | 0.163 |
| Total profiled gradient | 61.995 | 60.569 | 69.819 | 53.973 |

These are separate instrumented calls, not the unprofiled headline timings.
Each raw row's accounting is validated; marginal phase medians need not sum
exactly to the median total. Fine local/core timers are nested inside facet
execution and must not be added as independent elapsed phases.
Raw files retain builder/kernel/total samples for all three methods, PLS phase
profiles, and coarse/detailed RK profiles for every configuration.

Main candidate startup observations; each row is the median of two seeds'
single observations, not a dedicated construction benchmark:

| Input | Loading (ms) | Native construction (ms) |
| --- | --- | --- |
| 2D terrain, n=64 | 3.564 | 10.128 |
| 3D volume, n=16 | 9.920 | 51.577 |
| 3D volume, n=32 | 91.317 | 511.457 |
| 4D grid, n=4 | 1.170 | 9.793 |
| 5D grid, n=3 | 2.577 | 45.229 |
| 6D grid, n=2 | 0.652 | 19.010 |
| 7D grid, n=2 | 4.069 | 294.899 |

## Verification and limitations

Complete ordered old/new topology and gradient dumps agree before timing.
Every timed/profiled call matches its method's validated sequential reference.
Critical counts and Euler characteristics are unchanged on all 14 inputs.
All three methods agree on 13 inputs. The existing 4D/seed0 exception remains:
F-Max `[19, 46, 43, 16, 1]`, PLS/RK `[19, 46, 42, 15, 1]`.
The 7D vectors are `[5, 10, 7, 1, 0, 0, 0, 0]` and
`[3, 2, 0, 0, 0, 0, 0, 0]` for seeds 0 and 2.

Normal C++ and ASan/UBSan suites pass. The cached oracle tests cover sparse
vertex IDs, non-numerical bucket order, shared facets, plateau/tied/injective
filtrations, large-to-packed/graph scratch reuse, and one/two/four/eight workers
through dimension 9. The native Python suite has 164 passes and four optional
resident-benchmark skips; fallback runs 168 tests with 14 skips.
Eleven native benchmark-protocol tests and historical-header compatibility pass.
TTK was not supplied or rerun.

The audit checks phase nesting, hashes, identities, Euler characteristics,
balanced orders, summaries/intervals, and reversed schedules. It requires
unchanged sparse cells/entries/duplicate hits and unchanged local reduction
work. Indexed traversal visits must equal newly discovered faces plus
duplicate same-level hits. All diagnostics compile out of ordinary timed calls.

Validation assessment: **share with caveats**. The sequential 7D gain is
supported; parallel gain, low-dimensional behavior, memory impact, and
generalization beyond this small synthetic corpus remain limitations.

## Reproduction and raw evidence

From the repository root, using fresh output paths:

```sh
OMP_WAIT_POLICY=PASSIVE LC_ALL=C python3 tools/benchmark_simplicial_gradients.py \
  --baseline a066542 --candidate cb04b26 \
  --inputs ../rk-ab-inputs/terrain-n64-seed0.txt ../rk-ab-inputs/terrain-n64-seed2.txt \
    ../rk-ab-inputs/volume-n16-seed0.txt ../rk-ab-inputs/volume-n16-seed2.txt \
    ../rk-ab-inputs/volume-n32-seed0.txt ../rk-ab-inputs/volume-n32-seed2.txt \
  --grids 4:4 5:3 6:2 7:2 --seeds 0 2 --input-dir ../higher-dimensional-inputs \
  --workers 1 8 --blocks 6 --repeats 2 --warmups 2 --profiles 3 --rk-profiles 3 \
  --memory-repeats 0 --output ../rk-boundary-index-ab-main.json
```

Repeat with `--workers 8 1 --reverse-inputs` and a fresh confirmation path.

```sh
LC_ALL=C python3 tools/validate_simplicial_gradient_ab.py \
  ../rk-boundary-index-ab-main.json ../rk-boundary-index-ab-confirmation.json \
  --reversed-confirmation
```

Use the frozen study code when auditing historical helper hashes.
The initial `rk-closure-breakdown-baseline.json` is a clean diagnostic study;
`rk-boundary-index-profile-candidate.json` is the earlier exploratory candidate
profile with recorded dirty source status, not the headline evidence.
The raw files remain local; tools and this report are tracked.
All older study artifacts remain unchanged.

SHA-256 identities:

- `rk-boundary-index-ab-main.json`: `5fa4dba004ac366c282b18cadb7e6fea6b454d265920dba3737c9f3a7e728d2f`.
- `rk-boundary-index-ab-confirmation.json`: `51a330fe625143b6606709759f445575bbefb47076a31984ca62bf8f181101ed`.
- `rk-closure-breakdown-baseline.json`: `4064dedac42a757b4ca8d3ee499614a628257ff5c9f3224d5c0159b33822620c`.
- `rk-boundary-index-profile-candidate.json`: `38d4a7f9aa7e703e8320a3339690fc6359b4febe163ffb2f54feffc26171c3db`.
- Measured candidate headers: `12cb8cedcd884891a7af343c25b63ce3f9d5158ac563e4c081fd7eeecbac69aa`.

