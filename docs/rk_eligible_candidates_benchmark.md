# RK: ordered eligible-candidate scans

Historical snapshot: the subsequent [boundary-index study](rk_boundary_index_benchmark.md)
profiles and optimizes the remaining closure-preparation work with fresh evidence.

**The retained implementation reduces total sequential RK time by about 10.5%
on the tested 7D seed-2 input in both final studies.** It removes approximately
90% of repeated sparse candidate visits in 7D without changing coface queries
or the ordered gradient. Parallel medians improve, but a repeatable parallel
speedup is not established by the intervals. Seed-0 total-time gains are also
less conclusive.

Only RK's sparse local candidate scan changes. The previous
[ordered membership lookup](rk_ordered_lookup_benchmark.md), full cell,
on-demand closures, scheduling, pairing, and replay stay unchanged. PLS, F-Max,
native construction, and persistence are untouched. TTK is not rerun; the
[PLS key-storage study](pls_key_arena_benchmark.md) remains the latest direct
TTK comparison on its documented earlier RK revision.

## Implementation and correctness argument

Facet incidence is fixed during a round's isolated local reductions. A simplex
protected by more than one active facet cannot become eligible during that
local computation. For uncached cells larger than 16 simplices, build a list
of unprotected **original cell indices**, in the same order, once per facet.
Rescan this shorter list after each pairing, checking the existing removal flags.

The **full cell remains present for coface membership**: removing protected
simplices from that cell could falsely make a coface unique. Binary searches
still use canonical level-bucket ranks, not numerical simplex IDs. No query
or pairing decision is changed. The list is rebuilt for each facet computation,
so protection changes in later rounds are respected.

Small cells, externally cached cells, and packed kernels retain their prior
candidate traversal. External cache order need not be sorted. List storage is
constructed lazily, only for filtered cells; up to 16 indices use inline storage,
with vector overflow beyond that. There is no new dimension bound, shared
mutable candidate list, persistent cache, or scheduling change.

Diagnostic counts, identical across final repetitions and worker counts:

| Input | Sparse visits before | Sparse visits after | Protected visits remaining |
| --- | --- | --- | --- |
| 5D, seed 0 | 714,958 | 415,388 | 329,286 |
| 5D, seed 2 | 710,400 | 422,738 | 334,945 |
| 7D, seed 0 | 6,296,854 | 629,406 | 221,630 |
| 7D, seed 2 | 6,284,263 | 626,251 | 143,949 |

These are **repeated candidate-loop visits**, not all operations. The new
one-time filter scan is not included in that counter, but its time is included
in gradient timing. Protected visits remaining come from unchanged small-cell
scans. Removed-candidate visits, scan passes, coface visits, membership searches,
and membership comparisons are unchanged. A 90% visit reduction is not a
90% time reduction.

## Two final frozen studies

Baseline: `bfcf7d94bb1199ab20e87aefbd54bdec65616c60`.
Retained candidate: `3603b9f7e09d7b785dfdf0adb572b8d4bb90834b`.

Both studies ran on 10 September 2026 on the Apple M1 Max, native ARM64,
Apple Clang 15, C++17, `-O3 -DNDEBUG -pthread`. Main ran 15:42:28–15:44:10 UTC;
confirmation ran 15:44:18–15:46:01 UTC. Each used a clean worktree and the same
14 inputs: terrain n=64, volumes n=16/32, and 4D/5D/6D/7D grids of side
4/3/2/2, seeds 0 and 2, with one and eight workers.

Each configuration has two warmups and six paired blocks of two repetitions.
Version order is balanced, all six method orders occur equally, and confirmation
reverses input and worker order. Three separate PLS profiles and three coarse
plus three detailed RK profiles are retained per version/configuration.
No builds or tests ran during measurement. This is an interactive workstation,
not an isolated benchmarking host; raw outliers are retained.

Ratios are medians of paired block ratios, not ratios of marginal time medians.
The paired-block bootstrap 95% intervals describe within-session variability,
not independent-machine uncertainty.

| 7D seed | Workers | Main after/before [95% interval] | Confirmation after/before [95% interval] |
| --- | --- | --- | --- |
| 0 | 1 | 0.935 [0.886, 1.078] | 0.979 [0.960, 0.997] |
| 2 | 1 | 0.895 [0.879, 0.905] | 0.894 [0.855, 0.905] |
| 0 | 8 | 0.971 [0.868, 1.018] | 0.949 [0.916, 1.008] |
| 2 | 8 | 0.877 [0.841, 0.908] | 0.902 [0.876, 1.084] |

Sequential seed 2 improves by 10.5% and 10.6%, with both intervals below one.
Seed 0 has lower medians but only confirmation resolves a sequential gain.
At eight workers, seed 2 has 12.3% and 9.8% lower paired medians, but the
confirmation interval includes one. Both seed-0 parallel intervals include one.
Do not present these as a confirmed general parallel acceleration.

### Controls and the initial candidate

All final 2D/3D/4D/5D intervals include one. In 6D, sequential seed 2 favors the
candidate in both studies, but its reductions vary from 12.2% to 1.7%;
the remaining 6D cases are inconclusive. No final interval is wholly above one.
This is not proof of no regression: sequential volume32/seed2 has ratios
1.051 [0.871, 1.087] and 1.054 [0.972, 1.133], and remains a performance risk.

The first candidate, `febd89b3e25c0948e1f7460b441110b2f50bfc3a`, constructed
inline candidate storage even when filtering was disabled. Two complete studies
of that version showed a repeated sequential volume32/seed0 slowdown:
1.038 [1.010, 1.054] and 1.051 [1.030, 1.121].
That prompted the lazy-storage revision, not removal of those observations.

The final version gives 1.020 [0.992, 1.095] and 0.998 [0.987, 1.011] on that
control. This does **not** prove storage initialization caused the first
slowdown: compiler layout and session variability can also affect timings.
The initial candidate's stronger 7D timings are not substituted for the
retained version's less conclusive results. All four raw files remain intact.

## Current comparison with PLS and F-Max

Final eight-worker configurations below; F-Max is still sequential.

| Study | 7D seed | RK before (ms) | RK after (ms) | PLS (ms) | F-Max (ms) |
| --- | --- | --- | --- | --- | --- |
| Main | 0 | 20.739 | 20.186 | 22.419 | 31.267 |
| Main | 2 | 20.263 | 17.679 | 18.053 | 25.199 |
| Confirmation | 0 | 20.432 | 19.642 | 21.092 | 29.429 |
| Confirmation | 2 | 18.610 | 17.092 | 16.589 | 24.692 |

For seed 0, paired RK/PLS ratios are 0.883 [0.727, 0.978] and
0.902 [0.893, 0.969], favoring RK in both sessions. For seed 2 they are
0.978 [0.934, 1.055] and 1.004 [0.952, 1.154]: no resolved winner.
Eight-worker RK beats sequential F-Max on both 7D seeds in both studies.
With one worker, F-Max remains roughly 2.5 times faster than RK in 7D.
These two synthetic grids do not establish which algorithm is fastest generally,
nor should ratios from different optimization studies be multiplied.

## Timing boundary and phases

Start from the same finalized common native complex in memory. Include a fresh
builder, all algorithm-specific preparation, worker startup, local work, replay,
and natural in-method cleanup. Keep the returned gradient and builder alive at
the stop boundary. Loading and native construction are reported separately;
there is no persistence computation or untimed algorithm-specific cache.

List construction and any overflow allocation are inside the local-reduction
timer. List destruction occurs after that inner timer, but remains inside
facet execution and the total gradient measurement. Detailed counters compile
out of ordinary performance calls.

Separate sequential local-reduction diagnostic medians:

| Study | 7D seed | Local reduction before (ms) | Local reduction after (ms) |
| --- | --- | --- | --- |
| Main | 0 | 23.260 | 17.218 |
| Main | 2 | 25.851 | 17.142 |
| Confirmation | 0 | 23.955 | 15.996 |
| Confirmation | 2 | 27.146 | 16.670 |

Coarse RK diagnostic medians for the clearest total-time result:

| RK phase, 7D seed 2, one worker | Main before (ms) | Main after (ms) | Confirmation before (ms) | Confirmation after (ms) |
| --- | --- | --- | --- | --- |
| Builder | 0.149 | 0.181 | 0.159 | 0.176 |
| Workspace/pool setup | 0.164 | 0.182 | 0.195 | 0.165 |
| Level processing | 70.164 | 62.052 | 73.542 | 61.560 |
| Replay | 0.323 | 0.322 | 0.364 | 0.325 |
| Unattributed remainder | 0.145 | 0.148 | 0.143 | 0.139 |
| Total profiled gradient | 70.909 | 62.870 | 74.382 | 62.436 |

These are separate instrumented calls, not the unprofiled headline samples.
Per-row accounting is validated; marginal phase medians need not sum exactly
to the median total. Fine local/core timers are nested within facet/level
execution, and parallel totals are cumulative across tasks, not additional
wall-clock phases. Raw files retain builder/kernel/total samples for every
method and all PLS/RK phase profiles. No peak-memory conclusion is drawn.

Main candidate startup observations, each row the median of two seeds'
single observations, not a dedicated construction study:

| Input | Loading (ms) | Native construction (ms) |
| --- | --- | --- |
| 2D terrain, n=64 | 3.463 | 10.153 |
| 3D volume, n=16 | 13.109 | 51.649 |
| 3D volume, n=32 | 92.321 | 518.142 |
| 4D grid, n=4 | 1.204 | 9.892 |
| 5D grid, n=3 | 2.521 | 44.864 |
| 6D grid, n=2 | 0.709 | 19.526 |
| 7D grid, n=2 | 4.162 | 294.542 |

## Verification and evidence quality

Complete ordered old/new topology and gradient dumps agree before timing.
Every timed/profiled call matches its method's validated sequential reference.
Critical counts and Euler characteristics are unchanged on all 14 inputs.
All three methods agree on 13 inputs. The existing 4D/seed0 difference remains:
F-Max `[19, 46, 43, 16, 1]`, PLS/RK `[19, 46, 42, 15, 1]`.
The 7D vectors remain `[5, 10, 7, 1, 0, 0, 0, 0]` and
`[3, 2, 0, 0, 0, 0, 0, 0]` for seeds 0 and 2.

Normal C++ and ASan/UBSan suites pass for the retained code. The independent
cached full-scan oracle covers shared facets, sparse vertex IDs, non-numerical
bucket order, plateau/tied/injective filtrations, packed/graph controls,
and one/two/four/eight workers through dimension 9. Tests explicitly require
the visit decrease to equal the omitted protected visits, with unchanged
removed visits and coface queries. Native Python: 163 passed, four optional
resident-benchmark skips; fallback: 167 tests, 14 skips. No TTK adapter was
supplied in this pass.

The raw-evidence audit checks hashes, complete identities, phase accounting,
balanced orders, summaries/intervals and reversed schedules. Its opt-in
`--allow-protected-scan-elision` permits only the precisely accounted protected
visit reduction; all other local work remains invariant. Tests reject changes
to those invariants. The default stricter audit remains available for prior
studies. Validation assessment: **share with caveats**, especially the uncertain
parallel gains, positive volume-control medians, and small synthetic corpus.

## Reproduction

From the repository root, with fresh output paths:

```sh
OMP_WAIT_POLICY=PASSIVE LC_ALL=C python3 tools/benchmark_simplicial_gradients.py \
  --baseline bfcf7d9 --candidate 3603b9f \
  --inputs ../rk-ab-inputs/terrain-n64-seed0.txt ../rk-ab-inputs/terrain-n64-seed2.txt \
    ../rk-ab-inputs/volume-n16-seed0.txt ../rk-ab-inputs/volume-n16-seed2.txt \
    ../rk-ab-inputs/volume-n32-seed0.txt ../rk-ab-inputs/volume-n32-seed2.txt \
  --grids 4:4 5:3 6:2 7:2 --seeds 0 2 --input-dir ../higher-dimensional-inputs \
  --workers 1 8 --blocks 6 --repeats 2 --warmups 2 --profiles 3 --rk-profiles 3 \
  --memory-repeats 0 --output ../rk-eligible-lazy-ab-main.json
```

Repeat with `--workers 8 1 --reverse-inputs` and a fresh confirmation output.

```sh
LC_ALL=C python3 tools/validate_simplicial_gradient_ab.py \
  ../rk-eligible-lazy-ab-main.json ../rk-eligible-lazy-ab-confirmation.json \
  --reversed-confirmation --allow-protected-scan-elision
```

Use the frozen study code when auditing against recorded helper hashes.
Both initial studies use the same protocol with candidate `febd89b` and
`rk-eligible-candidates-ab-main.json` / `rk-eligible-candidates-ab-confirmation.json`.
Raw artifacts remain local; reproduction tools and this report are tracked.

SHA-256 identities:

- Final main: `7c62ea9d4fe8a78057d9116b2f29b0b60b21c77b2e876a5be37cc9dde2ee4797`.
- Final confirmation: `ba930d9d82fda87d2a96f669d507b06b8272c427a2dffc8eaa657bd95cda7358`.
- Final headers: `bf0a02634a602177ee8b4daefe5edc502dabf56b11e9656f61deb5ebbb5e3a64`.
- Initial main: `f61466e07ec4ca768fb2cffad004a8d5d3c1b56b1b72bb28e18a08d75a4dc52f`.
- Initial confirmation: `1a2436200d69fb842b21aee9804ad2a81a49f49bb8d4be5bd91b4300ad5286d9`.
