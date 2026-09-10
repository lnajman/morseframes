# RK: ordered lookup inside large local cells

Historical snapshot: the subsequent [eligible-candidate study](rk_eligible_candidates_benchmark.md)
tests the protected-scan target identified below, with fresh timing and caveats.

**Ordered coface lookup reduces sequential RK time by 7–13% on the tested 7D
inputs across two clean A/B studies.** It removes about 88% of the measured
membership comparisons without changing the gradient. Eight-worker results
are less consistent: both studies favor the optimization on seed 2, but neither
resolves a gain on seed 0. No consistent benefit is established in 5D.

This pass changes only RK's membership lookup inside local reduction. The
[on-demand closure optimization](rk_lazy_closure_benchmark.md), candidate scan
order, pair selection, scheduling, and replay are unchanged. PLS, F-Max,
construction, and persistence are not modified. TTK is not rerun; the
[PLS key-storage study](pls_key_arena_benchmark.md) remains the latest direct
TTK comparison, on its documented earlier RK revision.

## Diagnosis and one bounded change

New diagnostic counters separate membership comparisons in cells larger than
16 simplices, sparse candidate visits, scan passes, and visits skipped because
the candidate is removed or protected. Removed/protected categories are mutually
exclusive: removal is checked first. Packed-mask operations are not counted as
sparse searches. The counters compile out of ordinary performance calls.

The initial profile finds 23.4 million and 33.8 million linear-search
comparisons on the 7D seeds. Approximately 98–99% occur in cells larger than
16 simplices. The retained candidate uses binary search for these larger
workspace cells, comparing **level-bucket ranks, not numerical simplex IDs**.
It uses the existing cell and rank array, with no additional lookup table or
per-cell allocation. Small cells and externally cached cells retain linear
lookup; arbitrary external cache order is not assumed sorted. Packed kernels
retain their existing path. There is no new dimension bound.

The local cell is immutable during reduction. Local removals use the existing
separate flags, and unsuccessful searches still return the existing sentinel.
The first eligible candidate and unique-coface selection are unchanged, so
the optimization changes search cost, not matching decisions.

Exact diagnostic comparison counts, identical across repetitions and worker
counts within each version:

| Input | Membership searches | Comparisons before | Comparisons after |
| --- | ---: | ---: | ---: |
| 5D, seed 0 | 47,529 | 879,759 | 443,629 |
| 5D, seed 2 | 47,401 | 889,307 | 438,271 |
| 7D, seed 0 | 391,528 | 23,438,064 | 2,960,541 |
| 7D, seed 2 | 540,301 | 33,820,481 | 4,068,428 |

The 7D reductions are 87.4% and 88.0%. A comparison is not a CPU cycle:
binary search adds rank-array accesses and branches, and other work remains.
The smaller 5D comparison counts do not establish a time improvement.

Sequential local-reduction diagnostic medians, milliseconds:

| Study | 7D seed | Before | After |
| --- | ---: | ---: | ---: |
| Main | 0 | 29.411 | 24.396 |
| Main | 2 | 33.457 | 25.718 |
| Confirmation | 0 | 29.374 | 23.744 |
| Confirmation | 2 | 33.126 | 25.910 |

These separate instrumented calls support a reduction in local search cost.
They are not the headline timings. Detailed local/core timers are nested in
facet execution, and parallel diagnostic totals sum across tasks; they must
not be added as independent elapsed phases. No bandwidth or peak-memory result
is inferred from the comparison counts.

## Frozen before/after comparison

The baseline is `577ba2ccd595bf74268e379693941689f07c78b8`, which contains the
new diagnostics but retains linear membership lookup. The candidate is
`d364ab0fecd634349835e301222dd87ecf18b420`. Both use the same current benchmark
driver; new counters are optional so it can still compile against historical
header snapshots. Compatibility was checked against `ff32e44`.

Both studies ran on 10 September 2026 on the Apple M1 Max, native ARM64,
Apple Clang 15, C++17, `-O3 -DNDEBUG -pthread`, with clean worktrees. Each covers
14 inputs: terrain n=64, volumes n=16/32, and 4D/5D/6D/7D grids of side
4/3/2/2, each with seeds 0 and 2, at one and eight workers.

Each configuration uses two warmups and six paired blocks of two repetitions.
Version order is balanced and all six method orders occur equally. Confirmation
reverses input and worker-count order. There are three separate PLS profiles
and three coarse plus three detailed RK profiles per version/configuration.
No compilation or tests run during timing. Raw outliers are retained.

Timing starts with the same finalized native complex resident. It includes a
fresh builder and the full gradient call: algorithm-specific preparation,
thread startup, local work, replay, and in-method cleanup. Returned outputs and
builders remain alive at the stop boundary. Loading/native construction are
separate; no persistence computation or untimed algorithm-specific cache is
introduced. Raw samples retain builder, kernel, and total gradient times.

The ratios below are medians of paired block ratios, not quotients of marginal
time medians. The 95% paired-block bootstrap intervals describe within-session
variation on this machine, not independent-machine uncertainty.

| 7D seed | Workers | Main after/before [interval] | Confirmation after/before [interval] |
| --- | ---: | --- | --- |
| 0 | 1 | 0.911 [0.888, 0.954] | 0.926 [0.894, 0.976] |
| 2 | 1 | 0.870 [0.838, 0.924] | 0.892 [0.848, 0.903] |
| 0 | 8 | 0.921 [0.828, 1.014] | 0.897 [0.829, 1.023] |
| 2 | 8 | 0.863 [0.840, 0.966] | 0.962 [0.840, 0.995] |

All four sequential intervals favor the candidate, with paired reductions of
7.4–13.0%. Eight-worker seed-2 reductions vary substantially: 13.7% in the main
study and 3.8% in confirmation, with both intervals below one. Seed 0 has lower
paired medians in both studies, but both intervals overlap one. It must remain
an inconclusive parallel result, not a confirmed speedup.

### Controls and deviations

5D has no reproducible time benefit at either worker count. Main seed 0 with
eight workers is slower: ratio 1.132 [1.011, 1.198]. Confirmation gives
0.937 [0.860, 1.094]. All other 5D intervals include one. This slowdown is
retained, rather than discarded because it does not reproduce.

Other intervals favoring a slowdown are main volume32/seed2 with one worker,
1.041 [1.003, 1.052], and confirmation 6D/seed2 with one worker,
1.082 [1.042, 1.149]. Neither has a slowdown interval wholly above one in the
other study. Small 2D/4D timings also fluctuate, including positive median
changes on unchanged packed paths. There is no configuration with a slowdown
interval wholly above one in both studies; this is not proof of no regression
or a claim that the optimization helps every input.

## Current comparison with PLS and F-Max

Eight-worker configurations, milliseconds. F-Max remains sequential in these
same-run comparisons.

| Study | 7D seed | RK before | RK after | PLS | F-Max |
| --- | ---: | ---: | ---: | ---: | ---: |
| Main | 0 | 22.351 | 21.022 | 20.700 | 27.641 |
| Main | 2 | 20.669 | 18.593 | 16.577 | 25.002 |
| Confirmation | 0 | 22.370 | 20.019 | 22.001 | 27.666 |
| Confirmation | 2 | 20.031 | 19.015 | 16.891 | 24.731 |

RK/PLS seed-0 paired ratios are 1.009 [0.898, 1.045] and
0.956 [0.889, 1.022]: no resolved winner. Seed-2 ratios are
1.104 [1.035, 1.134] and 1.139 [1.100, 1.180]: PLS still leads. Eight-worker
RK beats sequential F-Max on both 7D seeds in both studies; sequential F-Max
remains substantially faster than sequential RK.

Do not multiply this study's ratios by the preceding closure study's ratios
to claim a measured cumulative gain: their timing sessions and baselines
differ. These synthetic grids are not a representative high-dimensional corpus.

## Construction, output quality, and verification

Main candidate loading/construction observations, milliseconds. Each row is
the median of two seed observations, not a dedicated construction benchmark.

| Input | Loading | Native construction |
| --- | ---: | ---: |
| 2D terrain, n=64 | 3.594 | 10.068 |
| 3D volume, n=16 | 10.020 | 50.480 |
| 3D volume, n=32 | 92.102 | 525.556 |
| 4D grid, n=4 | 1.229 | 9.830 |
| 5D grid, n=3 | 2.668 | 45.013 |
| 6D grid, n=2 | 0.682 | 19.072 |
| 7D grid, n=2 | 4.200 | 304.405 |

Complete ordered old/new complex and gradient dumps agree before timing;
every timed/profiled run matches its own validated sequential reference.
Critical counts and Euler characteristics are unchanged on all 14 inputs.
All three methods agree on 13 inputs. The existing 4D/seed0 exception is
F-Max `[19, 46, 43, 16, 1]` versus PLS/RK `[19, 46, 42, 15, 1]`.
The 7D count vectors remain `[5, 10, 7, 1, 0, 0, 0, 0]` and
`[3, 2, 0, 0, 0, 0, 0, 0]` for seeds 0 and 2. All vectors are retained in the
raw identities. Facet work, sparse scan passes, candidate visits, and membership
search counts are also unchanged; only comparison counts are expected to fall.

Normal C++ and ASan/UBSan suites pass. Tests compare with the retained
linear-search/cache oracle through dimension 9, using sparse vertex IDs,
non-numerical bucket order, shared facets, plateau/tied/injective filtrations,
graph controls, and one/two/four/eight workers. A large plateau checks fewer
membership comparisons without changed search or scan counts. The native
Python suite has 162 passes and four optional resident-benchmark skips (no TTK
adapter supplied); fallback runs 166 tests with 14 skips. Updated native
protocol tests and historical-header compatibility also pass.

The remaining local scans are a distinct next target: in 7D, 5,889,078 of
6,296,854 sparse candidate visits on seed 0 and 5,801,961 of 6,284,263 on seed 2
skip protected candidates. Those scans remain unchanged here. A future pass
can test preparing an ordered eligible-candidate list once per facet, counting
its preparation and retaining exact pairing order. Scheduling should remain a
separate experiment.

## Reproduction and retained evidence

From the repository root, with fresh output paths:

```sh
OMP_WAIT_POLICY=PASSIVE LC_ALL=C python3 tools/benchmark_simplicial_gradients.py \
  --baseline 577ba2c --candidate d364ab0 \
  --inputs ../rk-ab-inputs/terrain-n64-seed0.txt ../rk-ab-inputs/terrain-n64-seed2.txt \
    ../rk-ab-inputs/volume-n16-seed0.txt ../rk-ab-inputs/volume-n16-seed2.txt \
    ../rk-ab-inputs/volume-n32-seed0.txt ../rk-ab-inputs/volume-n32-seed2.txt \
  --grids 4:4 5:3 6:2 7:2 --seeds 0 2 --input-dir ../higher-dimensional-inputs \
  --workers 1 8 --blocks 6 --repeats 2 --warmups 2 --profiles 3 --rk-profiles 3 \
  --memory-repeats 0 --output ../rk-ordered-lookup-ab-main.json
```

Repeat with `--workers 8 1 --reverse-inputs` and a fresh confirmation output.
The automated audit recomputes summaries/intervals, checks phase accounting,
order balance, input/header/driver hashes, exact identities, unchanged work
counters, and the reversed schedule:

```sh
LC_ALL=C python3 tools/validate_simplicial_gradient_ab.py \
  ../rk-ordered-lookup-ab-main.json ../rk-ordered-lookup-ab-confirmation.json \
  --reversed-confirmation
```

Use the frozen study revision when auditing historical artifacts against their
recorded source hashes. The standalone RK profiler supports the same inputs
and one/eight-worker diagnostic passes without performance claims. Earlier raw
artifacts are unchanged; raw study files remain local, with tools and selected
results tracked. SHA-256 identities:

- Main: `dd13f8b482fd09391e4aad1c1faab8dfafa9495e6b1d6cd9f3e2bc89417e4d93`.
- Confirmation: `2af51d8cbccfb2da13733dce27e6c2f4663bae9571df491954be2c194f50fc2d`.
- Candidate headers: `9c9cea4e60150ddba2fc6fc6878f3ecffc5d81c13194ff92ec8113e817a12250`.
- Initial `rk-local-search-profile-baseline.json`: `0589b1afee107e0c7ed38fed612dd8d70b2534c78152161692f8bbe4f7c50f0b`.
- `rk-ordered-lookup-exploratory.json`: `c2b0c4dcea5b4372fd973f999e13402903654043082d446635c1d421eb9e7c18`.
