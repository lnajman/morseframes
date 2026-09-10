# ProcessLowerStars workspaces and higher-dimensional gradients

Revision `203fbb9` reduces the time of the native simplicial ProcessLowerStars
(PLS) implementation without changing its ordered gradient. Across two
controlled studies, all **48 paired configuration medians improve**; 47 of the
48 within-session intervals are wholly below one. Paired time reductions range
from **19–47% sequentially** and **14–30% with eight workers**.

The new dimension-independent comparison also qualifies the RK performance
claim: **F-Max is fastest among the sequential methods in dimensions 4–7**.
Eight-worker RK beats both alternatives clearly in the tested 4D/5D cases, but
does not establish a consistent advantage over sequential F-Max in 6D/7D.
On 7D seed 2, its comparison with parallel PLS is inconclusive in both studies.
These are small synthetic experiments, not a universal ranking.

## Bounded implementation change

Only the native PLS implementation in `include/morseframes/morse_sequence.hpp`
changes. The public API, vertex ordering, lexicographic Robins priority keys,
pair/critical selection, load-balancing policy and ordered event replay remain
unchanged. In particular, this is the existing in-house simplicial PLS/F-Max
variant, not a newly imported external reference implementation.

- Replace each star's hash table of local simplex indices with one immutable,
  dense direct-index array prepared inside the algorithm.
- Reuse classification, boundary-count/XOR and heap buffers across the stars
  assigned to a worker. Every worker owns its mutable buffers; the owner/index
  tables are read-only during local processing.
- Store star membership by vertex rank, allocating one membership vector per
  vertex rather than per simplex. Avoid a redundant vertex hash lookup during
  owner/key construction.

The heaps use the same comparator and insertion/pop semantics as before.
Events are still recorded and replayed, including on the sequential path;
callback order and callback prefixes are preserved. There is no native-complex
redesign, shared-constructor change, RK change or F-Max change in this pass.
Persistence callers using PLS can benefit, but persistence is not benchmarked.

## Correctness and experimental boundary

An independent test oracle rescans every remaining local boundary at each step,
without the production heaps, boundary counters/XORs or direct-index workspace.
Tests cover dimensions 0–7, three shuffled vertex functions, sparse vertex IDs,
non-pure complexes, repeated use of a builder, callbacks and workers 1/2/4/8.
Existing tests still reject tied vertex values and non-max-vertex extensions.

For all 12 benchmark inputs, full ordered complex and sequential-gradient dumps
match byte for byte between versions, for all three methods. Each sequential
reference is checked as a valid Morse sequence; every timed and diagnostic run
checks its fingerprint against its own method's sequential reference. Euler
characteristics and critical counts by dimension are also checked. Different
methods are **not** required to have identical sequences or critical counts.

Local ordinary C++ and ASan/UBSan tests pass. The native Python suite passes
all 156 tests, including the new native benchmark protocol checks; the fallback
suite runs 156 tests with 11 skips. Both nine-case RK correctness corpora pass
at workers 1/2/4. The benchmark itself never computes persistence.

The primary timer starts with the same finalized MorseFrames complex and vertex
function already resident. It includes a **fresh builder plus the full gradient
algorithm**: function-specific partitioning, keys, workspaces, executor startup,
local work and replay are not precomputed outside this timer. Builders and the
returned sequence remain alive when timing stops. Temporary allocations cleaned
up within the algorithm are included; subsequent output/builder destruction,
validation and protocol handling are excluded.

The common input arrays are loaded and used to construct the native complex
before that timer. Loading and native construction are reported separately,
as agreed. This is a prepared-complex gradient comparison, not a claim about
total latency from arrays. The [earlier RK/F-Max/TTK comparison](native_order_reuse_benchmark.md)
remains a separate historical snapshot: **TTK was not rerun here**, and no new
PLS-versus-TTK speed claim is made.

## Study design

Baseline: `9f4184ed588071fec03fbb375663e1a2d9dd2394`.
Candidate: `203fbb9b8ad683365573af3132db0c388e8d39e5`.
Both studies start with a clean worktree. Frozen header snapshots are compiled
with the same new native driver and `-std=c++17 -O3 -DNDEBUG -pthread`, using
Apple Clang 15 on the native ARM Apple M1 Max (10 cores).

Each study uses terrain n=64 and volume n=16 controls, plus Freudenthal grids
4D n=4, 5D n=3, 6D n=2 and 7D n=2. Each family has seeds 0/2 and workers 1/8:
12 inputs and 24 configurations. The high-dimensional inputs use deterministically
shuffled injective integer vertex ranks and the max-vertex simplex extension.
The generator limits face-generation attempts before enumerating cells.
Dimensions, mesh sizes and star-size distributions vary together: this is not
an isolated dimension-scaling experiment or a representative application corpus.

There are eight balanced old/new execution-order blocks, three repetitions
per version/block, and two warmups per worker setting. All six permutations
of the three methods occur equally often across the 24 measured repetitions.
The repeat reverses input and worker order. Separate runs collect three PLS
diagnostics and three fresh-process peak-memory samples per method/version/
configuration. No local build or test runs overlap measurement. Other desktop
applications remain active, and no timing sample is discarded.

Ratios use the median of candidate/baseline block-median ratios. The 95%
paired-block bootstrap intervals describe each session, not independent-machine
uncertainty. Main PLS ratios range from 0.550–0.723 sequentially and 0.723–0.860
with eight workers; repeat ranges are 0.532–0.807 and 0.701–0.827. The sole
inconclusive optimization interval is main 4D/seed 2/eight workers:
0.730–1.015, with paired median 0.767. All repeat intervals are below one.
Unchanged F-Max/RK controls fluctuate (paired new/old medians about 0.88–1.12
across these studies); their variation is not credited to this PLS change.

## Current gradient comparison

Main-study times in milliseconds: each entry is the median of the two
seed-specific medians. F-Max is sequential; the eight-worker setting applies
only to PLS/RK. Loading and construction are excluded from this table, but all
algorithm-specific preparation described above is included.

| Input | Simplices | F-Max, seq. | PLS, seq. | RK, seq. | PLS, 8 workers | RK, 8 workers |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2D terrain, n=64 | 24,067 | 1.351 | 5.529 | 1.057 | 4.505 | 0.699 |
| 3D volume, n=16 | 91,891 | 9.663 | 23.844 | 4.133 | 12.700 | 1.691 |
| 4D grid, n=4 | 15,307 | 1.381 | 4.908 | 2.053 | 2.394 | 0.793 |
| 5D grid, n=3 | 48,965 | 5.482 | 18.077 | 12.774 | 6.753 | 3.418 |
| 6D grid, n=2 | 18,731 | 1.581 | 6.464 | 5.344 | 2.521 | 1.790 |
| 7D grid, n=2 | 189,171 | 27.672 | 87.315 | 106.347 | 29.556 | 28.025 |

These pooled summaries are descriptive, not significance tests. Comparisons
use the F-Max observations collected in the *same worker-setting blocks*,
rather than mixing the separate one-worker and eight-worker sessions.
Eight-worker RK/PLS intervals are below one in 11/12 inputs in **both** studies;
the exception is 7D seed 2 (main ratio 1.021, interval 0.981–1.097; repeat
1.043, interval 0.990–1.080). Eight-worker RK/F-Max intervals are below one
on all eight 2D–5D inputs in both studies, but not on any of the four 6D/7D
inputs. RK is slower than F-Max on 6D seed 2 in both studies. The repeat also
has a wholly above-one RK/F-Max interval on 7D seed 2.

PLS is substantially improved but still slower than F-Max sequentially on every
input. The matching remains a useful independent baseline; simplicity alone
does not establish performance parity with an optimized external implementation.
Sequential RK is slower than sequential PLS on both 7D seeds. Optimizing RK
for this newly exposed regime is a more evidence-based next task than assuming
RK must already be the fastest method.

## PLS phase attribution

The instrumented runs are separate from the performance runs above. The
following main-study 7D times are milliseconds, taking the median of each
seed's three diagnostics and then the median across seeds.

| Version / workers | Builder | Setup | Local wall | Replay | Residual |
| --- | ---: | ---: | ---: | ---: | ---: |
| Previous / 1 | 1.263 | 15.556 | 95.465 | 0.846 | 7.579 |
| Optimized / 1 | 1.243 | 13.379 | 62.878 | 0.852 | 7.282 |
| Previous / 8 | 1.231 | 5.048 | 21.306 | 0.890 | 7.867 |
| Optimized / 8 | 1.271 | 4.184 | 15.119 | 0.827 | 7.423 |

The measured improvement is principally in local processing, consistent with
the removed hash work and reused buffers. Replay is unchanged and small here.
Residual is measured total minus builder/setup/local/replay, calculated **per
diagnostic run before aggregation**. It includes uninstrumented bookkeeping and
temporary-workspace cleanup; it must not be silently omitted or attributed
entirely to replay. Marginal medians need not sum to the median total. All
methods also retain builder/kernel/total timings in the raw performance data.

Large stars also limit how evenly work can be divided: the largest 7D stars
contain 27,134 and 18,056 simplices (14.34% and 9.54% of their respective
complexes). The 6D maxima are 19.09% and 11.55%, versus about 0.07% in the
2D controls. These are input descriptors, not a diagnosis of RK's bottleneck.

## Loading, construction and memory

Loading/construction below are the medians of the two candidate seed-specific
startup observations, not repeated construction benchmarks. Peak RSS uses
the median of three fresh processes per seed, then the median across seeds,
for eight-worker PLS. RSS includes input arrays, runtime, native complex,
constructor high-water mark and gradient workspaces; it is not isolated live
algorithm heap, and subtracting the constructor peak would not measure that.

| Input | Load (ms) | Native construction (ms) | Previous peak (MiB) | Optimized peak (MiB) |
| --- | ---: | ---: | ---: | ---: |
| 2D terrain, n=64 | 3.823 | 10.382 | 10.92 | 10.60 |
| 3D volume, n=16 | 9.941 | 54.596 | 35.31 | 33.93 |
| 4D grid, n=4 | 1.278 | 10.023 | 7.84 | 7.41 |
| 5D grid, n=3 | 2.836 | 44.935 | 21.35 | 20.24 |
| 6D grid, n=2 | 0.890 | 19.209 | 10.38 | 9.87 |
| 7D grid, n=2 | 4.335 | 302.932 | 87.97 | 81.31 |

No construction speedup is claimed. Per-process constructor and gradient peaks
for all three algorithms, both versions and worker settings remain in the raw
data. These small mesh results do not establish a universal memory reduction.

## Critical simplices

Counts are listed in dimension order, starting at dimension zero. PLS and RK
agree on every tested input. F-Max also agrees except on 4D seed 0, where it
creates one additional dimension-2 and one additional dimension-3 critical
simplex: `[19, 46, 43, 16, 1]`, versus `[19, 46, 42, 15, 1]` for PLS/RK.
Both results pass validation and have the same Euler characteristic.

| Input | Seed 0: PLS / RK | Seed 2: PLS / RK |
| --- | --- | --- |
| 2D terrain, n=64 | `[126, 224, 99]` | `[176, 334, 159]` |
| 3D volume, n=16 | `[21, 38, 29, 11]` | `[26, 54, 36, 7]` |
| 4D grid, n=4 | `[19, 46, 42, 15, 1]` | `[15, 38, 34, 10, 0]` |
| 5D grid, n=3 | `[15, 47, 48, 19, 4, 0]` | `[12, 36, 49, 26, 2, 0]` |
| 6D grid, n=2 | `[3, 5, 3, 0, 0, 0, 0]` | `[3, 2, 0, 0, 0, 0, 0]` |
| 7D grid, n=2 | `[5, 10, 7, 1, 0, 0, 0, 0]` | `[3, 2, 0, 0, 0, 0, 0, 0]` |

This observed agreement between PLS/RK is not asserted as a theorem for arbitrary
inputs. The benchmark intentionally accepts valid differences between methods.

## Reproduction and evidence

Run from the repository root; choose unused output paths. The controls are the
existing deterministic terrain/volume inputs from the earlier RK studies.

```sh
python3 tools/benchmark_simplicial_gradients.py \
  --baseline 9f4184ed588071fec03fbb375663e1a2d9dd2394 \
  --candidate 203fbb9b8ad683365573af3132db0c388e8d39e5 \
  --inputs ../rk-ab-inputs/terrain-n64-seed0.txt \
    ../rk-ab-inputs/terrain-n64-seed2.txt \
    ../rk-ab-inputs/volume-n16-seed0.txt \
    ../rk-ab-inputs/volume-n16-seed2.txt \
  --grids 4:4 5:3 6:2 7:2 --input-dir ../higher-dimensional-inputs \
  --workers 1 8 --output ../simplicial-pls-main.json
```

Repeat with `--workers 8 1 --reverse-inputs` and
`--output ../simplicial-pls-confirmation.json`. A high-dimensional-only run can
omit `--inputs`; its grid arrays are generated deterministically by the script.
Optional native protocol tests use
`MORSEFRAMES_SIMPLICIAL_GRADIENT_BENCHMARK=/absolute/path/to/worker`.

The raw JSON files are local study artifacts, not committed to this repository.
They retain every sample, execution order, derived interval, phase diagnostic,
critical vector, memory observation, input/header/compiler provenance and
completion flag. Recorded SHA-256 hashes:

- Main: `30fa851876469a6f8249086badf604c41cd6c247c1338f2c4df6c9f572bda86a`.
- Repeat: `8f5c4a990b69800853649c4a93745fa201b6cdc2ff46a44a38cccde7280ffdfd`.
- Baseline headers: `f161815c05f2822aaacd821fbbc8fc2f0ad9753d2872f9b7694d69aca3bc30cf`.
- Candidate headers: `d2bed87fe582e769e8dc39cb3fc7d4c2bb92a0d40005a1d72b88c9afb725204f`.

A separate post-run audit recomputed all summaries/comparisons, checked
phase partitions and nonnegative diagnostic residuals, all 24 configurations,
balanced method/version orders, source/input hashes, exact-version identities,
critical-cell Euler characteristics and memory/step counts. It also verified
that the repeat uses the same inputs and header snapshots in reversed order.
This validation is why uncertain 6D/7D comparisons and the F-Max critical-count
difference are retained rather than collapsed into a single ranking.
