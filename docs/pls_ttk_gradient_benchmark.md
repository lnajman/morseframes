# Direct comparison of optimized PLS with TTK

**The optimized native ProcessLowerStars (PLS) is still substantially slower
than the pinned TTK classic backend on the tested 2D/3D meshes.** This is now a
direct, interleaved comparison, not a ratio assembled from separate historical
runs. TTK has the lower algorithm time in all **288 paired performance
repetitions**, across two studies and 24 configurations. All 24 within-session
PLS/TTK bootstrap intervals are wholly above one.

The PLS workspace optimization remains in place. This pass changes the
benchmark harness, not the production algorithms. PLS, RK and F-Max use the
same header content as the [higher-dimensional study](simplicial_gradient_benchmark.md).
The present results do not estimate another optimization speedup and do not
replace those 4D–7D experiments.

## Current comparison

Main-study algorithm times are milliseconds, taking the median of the two
seed-specific medians. F-Max remains sequential in every row; the worker
setting applies to PLS, RK and TTK. Loading and native construction are
separate. **Fresh builders and all algorithm-specific preparation are included.**

| Input | Workers | F-Max | RK | PLS | TTK | Paired PLS/TTK |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2D terrain, n=64 | 1 | 1.554 | 1.186 | 6.212 | 1.113 | 5.394 |
| 2D terrain, n=64 | 8 | 1.868 | 0.997 | 5.605 | 0.918 | 6.000 |
| 3D volume, n=16 | 1 | 11.073 | 4.850 | 26.645 | 11.863 | 2.292 |
| 3D volume, n=16 | 8 | 11.784 | 1.842 | 14.684 | 3.360 | 4.384 |
| 3D volume, n=32 | 1 | 149.458 | 55.420 | 297.518 | 109.344 | 2.712 |
| 3D volume, n=32 | 8 | 146.948 | 15.770 | 150.598 | 24.770 | 6.256 |

Each ratio is computed from the paired repetition ratios, then aggregated
across seeds; it is **not** obtained by dividing the two displayed time medians.
The confirmation PLS/TTK ratios, at one/eight workers respectively, are:

- Terrain n=64: **5.271 / 6.619**.
- Volume n=16: **2.252 / 4.129**.
- Volume n=32: **2.706 / 6.443**.

Individual configuration medians range from 2.275–6.434 in the main study and
2.176–7.153 in the repeat. Thus the improvement to our PLS did not bring it to
TTK's speed on these inputs. Its dimension-independent implementation remains
useful as the RK baseline beyond this TTK adapter's 1D–3D domain.

RK retains lower paired-median algorithm times than all three other methods on
every volume configuration in both studies. This does not extend to a universal
RK lead: the terrain aggregates favor TTK, and the separate higher-dimensional
study contains regimes where sequential F-Max is better.

## Timing contract

Every measured method starts from identical resident vertex values and maximal
cell arrays. Each call constructs fresh native objects; no function-specific
partition, priority keys, builder, triangulation or previous gradient is supplied
for free. The comparison then separates each **raw** total into:

| Method | Native construction, reported separately | Algorithm time, compared above |
| --- | --- | --- |
| F-Max / RK / PLS | Owning representation, max-vertex filtration extension, finalization and connectivity | Fresh builder plus complete gradient construction |
| TTK | Native object setup, cell representation and connectivity preconditioning | Vertex ordering plus `buildGradient`, including lower-star construction and matching |

For PLS, the algorithm interval includes owner/priority preparation, star
partitioning, dense indices, worker workspaces, executor startup, local processing
and ordered replay. TTK's lower stars are built inside its timed gradient call.
The classic backend's cache-bypass path is used on every call. TTK performs
ordering before connectivity preparation in its existing call order; its
algorithm time is therefore the sum of its ordering and gradient intervals.

The returned native gradient remains alive when the timer stops. MorseFrames
returns its full ordered Morse sequence; TTK returns its native cell pairings.
No conversion of TTK's output to an ordered Morse sequence is included. Post-result
object destruction, validation, input generation and persistence are excluded.
Temporary cleanup performed *inside* gradient construction remains included.
The process and library runtimes are warmed; this is not a cold-process benchmark.

Performance calls use only outer phase-boundary clocks, with no MorseFrames
internal profiling. Separate diagnostic calls collect nested details. Components
are summed per sample before computing distributions; subtracting independently
aggregated medians would give a different, invalid partition.

## PLS phase measurements

Main-study PLS diagnostic times below are milliseconds: median of four
diagnostics per seed, then median across seeds. They are separate observations
from the performance table and need not add to its headline medians.

| Input | Workers | Lower-star setup | Local processing | Replay | Unattributed |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2D terrain, n=64 | 1 | 2.437 | 2.080 | 0.091 | 1.680 |
| 2D terrain, n=64 | 8 | 1.992 | 1.043 | 0.097 | 1.762 |
| 3D volume, n=16 | 1 | 7.538 | 14.293 | 0.489 | 4.106 |
| 3D volume, n=16 | 8 | 4.390 | 4.210 | 0.520 | 4.673 |
| 3D volume, n=32 | 1 | 78.053 | 159.149 | 5.882 | 41.134 |
| 3D volume, n=32 | 8 | 51.815 | 42.298 | 6.766 | 46.953 |

These nested phases partition the gradient-call time, **not** the separate
builder phase. The largest-volume/eight-worker builder performance median is
6.947 ms. Unattributed time is the gradient duration minus measured nested
phases, calculated per diagnostic before aggregation. It includes uninstrumented
bookkeeping and temporary-workspace cleanup; no particular allocation or
scheduling cause has been isolated for all of it.

On the largest eight-worker volume, replay accounts for a median **4.34%** of
PLS diagnostic algorithm time in the main study and **4.02%** in the repeat.
Those fractions are calculated per diagnostic, including builder time in the
denominator. Replay alone cannot explain the PLS/TTK gap. Preparation, local
processing and the uninstrumented remainder all warrant attention if PLS is
optimized further; the table does not establish that eliminating any one stage
would achieve TTK parity.

TTK exposes no separate lower-star-construction and local-matching timers here.
For the largest eight-worker volume, its main performance medians are 2.398 ms
for vertex ordering and 21.935 ms for the combined gradient call. Do not
compare the PLS local-only diagnostic column to that combined TTK interval as
if they measured the same operation.

## Native construction and loading

The main construction medians for PLS/TTK are **12.998 / 0.586 ms** on the
eight-worker terrain, **64.355 / 5.169 ms** on volume n=16, and
**622.813 / 42.057 ms** on volume n=32. These are kept outside the primary
gradient comparison. F-Max, RK and PLS use identical construction code; any
differences between their measured construction times are not algorithmic
construction improvements.

The generated `resident_gradient_pls_ttk_966ecf6_construction_table.tex` reports
all four methods, both worker settings and shared loading. Loading is one
read/parse/validation per native invocation, summarized across seeds; it is
neither a per-algorithm cost nor a cold-cache disk benchmark. Full
resident-to-gradient totals and every outer phase remain in the raw JSON.
There are no new peak-memory measurements in this pass.

## Critical simplices and validation

All four methods have identical dimensionwise critical counts on all six
inputs, in both worker settings and both studies. Counts are ordered from
dimension zero upward:

| Input | Seed 0 | Seed 2 |
| --- | --- | --- |
| 2D terrain, n=64 | `[126, 224, 99]` | `[176, 334, 159]` |
| 3D volume, n=16 | `[21, 38, 29, 11]` | `[26, 54, 36, 7]` |
| 3D volume, n=32 | `[128, 455, 379, 51]` | `[152, 649, 570, 72]` |

Agreement is an observation, not an enforced equality between methods and not
a general theorem. The benchmark accepts valid differences and verifies the
reported agreement flag. It checks Euler characteristics independently.
The previously observed F-Max critical-count difference on a 4D input remains
in the higher-dimensional study.

Each MorseFrames reference is validated as a Morse sequence. Every timed and
diagnostic MorseFrames run compares all sequence fields with its own sequential
reference. Every TTK run compares its complete critical/upward/downward pairing
signature with its own one-worker reference. Different algorithms are not
required to return identical gradients.

All 162 native Python tests pass, including the rebuilt native adapter's
1D–3D PLS, phase-accounting, order-balance and tie-rejection checks. The fallback
suite runs 162 tests with 13 skips. Core C++ and ASan/UBSan tests, and the native
nine-case RK corpus at workers 1/2/4, also pass. Historical v1/v2 parsing is
retained; all three published `ba57c79` v2 tables reproduce byte for byte.

A separate post-run audit recomputes all summaries and comparisons from raw
samples, allowing only floating-point roundoff in derived sums. It checks every
configuration, balanced execution positions and adjacent-method pairs, phase
partitions, diagnostic residuals, critical counts, matching source/input/binary
hashes and the repeat's reversed input/worker order. No outliers are removed.

## Protocol and limitations

Both studies start clean at harness revision
`966ecf66c183504e10f42a428a467185b438414c`. They use the same native ARM binary
on the Apple M1 Max (10 cores), Apple Clang 15, C++17, `-O3 -DNDEBUG`, and OpenMP
with `OMP_WAIT_POLICY=PASSIVE`. TTK is an unmodified, clean checkout of
`f4ffd1a1049d0ccf6e8f3eb4f7c096a6cc251ba0`, using its explicit triangulation and
classic gradient backend. No local build or test overlaps timing.
The existing TTK release configuration also defines `TTK_ENABLE_KAMIKAZE`;
these are production-configuration implementation timings, not an experiment
equalizing every library's runtime safety checks.

Each study has six synthetic inputs (terrain n=64, volumes n=16/32, seeds 0/2),
workers 1/8, twelve performance repetitions, four separate diagnostics and two
warmups per mode. The four-round Williams schedule gives each method every
position once and each ordered adjacent-method pair once within a cycle. Three
complete cycles are measured. The repeat reverses input and worker order and
shifts the method-order cycle by two rows.

The 95% intervals resample paired repetition ratios within each configuration.
They describe these sessions, not independent-machine uncertainty or a universal
speed ratio; successive repetitions need not be statistically independent.
Other desktop applications remain active, and thread placement is not controlled
by the runner. In particular, the main terrain/eight-worker data contain large
outliers, including a paired PLS/TTK ratio of 61.48; the complete samples and
their distributions are retained. The direction of the PLS/TTK comparison is
consistent in every paired repetition, but the precise multiplier should not
be treated as portable across workloads or machines.

## Reproduction and artifacts

From the repository root, build the optional native benchmark against the pinned
TTK source using `tools/build_ttk_gradient_benchmark.sh` with target
`morseframes_resident_gradient_benchmark`, then run:

```sh
OMP_WAIT_POLICY=PASSIVE python3 tools/benchmark_resident_gradients.py \
  --benchmark ../work/ttk-benchmark/build-f4ffd1a1049d0ccf6e8f3eb4f7c096a6cc251ba0/morseframes_resident_gradient_benchmark \
  --include-pls --terrain-sizes 64 --volume-sizes 16 32 --seeds 0 2 \
  --workers 1 8 --repeats 12 --diagnostics 4 --warmups 2 \
  --input-dir ../rk-ab-inputs --output ../resident-pls-ttk-main.json
```

For the repeat, use `--workers 8 1 --order-offset 2 --reverse-inputs` and
`--output ../resident-pls-ttk-confirmation.json`. Existing output files are
never overwritten by the runner. The optional `--include-pls` mode uses the new
`resident-gradient-v3` / `resident-gradient-study-v3` formats; default
three-method invocations continue using v2 and the six-order schedule.

```sh
python3 tools/render_resident_gradients.py \
  --input ../resident-pls-ttk-main.json \
  --table-output docs/resident_gradient_pls_ttk_966ecf6_algorithm_table.tex \
  --phases-output docs/resident_gradient_pls_ttk_966ecf6_phases_table.tex \
  --construction-output docs/resident_gradient_pls_ttk_966ecf6_construction_table.tex
```

The repository tracks the reproducible driver, tests, renderer and selected
tables. Raw JSON remains in local study artifacts, with SHA-256 hashes:

- Main: `cc168b2c523bcedb9738804b4063c70922a595e8c84d12d150d69ff286c499b9`.
- Repeat: `bda563a076b14b85d45bb9e955b9a5359c08cc39f1822343cb7cd0d44a46ea26`.
- Native executable: `5968a419048de72a72734decfe62088ca5c4a5d37f481d4dd1609fa8055aa284`.
- Production headers: `d2bed87fe582e769e8dc39cb3fc7d4c2bb92a0d40005a1d72b88c9afb725204f`.

The production-header hash is unchanged from the PLS workspace study. This
direct comparison supersedes the informal cross-run PLS/TTK estimates, while
leaving all earlier raw studies intact.
