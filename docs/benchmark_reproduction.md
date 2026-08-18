# Benchmark Reproduction

This page explains how to regenerate the public benchmark artifacts in this
repository. It is meant for software reproducibility: manuscript text and
discussion notes live outside the public repository, in the private manuscript
workspace until a public preprint or published version exists.

Run commands from the repository root.

```sh
cd morseframes
```

The examples below write raw CSV, Markdown summaries, and diagnostic prose to
`../work/`. The public repository tracks the scripts and selected LaTeX table
fragments, but not the manuscript prose built from them.

Some rendering commands below intentionally write tracked files under `docs/`.
Use those commands when you want to refresh the public table fragments. If you
only want to test the workflow on a local machine, redirect the table outputs to
`../work/` or restore the tracked table fragments afterward.

## Output Policy

Tracked public artifacts:

- `docs/*_table.tex`: LaTeX table fragments used to report benchmark results.
- `tools/*.py`: benchmark, validation, and table-rendering scripts.
- `benchmarks/benchmark_gudhi_view.cpp`: native GUDHI-view benchmark.

Local or private artifacts:

- `../work/*.csv`, `../work/*.md`, `../work/*.json`: raw benchmark outputs and
  summaries.
- `docs/*_prose.tex`: generated prose fragments. These are ignored by Git and
  should be copied into the private notes repository only when needed.
- report PDFs and manuscript drafts: private-note material, not public package
  documentation.

## Quick Validation

These checks are the fastest way to confirm that the source tree is usable.

```sh
MORSEFRAMES_DISABLE_CPP_BACKEND=1 \
  python3 -m unittest discover -s python/tests -p "test_*.py"

PYTHONPATH=python python3 python/examples/quickstart.py
PYTHONPATH=python python3 python/examples/prime_field_tutorial.py --modulus 3
```

To include the native C++ backend, install the package in editable mode:

```sh
python3 -m pip install -e ".[dev]"
python3 -c "import morseframes as mf; print(mf.__version__, mf.cpp_backend_available())"
```

The C++ smoke tests are:

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

## Synthetic Morse vs Standard Benchmarks

The main synthetic runner is `tools/benchmark_persistence.py`. It can run one
strategy, or the default strategy portfolio with `--sequence-algorithm
portfolio`.

Small smoke run:

```sh
mkdir -p ../work
PYTHONPATH=python python3 tools/benchmark_persistence.py \
  --preset smoke \
  --sequence-algorithm portfolio \
  --format summary \
  --output ../work/benchmark_smoke_summary.txt
```

Regenerate the public synthetic scale table:

```sh
mkdir -p ../work
PYTHONPATH=python python3 -c "import morseframes as mf; print(mf.cpp_backend_available())"

PYTHONPATH=python python3 tools/benchmark_persistence.py \
  --families lower-star plateau rips \
  --sizes 48 \
  --seeds 0 1 2 \
  --repeats 3 \
  --sequence-algorithm portfolio \
  --validation-mode core \
  --format csv \
  --output ../work/synthetic_scale_size48_portfolio.csv

PYTHONPATH=python python3 tools/render_synthetic_scale_table.py \
  --input ../work/synthetic_scale_size48_portfolio.csv \
  --table-output docs/synthetic_scale_table.tex \
  --prose-output ../work/synthetic_scale_prose.tex
```

The table reports `Std/Morse`, so values above `1` mean the Morse pipeline is
faster than ordinary full-complex persistence for that row.

CSV and JSON rows also retain gradient quality, rather than reporting timing
alone. The fields include the total critical count,
`critical_simplices_by_dimension`, `critical_ratio`, `num_regular_pairs`, and
`sequence_seconds_per_eliminated_simplex`. Timing is split into sequence and
downstream persistence phases, with `total_seconds` recording their complete
end-to-end pipeline. This makes comparisons between strategies meaningful even
when a faster constructor produces a larger Morse complex.

The tracked synthetic table is a native-backed core-mode benchmark. Before
replacing it, make sure the backend check above prints `True`; otherwise the CSV
will contain `cpp_backend=False` rows and the timing will describe the
pure-Python fallback instead of the optimized C++ backend.

## ProcessLowerStars Scaling

The focused ProcessLowerStars runner constructs two controlled simplicial
families. The balanced family gives every anchor vertex the same triangle-fan
workload. The skewed family keeps the same number of simplices but concentrates
most of that work in one lower star. All vertex values are injective and every
simplex uses the exact max-vertex lower-star extension.

```sh
PYTHONPATH=python python3 tools/benchmark_process_lower_stars.py \
  --anchors 16 \
  --balanced-fans 8 32 128 \
  --light-fan 2 \
  --workers 1 2 4 8 \
  --repeats 5 \
  --warmups 1 \
  --format csv \
  --output ../work/process_lower_stars_scaling.csv
```

For each scale, the heavy fan is selected automatically so both families have
the same total fan size and therefore the same simplex count. For example,
`16 * 8 = 98 + 15 * 2`. This isolates load imbalance from input size. Each
parallel row is checked against the sequential algorithm
for the exact Morse-step sequence and against ordinary persistence for the
barcode. The output records estimated task loads, critical counts by dimension,
construction and downstream persistence times, speedup, and parallel
efficiency. A `cpp_backend=False` row is a correctness run of the sequential
fallback, not a parallel-performance measurement.

The first controlled run on an Apple M1 Max (10 CPU cores), using seven
measured repetitions after one warm-up, is shown below. At 12,304 simplices the
balanced case reaches 1.81x sequence speedup and 1.60x end-to-end speedup. The
matched skewed case reaches only 1.19x and 1.13x, respectively, while its
estimated eight-worker load ratio rises to 6.79. The two cases have identical
critical counts `(2048, 2032, 0)`, so this gap is attributable to scheduling
imbalance rather than a different Morse complex.

![ProcessLowerStars scaling](process_lower_stars_scaling.svg)

The exact paper-ready values are generated in
`docs/process_lower_stars_scaling_table.tex`. These measurements are a local
scheduler study, not yet the comparison with Robins' implementation; that
external benchmark remains a separate stage.

## Roadmap and External Data

The benchmark runner also has Roadmap and CAM-style families:

```text
cam-s4-rips
roadmap-rips
```

Roadmap datasets are cached under `../work/roadmap-data` by default. Missing
Roadmap files are not downloaded unless requested explicitly:

```sh
PYTHONPATH=python python3 tools/benchmark_persistence.py \
  --preset roadmap \
  --sequence-algorithm portfolio \
  --download-roadmap-data \
  --format csv \
  --output ../work/roadmap_portfolio.csv
```

Use this only when network access is acceptable.

## Native GUDHI-View Benchmark

The native GUDHI benchmark compares three in-process paths on the same
`Gudhi::Simplex_tree` input:

- `Direct`: MorseFrames through a read-only `Simplex_tree` view.
- `Import`: copy into the compact owning MorseFrames complex first.
- `GUDHI`: GUDHI persistent cohomology on the original `Simplex_tree`.

This benchmark is optional because it needs GUDHI and Boost headers. Configure
them explicitly when CMake cannot find them:

```sh
cmake -S . -B build-gudhi \
  -DMORSEFRAMES_GUDHI_INCLUDE_DIR=/path/to/gudhi/include \
  -DMORSEFRAMES_BOOST_INCLUDE_DIR=/path/to/boost/include

cmake --build build-gudhi --target morseframes_benchmark_gudhi_view
```

Quick run:

```sh
mkdir -p ../work
./build-gudhi/morseframes_benchmark_gudhi_view \
  --quick \
  --repeats 3 \
  > ../work/native_gudhi_view_quick.csv

PYTHONPATH=python python3 tools/render_native_gudhi_view_table.py \
  --input ../work/native_gudhi_view_quick.csv \
  --output docs/native_gudhi_view_quick_table.tex \
  --summary

PYTHONPATH=python python3 tools/render_native_gudhi_stage_profile.py \
  --input ../work/native_gudhi_view_quick.csv \
  --table-output docs/native_gudhi_stage_profile_quick_table.tex \
  --prose-output ../work/native_gudhi_stage_profile_quick_prose.tex \
  --summary
```

Default-size repeat run:

```sh
./build-gudhi/morseframes_benchmark_gudhi_view \
  --repeats 30 \
  > ../work/native_gudhi_view_default_r30.csv

PYTHONPATH=python python3 tools/render_native_gudhi_view_table.py \
  --input ../work/native_gudhi_view_default_r30.csv \
  --output docs/native_gudhi_view_default_r30_table.tex \
  --caption-title "Native \\texttt{Gudhi::Simplex\\_tree} default benchmark." \
  --label tab:native-gudhi-view-default-r30 \
  --summary
```

Larger lean run:

```sh
./build-gudhi/morseframes_benchmark_gudhi_view \
  --large \
  --lean \
  --repeats 30 \
  > ../work/native_gudhi_large_lean_r30.csv

PYTHONPATH=python python3 tools/render_native_gudhi_view_table.py \
  --input ../work/native_gudhi_large_lean_r30.csv \
  --output docs/native_gudhi_large_lean_r30_table.tex \
  --caption-title "Native \\texttt{Gudhi::Simplex\\_tree} larger lean benchmark." \
  --label tab:native-gudhi-large-lean-r30 \
  --summary
```

In these tables, `GUDHI/Direct < 1` means GUDHI is faster end-to-end, while
`GUDHI/Reducer > 1` means the Morse reducer kernel alone is faster than GUDHI
persistence after the Morse input has already been built.

## Prime-Field Overhead

Prime-field coefficient experiments are generated by
`tools/benchmark_prime_field_overhead.py`.

Quick local run:

```sh
mkdir -p ../work
PYTHONPATH=python python3 tools/benchmark_prime_field_overhead.py \
  --families lower-star plateau rips \
  --sizes 8 12 16 \
  --seeds 0 1 \
  --algorithms saturated f-max same-level-reduction \
  --primes 3 5 \
  --repeats 5 \
  --output-csv ../work/prime_field_overhead_quick.csv \
  --output-md ../work/prime_field_overhead_quick.md
```

Composite moduli are intentionally rejected by the barcode API; these reducers
work over fields `F_p`.

## Profile-Selection Validation

The profile-selection scripts compare cheap strategy-selection metrics against
measured portfolio timings. These runs are more expensive than the smoke tests.

Preview the commands without executing them:

```sh
PYTHONPATH=python python3 tools/run_fair_profile_validation.py \
  --validation-preset report \
  --dry-run
```

Regenerate the public validation table from fresh timings:

```sh
mkdir -p ../work
PYTHONPATH=python python3 tools/run_fair_profile_validation.py \
  --validation-preset report \
  --output-dir ../work \
  --table-output docs/profile_metric_fair_validation_table.tex \
  --prose-output ../work/profile_metric_fair_validation_prose.tex \
  --manifest-output ../work/fair_profile_validation_manifest.md
```

If CSVs already exist in `../work`, summaries can be regenerated without
rerunning timings:

```sh
PYTHONPATH=python python3 tools/run_fair_profile_validation.py \
  --validation-preset report \
  --output-dir ../work \
  --summaries-only \
  --table-output docs/profile_metric_fair_validation_table.tex \
  --prose-output ../work/profile_metric_fair_validation_prose.tex \
  --manifest-output ../work/fair_profile_validation_manifest.md
```

Selector decision and feature diagnostic tables are rendered from the validation
CSVs:

```sh
PYTHONPATH=python python3 tools/summarize_selector_decisions.py \
  --table-output ../work/profile_selector_decision_summary.txt \
  --csv-output ../work/profile_selector_decision_summary.csv \
  --latex-output docs/profile_selector_decision_summary_table.tex \
  --prose-output ../work/profile_selector_decision_summary_prose.tex

PYTHONPATH=python python3 tools/analyze_selector_features.py \
  --table-output ../work/selector_feature_diagnostic.txt \
  --csv-output ../work/selector_feature_diagnostic.csv \
  --json-output ../work/selector_feature_diagnostic.json \
  --latex-output docs/selector_feature_diagnostic_table.tex \
  --prose-output ../work/selector_feature_diagnostic_prose.tex
```

## Benchmark Summary Page

The visible compact tables in `docs/benchmark_summary.md` are generated from the
tracked LaTeX table fragments:

```sh
python3 tools/render_benchmark_summary.py
```

CI checks that this generated block is up to date:

```sh
python3 tools/render_benchmark_summary.py --check
```

## Before Committing Regenerated Results

Before committing regenerated table fragments, run:

```sh
git diff -- docs tools benchmarks
git diff --check
python3 tools/render_benchmark_summary.py --check
MORSEFRAMES_DISABLE_CPP_BACKEND=1 \
  python3 -m unittest discover -s python/tests -p "test_*.py"
```

Commit only public artifacts that are meant to be reproducible from this
repository. Keep manuscript text, discussion packages, generated prose, and PDFs
in the private notes repository.
