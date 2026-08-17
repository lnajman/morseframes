Reduction-kernel correctness corpus
===================================

The reduction-kernel correctness corpus is a small, fixed collection of
simplicial complexes used to test the sequential and parallel implementations.
It is intentionally separate from the performance benchmark inputs: its job is
to expose semantic or determinism regressions, not to produce timing results.

The manifest is stored in
``tests/data/reduction_kernel_corpus.json``. It includes plateaus, lower-star
filtrations, non-pure complexes, multiple connected components, and explicitly
filtered cycles and fillings.

Validation contract
-------------------

For every corpus case, ``tools/validate_reduction_kernel_corpus.py`` checks:

* the parallel implementation produces exactly the same ordered Morse steps as
  the sequential reduction-kernel implementation with 1, 2, and 4 workers;
* sequential and parallel Morse persistence agrees with ordinary
  boundary-matrix reduction after omitting zero-length intervals, which Morse
  cancellation may remove; and
* the corpus can be constructed and validated by both the native extension and
  the pure-Python fallback.

The native and fallback runs are separate CI steps. To reproduce the checks
locally after installing MorseFrames in editable mode, run:

.. code-block:: sh

   python tools/validate_reduction_kernel_corpus.py
   MORSEFRAMES_DISABLE_CPP_BACKEND=1 python tools/validate_reduction_kernel_corpus.py

Use ``--workers`` to test another fixed worker set, or ``--json`` to obtain a
machine-readable report. Each report includes a short digest of the ordered
step sequence; equality is checked on the full sequence, not on this digest.

Extending the corpus
--------------------

New cases should be named, minimal examples of a distinct correctness risk.
Inputs must be deterministic and independent of machine characteristics. A
``simplices`` case must list every face and use a monotone filtration;
``constant_facets`` and ``lower_star`` cases are closed under faces by the
builder. Random or large inputs belong in the benchmark datasets rather than
this manifest.
