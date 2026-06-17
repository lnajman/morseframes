# Phase 0 Conventions

This prototype follows these conventions for the current MorseFrames reference
and coreference pipelines.

## Input

- The input is an explicit finite simplicial complex.
- Every simplex must be inserted explicitly.
- The complex must be closed under faces.
- The filtration must be monotone: if `eta < sigma`, then `F(eta) <= F(sigma)`.
- Filtration values are compressed into integer levels after `finalize()`.

Face insertion from maximal simplices is intentionally not implemented yet, because generated
faces need a clear filtration policy. A lower-star helper should be added separately.

## Ordering

- The saturated `F`-sequence processes filtration levels in increasing order.
- Inside one level, available regular pairs are preferred over critical fillings.
- Ties are deterministic and follow the complex's level bucket order:
  dimension first, then lexicographic vertex order.
- Critical labels are assigned in the order critical fillings appear in the `F`-sequence.
- The latest pivot is therefore the largest critical label in an annotation.

## Intervals

- A reference persistence pair `(birth, death)` represents `[F(birth), F(death))`.
- The reported interval dimension is `dim(birth)`.
- Essential intervals are represented as `[F(birth), infinity)`.
- Zero-length intervals are retained in raw output and filtered by `off_diagonal_pairs()`.

## Coefficients

- Morse sequence construction is coefficient-independent.
- The reference-map Morse persistence pipeline supports `Z2` and prime fields
  `F_p` in the C++ core and Python API with
  `compute_morse_persistence(..., modulus=p)`.
- The coreference-map Morse persistence pipeline also supports prime fields
  `F_p` in the C++ core and Python API with
  `compute_morse_coreference_persistence(..., modulus=p)`.
- Ordinary full-complex persistence also supports prime fields `F_p` in the
  C++ core and Python API with `compute_standard_persistence(..., modulus=p)`.
- An annotation is a sorted vector of critical labels.
- Annotation addition is symmetric difference.
- A prime-field annotation is a sorted vector of `(critical_id, coefficient)`
  pairs.
- Composite `Z_n` coefficients are intentionally rejected for the barcode API;
  arbitrary rings need a separate algebraic design.
- The C++ prime-field Morse reducers use compact working-set tables for the
  persistence-reduction phase and inverse-indexed pivot updates, matching the
  structure of the optimized `Z2` reducer.

## Current Scope

Implemented:

- explicit complex finalization;
- boundary and coboundary construction;
- saturated `F`-sequence construction;
- full reference map for all simplices;
- Morse-reference persistence reduction;
- full coreference map for all simplices;
- Morse-coreference persistence reduction;
- ordinary full-complex `Z2` persistence for validation;
- ordinary full-complex prime-field persistence in the C++ core and Python API;
- reference-side Morse prime-field persistence in the C++ core and Python API;
- coreference-side Morse prime-field persistence in the C++ core and Python API;
- debug invariant checks;
- structural, annotation, and timing instrumentation;
- lazy inverse lists for pivot updates;
- reducer storage restricted to `W_boundary_plus` / `W_coboundary_plus`;
- compact working-set C++ prime-field Morse reducers;
- inverse-indexed prime-field pivot updates;
- tiny validation tests, including tetrahedron and lower-star examples.

Not implemented yet:

- lower-star or maximal-simplex input helpers.
