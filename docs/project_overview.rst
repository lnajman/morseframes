Overview
========

MorseFrames is a small research library for computing persistent homology from
Morse data. Its input is a filtered simplicial complex together with, or used
to construct, a Morse sequence compatible with the filtration. Its output is a
set of critical simplexes, reference and coreference data, Morse complexes, and
persistence diagrams.

The guiding idea is that the persistence computation should work on the Morse
frame directly. Instead of first materializing a reduced filtered Morse complex
and then running an ordinary persistence implementation on that new complex,
MorseFrames computes reference and coreference annotations and reduces those
annotations. This keeps the main reduction indexed by critical simplexes and by
the working sets actually needed for their boundaries or coboundaries.

What Is a Morse Sequence?
-------------------------

Informally, a Morse sequence is an ordered way of simplifying a filtered
complex while remembering enough information to recover its homology. At each
step, the algorithm either declares a simplex to be critical, or pairs a
simplex with one of its cofaces. A pair represents a local cancellation: the
two simplexes are treated as regular and do not become generators of the Morse
complex. A critical simplex is one that survives this cancellation process and
acts as a generator in the reduced description.

The word "sequence" is important. The order of the decisions records how the
simplification respects the filtration, including levels where many simplexes
have the same value. Different strategies may choose different regular pairs
and therefore produce different numbers of critical simplexes, but a valid
sequence preserves the persistence information. MorseFrames keeps this sequence
as a first-class object because it is useful on its own, not only as an
intermediate step before persistence.

Core Objects
------------

``FilteredComplex``
   Stores the filtered simplicial complex. The Python interface exposes
   membership checks, simplex records, boundary/coboundary access, and
   conversion helpers. The C++ core stores compact integer simplex identifiers
   and gives the algorithms stable face/coface queries.

``MorseSequence``
   Records the sequence of critical declarations and regular pairs. The
   current strategies include saturated, F-Min, F-Max, same-level reduction,
   plateau-greedy, and flooding variants. The strategy interface is deliberately
   explicit so that new sequence algorithms can be added without changing the
   persistence interface. Saturated, same-level, plateau-greedy, and the named
   flooding variants are filtration-monotone flooding constructions; F-Min and
   F-Max are global seed-and-expand ``F``-sequence builders.

``Reference map``
   Expresses the boundary information needed by critical simplexes in terms of
   critical simplexes. It is the main object used by Morse-reference
   persistence.

``Coreference map``
   The dual construction, using coboundary information. It is used by the
   Morse-coreference path and is useful for inspecting cocritical structure.

``PersistenceDiagram``
   Stores finite and essential intervals. The same public API is available for
   ordinary persistence, Morse-reference persistence, and Morse-coreference
   persistence, over ``Z_2`` and prime fields ``F_p``.

Why Plateau Data Matters
------------------------

The sequence builders operate on filtered complexes directly, including
complexes with plateaus. They do not require converting the input to a
simplex-wise lower-star filtration before a Morse sequence can be computed.
This is important for image-like or terrain-like data where many cells or
simplexes naturally share the same value.

Current Scope
-------------

The public repository contains:

* a header-only C++ core for filtered simplicial complexes;
* a Python interface with a pure-Python fallback and optional native backend;
* persistence over ``Z_2`` and prime fields ``F_p``;
* a simplex-tree-like builder for convenient construction;
* an experimental GUDHI-facing adapter;
* reproducible benchmark scripts and generated benchmark table fragments.

The library is still research code. Names and interfaces may evolve while the
paper and the GUDHI integration mature, but the current structure is already
usable for experiments and examples.
