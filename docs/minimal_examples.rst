Minimal Examples
================

The examples below are intentionally small and can be pasted into a Python
session after installing the package, or from the repository root with
``PYTHONPATH=python``.

One Edge
--------

This example builds two vertices and one edge. The edge kills one connected
component, so the finite barcode contains one zero-dimensional interval.

.. code-block:: python

   import morseframes as mf

   complex_ = mf.FilteredComplex.from_simplices([
       ([0], 0.0),
       ([1], 0.0),
       ([0, 1], 1.0),
   ])

   frame = mf.compute_morse_sequence_and_reference_map(
       complex_,
       algorithm="f-max",
   )
   diagram = mf.compute_morse_persistence(
       complex_,
       frame.sequence,
       frame.references,
   )

   print(frame.sequence.critical_simplices_as_simplices(complex_))
   print(diagram.finite_barcode())
   print(diagram.essential_barcode())

Prime-Field Coefficients
------------------------

Prime fields are selected with the ``modulus`` argument. MorseFrames validates
that the modulus is prime.

.. code-block:: python

   import morseframes as mf

   triangle = mf.FilteredComplex.from_simplices([
       ([0], 0.0),
       ([1], 0.0),
       ([2], 0.0),
       ([0, 1], 1.0),
       ([0, 2], 1.0),
       ([1, 2], 1.0),
       ([0, 1, 2], 2.0),
   ])

   diagram = mf.compute_morse_persistence(
       triangle,
       algorithm="f-min",
       modulus=3,
   )

   print(diagram.finite_barcode())
   print(diagram.essential_barcode())

Strategy Comparison
-------------------

Different Morse sequence strategies can produce different numbers of critical
simplexes while preserving the persistence barcode.

.. code-block:: python

   import morseframes as mf

   plateau = mf.FilteredComplex.from_simplices([
       ([0], 1.0),
       ([1], 1.0),
       ([2], 1.0),
       ([3], 2.0),
       ([0, 3], 2.0),
       ([1, 2], 2.0),
       ([1, 3], 2.0),
       ([2, 3], 2.0),
       ([1, 2, 3], 2.0),
   ])

   for algorithm in [
       "saturated",
       "f-min",
       "f-max",
       "same-level-reduction",
       "plateau-greedy",
   ]:
       sequence = mf.compute_morse_sequence(plateau, algorithm=algorithm)
       diagram = mf.compute_morse_persistence(
           plateau,
           sequence,
           mf.compute_reference_map(plateau, sequence),
       )
       print(
           algorithm,
           len(sequence.critical_simplices),
           diagram.finite_barcode(),
           diagram.essential_barcode(),
       )

One-Pass Reference Construction
-------------------------------

When the reference map is needed immediately after the sequence, use the fused
entry point. It computes the sequence and reference map together.

.. code-block:: python

   frame = mf.compute_morse_sequence_and_reference_map(
       plateau,
       algorithm="f-max",
   )
   diagram = mf.compute_morse_persistence(
       plateau,
       frame.sequence,
       frame.references,
   )
