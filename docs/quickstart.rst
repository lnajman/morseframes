Quick Start
===========

The Python package can be installed in editable mode from the repository root:

.. code-block:: sh

   python3 -m pip install -e .

For documentation-only or pure-Python experiments, the optional C++ backend can
be disabled:

.. code-block:: sh

   MORSEFRAMES_DISABLE_CPP_BACKEND=1 PYTHONPATH=python python3

The minimal workflow is to build a filtered complex, compute a Morse sequence,
construct a reference map, and reduce the resulting Morse frame:

.. code-block:: python

   import morseframes as mf

   complex_ = mf.FilteredComplex.from_simplices([
       ([0], 0.0),
       ([1], 0.0),
       ([2], 0.0),
       ([0, 1], 1.0),
       ([1, 2], 1.0),
       ([0, 2], 1.0),
       ([0, 1, 2], 2.0),
   ])

   sequence = mf.compute_morse_sequence(complex_, algorithm="f-max")
   references = mf.compute_reference_map(complex_, sequence)
   diagram = mf.compute_morse_persistence(complex_, sequence, references)

   print(diagram.finite_barcode())
   print(diagram.essential_barcode())

The sequence and reference map can also be computed in one pass:

.. code-block:: python

   frame = mf.compute_morse_sequence_and_reference_map(
       complex_,
       algorithm="f-max",
   )
   diagram = mf.compute_morse_persistence(
       complex_,
       frame.sequence,
       frame.reference_map,
   )

Prime-field coefficients are available by passing a prime modulus:

.. code-block:: python

   diagram = mf.compute_morse_persistence(
       complex_,
       algorithm="f-min",
       modulus=3,
   )

See :doc:`python_api`, :doc:`python_prime_field_tutorial`, and
:doc:`strategies` for the fuller interface.
