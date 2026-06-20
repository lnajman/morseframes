Inspecting Morse Data
=====================

The Python interface exposes the intermediate objects used by the Morse-frame
pipeline. This is useful when debugging a sequence strategy, preparing a plot,
or checking how a reference or coreference map represents the reduced
boundary/coboundary information.

The example below uses a small plateau complex. The saturated strategy leaves
five critical simplexes and creates two regular pairs, so the printed objects
show both critical generators and cancellations.

Complete Example
----------------

.. code-block:: python

   from pprint import pprint
   import morseframes as mf


   def by_dimension(simplexes):
       groups = {}
       for simplex in simplexes:
           groups.setdefault(len(simplex) - 1, []).append(simplex)
       return {
           dimension: tuple(values)
           for dimension, values in sorted(groups.items())
       }


   def nonempty_support(complex_, supports):
       return {
           simplex: support
           for simplex, support in zip(complex_.simplex_list(), supports)
           if support
       }


   complex_ = mf.FilteredComplex.from_simplices([
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

   frame = mf.compute_morse_sequence_and_reference_map(
       complex_,
       algorithm="saturated",
   )
   sequence = frame.sequence
   references = frame.references
   coreferences = mf.compute_coreference_map(
       complex_,
       sequence,
       algorithm="saturated",
   )
   morse_complex = mf.compute_morse_complex(complex_, sequence, references)
   diagram = mf.compute_morse_persistence(complex_, sequence, references)

   print("simplices in filtration order")
   pprint(complex_.filtration_list())

   print("\ncritical simplexes by dimension")
   pprint(by_dimension(sequence.critical_simplices_as_simplices(complex_)))

   print("\nregular pairs")
   pprint(sequence.pairs_as_simplices(complex_))

   print("\nsequence steps")
   pprint([
       (step.type, step.sigma, step.tau, step.level)
       for step in sequence.steps_as_simplices(complex_)
   ])

   print("\nreference support")
   pprint(nonempty_support(
       complex_,
       mf.reference_map_as_simplices(complex_, sequence, references),
   ))

   print("\ncoreference support")
   pprint(nonempty_support(
       complex_,
       mf.coreference_map_as_simplices(complex_, sequence, coreferences),
   ))

   print("\nMorse boundary support")
   pprint(morse_complex.boundaries_as_simplices(complex_))

   print("\npersistence")
   print("finite", diagram.finite_barcode())
   print("essential", diagram.essential_barcode())

Representative Output
---------------------

The first block confirms the filtered complex. The values are not
simplex-wise distinct: several simplexes share level ``2.0``.

.. code-block:: text

   simplices in filtration order
   (((0,), 1.0),
    ((1,), 1.0),
    ((2,), 1.0),
    ((3,), 2.0),
    ((0, 3), 2.0),
    ((1, 2), 2.0),
    ((1, 3), 2.0),
    ((2, 3), 2.0),
    ((1, 2, 3), 2.0))

The sequence has three critical vertices, two critical edges, and two regular
pairs. Each regular pair is a cancellable face/coface pair.

.. code-block:: text

   critical simplexes by dimension
   {0: ((0,), (1,), (2,)), 1: ((1, 2), (1, 3))}

   regular pairs
   (((3,), (0, 3)), ((2, 3), (1, 2, 3)))

   sequence steps
   [('critical', (0,), None, 0),
    ('critical', (1,), None, 0),
    ('critical', (2,), None, 0),
    ('regular_pair', (3,), (0, 3), 1),
    ('critical', (1, 2), None, 1),
    ('critical', (1, 3), None, 1),
    ('regular_pair', (2, 3), (1, 2, 3), 1)]

The reference map writes each simplex boundary contribution in terms of
critical simplexes. The coreference map is the corresponding coboundary-side
view. Empty annotations are omitted below only to keep the display compact.

.. code-block:: text

   reference support
   {(0,): ((0,),),
    (1,): ((1,),),
    (1, 2): ((1, 2),),
    (1, 3): ((1, 3),),
    (2,): ((2,),),
    (2, 3): ((1, 2), (1, 3)),
    (3,): ((0,),)}

   coreference support
   {(0,): ((0,),),
    (0, 3): ((1, 3),),
    (1,): ((1,),),
    (1, 2): ((1, 2),),
    (1, 3): ((1, 3),),
    (2,): ((2,),)}

The Morse complex boundary is indexed by critical simplexes. Here the two
critical edges kill two of the three initially critical connected components.

.. code-block:: text

   Morse boundary support
   {(0,): (),
    (1,): (),
    (1, 2): ((1,), (2,)),
    (1, 3): ((0,), (1,)),
    (2,): ()}

   persistence
   finite ((0, 1.0, 2.0), (0, 1.0, 2.0))
   essential ((0, 1.0),)

Useful Accessors
----------------

The most common inspection helpers are:

* ``complex_.simplex_list()`` for the complex as a list of simplexes;
* ``complex_.filtration_list()`` for simplexes in filtration order;
* ``simplex in complex_`` or ``complex_.contains(simplex)`` for membership;
* ``complex_.boundary_simplices(simplex)`` and
  ``complex_.coboundary_simplices(simplex)`` for local incidence;
* ``sequence.steps_as_simplices(complex_)`` for the full Morse sequence;
* ``sequence.critical_simplices_as_simplices(complex_)`` for critical
  simplexes;
* ``sequence.pairs_as_simplices(complex_)`` for regular pairs;
* ``mf.reference_map_as_simplices(complex_, sequence, references)`` for the
  reference map;
* ``mf.coreference_map_as_simplices(complex_, sequence, coreferences)`` for
  the coreference map;
* ``mf.compute_morse_complex(...).boundaries_as_simplices(complex_)`` for the
  Morse complex boundary;
* ``diagram.finite_barcode()`` and ``diagram.essential_barcode()`` for
  persistence intervals.
