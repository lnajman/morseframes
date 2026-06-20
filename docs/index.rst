MorseFrames
===========

MorseFrames is an experimental C++ and Python library for Morse sequences,
reference and coreference maps, Morse complexes, and Morse-based persistent
homology.

The documentation is split between user-facing Python material, implementation
notes for the C++ core and GUDHI adapter, and reproducibility instructions for
the benchmark tables.

The current public alpha can be installed from PyPI:

.. code-block:: sh

   python3 -m pip install "morseframes==0.1.0a1"

.. toctree::
   :maxdepth: 2
   :caption: User guide

   project_overview
   quickstart
   minimal_examples
   inspecting_morse_data
   python_api
   python_prime_field_tutorial
   strategies
   plotting_demo

.. toctree::
   :maxdepth: 2
   :caption: C++ and integration

   cpp_complex_view_api
   gudhi_contribution_design_note
   gudhi_upstream_patch_map

.. toctree::
   :maxdepth: 2
   :caption: Experiments and reproducibility

   benchmark_summary
   benchmark_reproduction
   phase0_conventions

.. toctree::
   :maxdepth: 2
   :caption: API reference

   api
