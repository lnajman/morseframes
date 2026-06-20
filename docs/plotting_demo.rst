Plotly Demo
===========

The repository includes a small interactive Plotly demonstration for a
triangulated square. It overlays Morse critical simplexes and Morse-complex
incidence lines on the surface, and includes controls for the persistence
threshold and surface opacity.

The demo can be regenerated from the repository root:

.. code-block:: sh

   python3 python/examples/plotly_morse_square_demo.py \
     --output docs/plotly_morse_square_demo.html

The generated HTML file is intentionally committed as a public, browser-ready
artifact:

:download:`Download the interactive demo <plotly_morse_square_demo.html>`
