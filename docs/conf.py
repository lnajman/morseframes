from __future__ import annotations

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

# Read the Docs should not need to compile the optional C++ backend just to
# build API documentation for the Python interface.
os.environ.setdefault("MORSEFRAMES_DISABLE_CPP_BACKEND", "1")

project = "MorseFrames"
author = "Laurent Najman"
copyright = "2026, Laurent Najman"

try:
    import morseframes

    release = morseframes.__version__
except Exception:
    release = "0.1.0"

version = release

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.napoleon",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
master_doc = "index"

exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "*_table.tex",
    "plotly_morse_square_demo.html",
]

html_theme = "sphinx_rtd_theme"
html_title = "MorseFrames documentation"
html_static_path = []

autodoc_typehints = "description"
autosectionlabel_prefix_document = True
myst_heading_anchors = 3
