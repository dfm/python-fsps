#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (
    division,
    print_function,
    absolute_import,
    unicode_literals,
)

import os
import sys

# Patch the path to get the local version.
d = os.path.dirname
sys.path.insert(0, d(d(os.path.abspath(__file__))))

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
]

source_suffix = ".rst"
master_doc = "index"
exclude_patterns = ["_build", "_themes"]
html_static_path = ["_static"]
html_theme_path = ["_themes"]
templates_path = ["_templates"]

# General information about the project.
project = "Python FSPS"
copyright = "2013-2021 Python-FSPS developers"
# version =
# release =

html_show_sphinx = True
html_show_sourcelink = False
html_use_smartypants = True
pygments_style = "sphinx"
sys.path.append(os.path.abspath("_themes"))
html_theme_path = ["_themes"]
html_theme = "dfm"
# Custom sidebar templates, maps document names to template names.
html_sidebars = {
    "index": ["sidebarintro.html", "searchbox.html"],
    "**": [
        "sidebarlogo.html",
        "localtoc.html",
        "relations.html",
        "searchbox.html",
    ],
}
