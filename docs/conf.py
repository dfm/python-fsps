#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from importlib.metadata import version as get_version
from pathlib import Path

if "SPS_HOME" not in os.environ:
    path = Path(__file__).absolute()
    sps_home = path.parent.parent / "src" / "fsps" / "libfsps"
    os.environ["SPS_HOME"] = str(sps_home)

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
copyright = "2013-2023 Python-FSPS developers"
version = get_version("fsps")

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
