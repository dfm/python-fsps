#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import os
import sys

extensions = ["sphinx.ext.autodoc", "sphinx.ext.intersphinx",
              "sphinx.ext.mathjax"]

source_suffix = ".rst"
master_doc = "index"
exclude_patterns = ["_build", "_themes"]
html_static_path = ["_static"]
html_theme_path = ["_themes"]
templates_path = ["_templates"]

# General information about the project.
project = "Python FSPS"
copyright = "2013, Dan Foreman-Mackey"
# version =
# release =

html_show_sphinx = True
html_show_sourcelink = False
html_use_smartypants = True
pygments_style = "sphinx"
# html_theme = "daftish"
# html_theme_options = {
#         "tagline": "Rapid exoplanet transit modeling in Python.",
#         "github": "https://github.com/dfm/bart",
#         "license_name": "MIT License",
#         "license_link": "https://raw.github.com/dfm/bart/master/LICENSE.rst",
#         "google_analytics": "UA-22909046-1",
#     }
# html_sidebars = {
#             "**": ["relations.html", "searchbox.html"]
#         }
# html_additional_pages = {}
