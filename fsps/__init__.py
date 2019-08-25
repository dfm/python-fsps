#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import os
import pygit2

__version__ = "0.3.1"

# Check to make sure that the required environment variable is present.
try:
    sps_home = os.environ["SPS_HOME"]
except KeyError:
    raise ImportError("You need to have the SPS_HOME environment variable")

# Check the git history to make sure the required FSPS updates are present
REQUIRED_GITHASH = 'eb4341d'

try:
    fsps_repo = pygit2.Repository(sps_home)
except pygit2.GitError:
    raise ImportError("Your SPS_HOME does not point to a repository under git version "
                      "control. FSPS is available on github at "
                      "https://github.com/cconroy20/fsps and should be cloned from there")

accepted = any(REQUIRED_GITHASH == commit.short_id
        for commit in fsps_repo.walk(fsps_repo.head.target))
if not accepted:
    raise ImportError("Your FSPS version does not have correct history. "
                      "Perhaps you need to pull the most up-to-date FSPS version. "
                      "Please make sure that you have the following commit "
                      "in your git history: {0}".format(REQUIRED_GITHASH))

# Only import the module if not run from the setup script.
try:
    __FSPS_SETUP__
except NameError:
    __FSPS_SETUP__ = False
if not __FSPS_SETUP__:
    __all__ = ["StellarPopulation", "find_filter", "get_filter",
               "list_filters"]
    from .fsps import StellarPopulation
    from .filters import find_filter, get_filter, list_filters
