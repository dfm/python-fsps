#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__version__ = "0.1.1"

import os
import re
import subprocess
import warnings

def run_command(cmd):
    """
    Open a child process, and return its exit status and stdout.

    """
    child = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE,
                             stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    out = [s for s in child.stdout]
    err = [s for s in child.stderr]
    w = child.wait()
    return os.WEXITSTATUS(w), out, err


# Check to make sure that the required environment variable is present.
try:
    ev = os.environ["SPS_HOME"]
except KeyError:
    raise ImportError("You need to have the SPS_HOME environment variable")

# Check the githash, and if there is none check the SVN version
cmd = 'cd {0}; git log --format="format:%h"'.format(ev)
stat, out, err = run_command(cmd)
accepted = (len(out) > 0) and (len(err) == 0)
if not accepted:
    warnings.warn("Your FSPS version is not under git version "
                  "control. FSPS is now available on github at "
                  "https://github.com/cconroy20/fsps")
    
    # Check the SVN revision number.
    ACCEPTED_FSPS_REVISIONS = [158, 160, 166, 168, 170]
    cmd = ["svnversion", ev]
    stat, out, err = run_command(" ".join(cmd))
    fsps_vers = int(re.match("^([0-9])+", out[0]).group(0))

    # Make sure you don't have some weird mixed version.
    accepted = ((fsps_vers in ACCEPTED_FSPS_REVISIONS) and
                (len(out[0].split(':')) == 1) and
                stat == 0)

    if not accepted:
        raise ImportError("Your FSPS revision, {0}, is not known to work with "
                        "this version of python-fsps. You can checkout an "
                        "accepted FSPS revision with "
                        "'svn update -r rev_number'. "
                        "The accepted FSPS rev_numbers are: "
                        "{1}".format(out[0].rstrip('\n'),
                                    ACCEPTED_FSPS_REVISIONS))
else:
    # Store the githash.  If any version checking is going to happen,
    # it should happen here
    fsps_vers = out[0]
    
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
