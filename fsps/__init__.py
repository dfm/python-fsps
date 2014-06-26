#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = []

import os
import subprocess
import re

try:
    ev = os.environ["SPS_HOME"]
except KeyError:
    raise ImportError("You need to have the SPS_HOME environment variable")

ACCEPTED_FSPS_REVISIONS = [140, 143]

def run_command(cmd):
    """Open a child process, and return its exit status and stdout"""
    child = subprocess.Popen(cmd, shell =True, stderr = subprocess.PIPE, 
                             stdin=subprocess.PIPE, stdout = subprocess.PIPE)
    out = [s for s in child.stdout]
    w = child.wait()
    return os.WEXITSTATUS(w), out


cmd = ['svnversion', ev]
stat, out = run_command(' '.join(cmd))
fsps_vers = int(re.match("^([0-9])+", out[0]).group(0)) #pull out only the beginning numeric characters

accepted = ((fsps_vers in ACCEPTED_FSPS_REVISIONS) and
            (len(out[0].split(':')) == 1) and
            stat == 0) #make sure you don't have some weird mixed version
if not accepted:
    raise ImportError("Your FSPS revision, {0}, is not known to work with this " \
                      "version of python-fsps. You can checkout an accepted FSPS revision "\
                      "with 'svn update -r rev_number'. The accepted FSPS rev_numbers are: " \
                      "{1}".format(out[0].rstrip('\n'), ACCEPTED_FSPS_REVISIONS))
    

from .fsps import StellarPopulation  # NOQA
from .filters import find_filter  # NOQA
