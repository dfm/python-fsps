#!/usr/bin/env python3

import os

import numpy.f2py

try:
    incl_dir = numpy.f2py.get_include()
except AttributeError:
    incl_dir = os.path.join(os.path.dirname(numpy.f2py.__file__), "src")

print(os.path.normpath(incl_dir))
