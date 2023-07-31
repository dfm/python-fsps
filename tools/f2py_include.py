#!/usr/bin/env python3

import numpy.f2py

try:
    print(numpy.f2py.get_include())
except AttributeError:
    import os

    print(os.path.join(os.path.dirname(numpy.f2py.__file__), "src"))
