# -*- coding: utf-8 -*-

__all__ = [
    "__version__",
    "StellarPopulation",
    "find_filter",
    "get_filter",
    "list_filters",
]

import os
import sys

here = os.path.dirname(os.path.abspath(__file__))

# A hack to handle extra dll libraries for some environments on Windows; ref:
# https://github.com/numpy/numpy/blob/01c9bfe4f48d23ec2d2db50ffc58d6e5e42cbe93/numpy/distutils/misc_util.py#L2309
extra_dll_dir = os.path.join(here, ".libs")
if sys.platform == "win32" and os.path.isdir(extra_dll_dir):
    if sys.version_info >= (3, 8):
        os.add_dll_directory(extra_dll_dir)
    else:
        os.environ.setdefault("PATH", "")
        os.environ["PATH"] += os.pathsep + extra_dll_dir

# Check the that SPS_HOME variable is set properly
from .sps_home import check_sps_home  # noqa

check_sps_home()

# Finally actually import the functions
try:
    from .fsps_version import version as __version__  # noqa
except ImportError:
    __version__ = "develop"
from .fsps import StellarPopulation  # noqa
from .filters import find_filter, get_filter, list_filters  # noqa
