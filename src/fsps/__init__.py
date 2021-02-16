# -*- coding: utf-8 -*-

__all__ = [
    "__version__",
    "StellarPopulation",
    "find_filter",
    "get_filter",
    "list_filters",
]

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
