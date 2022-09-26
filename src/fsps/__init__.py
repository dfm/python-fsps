# Check the that SPS_HOME variable is set properly
from fsps.sps_home import check_sps_home

check_sps_home()

del check_sps_home
# End check

from fsps.fsps_version import version as __version__
from fsps.fsps import StellarPopulation as StellarPopulation
from fsps.filters import (
    find_filter as find_filter,
    get_filter as get_filter,
    list_filters as list_filters,
)
