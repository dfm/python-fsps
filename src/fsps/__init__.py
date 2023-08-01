# Check the that SPS_HOME variable is set properly
from fsps.sps_home import check_sps_home

check_sps_home()

del check_sps_home
# End check

from fsps.filters import find_filter as find_filter
from fsps.filters import get_filter as get_filter
from fsps.filters import list_filters as list_filters
from fsps.fsps import StellarPopulation as StellarPopulation
from fsps.fsps_version import __version__ as __version__
