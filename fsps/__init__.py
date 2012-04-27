import os

try:
    ev = os.environ["SPS_HOME"]
except KeyError:
    raise ImportError("You need to have the SPS_HOME environment variable")

from fsps import *
