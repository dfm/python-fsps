import os

try:
    ev = os.environ["SPS_HOME"]
except KeyError:
    raise ImportError("You need to have the SPS_HOME environment variable")

expected = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fsps")
if ev[-1] == "/":
    ev = ev[:-1]
if ev != expected:
    print "Warning: $SPS_HOME is set to %s. Expected %s."%(ev, expected)

from fsps import *

