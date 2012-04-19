import os

try:
    ev = os.environ["SPS_HOME"]
except KeyError:
    raise ImportError("You need to have the SPS_HOME environment variable")

expected = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fsps")
if ev != expected:
    print "Warning: $SPS_HOME is set to %s. Expected %s."%(ev, expected)

from _fsps import fsps

# Because I don't know Fortran all that well, this is how I hack things to
# make sure that my magic numbers are right.
ntfull, nspec, nbands = fsps.get_dims()
assert ntfull == fsps.ntfull2 and nspec == fsps.nspec2 \
        and nbands == fsps.nbands2, \
        "There's a problem with the dimensions in fsps.f95"

