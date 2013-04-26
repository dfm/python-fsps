#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import pyfits
import numpy as np
import matplotlib.pyplot as pl
import fsps

# Measurements of cluster parameters.
tage = 10. ** (8.04 - 9)
logmass = 4.09
dist_mod = 24.5

# Set up the stellar population model.
sp = fsps.StellarPopulation(imf_type=2, dust_type=1, mwr=3.1, dust2=0.3)

# The measured magnitudes from the literature.
data = {"wfc3_f160w": 16.386,
        "wfc3_f275w": 17.398,
        "wfc_acs_f814w": 17.155}

# There are also a few filters that we have data for but aren't included in
# the standard FSPS install:
#   "F110W": 16.877,
#   "F336W": 17.349,
#   "F475W": 17.762,

# Load the observed spectrum.
f = pyfits.open("DAO69.fits")
obs_spec = np.array(f[0].data, dtype=float)
f.close()

obs_spec /= 5e-20

# The observed wavelength grid in the data is magically this:
obs_lambda = np.arange(0, 4540) * 1.2 + 3700

# Compute the model magnitudes.
for b, v in data.iteritems():
    print(b, v, sp.get_mags(zmet=20, tage=tage, band=b) - 2.5 * logmass
          + dist_mod)

# Compute the model spectrum in ``L_sun / A``.
lam, spec = sp.get_spectrum(zmet=20, tage=tage, peraa=True)
spec *= 3.839e33 * 10. ** (logmass - dist_mod / 2.5)

f = 1.0  # obs_spec[0, obs_lambda < 5000.][-1] / spec[lam < 5000.][-1]
print(obs_spec[0, obs_lambda < 5000.][-1] / spec[lam < 5000.][-1])

pl.loglog(obs_lambda, obs_spec[0], "k")
pl.loglog(lam, spec * f, "r")

pl.xlim(3700, 7000)
# pl.ylim(10 ** 3, 10 ** 4.5)
pl.savefig("spectrum.png")
