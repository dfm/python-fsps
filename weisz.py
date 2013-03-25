#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import time
import fsps
import pyfits
import numpy as np
import matplotlib.pyplot as pl

strt = time.time()
sp = fsps.StellarPopulation(imf_type=2, dust_type=1, mwr=3.1, dust2=0.3)
print("Loading took {} seconds".format(time.time() - strt))

# "nicmos_f110w": 16.877,
# "F336W": 17.349,
# "F475W": 17.762,

data = {"wfc3_f160w": 16.386,
        "wfc3_f275w": 17.398,
        "wfc_acs_f814w": 17.155}

tage = 10. ** (8.04 - 9)
logmass = 4.09

strt = time.time()
for b, v in data.iteritems():
    # print(b, v, sp.get_mags(zmet=20, tage=tage, band=b) - 2.5 * logmass + 24.5)
    sp.get_mags(zmet=20, tage=tage, band=b) - 2.5 * logmass + 24.5
print("Computing mags took {} seconds".format(time.time() - strt))

# Plot observations.
f = pyfits.open("DAO69.fits")
obs_spec = np.array(f[0].data, dtype=float)
f.close()
obs_lambda = np.arange(0, 4540) * 1.2 + 3700
pl.loglog(obs_lambda, obs_spec[0], "k")

# Plot model.
# factor = 3e18 * 5e-20
strt = time.time()
lam, spec = sp.get_spectrum(zmet=20, tage=tage)
print("Computing spectrum took {} seconds".format(time.time() - strt))

spec /= lam ** 2
f = obs_spec[0, obs_lambda < 5000.][-1] / spec[lam < 5000.][-1]
pl.loglog(lam, spec * f, "r")

pl.xlim(3700, 9000)
pl.ylim(10 ** 3, 10 ** 4.5)
pl.savefig("spectrum.png")
