#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

from multiprocessing import Pool
import numpy as np
from numpy.testing import assert_allclose
from .fsps import StellarPopulation


pop = StellarPopulation(imf_type=2, imf3=2.3)


def test_imf3():
    pop.params["imf3"] = 2.3
    model1 = pop.get_spectrum(tage=0.2)[1]
    pop.params["imf3"] = 8.3
    assert pop.params.dirtiness == 2
    model2 = pop.get_spectrum(tage=0.2)[1]
    original = np.array(model1 - model2)

    pop.params["imf3"] = 2.3
    assert pop.params.dirtiness == 2
    model1 = pop.get_spectrum(tage=0.2)[1]
    pop.params["imf3"] = 8.3
    assert pop.params.dirtiness == 2
    model2 = pop.get_spectrum(tage=0.2)[1]

    assert_allclose(original, model1 - model2)


def _get_model(theta):
    pop.params["imf3"] = theta
    assert pop.params.dirtiness == 2
    return pop.get_spectrum(tage=0.2)[1]


def test_imf3_multiprocessing():
    pool = Pool()
    thetas = np.linspace(2.3, 8.3, 4)

    single = map(_get_model, thetas)
    multi = pool.map(_get_model, thetas)

    assert_allclose(single, multi)


def test_get_mags():
    fuv1 = pop.get_mags(bands=["galex_fuv"])[:, 0]
    mags = pop.get_mags()
    fuv2 = mags[:, 61]
    fuv3 = mags[:, 62]
    assert np.all(fuv1 == fuv2)
    assert np.all(fuv1 != fuv3)


def test_smoothspec():
    wave, spec = pop.get_spectrum(tage=1, peraa=True)
    spec2 = pop.smoothspec(wave, spec, 160., minw=1e3, maxw=1e4)
    assert (spec-spec2 == 0.).sum() > 0


def test_ztinterp():
    wave, s2 = pop.get_spectrum(tage=1, zmet=2, peraa=True)
    wave, s3 = pop.get_spectrum(tage=1, zmet=3, peraa=True)
    s, m, l = pop.ztinterp(-0.5, 1., peraa=True)

    optical = ((wave > 2000) & (wave < 5e3))
    assert s[optical].sum() < s2[optical].sum()
    assert s[optical].sum() > s3[optical].sum()
    assert m < 1


def test_exposedspec():
    allspec = pop.all_ssp_spec(peraa=True)
    assert len(allspec.shape) == 3
    assert allspec.shape[0] == len(pop.wavelengths)
