# -*- coding: utf-8 -*-

from __future__ import division, print_function

from multiprocessing import Pool
import numpy as np
from numpy.testing import assert_allclose
from .fsps import StellarPopulation


pop = StellarPopulation(zcontinuous=1)
default_params = dict([(k, pop.params[k]) for k in pop.params.all_params])


def _reset_default_params():
    for k in pop.params.all_params:
        pop.params[k] = default_params[k]


def test_imf3():
    _reset_default_params()
    pop.params["imf_type"] = 2
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


# def test_imf3_multiprocessing():
#     pool = Pool()
#     thetas = np.linspace(2.3, 8.3, 4)

#     single = map(_get_model, thetas)
#     multi = pool.map(_get_model, thetas)

#     assert_allclose(single, multi)


def test_get_mags():
    fuv1 = pop.get_mags(bands=["galex_fuv"])[:, 0]
    mags = pop.get_mags()
    fuv2 = mags[:, 61]
    fuv3 = mags[:, 62]
    assert np.all(fuv1 == fuv2)
    assert np.all(fuv1 != fuv3)


def test_smoothspec():
    _reset_default_params()
    wave, spec = pop.get_spectrum(tage=1, peraa=True)
    spec2 = pop.smoothspec(wave, spec, 160., minw=1e3, maxw=1e4)
    assert (spec-spec2 == 0.).sum() > 0


def test_ssp():
    _reset_default_params()
    pop.params["sfh"] = 0
    wave, spec = pop.get_spectrum(tage=1, peraa=True)
    assert (wave[0] > 0) & (wave[0] < wave[-1])
    assert (wave[-1] > 1e6) & (wave[-1] < 1e10)
    Mv = 4.62 # AB absolute magnitude for a Zsol 1Gyr old SSP
    mag = pop.get_mags(tage=1, bands=["v"])
    assert np.all( abs(mag - Mv) < 1.0)
    assert np.all( (pop.stellar_mass < 1.0) & (pop.stellar_mass >0))


def test_csp():
    _reset_default_params()
    pop.params["sfh"] = 1
    pop.params["tau"] = 1.0
    wave, spec = pop.get_spectrum(tage=1.0)
    pop.params["tau"] = 3.0
    assert pop.params.dirtiness == 1
    pop.params["sfh"] = 0
    pop.params["tau"] = 1.0


def test_redshift():
    _reset_default_params()
    pop.params["sfh"] = 0
    pop.params["zred"] = 0.0
    v1 = pop.get_mags(redshift=1.0, tage=1.0, bands=["v"])
    v2 = pop.get_mags(redshift=1.0, tage=1.0, bands=["v"])
    assert np.all(v1 == v2)

    pop.params["zred"] = 1.0
    v3 = pop.get_mags(redshift=0.0, tage=1.0, bands=["v"])
    v4 = pop.get_mags(redshift=0.0, tage=1.0, bands=["v"])
    assert np.all(v3 == v4)

    # The following fails for now, because of how redshifting and filter projection is
    # delegated in and accessed from FSPS.  The difference will be dist. mod. - 2.5*log(1+zred)

    # assert np.all(v3 == v1)

    v5 = pop.get_mags(redshift=1.0, tage=1.0, bands=["v"])
    assert np.all(v5 != v4)
    
