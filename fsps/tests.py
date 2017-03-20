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

def test_nebemlineinspec():
    _reset_default_params()
    pop.params["sfh"] = 0
    pop.params['nebemlineinspec'] = False
    pop.params['add_neb_emission'] = True
    wave, spec_neboff = pop.get_spectrum()
    pop.params.dirtiness = 2
    pop.params['nebemlineinspec'] = True
    wave, spec_nebon = pop.get_spectrum()
    assert (spec_nebon-spec_neboff).sum() > 0
    ha_idx = (wave > 6556) & (wave < 6573)
    assert (spec_nebon-spec_neboff)[ha_idx].sum() > 0

def test_ssp():
    _reset_default_params()
    pop.params["sfh"] = 0
    wave, spec = pop.get_spectrum(tage=1, peraa=True)
    assert (wave[0] > 0) & (wave[0] < wave[-1])
    assert (wave[-1] > 1e6) & (wave[-1] < 1e10)
    Mv = 4.62  # AB absolute magnitude for a Zsol 1Gyr old SSP
    mag = pop.get_mags(tage=1, bands=["v"])
    assert np.all(abs(mag - Mv) < 1.0)
    assert np.all((pop.stellar_mass < 1.0) & (pop.stellar_mass > 0))


def test_csp():
    _reset_default_params()
    pop.params["sfh"] = 1
    pop.params["tau"] = 1.0
    wave, spec = pop.get_spectrum(tage=1.0)
    pop.params["tau"] = 3.0
    assert pop.params.dirtiness == 1


def test_redshift():
    _reset_default_params()
    pop.params["sfh"] = 0
    pop.params["zred"] = 0.0
    pop.params["add_igm_absorption"] = False
    v1 = pop.get_mags(redshift=1.0, tage=1.0, bands=["v"])
    v2 = pop.get_mags(redshift=1.0, tage=1.0, bands=["v"])
    assert np.all(v1 == v2)

    pop.params["zred"] = 1.0
    v3 = pop.get_mags(redshift=None, tage=1.0, bands=["v"])
    v4 = pop.get_mags(redshift=None, tage=1.0, bands=["v"])
    v5 = pop.get_mags(redshift=0.0, tage=1.0, bands=["v"])
    assert np.all(v3 == v4)

    assert np.all(v3 == v1)
    assert np.all(v5 != v4)


def test_tabular():
    _reset_default_params()

    import os
    fn = os.path.join(os.environ['SPS_HOME'], 'data/sfh.dat')
    age, sfr, z = np.genfromtxt(fn, unpack=True, skip_header=0)

    pop.params['sfh'] = 3
    pop.set_tabular_sfh(age, sfr)
    w, spec = pop.get_spectrum(tage=0)
    pop.set_tabular_sfh(age, sfr)
    assert pop.params.dirty
    w, spec = pop.get_spectrum(tage=0)
    assert spec.shape[0] == len(pop.ssp_ages)
    assert pop.params['sfh'] == 3
    w, spec_last = pop.get_spectrum(tage=-99)
    assert spec_last.ndim == 1
    w, spec = pop.get_spectrum(tage=age.max())
    assert np.allclose(spec / spec_last - 1., 0.)
    pop.params['logzsol'] = -1
    w, spec_lowz = pop.get_spectrum(tage=age.max())
    assert not np.allclose(spec / spec_lowz - 1., 0.)


def test_mformed():
    _reset_default_params()
    pop.params['sfh'] = 1
    pop.params['const'] = 0.5
    w, s = pop.get_spectrum(tage=0)
    assert pop.formed_mass[-1] == 1
    assert pop.formed_mass[50] < 1.0
    assert pop.formed_mass[50] > 0.0
    w, s = pop.get_spectrum(tage=0)
    assert pop.formed_mass[-1] == 1.0


def test_light_ages():
    _reset_default_params()
    tmax = 5.0
    pop.params['sfh'] = 1
    pop.params['const'] = 0.5
    w, spec = pop.get_spectrum(tage=tmax)
    mstar = pop.stellar_mass
    lbol = pop.log_lbol
    pop.params['compute_light_ages'] = True
    w, light_age = pop.get_spectrum(tage=tmax)
    assert np.all(np.abs(np.log10(spec / light_age)) > 1)
    # make sure fuv really from young stars
    fuv = (w > 1220) & (w < 2000)
    assert (light_age[fuv]).max() < 0.1
    assert (light_age[fuv]).max() > 1e-5
    assert pop.log_lbol != lbol
    assert pop.stellar_mass != mstar
    assert pop.stellar_mass < tmax
    # luminosity weighted age always less than mass-weighted age
    # assert pop.log_lbol < pop.stellar_mass


def test_libraries():
    _reset_default_params()
    ilib, splib = pop.libraries
    assert ilib == pop.isoc_library
    assert splib == pop.spec_library


def test_smooth_lsf():
    _reset_default_params()
    tmax = 1.0
    wave_lsf = np.arange(4000, 7000., 10)
    x = (wave_lsf - 5500) / 1500.
    # a quadratic lsf dependence that goes from ~50 to ~100 km/s
    sigma_lsf = 50 * (1.0 + 0.4 * x + 0.6 * x**2) 
    w, spec = pop.get_spectrum(tage=tmax)
    pop.params['smooth_lsf'] = True
    assert pop.params.dirtiness == 2
    pop.set_lsf(wave_lsf, sigma_lsf)
    w, smspec = pop.get_spectrum(tage=tmax)
    hi = w > 7100
    sm = (w < 7000) & (w > 3000)
    assert np.allclose(spec[hi]/smspec[hi] - 1., 0.)
    assert not np.allclose(spec[sm] / smspec[sm] - 1., 0.)
    pop.set_lsf(wave_lsf, sigma_lsf * 2)
    assert pop.params.dirtiness == 2
    
# Requires scipy
# def test_sfr_avg():

#    _reset_default_params()
#    tmax = 5.0
#    pop.params['sfh'] = 1.0
#    pop.params['const'] = 0.5
#    w, spec = pop.get_spectrum(tage=0)
#    sfr6 = pop.sfr_avg(dt=1e-3)
#    dsfr = np.log10(pop.sfr/pop.sfr6)
#    good = pop.log_age > 6
#    assert np.all(np.abs(dsfr[good]) < 1e-2)
