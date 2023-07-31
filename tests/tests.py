# -*- coding: utf-8 -*-

import os
import sys

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fsps import StellarPopulation, filters

skip_slow_tests = pytest.mark.skipif(
    (sys.platform.startswith("win") or sys.platform.startswith("darwin"))
    and "CI" in os.environ,
    reason="Slow tests only run on Linux CI",
)


@pytest.fixture(scope="module")
def pop_and_params():
    pop = StellarPopulation(zcontinuous=1)
    params = dict([(k, pop.params[k]) for k in pop.params.all_params])
    return pop, params


def _reset_default_params(pop, params):
    pop._zcontinuous = 1
    for k in pop.params.all_params:
        pop.params[k] = params[k]


@skip_slow_tests
def test_isochrones(pop_and_params):
    """Just test that `isochrones()` method runs"""

    # recomputes SSPs

    pop, params = pop_and_params
    _reset_default_params(pop, params)
    pop.params["imf_type"] = 0
    pop.isochrones()


@skip_slow_tests
def test_imf3(pop_and_params):
    """Make sure that changing the (upper) imf changes the parameter dirtiness
    and also the SSP spectrum"""

    # recomputes SSPs

    pop, params = pop_and_params
    _reset_default_params(pop, params)
    pop.params["imf_type"] = 2
    pop.params["imf3"] = 2.3
    w, model1 = pop.get_spectrum(tage=0.2)

    # check that changing the IMF does something
    pop.params["imf3"] = 8.3
    assert pop.params.dirtiness == 2
    w, model2 = pop.get_spectrum(tage=0.2)
    assert not np.allclose(model1 / model2 - 1.0, 0.0)

    # Do we *really* need to do this second check?
    pop.params["imf3"] = 2.3
    assert pop.params.dirtiness == 2
    w, model1b = pop.get_spectrum(tage=0.2)
    assert pop.params.dirtiness == 0

    assert_allclose(model1 / model1b - 1.0, 0.0)


def test_param_checks(pop_and_params):
    # recomputes SSPs

    pop, params = pop_and_params
    _reset_default_params(pop, params)
    pop.params["sfh"] = 1
    pop.params["tage"] = 2
    pop.params["sf_start"] = 0.5
    # this should never raise an error:
    w, s = pop.get_spectrum(tage=pop.params["tage"])
    # this used to raise an assertion error in the setter:
    pop.params["sf_start"] = pop.params["tage"] + 0.1
    # this also used to raise an assertion error in the setter:
    pop.params["imf_type"] = 8
    # fix the invalid IMF but leave the invalid sf_start > tage
    pop.params["imf_type"] = 1
    with pytest.raises(AssertionError):
        w, s = pop.get_spectrum(tage=pop.params["tage"])
    pop.params["tage"] = 1.0
    pop.params["sf_start"] = 0.1
    w, s = pop.get_spectrum(tage=pop.params["tage"])


def test_smooth_lsf(pop_and_params):
    # recomputes SSPs

    pop, params = pop_and_params
    _reset_default_params(pop, params)
    tmax = 1.0
    wave_lsf = np.arange(4000, 7000.0, 10)
    x = (wave_lsf - 5500) / 1500.0
    # a quadratic lsf dependence that goes from ~50 to ~100 km/s
    sigma_lsf = 50 * (1.0 + 0.4 * x + 0.6 * x**2)
    w, spec = pop.get_spectrum(tage=tmax)
    pop.params["smooth_lsf"] = True
    assert pop.params.dirtiness == 2
    pop.set_lsf(wave_lsf, sigma_lsf)
    w, smspec = pop.get_spectrum(tage=tmax)
    hi = w > 7100
    sm = (w < 7000) & (w > 3000)
    assert np.allclose(spec[hi] / smspec[hi] - 1.0, 0.0)
    assert not np.allclose(spec[sm] / smspec[sm] - 1.0, 0.0)
    pop.set_lsf(wave_lsf, sigma_lsf * 2)
    assert pop.params.dirtiness == 2


@skip_slow_tests
def test_tabular(pop_and_params):
    """Test that you get the right shape spectral arrays for tabular SFHs, that
    the parameter dirtiness is appropriately managed for changing tabular SFH,
    and that multi-metallicity SFH work."""

    # uses default SSPs, but makes them for every metallicity

    pop, params = pop_and_params
    _reset_default_params(pop, params)

    import os

    fn = os.path.join(os.environ["SPS_HOME"], "data/sfh.dat")
    age, sfr, z = np.genfromtxt(fn, unpack=True, skip_header=0)

    # Mono-metallicity
    pop.params["sfh"] = 3
    pop.set_tabular_sfh(age, sfr)
    w, spec = pop.get_spectrum(tage=0)
    pop.set_tabular_sfh(age, sfr)
    assert pop.params.dirty
    w, spec = pop.get_spectrum(tage=0)
    assert spec.shape[0] == len(pop.ssp_ages)
    assert pop.params["sfh"] == 3
    w, spec_last = pop.get_spectrum(tage=-99)
    assert spec_last.ndim == 1
    w, spec = pop.get_spectrum(tage=age.max())
    assert np.allclose(spec / spec_last - 1.0, 0.0)
    pop.params["logzsol"] = -1
    w, spec_lowz = pop.get_spectrum(tage=age.max())
    assert not np.allclose(spec / spec_lowz - 1.0, 0.0)

    # test the formed mass for single age
    assert np.allclose(np.trapz(sfr, age) * 1e9, pop.formed_mass)

    # Multi-metallicity
    pop._zcontinuous = 3
    pop.set_tabular_sfh(age, sfr, z)
    w, spec_multiz = pop.get_spectrum(tage=age.max())
    assert not np.allclose(spec_lowz / spec_multiz - 1.0, 0.0)

    pop._zcontinuous = 1
    pop.set_tabular_sfh(age, sfr)
    # get mass weighted metallicity
    mbin = np.gradient(age) * sfr
    mwz = (z * mbin).sum() / mbin.sum()
    pop.params["logzsol"] = np.log10(mwz / pop.solar_metallicity)
    w, spec_onez = pop.get_spectrum(tage=age.max())
    assert not np.allclose(spec_onez / spec_multiz - 1.0, 0.0)


def test_get_mags(pop_and_params):
    """Basic test for supplying filter names to get_mags"""

    # uses default SSPs

    pop, params = pop_and_params
    _reset_default_params(pop, params)
    fuv1 = pop.get_mags(bands=["galex_fuv"])[:, 0]
    mags = pop.get_mags()
    fuv2 = mags[:, 61]  # this should be galex_FUV
    fuv3 = mags[:, 62]  # this should *not* be galex_FUV
    assert np.all(fuv1 == fuv2)
    assert np.all(fuv1 != fuv3)


def test_ssp(pop_and_params):
    """Test that you get a sensible wavelength array, and that you get a
    sensible V-band magnitude for a 1-Gyr SSP"""

    # uses default SSPs

    pop, params = pop_and_params
    _reset_default_params(pop, params)
    pop.params["sfh"] = 0
    wave, spec = pop.get_spectrum(tage=1, peraa=True)
    assert (wave[0] > 0) & (wave[0] < wave[-1]) & (wave[0] < 912.0)
    assert (wave[-1] > 1e6) & (wave[-1] < 1e10)
    Mv = 4.62  # AB absolute magnitude for a Zsol 1Gyr old SSP
    # This also tests get_mags
    mag = pop.get_mags(tage=1, bands=["v"])
    assert np.all(abs(mag - Mv) < 1.0)
    assert np.all((pop.stellar_mass < 1.0) & (pop.stellar_mass > 0))
    assert pop.params.dirtiness == 0


def test_libraries(pop_and_params):
    """This does not require or build clean SSPs"""

    # uses default SSPs

    pop, params = pop_and_params
    _reset_default_params(pop, params)
    ilib, splib, dlib = pop.libraries
    assert ilib == pop.isoc_library
    assert splib == pop.spec_library
    assert dlib == pop.duste_library


def test_filters():
    """Test all the filters got transmission data loaded."""

    # uses default SSPs

    flist = filters.list_filters()
    # force trans cache to load
    filters.FILTERS[flist[0]]._load_transmission_cache()
    for f in flist:
        assert f in filters.TRANS_CACHE, "transmission not loaded for {}".format(f)


def test_csp_dirtiness(pop_and_params):
    """Make sure that changing CSP parameters increases dirtiness to 1"""
    # uses default SSPs

    pop, params = pop_and_params
    _reset_default_params(pop, params)
    pop.params["sfh"] = 1
    pop.params["tau"] = 1.0
    wave, spec = pop.get_spectrum(tage=1.0)
    assert pop.params.dirtiness == 0
    pop.params["tau"] = 3.0
    assert pop.params.dirtiness == 1


def test_redshift(pop_and_params):
    """Test redshifting, make sure that
    1. redshifting does not persist in cached arrays
    2. specifying redshift via get_mags keyword or param key are consistent
    """

    # uses default SSPs

    pop, params = pop_and_params
    _reset_default_params(pop, params)
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


def test_dust3(pop_and_params):
    """Make sure nebular lines are actually added."""

    # uses default SSPs

    pop, params = pop_and_params
    _reset_default_params(pop, params)
    pop.params["sfh"] = 4
    pop.params["dust_type"] = 4
    pop.params["tau"] = 5.0
    pop.params["dust1"] = 0
    pop.params["dust2"] = 0.5

    # make sure dust3 affects the population when there are old stars
    pop.params["dust3"] = 0.0
    mag1 = pop.get_mags(tage=1.0, bands=["u", "b"])
    pop.params["dust3"] = 1.0
    mag2 = pop.get_mags(tage=1.0, bands=["u", "b"])
    assert np.all(mag2 > mag1)

    # make sure the dust3 isn't affecting young populations
    pop.params["dust_tesc"] = 8
    pop.params["dust3"] = 0.0
    mag3 = pop.get_mags(tage=0.05, bands=["u", "b"])
    pop.params["dust3"] = 1.0
    mag4 = pop.get_mags(tage=0.05, bands=["u", "b"])
    assert_allclose(mag3, mag4)


def test_nebemlineinspec(pop_and_params):
    """Make sure nebular lines are actually added."""

    # uses default SSPs

    pop, params = pop_and_params
    _reset_default_params(pop, params)
    pop.params["sfh"] = 4
    pop.params["tau"] = 5.0
    pop.params["add_neb_emission"] = True
    pop.params["nebemlineinspec"] = False
    wave, spec_neboff = pop.get_spectrum(tage=1.0)
    pop.params["nebemlineinspec"] = True
    wave, spec_nebon = pop.get_spectrum(tage=1.0)
    assert (spec_nebon - spec_neboff).sum() > 0
    assert np.all(np.isfinite(pop.emline_luminosity))
    assert np.all(np.isfinite(pop.emline_wavelengths))
    ha_idx = (wave > 6556) & (wave < 6573)
    assert (spec_nebon - spec_neboff)[ha_idx].sum() > 0


def test_mformed(pop_and_params):
    """Make sure formed mass integrates to 1 for parameteric SFH"""
    # uses default SSPs

    pop, params = pop_and_params
    _reset_default_params(pop, params)
    pop.params["sfh"] = 1
    pop.params["const"] = 0.5
    w, s = pop.get_spectrum(tage=0)
    assert pop.formed_mass[-1] == 1
    assert pop.formed_mass[50] < 1.0
    assert pop.formed_mass[50] > 0.0
    w, s = pop.get_spectrum(tage=0)
    assert pop.formed_mass[-1] == 1.0


def test_light_ages(pop_and_params):
    """Make sure light weighting works, and gives sensible answers for the
    light-weighted age in the FUV and mass-weighted age stored in
    `stellar_mass`"""
    # uses default SSPs

    pop, params = pop_and_params
    _reset_default_params(pop, params)
    tmax = 5.0
    pop.params["sfh"] = 1
    pop.params["const"] = 0.5
    w, spec = pop.get_spectrum(tage=tmax)
    mstar = pop.stellar_mass
    lbol = pop.log_lbol
    pop.params["compute_light_ages"] = True
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


def test_smoothspec(pop_and_params):
    # FIXME: This is not very stringent

    # uses default SSPs

    pop, params = pop_and_params
    _reset_default_params(pop, params)
    wave, spec = pop.get_spectrum(tage=1, peraa=True)
    spec2 = pop.smoothspec(wave, spec, 160.0, minw=1e3, maxw=1e4)
    assert (spec - spec2 == 0.0).sum() > 0


@skip_slow_tests
def test_ssp_weights(pop_and_params):
    """Check that weights dotted into ssp is the same as the returned spectrum
    when there's no dust or emission lines and zcontinuous=0"""

    # uses default SSPs

    pop, params = pop_and_params
    _reset_default_params(pop, params)

    import os

    fn = os.path.join(os.environ["SPS_HOME"], "data/sfh.dat")
    age, sfr, z = np.genfromtxt(fn, unpack=True, skip_header=0)
    pop.params["sfh"] = 3
    pop.set_tabular_sfh(age, sfr)
    zind = -3
    pop.params["logzsol"] = np.log10(pop.zlegend[zind] / pop.solar_metallicity)

    wave, spec = pop.get_spectrum(tage=age.max())
    mstar = pop.stellar_mass
    wght = pop._ssp_weights
    ssp, smass, slbol = pop._all_ssp_spec()

    assert np.allclose((smass[:, zind] * wght[:, 0]).sum(), mstar)


# Requires scipy
# def test_sfr_avg():

#    _reset_default_params()
#    pop.params['sfh'] = 1.0
#    pop.params['const'] = 0.5
#    w, spec = pop.get_spectrum(tage=0)
#    sfr6 = pop.sfr_avg(dt=1e-3)
#    dsfr = np.log10(pop.sfr/pop.sfr6)
#    good = pop.ssp_age > 6
#    assert np.all(np.abs(dsfr[good]) < 1e-2)


# def test_imf3_multiprocessing():
#     from multiprocessing import Pool
#     pool = Pool()
#     thetas = np.linspace(2.3, 8.3, 4)
#     single = map(_get_model, thetas)
#     multi = pool.map(_get_model, thetas)
#     assert_allclose(single, multi)
