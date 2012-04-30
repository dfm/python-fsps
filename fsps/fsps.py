__all__ = ["bands", "ntfull", "ssps", "compute", "get_stats"]

import os
import numpy as np
from _fsps import driver

bands = ["v", "u", "deep_b", "deep_r", "deep_i", "twomass_j", "twomass_h",
    "twomass_k", "sdss_u", "sdss_g", "sdss_r", "sdss_i", "sdss_z",
    "wfc_f435w", "wfc_f606w", "wfc_f775w", "wfc_f814w", "wfc_f850lp",
    "irac1", "irac2", "irac3", "irac4", "nicmos_f110w", "nicmos_f160w",
    "fors_v", "fors_r", "galex_nuv", "galex_fuv", "wfcam_z", "wfcam_y",
    "wfcam_j", "wfcam_h", "wfcam_k", "b", "r", "i", "b3", "wfpc2_f555w",
    "wfpc2_f814w", "wfc3_f275w"]

# Load the wavelength values for the spectra.
_fn = os.path.join(os.environ["SPS_HOME"], "SPECTRA",
                    "BaSeL3.1", "basel.lambda")
wavelengths = np.array([l for l in open(_fn)], dtype=float)

# Set up the FSPS module.
driver.setup()

# Get some of the parameters of the module.
ntfull = driver.get_ntfull()
nspec = driver.get_nspec()
nbands = driver.get_nbands()

# Get the wavelength bins.
wavelengths = driver.get_lambda(nspec)


def ssps(imf_type=0, imf1=1.3, imf2=2.3, imf3=2.3, vdmc=0.08, mdave=0.5,
        dell=0.0, delt=0.0, sbss=0.0, fbhb=0.0, pagb=1.0):
    # First check to make sure that the inputs are sane.
    assert imf_type in range(6)

    driver.ssps(imf_type, imf1, imf2, imf3, vdmc, mdave, dell, delt, sbss,
            fbhb, pagb)


def compute(zmet, dust_type=0, sfh=0, tau=1.0, const=0.0, fburst=0.0,
        tburst=0.0, dust_tesc=7.0, dust1=0.0, dust2=0.0, dust_clumps=-99.0,
        frac_no_dust=0.0, dust_index=-0.7, mwr=3.1, wgp1=1, wgp2=1, wgp3=0,
        tage=0.0):
    assert zmet in range(1, 23)
    assert sfh in range(5)
    assert 0.1 < tau < 100
    assert 0 <= const <= 1
    assert 0 <= fburst <= 1
    assert wgp1 in range(1, 19)
    assert wgp2 in range(1, 7)
    assert wgp3 in [0, 1]

    driver.compute(dust_type, zmet, sfh, tau, const, fburst, tburst,
            dust_tesc, dust1, dust2, dust_clumps, frac_no_dust, dust_index,
            mwr, wgp1, wgp2, wgp3, tage)


def get_stats():
    stats = np.vstack(driver.get_stats(ntfull)).T
    dt = [("age", float), ("mass", float), ("lbol", float), ("sfr", float),
            ("mdust", float)]
    return np.array([tuple(d) for d in stats], dtype=dt)


def get_spec():
    return wavelengths, driver.get_spec(nspec, ntfull)


def get_isochrones(zmet):
    assert zmet in range(1, 23)
    n_age, n_mass = driver.get_isochrone_dimensions()
    [driver.get_isochrone(zmet, t,
        driver.get_nmass_isochrone(zmet, t), nbands) for t in range(n_age)]

if __name__ == "__main__":
    import time

    strt = time.time()
    ssps()
    print("Calculating the SSPs took {0:.4f} seconds."
            .format(time.time() - strt))
    strt = time.time()
    compute(1)
    print("Calculating the CSP took {0:.4f} seconds."
            .format(time.time() - strt))

    get_isochrones(5)
