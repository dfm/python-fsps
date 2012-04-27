__all__ = ["compute"]

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
fn = os.path.join(os.environ["SPS_HOME"], "SPECTRA",
                    "BaSeL3.1", "basel.lambda")
wavelengths = np.array([l for l in open(fn)], dtype=float)


def compute(imf_type=0, zmet=1, sfh=0):
    assert imf_type in range(6)
    assert zmet in range(1, 23)
    assert sfh in range(5)

    driver.setup()
    driver.compute(imf_type, zmet, sfh)


def get_stats():
    nt = driver.get_ntfull()
    stats = np.vstack(driver.get_stats(nt)).T
    dt = [("age", float), ("mass", float), ("lbol", float), ("sfr", float),
            ("mdust", float)]
    return np.array([tuple(d) for d in stats], dtype=dt)


# def get_mag(band):
#     return fsps.mags[bands.index(band)]


# def get_spec():
#     return wavelengths, fsps.spec

if __name__ == "__main__":
    compute()
    print get_stats()
