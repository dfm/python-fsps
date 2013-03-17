#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["StellarPopulation"]

from ._fsps import driver


# Hard-set FSPS parameters.
NZ = driver.get_nz()
NTFULL = driver.get_ntfull()
NSPEC = driver.get_nspec()
NBANDS = driver.get_nbands()
LAMBDA_GRID = driver.get_lambda(NSPEC)
NAGE, NMASS = driver.get_isochrone_dimensions()


def _get_nmass_isochrone(z, t):
    return driver.get_nmass_isochrone(z, t)


class StellarPopulation(object):
    """

    """

    def __init__(self, **kwargs):
        # Before the first time we interact with the FSPS driver, we need to
        # run the ``setup`` method.
        if not driver.is_setup:
            driver.setup()

        # Set up the parameters to their default values.
        self.params = dict(
            dust_type=0,
            imf_type=2,
            compute_vega_mags=0,
            redshift_colors=0,
            pagb=1.0,
            dell=0.0,
            delt=0.0,
            fbhb=0.0,
            sbss=0.0,
            tau=1.0,
            const=0.0,
            tage=0.0,
            fburst=0.0,
            tburst=11.0,
            dust1=0.0,
            dust2=0.0,
            logzsol=-0.2,
            zred=0.0,
            pmetals=0.02,
            imf1=1.3,
            imf2=2.3,
            imf3=2.3,
            vdmc=0.08,
            dust_clumps=-99.,
            frac_nodust=0.0,
            dust_index=-0.7,
            dust_tesc=7.0,
            frac_obrun=0.0,
            uvb=1.0,
            mwr=3.1,
            redgb=1.0,
            dust1_index=-1.0,
            mdave=0.5,
            sf_start=0.0,
            sf_trunc=0.0,
            sf_theta=0.0,
            duste_gamma=0.01,
            duste_umin=1.0,
            duste_qpah=3.5,
            fcstar=1.0,
            masscut=150.0,
            zmet=1,
            sfh=0,
            wgp1=1,
            wgp2=1,
            wgp3=1,
            evtype=-1
        )
        self._dirty = True
        for k, v in self.params.iteritems():
            self.params[k] = kwargs.pop(k, v)

        if len(kwargs):
            raise TypeError("__init__() got an unexpected keyword argument "
                            "'{}'".format(kwargs.keys()[0]))

    def _update_params(self):
        keys = ["dust_type", "imf_type", "compute_vega_mags",
                "redshift_colors", "zmet", "sfh", "wgp1", "wgp2", "wgp3",
                "evtype", "pagb", "dell", "delt", "fbhb", "sbss", "tau",
                "const", "tage", "fburst", "tburst", "dust1", "dust2",
                "logzsol", "zred", "pmetals", "imf1", "imf2", "imf3", "vdmc",
                "dust_clumps", "frac_nodust", "dust_index", "dust_tesc",
                "frac_obrun", "uvb", "mwr", "redgb", "dust1_index", "mdave",
                "sf_start", "sf_trunc", "sf_theta", "duste_gamma",
                "duste_umin", "duste_qpah", "fcstar", "masscut"]
        driver.set_params(*[self.params[k] for k in keys])
        self._dirty = False

    def compute_ssp(self, zi=None):
        if self._dirty:
            self._update_params()

        if zi is None:
            driver.ssps()
        else:
            assert 0 <= zi < NZ
            driver.ssp(int(zi) + 1)

    def compute_csp(self, zmet):
        if self._dirty:
            self._update_params()

        driver.compute(int(zmet) + 1)

    @property
    def wavelengths(self):
        return LAMBDA_GRID
