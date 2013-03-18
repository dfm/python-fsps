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
    This is the main interface to use when interacting with FSPS from Python.
    Most of the Fortran API is exposed through Python hooks with various
    features added for user friendliness. When initializing, you can set any
    of the parameters of the system using keyword arguments. Below, you'll
    find a list of the options that you can include (with the comments taken
    directly from the `FSPS docs
    <http://www.ucolick.org/~cconroy/FSPS_files/MANUAL.pdf>`_. To change
    these values later, use the ``params`` property—which is ``dict``-like.
    For example:

    ::

        sp = StellarPopulation(imf_type=2)
        sp.params["imf_type"] = 1

    :param dust_type: (default: 0)
        Common variable deﬁning the extinction curve for dust around old
        stars:

        * 0: power law with index dust index set by ``dust_index``.
        * 1: Milky Way extinction law (with :math:`R = A_V /E(B - V)` value
          ``mwr``) parameterized by Cardelli et al. (1989).
        * 2: Calzetti et al. (2000) attenuation curve. Note that if this
          value is set then the dust attenuation is applied to all starlight
          equally (not split by age), and therefore the only relevant
          parameter is ``dust2``, which sets the overall normalization.
        * 3: allows the user to access a variety of attenuation curve models
          from Witt & Gordon (2000) using the parameters ``wgp1`` and
          ``wgp2``.

    :param imf_type: (default: 2)
        Common variable defining the IMF type:

        * 0: Salpeter (1955)
        * 1: Chabrier (2003)
        * 2: Kroupa (2001)
        * 3: van Dokkum (2008)
        * 4: Dave (2008)
        * 5: tabulated piece-wise power law IMF, specified in ``imf.dat``
          file located in the data directory

    :param compute_vega_mags: (default: True)
        A switch that sets the zero points of the magnitude system: ``True``
        uses Vega magnitudes versus AB magnitudes.

    :param redshift_colors: (default: False)

        * ``False``: Magnitudes are computed at a fixed redshift specified
          by ``zred``.
        * ``True``: Magnitudes are computed at a redshift that corresponds
          to the age of the output SSP/CSP (assuming a redshift–age relation
          appropriate for a WMAP5 cosmology). This switch is useful if
          the user wants to compute the evolution in observed colors of a
          SSP/CSP.

    :param pagb: (default: 1.0)
        Weight given to the post–AGB phase. A value of 0.0 turns off post-AGB
        stars; a value of 1.0 implies that the Vassiliadis & Wood (1994)
        tracks are implemented as–is.

    :param dell: (default: 0.0)
        Shift in :math:`\log L_\mathrm{bol}` of the TP-AGB isochrones. Note
        that the meaning of this parameter and the one below has changed to
        reflect the updated calibrations presented in Conroy & Gunn (2009).
        That is, these parameters now refer to a modification about the
        calibrations presented in that paper.

    :param delt: (default: 0.0)
        Shift in :math:`\log T_\mathrm{eff}` of the TP-AGB isochrones.

    :param fbhb: (default: 0.0)
        Fraction of horizontal branch stars that are blue. The blue HB stars
        are uniformly spread in :math:`\log T_\mathrm{eff}` to `10^4` K. See
        Conroy et al. (2009a) for details and a plausible range.

    :param sbss: (default: 0.0)
        Specific frequency of blue straggler stars. See Conroy et al. (2009a)
        for details and a plausible range.

    :param tau: (default: 1.0)
        Defines e-folding time for the SFH, in Gyr. Only used if ``sfh=1`` or
        ``sfh=4``. The range is :math:`0.1 < \\tau < 10^2`.

    :param const: (default: 0.0)
        Defines the constant component of the SFH. This quantity is defined
        as the fraction of mass formed in a constant mode of SF; the range
        is therefore :math:`0 \le C \le 1`. Only used if ``sfh=1`` or
        ``sf=4``.

    :param tage: (default: 0.0)
        If set to a non-zero value, the
        :func:`fsps.StellarPopulation.compute_csp` method will compute the
        spectra and magnitudes only at this age, and will therefore only
        output one age result. The units are Gyr. (The default is to compute
        and return results from :math:`t \\approx 0` to the maximum age in
        the isochrones).

    :param fburst: (default: 0.0)
        Deﬁnes the fraction of mass formed in an instantaneous burst of star
        formation. Only used if ``sfh=1`` or ``sfh=4``.

    :param tburst: (default: 11.0)
        Defines the age of the Universe when the burst occurs. If
        ``tburst > tage`` then there is no burst. Only used if ``sfh=1`` or
        ``sfh=4``.

    :param dust1: (default: 0.0)
        Dust parameter describing the attenuation of young stellar light,
        i.e. where ``t <= dust_tesc`` (for details, see Conroy et al. 2009a).

    :param dust2: (default: 0.0)
        Dust parameter describing the attenuation of old stellar light,
        i.e. where ``t > dust_tesc`` (for details, see Conroy et al. 2009a).

    :param logzsol: (default: -0.2)
        Undocumented.

    """

    def __init__(self, **kwargs):
        # Before the first time we interact with the FSPS driver, we need to
        # run the ``setup`` method.
        if not driver.is_setup:
            driver.setup()

        # Set up the parameters to their default values.
        self.params = ParameterSet(
            dust_type=0,
            imf_type=2,
            compute_vega_mags=True,
            redshift_colors=False,
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
        self.params.dirty = False

    def compute_ssp(self, zi=None):
        if self.params.dirty:
            self._update_params()

        if zi is None:
            driver.ssps()
        else:
            assert 0 <= zi < NZ
            driver.ssp(int(zi) + 1)

    def compute_csp(self, zmet):
        if self.params.dirty:
            self._update_params()

        driver.compute(int(zmet) + 1)

    @property
    def wavelengths(self):
        return LAMBDA_GRID


class ParameterSet(object):

    def __init__(self, **kwargs):
        self.dirty = True
        self._params = kwargs
        self.iteritems = self._params.iteritems

    def __getitem__(self, k):
        return self._params[k]

    def __setitem__(self, k, v):
        self.dirty = True
        self._params[k] = v
