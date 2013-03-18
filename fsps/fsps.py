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

    :param zred: (default: 0.0)
        Redshift. If this value is non-zero and if ``redshift_colors=1``,
        the magnitudes will be computed for the spectrum placed at redshift
        ``zred``.

    :param pmetals: (default: 0.02)
        Undocumented.

    :param imf1: (default: 1.3)
        Logarithmic slope of the IMF over the range :math:`0.08 < M < 0.5
        M_\odot`. Only used if ``imf_type=2``.

    :param imf2: (default: 2.3)
        Logarithmic slope of the IMF over the range :math:`0.5 < M < 1
        M_\odot`. Only used if ``imf_type=2``.

    :param imf3: (default: 2.3)
        Logarithmic slope of the IMF over the range :math:`1.0 < M < 100
        M_\odot`. Only used if ``imf_type=2``.

    :param vdmc: (default: 0.08)
        IMF parameter defined in van Dokkum (2008). Only used if
        ``imf_type=3``.

    :param dust_clumps: (default: -99.)
        Dust parameter describing the dispersion of a Gaussian PDF density
        distribution for the old dust. Setting this value to -99.0 sets the
        distribution to a uniform screen. See Conroy et al. (2009b) for
        details.

    :param frac_nodust: (default: 0.0)
        Fraction of starlight that is not attenuated by the diffuse dust
        component (i.e. that is not affected by ``dust2``).

    :param dust_index: (default: -0.7)
        Power law index of the attenuation curve. Only used when
        ``dust_type=0``.

    :param dust_tesc: (default: 7.0)
        Stars younger than ``dust_tesc`` are attenuated by both ``dust1`` and
        ``dust2``, while stars older are attenuated by ``dust2`` only. Units
        are :math:`\\log (\\mathrm{yrs})`.

    :param frac_obrun: (default: 0.0)
        Undocumented.

    :param uvb: (default: 1.0)
        Parameter characterizing the strength of the 2175A extinction feature
        with respect to the standard Cardelli et al. determination for the
        MW. Only used when ``dust_type=1``.

    :param mwr: (default: 3.1)
        The ratio of total to selective absorption which characterizes the MW
        extinction curve: :math:`R = A_V /E(B - V)`. Only used when
        ``dust_type=1``.

    :param redgb: (default: 1.0)
        Undocumented.

    :param dust1_index: (default: -1.0)
        Undocumented.

    :param mdave: (default: 0.5)
        IMF parameter defined in Dave (2008). Only used if ``imf_type=4``.

    :param sf_start: (default: 0.0)
        Start time of the SFH, in Gyr.

    :param sf_trunc: (default: 0.0)
        Undocumented.

    :param sf_theta: (default: 0.0)
        Undocumented.

    :param duste_gamma: (default: 0.01)
        Parameter of the Draine & Li (2007) dust emission model. Specifies
        the relative contribution of dust heated at a radiation field
        strength of :math:`U_\mathrm{min}` and dust heated at
        :math:`U_\mathrm{min} < U \le U_\mathrm{max}`. Allowable range is 0.0
        – 1.0.

    :param duste_umin: (default: 1.0)
        Parameter of the Draine & Li (2007) dust emission model. Specifies
        the minimum radiation field strength in units of the MW value. Valid
        range is 0.1 – 25.0.

    :param duste_qpah: (default: 3.5)
        Parameter of the Draine & Li (2007) dust emission model. Specifies
        the grain size distribution through the fraction of grain mass in
        PAHs. This parameter has units of % and a valid range of 0.0 − 10.0.

    :param fcstar: (default: 1.0)
        Fraction of stars that the Padova isochrones identify as Carbon stars
        that FSPS assigns to a Carbon star spectrum. Set this to 0.0 if for
        example the users wishes to turn all Carbon stars into regular M-type
        stars.

    :param masscut: (default: 150.0)
        Undocumented.

    :param zmet: (default: 1)
        The metallicity is specified as an integer ranging between 1 and 22.

    :param sfh: (default: 0)
        Defines the type of star formation history, normalized such that one
        solar mass of stars is formed over the full SFH. Default value is 0.

        * 0: Compute an SSP
        * 1: Compute a five parameter SFH (see below).
        * 2: Compute a tabulated SFH defined in a file called ``sfh.dat``
          that must reside in the data directory. The file must contain three
          rows. The first column is time since the Big Bang in Gyr, the
          second is the SFR in units of solar masses per year, the third is
          the absolute metallicity. An example is provided in the data
          directory. The time grid in this file can be arbitrary (so long as
          the units are correct), but it is up to the user to ensure that the
          tabulated SFH is well-sampled so that the outputs are stable.
          Obviously, highly oscillatory data require dense sampling.
        * 4: Delayed tau-model. This is the same as option 1 except that the
          tau-model component takes the form :math:`t\,e^{−t/\\tau}`.

    :param wgp1: (default: 1)
        Integer specifying the optical depth in the Witt & Gordon (2000)
        models. Values range from 1 − 18, corresponding to optical depths of
        0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50,
        5.00, 5.50, 6.00, 7.00, 8.00, 9.00, 10.0. Note that these optical
        depths are defined differently from the optical depths defined by
        the parameters ``dust1`` and ``dust2``. See Witt & Gordon (2000)
        for details.

    :param wgp2: (default: 1)
        Integer specifying the type of large-scale geometry and extinction
        curve. Values range from 1-6, corresponding to MW+dusty, MW+shell,
        MW+cloudy, SMC+dusty, SMC+shell, SMC+cloudy. Dusty, shell, and cloudy
        specify the geometry and are described in Witt & Gordon (2000).

    :param wgp3: (default: 1)
        Integer specifying the local geometry for the Witt & Gordon (2000)
        dust models. A value of 0 corresponds to a homogeneous distribution,
        and a value of 1 corresponds to a clumpy distribution. See Witt &
        Gordon (2000) for details.

    :param evtype: (default: -1)
        Undocumented.

    """

    def __init__(self, compute_vega_mags=True, redshift_colors=False,
                 **kwargs):
        # Set up the parameters to their default values.
        self.params = ParameterSet(
            dust_type=0,
            imf_type=2,
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

        # Parse any input options.
        for k, v in self.params.iteritems():
            self.params[k] = kwargs.pop(k, v)

        # Make sure that we didn't get any unknown options.
        if len(kwargs):
            raise TypeError("__init__() got an unexpected keyword argument "
                            "'{}'".format(kwargs.keys()[0]))

        # Before the first time we interact with the FSPS driver, we need to
        # run the ``setup`` method.
        if not driver.is_setup:
            driver.setup(compute_vega_mags, redshift_colors)

    def _update_params(self):
        if self.params.dirtiness == 2:
            driver.set_ssp_params(*[self.params[k]
                                    for k in self.params.ssp_params])
        if self.params.dirtiness >= 1:
            driver.set_csp_params(*[self.params[k]
                                    for k in self.params.csp_params])
        self.params.dirtiness = 0

    def _compute_csp(self):
        self._update_params()
        driver.compute()

    @property
    def spectrum(self):
        """
        A grid (in age) of the spectra for the current CSP. The shape of
        the resulting object is ``(NAGE, NSPEC)``.

        """
        if self.params.dirty:
            self._compute_csp()
        return driver.get_spec(NSPEC, NAGE)

    @property
    def wavelengths(self):
        """
        The wavelength scale for the computed spectra.

        """
        return LAMBDA_GRID


class ParameterSet(object):

    ssp_params = ["imf_type", "imf1", "imf2", "imf3", "vdmc", "mdave",
                  "dell", "delt", "sbss", "fbhb", "pagb"]

    csp_params = ["dust_type",
                  "zmet", "sfh", "wgp1", "wgp2", "wgp3", "evtype", "tau",
                  "const", "tage", "fburst", "tburst", "dust1", "dust2",
                  "logzsol", "zred", "pmetals", "dust_clumps", "frac_nodust",
                  "dust_index", "dust_tesc", "frac_obrun", "uvb", "mwr",
                  "redgb", "dust1_index", "sf_start", "sf_trunc", "sf_theta",
                  "duste_gamma", "duste_umin", "duste_qpah", "fcstar",
                  "masscut"]

    @property
    def all_params(self):
        return self.ssp_params + self.csp_params

    @property
    def dirty(self):
        return self.dirtiness > 0

    def __init__(self, **kwargs):
        self.dirtiness = 2
        self._params = kwargs
        self.iteritems = self._params.iteritems

    def __getitem__(self, k):
        return self._params[k]

    def __setitem__(self, k, v):
        if k in self.ssp_params:
            self.dirtiness = 2
        elif k in self.csp_params:
            self.dirtiness = 1
        else:
            raise KeyError("Unrecognized parameter {}".format(k))

        self._params[k] = v
