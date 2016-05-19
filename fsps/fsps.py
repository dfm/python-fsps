#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["StellarPopulation"]

import os
import numpy as np
from ._fsps import driver
from .filters import FILTERS


class StellarPopulation(object):
    """
    This is the main interface to use when interacting with FSPS from Python.
    Most of the Fortran API is exposed through Python hooks with various
    features added for user friendliness. When initializing, you can set any
    of the parameters of the system using keyword arguments. Below, you'll
    find a list of the options that you can include (with the comments taken
    directly from the `FSPS docs
    <https://github.com/cconroy20/fsps/blob/master/doc/MANUAL.pdf>`_. To change
    these values later, use the ``params`` property—which is ``dict``-like.
    For example:

    ::

        sp = StellarPopulation(imf_type=2)
        sp.params["imf_type"] = 1

    :param compute_vega_mags: (default: False)
        A switch that sets the zero points of the magnitude system: ``True``
        uses Vega magnitudes versus AB magnitudes.

    :param vactoair_flag: (default: False)
        If ``True``, output wavelengths in air (rather than vac).

    :param zcontinuous: (default: 0)
        Flag specifying how interpolation in metallicity is performed before
        computing composite models:

        * 0: No interpolation, use the metallicity index specified by ``zmet``.
        * 1: The SSPs are interpolated to the value of ``logzsol`` before the
          spectra and magnitudes are computed, and the value of ``zmet`` is
          ignored.
        * 2: The SSPs are convolved with a metallicity distribution function
          specified by the ``logzsol`` and ``pmetals`` parameters. The value of
          ``zmet`` is ignored.

    :param redshift_colors: (default: False)
        Flag specifying how to compute magnitudes:

        * ``False``: Magnitudes are computed at a fixed redshift specified by
          ``zred``.
        * ``True``: Magnitudes are computed at a redshift that corresponds to
          the age of the output SSP/CSP (assuming a redshift–age relation
          appropriate for a WMAP5 cosmology). This switch is useful if the user
          wants to compute the evolution in observed colors of a SSP/CSP.

    :param smooth_velocity: (default: True)
        Switch to choose smoothing in velocity space (``True``) or wavelength
        space.

    :param add_stellar_remnants: (default: True)
        Switch to add stellar remnants in the stellar mass computation.

    :param add_igm_absorption: (default: False)
        Switch to include IGM absorption via Madau (1995).  The ``zred``
        parameter must be non-zero for this switch to have any effect. The
        optical depth can be scaled using the ``igm_factor`` parameter.

    :param add_neb_emission: (default: False)
        Switch to turn on/off a Cloudy-based nebular emission model.

    :param add_neb_continuum: (default: True)
        Switch to turn on/off the nebular continuum component (automatically
        turned off if ``add_neb_emission`` is ``False``).

    :param cloudy_dust: (default: False)
        Switch to include dust in the Cloudy tables.

    :param add_dust_emission: (default: True)
        Switch to turn on/off the Draine & Li 2007 dust emission model.

    :param add_agb_dust_model: (default: True)
        Switch to turn on/off the AGB circumstellar dust model. NB: The AGB
        dust emission is scaled by the parameter `agb_dust`.

    :param tpagb_norm_type: (default: 2)
        Flag specifying TP-AGB normalization scheme:

        * 0: default Padova 2007 isochrones
        * 1: Conroy & Gunn 2010 normalization
        * 2: Villaume, Conroy, Johnson 2015 normalization

    :param dust_type: (default: 0)
        Common variable deﬁning the extinction curve for dust around old stars:

        * 0: power law with index dust index set by ``dust_index``.
        * 1: Milky Way extinction law (with the :math:`R = A_V /E(B - V)` value
          given by ``mwr``) parameterized by Cardelli et al. (1989).
        * 2: Calzetti et al. (2000) attenuation curve. Note that if this value
          is set then the dust attenuation is applied to all starlight equally
          (not split by age), and therefore the only relevant parameter is
          ``dust2``, which sets the overall normalization.
        * 3: allows the user to access a variety of attenuation curve models
          from Witt & Gordon (2000) using the parameters ``wgp1`` and ``wgp2``.
        * 4: Kriek & Conroy 2013 attenuation curve.  In this model the slope of
          the curve, set by the parameter ``dust_index``, is linked to the
          strength of the UV bump.

    :param imf_type: (default: 2)
        Common variable defining the IMF type:

        * 0: Salpeter (1955)
        * 1: Chabrier (2003)
        * 2: Kroupa (2001)
        * 3: van Dokkum (2008)
        * 4: Dave (2008)
        * 5: tabulated piece-wise power law IMF, specified in ``imf.dat`` file
          located in the data directory

    :param pagb: (default: 1.0)
        Weight given to the post–AGB phase. A value of 0.0 turns off post-AGB
        stars; a value of 1.0 implies that the Vassiliadis & Wood (1994) tracks
        are implemented as–is.

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
        are uniformly spread in :math:`\log T_\mathrm{eff}` to :math:`10^4`
        K. See Conroy et al. (2009a) for details and a plausible range.

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
        ``sfh=4``.

    :param tage: (default: 0.0)
        If set to a non-zero value, the
        :func:`fsps.StellarPopulation.compute_csp` method will compute the
        spectra and magnitudes only at this age, and will therefore only output
        one age result. The units are Gyr. (The default is to compute and
        return results from :math:`t \\approx 0` to the maximum age in the
        isochrones).

    :param fburst: (default: 0.0)
        Deﬁnes the fraction of mass formed in an instantaneous burst of star
        formation. Only used if ``sfh=1`` or ``sfh=4``.

    :param tburst: (default: 11.0)
        Defines the age of the Universe when the burst occurs. If ``tburst >
        tage`` then there is no burst. Only used if ``sfh=1`` or ``sfh=4``.

    :param dust1: (default: 0.0)
        Dust parameter describing the attenuation of young stellar light,
        i.e. where ``t <= dust_tesc`` (for details, see Conroy et al. 2009a).

    :param dust2: (default: 0.0)
        Dust parameter describing the attenuation of old stellar light,
        i.e. where ``t > dust_tesc`` (for details, see Conroy et al. 2009a).

    :param logzsol: (default: 0.0)
        Parameter describing the metallicity, given in units of :math:`\log
        (Z/Z_\odot)`.  Only used if ``zcontinuous > 0``.

    :param zred: (default: 0.0)
        Redshift. If this value is non-zero and if ``redshift_colors=1``, the
        magnitudes will be computed for the spectrum placed at redshift
        ``zred``.

    :param pmetals: (default: 2.0)
       The power for the metallicty distribution function.  The MDF is given by
       :math:`(Z \\, e^{{-Z}})^{{\mathrm{pmetals}}}` where :math:`Z =
       z/(z_\\odot \\, 10^{{\mathrm{logzsol}}})` and z is the metallicity in
       linear units (i.e., :math:`z_\odot = 0.019`).  Using a negative value
       will result in smoothing of the SSPs by a three-point triangular kernel
       before linear interpolation (in :math:`\log Z`) to the requested
       metallicity.

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
        Fraction of the young stars (age < dust_tesc) that are not attenuated
        by ``dust1``, representing runaway OB stars.  These stars are still
        attenuated by ``dust2``.

    :param uvb: (default: 1.0)
        Parameter characterizing the strength of the 2175A extinction feature
        with respect to the standard Cardelli et al. determination for the
        MW. Only used when ``dust_type=1``.

    :param mwr: (default: 3.1)
        The ratio of total to selective absorption which characterizes the MW
        extinction curve: :math:`R = A_V /E(B - V)`. Only used when
        ``dust_type=1``.

    :param redgb: (default: 1.0)
        Weighting of red  giant branch.

    :param dust1_index: (default: -1.0)
        Power law index of the attenuation curve affecting stars younger than
        dust_tesc corresponding to ``dust1``. Only used when ``dust_type=0``.

    :param mdave: (default: 0.5)
        IMF parameter defined in Dave (2008). Only used if ``imf_type=4``.

    :param sf_start: (default: 0.0)
        Start time of the SFH, in Gyr.

    :param sf_trunc: (default: 0.0)
        Undocumented.

    :param sf_slope: (default: 0.0)
        Undocumented.

    :param duste_gamma: (default: 0.01)
        Parameter of the Draine & Li (2007) dust emission model. Specifies the
        relative contribution of dust heated at a radiation field strength of
        :math:`U_\mathrm{min}` and dust heated at :math:`U_\mathrm{min} < U \le
        U_\mathrm{max}`. Allowable range is 0.0 – 1.0.

    :param duste_umin: (default: 1.0)
        Parameter of the Draine & Li (2007) dust emission model. Specifies the
        minimum radiation field strength in units of the MW value. Valid range
        is 0.1 – 25.0.

    :param duste_qpah: (default: 3.5)
        Parameter of the Draine & Li (2007) dust emission model. Specifies the
        grain size distribution through the fraction of grain mass in
        PAHs. This parameter has units of % and a valid range of 0.0 − 10.0.

    :param fcstar: (default: 1.0)
        Fraction of stars that the Padova isochrones identify as Carbon stars
        that FSPS assigns to a Carbon star spectrum. Set this to 0.0 if for
        example the users wishes to turn all Carbon stars into regular M-type
        stars.

    :param masscut: (default: 150.0)
        Truncate the IMF above this value.

    :param zmet: (default: 1)
        The metallicity is specified as an integer ranging between 1 and nz. If
        ``zcontinuous > 0`` then this parameter is ignored.

    :param sfh: (default: 0)
        Defines the type of star formation history, normalized such that one
        solar mass of stars is formed over the full SFH. Default value is 0.

        * 0: Compute an SSP
        * 1: Compute a five parameter SFH (see below).
        * 2: Compute a tabulated SFH defined in a file called ``sfh.dat`` that
          must reside in the data directory. The file must contain three
          rows. The first column is time since the Big Bang in Gyr, the second
          is the SFR in units of solar masses per year, the third is the
          absolute metallicity. An example is provided in the data
          directory. The time grid in this file can be arbitrary (so long as
          the units are correct), but it is up to the user to ensure that the
          tabulated SFH is well-sampled so that the outputs are stable.
          Obviously, highly oscillatory data require dense sampling.  This
          option is not supported in Python-FSPS.
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
        Compute SSPs for only the given evolutionary type.

    :param sigma_smooth: (default: 0.0)
        If smooth_velocity is True, this gives the velocity dispersion in km/s.
        Otherwise, it gives the width of the gaussian wavelength smoothing in
        Angstroms.  These widths are in terms of :math:`\sigma`, *not* FWHM.

    :param agb_dust: (default: 1.0)
        Scales the circumstellar AGB dust emission.

    :param min_wave_smooth: (default: 1e3)
        Minimum wavelength to consider when smoothing the spectrum.

    :param max_wave_smooth: (default: 1e4)
        Maximum wavelength to consider when smoothing the spectrum.

    :param gas_logu: (default: -2)
        Log of the gas ionization parameter, for determining the nebular
        emission.

    :param gas_logz: (default: 0.0)
        Log of the gas-phase metallicity, for determining the nebular emission.
        In units of :math:`\log (Z/Z_\odot)`.

    :param igm_factor: (default: 1.0)
        Factor used to scale the IGM optical depth.

    :param fagn: (default: 0.0)
        The total luminosity of the AGN, expressed as a fraction of the
        bolometric stellar luminosity (so it can be greater than 1). The shape
        of the AGN SED is from the Nenkova et al. 2008 templates.

    :param agn_tau: (default: 10)
        Optical depth of the AGN dust torus, which affects the shape of the AGN
        SED.  Outside the range (5, 150) the AGN SED is an
        extrapolation.
    """

    def __init__(self, compute_vega_mags=False, zcontinuous=0,
                 **kwargs):

        # Set up the parameters to their default values.
        self.params = ParameterSet(
            smooth_velocity=True,
            vactoair_flag=False,
            redshift_colors=False,
            dust_type=0,
            add_dust_emission=True,
            add_agb_dust_model=True,
            add_neb_emission=False,
            add_neb_continuum=True,
            cloudy_dust=False,
            add_igm_absorption=False,
            add_stellar_remnants=True,
            imf_type=2,
            tpagb_norm_type=2,
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
            logzsol=0.0,
            zred=0.0,
            pmetals=2.0,
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
            sf_slope=0.0,
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
            evtype=-1,
            sigma_smooth=0.0,
            agb_dust=1.0,
            min_wave_smooth=1e3,
            max_wave_smooth=1e4,
            gas_logu=-2,
            gas_logz=0.0,
            igm_factor=1.0,
            fagn=0.0,
            agn_tau=10.0
        )

        # Parse any input options.
        for k, v in self.params.iteritems():
            self.params[k] = kwargs.pop(k, v)

        # Make sure that we didn't get any unknown options.
        if len(kwargs):
            raise TypeError("__init__() got an unexpected keyword argument "
                            "'{0}'".format(kwargs.keys()[0]))

        # Before the first time we interact with the FSPS driver, we need to
        # run the ``setup`` method.
        if not driver.is_setup:
            driver.setup(compute_vega_mags)
        else:
            cvms = driver.get_setup_vars()
            assert compute_vega_mags == bool(cvms)

        self._zcontinuous = zcontinuous
        # Caching.
        self._wavelengths = None
        self._zlegend = None
        self._ssp_ages = None
        self._stats = None

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

        NSPEC = driver.get_nspec()
        NTFULL = driver.get_ntfull()
        driver.compute_zdep(NSPEC, NTFULL, self._zcontinuous)

        self._stats = None

    def get_spectrum(self, zmet=None, tage=0.0, peraa=False):
        """
        Return spectra for the current CSP.

        :param zmet: (default: None)
            The (integer) index of the metallicity to use. By default, use
            the current value of ``self.params["zmet"]``.

        :param tage: (default: 0.0)
            The age of the stellar population for which to obtain a spectrum. By
            default, this will compute a grid of ages from :math:`t \\approx 0`
            to the maximum age in the isochrones.

        :param peraa: (default: False)
            If ``True``, return the spectrum in :math:`L_\odot/\AA`.
            Otherwise, return the spectrum in the FSPS standard
            :math:`L_\odot/\mathrm{Hz}`.

        :returns wavelengths:
            The wavelength grid in Angstroms.

        :returns spectrum:
            The spectrum in :math:`L_\odot/\mathrm{Hz}` or :math:`L_\odot/\AA`.
            If an age was provided by the ``tage`` parameter then the result
            is a 1D array with ``NSPEC`` values. Otherwise, it is a 2D array
            with shape ``(NTFULL, NSPEC)``.
        """
        self.params["tage"] = tage
        if zmet is not None:
            self.params["zmet"] = zmet

        if self.params.dirty:
            self._compute_csp()

        wavegrid = self.wavelengths
        if peraa:
            factor = 3e18 / wavegrid ** 2

        else:
            factor = np.ones_like(wavegrid)

        NSPEC = driver.get_nspec()
        if tage > 0.0:
            return wavegrid, driver.get_spec(NSPEC, 1)[0] * factor

        NTFULL = driver.get_ntfull()
        return wavegrid, driver.get_spec(NSPEC, NTFULL) * factor[None, :]

    def get_mags(self, zmet=None, tage=0.0, redshift=0.0, bands=None):
        """
        Get the magnitude of the CSP.

        :param zmet: (default: None)
            The (integer) index of the metallicity to use. By default, use
            the current value of ``self.params["zmet"]``.

        :param tage: (default: 0.0)
            The age of the stellar population. By default, this will compute
            a grid of ages from :math:`t \approx 0` to the maximum age in the
            isochrones.

        :param redshift: (default: 0.0)
            Optionally redshift the spectrum first.

        :param bands: (default: None)
            The names of the filters that you would like to compute the
            magnitude for. This should correspond to the result of
            :func:`fsps.find_filter`.

        :returns mags:
            The magnitude grid. If an age was was provided by the ``tage``
            parameter then the result is a 1D array with ``NBANDS`` values.
            Otherwise, it is a 2D array with shape ``(NTFULL, NBANDS)``. If a
            particular set of bands was requested then this return value will
            be properly compressed along that axis, ordered according to the
            ``bands`` argument.
        """
        self.params["tage"] = tage
        if zmet is not None:
            self.params["zmet"] = zmet

        if self.params.dirty:
            self._compute_csp()

        if tage > 0.0:
            NTFULL = 1
        else:
            NTFULL = driver.get_ntfull()
        NBANDS = driver.get_nbands()
        NSPEC = driver.get_nspec()
        band_array = np.ones(NBANDS, dtype=bool)
        if bands is not None:
            user_sorted_inds = np.array([FILTERS[band.lower()].index
                                         for band in bands])
            band_array[np.array([i not in user_sorted_inds
                                 for i in range(NBANDS)],
                                dtype=bool)] = False

        inds = np.array(band_array, dtype=int)
        mags = driver.get_mags(NSPEC, NTFULL, redshift, inds)

        if tage > 0.0:
            if bands is not None:
                return mags[0, user_sorted_inds]
            else:
                return mags[0, :]
        else:
            if bands is not None:
                return mags[:, user_sorted_inds]
            else:
                return mags

    def ztinterp(self, zpos, tpos, peraa=False):
        """
        Return an SSP spectrum, mass, and luminosity interpolated to a target
        metallicity and age.  This effectively wraps the ZTINTERP subroutine.
        Only the SSPs bracketing a given metallicity will be regenerated, if
        parameters are dirty.

        :param zpos:
            The metallicity, in units of :math:`\log(Z/Z_\odot)`

        :param tpos:
            The desired age, in Gyr.

        :param peraa: (default: False)
            If true, return spectra in units of :math:`L_\odot/\AA`, otherwise
            :math:`L_\odot/\mathrm{Hz}`

        :returns spec:
            The SSP spectrum, interpolated to zpos and tpos.

        :returns mass:
            The stellar mass of the SSP at tpos.

        :returns lbol:
            The bolometric luminosity of the returned SSP.
        """
        if self.params.dirtiness == 2:
            self._update_params()

        NSPEC = driver.get_nspec()
        spec, mass, lbol = np.zeros(NSPEC), np.zeros(1), np.zeros(1)
        logt_yrs = np.log10(tpos * 1e9)
        driver.interp_ssp(zpos, logt_yrs, spec, mass, lbol)

        if peraa:
            wavegrid = self.wavelengths
            factor = 3e18 / wavegrid ** 2
            spec *= factor

        return spec, mass, lbol

    def _all_ssp_spec(self, update=True, peraa=False):
        """
        Return the contents of the ssp_spec_zz array.

        :param update: (default: True)
            If True, forces an update of the SSPs if the ssp parameters have
            changed. Otherwise simply dumps the current contents of the
            ``ssp_spec_zz`` array.

        :param peraa: (default: False)
            If true, return spectra in units of :math:`L_\odot/\AA`, otherwise
            :math:`L_\odot/\mathrm{Hz}`

        :returns spec:
            The spectra of the SSPs, having shape (nspec, ntfull, nz).

        :returns mass:
            The mass of the SSPs, having shape (ntfull, nz).

        :returns lbol:
            The bolometric luminosity of the SSPs, having shape (ntfull, nz).
        """

        if (self.params.dirtiness == 2) and update:
            self._update_params()

        NSPEC = driver.get_nspec()
        NTFULL = driver.get_ntfull()
        NZ = driver.get_nz()
        spec = np.zeros([NSPEC, NTFULL, NZ], order='F')
        mass = np.zeros([NTFULL, NZ], order='F')
        lbol = np.zeros([NTFULL, NZ], order='F')
        driver.get_ssp_spec(spec, mass, lbol)

        if peraa:
            wavegrid = self.wavelengths
            factor = 3e18 / wavegrid ** 2
            spec *= factor[:, None, None]

        return spec, mass, lbol

    def get_stellar_spectrum(self, mact, logt, lbol, logg, phase, comp,
                             mdot=0, weight=1, zmet=None, peraa=True):
        """
        Get the spectrum of a star with a given set of physical parameters.
        This uses the metallicity given by the current value of ``zmet``.

        :param mact:
            Actual stellar mass (after taking into account mass loss).  Used to
            calculate surface gravity.

        :param logt:
            The log of the effective temperature.

        :param lbol:
            Stellar luminosity, in units of :math:`L_\odot`

        :param logg:
            Log of the surface gravity g.  Note that this variable is actually
            ignored, and logg is calculated internally using ``mact``,
            ``lbol``, and ``logt``.

        :param phase:
            The evolutionary phase, 0 through 6.

        :param comp:
            Composition, in terms of C/O ratio.  Only used for AGB stars
            (phase=5), where the division between carbon and oxyygen rich stars
            is :math:`C/O = 1`.

        :param mdot:
            The log of the mass loss rate.

        :param weight:
            The IMF weight

        :returns outspec:
            The spectrum of the star, in :math:`L_\odot/\mathrm{Hz}`
        """
        if zmet is not None:
            self.params["zmet"] = zmet

        if self.params.dirty:
            self._update_params()

        NSPEC = driver.get_nspec()
        outspec = np.zeros(NSPEC)
        driver.stellar_spectrum(mact, logt, lbol, logg, phase,
                                comp, mdot, weight, outspec)
        if peraa:
            wavegrid = self.wavelengths
            factor = 3e18 / wavegrid ** 2
            outspec *= factor

        return outspec

    def isochrones(self, outfile='pyfsps_tmp'):
        """
        Write the isochrone data (age, mass, weights, phases, magnitudes, etc.)
        to a .cmd file, then read it into a huge numpy array.

        :param outfile: (default: 'pyfsps_tmp')
            The file root name of the .cmd file, which will be placed in the
            $SPS_HOME/OUTPUTS/ directory

        :returns dat:
            A huge numpy array containing information about every isochrone
            point for the current metallicity.

        :returns header:
            A list of the column names pulled from the header of the .cmd file.
            The number of column names will not match the number of columns in
            dat since there are as many ``mags`` as filters.  Use
            ``fsps.list_filters()`` to get the ordered filter name list.
        """
        if self.params.dirty:
            self._compute_csp()

        from . import ev, list_filters
        absfile = os.path.join(ev, 'OUTPUTS', outfile+'.cmd')
        driver.write_isoc(outfile)

        with open(absfile, 'r') as f:
            # drop the comment hash and mags field
            header = f.readline().split()[1:-1]
        header += list_filters()
        cmd_data = np.loadtxt(absfile, comments='#',
                              dtype=np.dtype([(n, np.float) for n in header]))
        return cmd_data

    def smoothspec(self, wave, spec, sigma, minw=None, maxw=None):
        """
        Smooth a spectrum by a gaussian with standard deviation given by sigma.
        Whether the smoothing is in velocity space or in wavelength space
        depends on the value of the value of smooth_velocity.

        :param wave:
            The input wavelength grid.

        :param spec:
            The input spectrum.

        :param sigma:
            The standard deviation of the gaussian broadening function.

        :param minw:
            Optionally set the minimum wavelength to consider when
            broadening.

        :param maxw:
            Optionally set the maximum wavelength to consider when
            broadening.

        :returns outspec:
            The smoothed spectrum, on the same wavelength grid as the input.
        """
        if maxw is None:
            maxw = np.max(wave)
        if minw is None:
            minw = np.min(wave)
        assert len(wave) == len(spec)
        outspec = np.array(spec)
        driver.smooth_spectrum(wave, outspec, sigma, minw, maxw)
        return outspec

    def filter_data(self):
        """Return effective wavelengths, and vega and solar magnitudes
        of all filters.
        """
        NBANDS = driver.get_nbands()
        lambda_eff, magvega, magsun = driver.get_filter_data(NBANDS)
        return lambda_eff, magvega, magsun

    @property
    def wavelengths(self):
        """
        The wavelength scale for the computed spectra, in :math:`\AA`

        """
        if self._wavelengths is None:
            NSPEC = driver.get_nspec()
            self._wavelengths = driver.get_lambda(NSPEC)
        return self._wavelengths

    @property
    def zlegend(self):
        """
        The available metallicities.

        """
        if self._zlegend is None:
            NZ = driver.get_nz()
            self._zlegend = driver.get_zlegend(NZ)
        return self._zlegend

    @property
    def ssp_ages(self):
        """
        The age grid of the SSPs, in log(years), used by FSPS.

        """
        if self._ssp_ages is None:
            NTFULL = driver.get_ntfull()
            self._ssp_ages = driver.get_timefull(NTFULL)
        return self._ssp_ages

    @property
    def log_age(self):
        """log10(age/yr)."""
        return self._stat(0)

    @property
    def stellar_mass(self):
        """Stellar mass in solar masses (including remants if the FSPS
        parameter `add_stellar_remants=1`).
        """
        return self._stat(1)

    @property
    def log_lbol(self):
        """log(bolometric luminosity / :math:`L_\odot/`)."""
        return self._stat(2)

    @property
    def sfr(self):
        """Star formation rate (solar masses per year)."""
        return self._stat(3)

    @property
    def dust_mass(self):
        """Dust mass, in solar masses."""
        return self._stat(4)

    def _get_grid_stats(self):
        if self.params.dirty:
            self._compute_csp()

        if self._stats is None:
            self._stats = driver.get_stats(driver.get_ntfull())

        return self._stats

    def _stat(self, k):
        stats = self._get_grid_stats()
        if self.params["tage"] > 0:
            return stats[k][0]
        return stats[k]


class ParameterSet(object):

    ssp_params = ["imf_type", "imf1", "imf2", "imf3", "vdmc", "mdave",
                  "dell", "delt", "sbss", "fbhb", "pagb",
                  "add_stellar_remnants", "tpagb_norm_type",
                  "add_agb_dust_model", "agb_dust", "redgb", "masscut",
                  "fcstar", "evtype"]

    csp_params = ["smooth_velocity", "vactoair_flag", "redshift_colors",
                  "dust_type", "add_dust_emission", "add_neb_emission",
                  "add_neb_continuum", "cloudy_dust", "add_igm_absorption",
                  "zmet", "sfh", "wgp1", "wgp2", "wgp3",
                  "tau", "const", "tage", "fburst", "tburst",
                  "dust1", "dust2", "logzsol", "zred", "pmetals",
                  "dust_clumps", "frac_nodust", "dust_index", "dust_tesc",
                  "frac_obrun", "uvb", "mwr", "dust1_index",
                  "sf_start", "sf_trunc", "sf_slope", "duste_gamma",
                  "duste_umin", "duste_qpah", "sigma_smooth",
                  "min_wave_smooth", "max_wave_smooth", "gas_logu",
                  "gas_logz", "igm_factor", "fagn", "agn_tau"]

    @property
    def all_params(self):
        return self.ssp_params + self.csp_params

    @property
    def dirty(self):
        return self.dirtiness > 0

    def __init__(self, **kwargs):
        self.dirtiness = 2
        self._params = kwargs
        try:
            self.iteritems = self._params.iteritems
        except AttributeError:
            self.iteritems = self._params.items

    def check_params(self):
        NZ = driver.get_nz()
        assert self._params["zmet"] in range(1, NZ + 1), \
            "zmet={0} out of range [1, {1}]".format(self._params["zmet"], NZ)
        assert self._params["dust_type"] in range(5), \
            "dust_type={0} out of range [0, 4]".format(
                self._params["dust_type"])
        assert self._params["imf_type"] in range(6), \
            "imf_type={0} out of range [0, 5]".format(self._params["imf_type"])

    def __getitem__(self, k):
        return self._params[k]

    def __setitem__(self, k, v):
        original = self._params[k]
        is_changed = original != v

        if is_changed:
            if k in self.ssp_params:
                self.dirtiness = 2
            elif k in self.csp_params:
                self.dirtiness = max(1, self.dirtiness)

            self._params[k] = v
            self.check_params()
