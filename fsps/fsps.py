#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["StellarPopulation", "find_filter"]

import numpy as np
from ._fsps import driver


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

    :param smooth_velocity: (default: True)
        Switch to choose smoothing in velocity space (``True``) or
        wavelength space.

    :param add_stellar_remanants: (default: True)
        Switch to add stellar remnnants in the stellar mass
        computation.

    :param add_neb_emission: (default: False)
        Switch to turn on/off a Cloudy-based nebular emission
        model. NB: this feature is currently under development, do not
        use!

    :param add_dust_emission: (default: True)
        Switch to turn on/off the Draine & Li 2007 dust emission
        model.

    :param add_agb_dust_model: (default: False)
        Switch to turn on/off the AGB circumstellar dust model. NB:
        this feature is currently under development, do not use!  If
        you do use it, note that the AGB dust emission is scaled by
        the parameter `agb_dust`, which defaults to 0!
   
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
        Weighting of red  giant branch.

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
        Truncate the IMF above this value.

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
        Compute SSPs for only the given evolutionary type

    :param sigma_smooth: (default: 0.0)
        If smooth_velocity is True, this gives the velocity dispersion in
        km/s.  Otherwise, it gives the width of the gaussian wavelength
        smoothing in Angstroms.

    :param agb_dust: (default: 1.0)
        Scales the circumstellar AGB dust emission.

    :param min_wave_smooth: (default: 1e3)
        Undocumented.

    :param max_wave_smooth: (default: 1e4)
        Undocumented.

    """

    def __init__(self, compute_vega_mags=True, redshift_colors=False,
                 smooth_velocity=True, add_stellar_remnants =True,
                 add_dust_emission = True, add_agb_dust_model = False,
                 add_neb_emission = False, **kwargs):
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
            evtype=-1,
            sigma_smooth=0.0,
            agb_dust=1.0,
            min_wave_smooth=1e3,
            max_wave_smooth=1e4,
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
            driver.setup(compute_vega_mags, redshift_colors, smooth_velocity,
                         add_stellar_remnants, add_neb_emission,
                         add_dust_emission, add_agb_dust_model)

        else:
            cvms, rcolors, svel, asr, ane, ade, agbd = driver.get_setup_vars()
            assert compute_vega_mags == bool(cvms)
            assert redshift_colors == bool(rcolors)
            assert smooth_velocity == bool(svel)
            assert add_stellar_remnants == bool(asr)
            assert add_neb_emission == bool(ane)
            assert add_dust_emission == bool(ade)
            assert add_agb_dust_model == bool(agbd)

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
        driver.compute()
        self._stats = None

    def get_spectrum(self, zmet=None, tage=0.0, peraa=False):
        """
        A grid (in age) of the spectra for the current CSP.

        :param zmet: (default: None)
            The (integer) index of the metallicity to use. By default, use
            the current value of ``self.params["zmet"]``.

        :param tage: (default: 0.0)
            The age of the stellar population. By default, this will compute
            a grid of ages from :math:`t \approx 0` to the maximum age in the
            isochrones.

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
        NTFULL = driver.get_ntfull()
        if tage > 0.0:
            return wavegrid, driver.get_spec(NSPEC, NTFULL)[0] * factor

        return wavegrid, driver.get_spec(NSPEC, NTFULL) * factor[None, :]

    @property
    def wavelengths(self):
        """
        The wavelength scale for the computed spectra.

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
            Otherwise, it is a 2D array with shape ``(NTFULL, NBANDS)``. If
            a particular band was requested then this return value will be
            properly compressed along that axis.

        """
        self.params["tage"] = tage
        if zmet is not None:
            self.params["zmet"] = zmet

        if self.params.dirty:
            self._compute_csp()

        NTFULL = driver.get_ntfull()
        NBANDS = driver.get_nbands()

        band_array = np.ones(NBANDS, dtype=bool)
        if bands is not None:
            inds = [FILTERS[band.lower()].index for band in bands]
            band_array[np.array([i not in inds for i in range(NBANDS)],
                                dtype=bool)] = False

        inds = np.array(band_array, dtype=int)
        mags = driver.get_mags(NTFULL, redshift, inds)

        if tage > 0.0:
            return mags[0, band_array]
        return mags[:, band_array]

    
    def ztinterp(self, zpos, tpos, peraa = False):
        """Return an SSP spectrum, mass, and luminosity interpolated
        to a target metallicity and age.  This effectively wraps the
        ZTINTERP subroutine.  Only the SSPs bracketing a given
        metallicity will be regenerated, if parameters are dirty.

        :param zpos:
            The metallicity, in units of log(Z/Z_sun)
            
        :param tpos:
            The desired age, in Gyr.

        :param peraa: (default: False)
            If true, return spectra in units of L_sun/AA, otherwise
            L_sun/Hz
            
        :returns spec:
            The SSP spectrum, interpolated to zpos and tpos.
            
        :returns mass:
            The stellar mass of the SSP at tpos.
            
        :returns lbol:
            The bolometric luminosity of the returned SSP.
        """
        
        if self.params.dirty:
            self._update_params()

        NSPEC = driver.get_nspec()
        spec, lbol, mass = np.zeros(NSPEC), np.zeros(1), np.zeros(1)
        logt_yrs = np.log10(tpos * 1e9)
        driver.interp_ssp(zpos,logt_yrs,spec,mass,lbol)

        if peraa:
            wavegrid = self.wavelengths
            factor = 3e18 / wavegrid ** 2
            spec *= factor
        
        return spec, mass, lbol
    
    def all_ssp_spec(self, update = True, peraa = False):
        """Return the contents of the ssp_spec_zz array.

        :param update: (default: True)
            If True, forces an update of the SSPs if the ssp
            parameters have changed. Otherwise simply dumps the
            current contents of the ssp_spec_zz array.
            
        :param peraa: (default: False)
            If true, return spectra in units of L_sun/AA, otherwise
            L_sun/Hz

        :returns specarray:
            The spectra of the SSPs, having shape (nspec, ntfull, nz).
            
        """

        if self.params.dirty and update:
            self._update_params()

        NSPEC = driver.get_nspec()
        NTFULL = driver.get_ntfull()
        NZ = driver.get_nz()
        specarray = driver.get_ssp_spec(NSPEC, NTFULL, NZ)
                
        if peraa:
            wavegrid = self.wavelengths
            factor = 3e18 / wavegrid ** 2
            specarray *= factor[:,None,None]
            
        return specarray
    
    def smoothspec(self, wave, spec, sigma, minw = None, maxw = None):
        """Smooth a spectrum by a gaussian with standard deviation
        given by sigma.  Whether the smoothing is in velocity space or
        in wavelength space depends on the value of the value of
        smooth_velocity.

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
        ns = len(wave)
        if maxw is None:
            maxw = np.max(wave)
        if minw is None:
            minw = np.min(wave)
        assert len(wave) == len(spec)
        outspec = np.array(spec)
        driver.smooth_spectrum(wave, outspec, sigma, minw, maxw)
        return outspec

    @property
    def log_age(self):
        """log10(age/yr)."""
        return self._stat(0)

    @property
    def stellar_mass(self):
        """Stellar mass (including remants if the FSPS parameters
        `add_stellar_remants=1`) in solar masses.
        """
        return self._stat(1)

    @property
    def log_lbol(self):
        """log(bolometric luminosity / L_solar)."""
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
                  "dell", "delt", "sbss", "fbhb", "pagb", "agb_dust",
                  "redgb", "masscut", "fcstar", "evtype"]

    csp_params = ["dust_type", "zmet", "sfh", "wgp1", "wgp2", "wgp3",
                  "tau", "const", "tage", "fburst", "tburst",
                  "dust1", "dust2", "logzsol", "zred", "pmetals",
                  "dust_clumps", "frac_nodust", "dust_index", "dust_tesc",
                  "frac_obrun", "uvb", "mwr", "dust1_index",
                  "sf_start", "sf_trunc", "sf_theta", "duste_gamma",
                  "duste_umin", "duste_qpah","sigma_smooth",
                  "min_wave_smooth", "max_wave_smooth"]

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

    def check_params(self):
        NZ = driver.get_nz()
        assert self._params["zmet"] in range(1, NZ + 1), \
            "zmet={0} out of range [1, {1}]".format(self._params["zmet"], NZ)
        assert self._params["dust_type"] in range(4), \
            "dust_type={0} out of range [0, 3]".format(
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


class Filter(object):

    def __init__(self, index, name, fullname):
        self.index = index - 1
        self.name = name.lower()
        self.fullname = fullname

    def __str__(self):
        return "<Filter({0})>".format(self.name)

    def __repr__(self):
        return "<Filter({0})>".format(self.name)


FILTERS = [(1, "V", "Johnson V (from Bessell 1990 via M. Blanton) - this "
            "defines V=0 for the Vega system"),
           (2, "U", "Johnson U (from Bessell 1990 via M. Blanton)"),
           (3, "B", "Johnson B (from Bessell 1990 via M. Blanton)"),
           (4, "Buser_B2", "Johnson B (from BC03. This is the B2 filter from "
            "Buser)"),
           (5, "Cousins_R", "Cousins R (from Bessell 1990 via M. Blanton)"),
           (6, "Cousins_I", "Cousins I (from Bessell 1990 via M. Blanton)"),
           (7, "CFHT_B", "CFHT B-band (from Blanton's kcorrect)"),
           (8, "CFHT_R", "CFHT R-band (from Blanton's kcorrect)"),
           (9, "CFHT_I", "CFHT I-band (from Blanton's kcorrect)"),
           (10, "2MASS_J", "2MASS J filter (total response w/atm)"),
           (11, "2MASS_H", "2MASS H filter (total response w/atm))"),
           (12, "2MASS_Ks", "2MASS Ks filter (total response w/atm)"),
           (13, "SDSS_u", "SDSS Camera u Response Function, airmass = 1.3 "
                "(June 2001)"),
           (14, "SDSS_g", "SDSS Camera g Response Function, airmass = 1.3 "
                "(June 2001)"),
           (15, "SDSS_r", "SDSS Camera r Response Function, airmass = 1.3 "
                "(June 2001)"),
           (16, "SDSS_i", "SDSS Camera i Response Function, airmass = 1.3 "
                "(June 2001)"),
           (17, "SDSS_z", "SDSS Camera z Response Function, airmass = 1.3 "
                "(June 2001)"),
           (18, "WFPC2_F255W", "HST WFPC2 F255W "
               "(http://acs.pha.jhu.edu/instrument/photometry/)"),
           (19, "WFPC2_F300W", "HST WFPC2 F300W "
               "(http://acs.pha.jhu.edu/instrument/photometry/)"),
           (20, "WFPC2_F336W", "HST WFPC2 F336W "
               "(http://acs.pha.jhu.edu/instrument/photometry/)"),
           (21, "WFPC2_F439W", "HST WFPC2 F439W "
               "(http://acs.pha.jhu.edu/instrument/photometry/)"),
           (22, "WFPC2_F450W", "HST WFPC2 F450W "
               "(http://acs.pha.jhu.edu/instrument/photometry/)"),
           (23, "WFPC2_F555W", "HST WFPC2 F555W "
               "(http://acs.pha.jhu.edu/instrument/photometry/)"),
           (24, "WFPC2_F606W", "HST WFPC2 F606W "
               "(http://acs.pha.jhu.edu/instrument/photometry/)"),
           (25, "WFPC2_F814W", "HST WFPC2 F814W "
               "(http://acs.pha.jhu.edu/instrument/photometry/)"),
           (26, "WFPC2_F850LP", "HST WFPC2 F850LP "
               "(http://acs.pha.jhu.edu/instrument/photometry/)"),
           (27, "WFC_ACS_F435W", "HST WFC ACS F435W "
               "(http://acs.pha.jhu.edu/instrument/photometry/)"),
           (28, "WFC_ACS_F475W", "HST ACS F475W "
               "(http://acs.pha.jhu.edu/instrument/photometry/)"),
           (29, "WFC_ACS_F555W", "HST ACS F555W "
                "(http://acs.pha.jhu.edu/instrument/photometry/"),
           (30, "WFC_ACS_F606W", "HST WFC ACS F606W "
                "(http://acs.pha.jhu.edu/instrument/photometry/)"),
           (31, "WFC_ACS_F625W", "HST ACS F625W "
                "(http://acs.pha.jhu.edu/instrument/photometry/)"),
           (32, "WFC_ACS_F775W", "HST WFC ACS F775W "
                "(http://acs.pha.jhu.edu/instrument/photometry/)"),
           (33, "WFC_ACS_F814W", "HST WFC ACS F814W "
                "(http://acs.pha.jhu.edu/instrument/photometry/)"),
           (34, "WFC_ACS_F850LP", "WFC ACS F850LP "
                "(http://acs.pha.jhu.edu/instrument/photometry/)"),
           (35, "WFC3_UVIS_F218W", "HST WFC3 UVIS F218W "
                "(http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/) Chip #1"),
           (36, "WFC3_UVIS_F225W", "HST WFC3 UVIS F225W "
                "(http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/) Chip #1"),
           (37, "WFC3_UVIS_F275W", "HST WFC3 UVIS F275W "
                "(http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/) Chip #1"),
           (38, "WFC3_UVIS_F336W", "HST WFC3 UVIS F336W "
                "(http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/) Chip #1"),
           (39, "WFC3_UVIS_F390W", "HST WFC3 UVIS F390W "
                "(http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/) Chip #1"),
           (40, "WFC3_UVIS_F438W", "HST WFC3 UVIS F438W "
                "(http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/) Chip #1"),
           (41, "WFC3_UVIS_F475W", "HST WFC3 UVIS F475W "
                "(http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/) Chip #1"),
           (42, "WFC3_UVIS_F555W", "HST WFC3 UVIS F555W "
                "(http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/) Chip #1"),
           (43, "WFC3_UVIS_F606W", "HST WFC3 UVIS F606W "
                "(http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/) Chip #1"),
           (44, "WFC3_UVIS_F775W", "HST WFC3 UVIS F775W "
                "(http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/) Chip #1"),
           (45, "WFC3_UVIS_F814W", "HST WFC3 UVIS F814W "
                "(http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/) Chip #1"),
           (46, "WFC3_UVIS_F850LP", "HST WFC3 UVIS F850LP "
                "(http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/) Chip #1"),
           (47, "WFC3_IR_F098M", "HST WFC3 IR F098M "
                "(http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/)"),
           (48, "WFC3_IR_F105W", "HST WFC3 IR F105W "
                "(http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/)"),
           (49, "WFC3_IR_F110W", "HST WFC3 IR F110W "
                "(http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/)"),
           (50, "WFC3_IR_F125W", "HST WFC3 IR F125W "
                "(http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/)"),
           (51, "WFC3_IR_F140W", "HST WFC3 IR F140W "
                "(http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/)"),
           (52, "WFC3_IR_F160W", "HST WFC3 IR F160W "
                "(http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/)"),
           (53, "IRAC_1", "Spitzer IRAC Channel 1 (3.6um)"),
           (54, "IRAC_2", "Spitzer IRAC Channel 2 (4.5um)"),
           (55, "IRAC_3", "Spitzer IRAC Channel 3 (5.8um)"),
           (56, "IRAC_4", "Spitzer IRAC Channel 4 (8.0um)"),
           (57, "ISAAC_Ks", "ISAAC Ks"),
           (58, "FORS_V", "FORS V"),
           (59, "FORS_R", "FORS R"),
           (60, "NICMOS_F110W", "HST NICMOS F110W"),
           (61, "NICMOS_F160W", "HST NICMOS F160W"),
           (62, "GALEX_FUV", "GALEX FUV"),
           (63, "GALEX_NUV", "GALEX NUV"),
           (64, "DES_g", "DES g (from Huan Lin, for DES camera)"),
           (65, "DES_r", "DES r (from Huan Lin, for DES camera)"),
           (66, "DES_i", "DES i (from Huan Lin, for DES camera)"),
           (67, "DES_z", "DES z (from Huan Lin, for DES camera)"),
           (68, "DES_Y", "DES Y (from Huan Lin, for DES camera)"),
           (69, "WFCAM_Z", "WFCAM Z (from Hewett et al. 2006, via A. Smith)"),
           (70, "WFCAM_Y", "WFCAM Y (from Hewett et al. 2006, via A. Smith)"),
           (71, "WFCAM_J", "WFCAM J (from Hewett et al. 2006, via A. Smith)"),
           (72, "WFCAM_H", "WFCAM H (from Hewett et al. 2006, via A. Smith)"),
           (73, "WFCAM_K", "WFCAM K (from Hewett et al. 2006, via A. Smith)"),
           (74, "Steidel_Un", "Steidel Un (via A. Shapley; see Steidel et al. "
                "2003)"),
           (75, "Steidel_G", "Steidel G  (via A. Shapley; see Steidel et al. "
                "2003)"),
           (76, "Steidel_Rs", "Steidel Rs (via A. Shapley; see Steidel et al. "
                "2003)"),
           (77, "Steidel_I", "Steidel I  (via A. Shapley; see Steidel et al. "
                "2003)"),
           (78, "MegaCam_u", "CFHT MegaCam u* "
                "(http://cadcwww.dao.nrc.ca/megapipe/docs/filters.html, "
                "Dec 2010)"),
           (79, "MegaCam_g", "CFHT MegaCam g' "
                "(http://cadcwww.dao.nrc.ca/megapipe/docs/filters.html)"),
           (80, "MegaCam_r", "CFHT MegaCam r' "
                "(http://cadcwww.dao.nrc.ca/megapipe/docs/filters.html)"),
           (81, "MegaCam_i", "CFHT MegaCam i' "
                "(http://cadcwww.dao.nrc.ca/megapipe/docs/filters.html)"),
           (82, "MegaCam_z", "CFHT MegaCam z' "
                "(http://cadcwww.dao.nrc.ca/megapipe/docs/filters.html)"),
           (83, "WISE_W1", "3.4um WISE W1 "
                "(http://www.astro.ucla.edu/~wright/WISE/passbands.html)"),
           (84, "WISE_W2", "4.6um WISE W2 "
                "(http://www.astro.ucla.edu/~wright/WISE/passbands.html)"),
           (85, "WISE_W3", "12um WISE W3 "
                "(http://www.astro.ucla.edu/~wright/WISE/passbands.html)"),
           (86, "WISE_W4", "22um WISE W4 22um "
                "(http://www.astro.ucla.edu/~wright/WISE/passbands.html)"),
           (87, "UVOT_W2", "UVOT W2 (from Erik Hoversten, 2011)"),
           (88, "UVOT_M2", "UVOT M2 (from Erik Hoversten, 2011)"),
           (89, "UVOT_W1", "UVOT W1 (from Erik Hoversten, 2011)"),
           (90, "MIPS_24", "Spitzer MIPS 24um"),
           (91, "MIPS_70", "Spitzer MIPS 70um"),
           (92, "MIPS_160", "Spitzer MIPS 160um"),
           (93, "SCUBA_450WB", "JCMT SCUBA 450WB "
            "(www.jach.hawaii.edu/JCMT/continuum/background/background.html)"),
           (94, "SCUBA_850WB", "JCMT SCUBA 850WB"),
           (95, "PACS_70", "Herschel PACS 70um"),
           (96, "PACS_100", "Herschel PACS 100um"),
           (97, "PACS_160", "Herschel PACS 160um"),
           (98, "SPIRE_250", "Herschel SPIRE 250um"),
           (99, "SPIRE_350", "Herschel SPIRE 350um"),
           (100, "SPIRE_500", "Herschel SPIRE 500um"),
           (101, "IRAS_12", "IRAS 12um"),
           (102, "IRAS_25", "IRAS 25um"),
           (103, "IRAS_60", "IRAS 60um"),
           (104, "IRAS_100", "IRAS 100um"),
           (105, "Bessell_L", "Bessell & Brett (1988) L band"),
           (106, "Bessell_LP", "Bessell & Brett (1988) L' band"),
           (107, "Bessell_M", "Bessell & Brett (1988) M band"),
           (108, "Stromgren_u", "Stromgren u (Bessell 2011)"),
           (109, "Stromgren_v", "Stromgren v (Bessell 2011)"),
           (110, "Stromgren_b", "Stromgren b (Bessell 2011)"),
           (111, "Stromgren_y", "Stromgren y (Bessell 2011)"),
           (112, "1500A", "Idealized 1500A bandpass with 15% bandwidth, "
            "FWHM = 225A from M. Dickinson"),
           (113, "2300A", "Idealized 2300A bandpass with 15% bandwidth, "
            "FWHM = 345A from M. Dickinson"),
           (114, "2800A", "Idealized 2800A bandpass with 15% bandwidth, "
            "FWHM = 420A from M. Dickinson"),
           (115, "JWST_F070W",
            "(http://www.stsci.edu/jwst/instruments/nircam/instrumentdesign/filters/)"),
           (116, "JWST_F090W",
            "(http://www.stsci.edu/jwst/instruments/nircam/instrumentdesign/filters/)"),
           (117, "JWST_F115W",
            "(http://www.stsci.edu/jwst/instruments/nircam/instrumentdesign/filters/)"),
           (118, "JWST_F150W",
            "(http://www.stsci.edu/jwst/instruments/nircam/instrumentdesign/filters/)"),
           (119, "JWST_F200W",
            "(http://www.stsci.edu/jwst/instruments/nircam/instrumentdesign/filters/)"),
           (120, "JWST_F277W",
            "(http://www.stsci.edu/jwst/instruments/nircam/instrumentdesign/filters/)"),
           (121, "JWST_F356W",
            "(http://www.stsci.edu/jwst/instruments/nircam/instrumentdesign/filters/)"),
           (122, "JWST_F444W",
            "(http://www.stsci.edu/jwst/instruments/nircam/instrumentdesign/filters/)")]

FILTERS = dict([(f[1].lower(), Filter(*f)) for f in FILTERS])


def find_filter(band):
    """
    Find the FSPS name for a filter.

    Usage:

    ::

        >>> import fsps
        >>> fsps.find_filter("F555W")
        ['wfpc2_f555w', 'wfc_acs_f555w']

    :param band:
        Something like the name of the band.

    """
    b = band.lower()
    possible = []
    for k in FILTERS.keys():
        if b in k:
            possible.append(k)
    return possible
