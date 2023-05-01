Modelling Stellar Populations
=============================

:class:`fsps.StellarPopulation` drives FSPS and is generally the only API needed.
:func:`fsps.find_filter` can also be used interactively to search for a filter (see also :doc:`filters`).

Example
-------

Lets start by initializing a simple stellar population
with solar metallicity
and some dust with a Calzetti et al (2000) extinction curve::

   >>> import fsps
   >>> sp = fsps.StellarPopulation(compute_vega_mags=False, zcontinuous=1,
                                   sfh=0, logzsol=0.0, dust_type=2, dust2=0.2)
   >>> sp.libraries
   ('mist', 'miles', 'DL07')

The last line indicates that we are using the MIST isochrones and MILES spectral
library. These can be changed only by reinstalling python-FSPS with appropriate
flags.

Let's get the AB magnitudes in SDSS bands, for a simple stellar population (SSP)
that is 13.7 Gyr old::

   >>> sdss_bands = fsps.find_filter('sdss')
   >>> print(sdss_bands)
   ['sdss_u', 'sdss_g', 'sdss_i', 'sdss_r', 'sdss_z']
   >>> sp.get_mags(tage=13.7, bands=sdss_bands)
   array([ 9.85484694,  7.8785663 ,  6.56580511,  6.99828574,  6.16779663])

Now we can change a parameter (say, lower the metallicity) and get a new set of magnitudes::

   >>> sp.params['logzsol'] = -1
   >>> sp.get_mags(tage=13.7, bands=sdss_bands)
   array([ 8.5626572 ,  7.07918435,  6.05304881,  6.38592117,  5.84199   ])

We can also get the spectrum, here in units of :math:`L_\odot/\mathrm{Hz}`,
as well as the total stellar mass formed by 13.7 Gyr
and the surviving stellar mass at 13.7 Gyr (both in units of :math:`M_\odot`)::

   >>> wave, spec = sp.get_spectrum(tage=13.7)
   >>> sp.formed_mass
   1.0
   >>> sp.stellar_mass
   0.57809244144339011

It is highly recommended that only one instance of :class:`fsps.StellarPopulation`
be used in a given program.

Example using nebular emission
------------------------------

We initialize a simple stellar population and set the flag to
include nebular emission::

   >>> sp = fsps.StellarPopulation(zcontinuous=1,
                                   add_neb_emission=1)

We can change the stellar metallicity, the gas phase metallicity, the
gas ionization parameter, and then return the total spectrum at 1 Myr::

   >>> sp.params['logzsol'] = -1.0
   >>> sp.params['gas_logz'] = -1.0
   >>> sp.params['gas_logu'] = -2.5
   >>> wave, spec = sp.get_spectrum(tage=0.001)

Note: for the nebular model to be fully self-consistent, the gas phase
metallicity and the stellar metallicity should be set to the same value.
This effectively adds the emission spectrum to the same stellar spectrum
that was used as the ionizing spectrum in Cloudy.
If users choose to vary the gas phase metallicity at constant stellar
metallicity, expect non-hydrogenic emission lines to be accurate within 1-15%.

Emission line wavelengths and line luminosities can be accessed through
the stellar population object::

   >>> sp.emline_wavelengths
   >>> sp.emline_luminosity


API Reference
-------------

.. autoclass:: fsps.StellarPopulation
   :inherited-members:
