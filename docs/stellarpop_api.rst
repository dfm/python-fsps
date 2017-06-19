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
   ('pdva', 'miles')

The last line indicates that we are using the Padova isochrones and MILES spectral library.
These can be changed only within FSPS itself (which must recompiled and python-FSPS reinstalled for the changes to take effect.)

Let's get the AB magnitudes in SDSS bands, for an SSP that is 13.7 Gyr old::

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

It is highly recommended that only one instance of :class:`fsps.StellarPopulation` be used in a given program.

API Reference
-------------

.. autoclass:: fsps.StellarPopulation
   :inherited-members:
