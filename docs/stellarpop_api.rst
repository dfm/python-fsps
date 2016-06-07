Modelling Stellar Populations
=============================

:class:`fsps.StellarPopulation` drives FSPS and is generally the only API needed. 
:func:`fsps.find_filter` can also be used interactively to search for a filter (see also :doc:`filters`).

Example
-------

Lets start by initializing a simple stellar population, born with solar metallicity 6 Gyr after the Big Bang, and some dust with a Calzetti et al (2000) extinction curve::

   >>> import fsps
   >>> sp = fsps.StellarPopulation(compute_vega_mags=False, sfh=0, zmet=20,
                                   dust_type=2, dust2=0.2, sf_start=6.)

Lets get the AB magnitudes in SDSS bands, when the SSP has evolved to 13.7 Gyr after the Big Bang::

   >>> sdss_bands = fsps.find_filter('sdss')
   >>> print(sdss_bands)
   ['sdss_u', 'sdss_g', 'sdss_i', 'sdss_r', 'sdss_z']
   >>> sp.get_mags(tage=13.7, bands=sdss_bands)
   array([ 9.83712169,  7.89895578,  6.58329011,  6.9877158 ,  6.21844951])

Now we can change a parameter (say, lower the metallicity) and get a new set of magnitudes::

   >>> sp.params['zmet'] = 10
   >>> sp.get_mags(tage=13.7, bands=sdss_bands)
   array([ 8.56374292,  7.12518721,  6.14903857,  6.43916996,  5.9635293 ])

We can also get the spectrum, here in units of :math:`L_\odot/\mathrm{Hz}`::

   >>> wave, spec = sp.get_spectrum(tage=13.7)

 It is highly recemmended that only one instance of :class:`fsps.StellarPopulation` be used in a given program.

API Reference
-------------

.. autoclass:: fsps.StellarPopulation
   :inherited-members:
