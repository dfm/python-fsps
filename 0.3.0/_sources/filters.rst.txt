Table of FSPS Filters
=====================

For convenience, we provide provide a table of filters supported by FSPS.
These filters can be referred to by providing their *Name* string, which is case-insensitive.
For example, to model an SED in HST/ACS F660W and F814W bandpasses::

    import fsps
    sp = fsps.StellarPopulation()
    sp.get_mags(bands=['WFC_ACS_F606W', 'WFC_ACS_F814W'])


.. include:: filter_table.rst
