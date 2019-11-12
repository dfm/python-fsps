Working with Filters
====================

Python-FSPS provides simple access to FSPS's filter set.
To find what filters are available, use :func:`fsps.find_filter` or :func:`fsps.list_filters` (see also :doc:`filters`).
:func:`fsps.get_filter` is used to get a :class:`fsps.filters.Filter` instance, which in turn provides properties such as the filter's transmission curve, Solar absolute magnitude, and effective wavelength.

Example
-------

::

   >>> import fsps
   >>> fsps.find_filter('sdss')
   ['sdss_u', 'sdss_g', 'sdss_i', 'sdss_r', 'sdss_z']
   >>> g = fsps.get_filter('sdss_g')
   >>> g.msun_ab
   5.12
   >>> g.lambda_eff
   4702.5


API Reference
-------------

.. autofunction:: fsps.find_filter

.. autofunction:: fsps.list_filters

.. autofunction:: fsps.get_filter

.. autoclass:: fsps.filters.Filter
   :inherited-members:
