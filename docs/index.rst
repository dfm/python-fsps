Python FSPS
===========

These are a set of Python binding's to `Charlie Conroy's Flexible Stellar
Population Synthesis (FSPS) Fortran library
<http://www.ucolick.org/~cconroy/FSPS.html>`_.


Installation
------------

To install, first you need to follow the directions on `the FSPS page
<http://www.ucolick.org/~cconroy/FSPS.html>`_ to checkout and compile
the FSPS package. In particular, these bindings rely on the value of the
``SPS_HOME`` environment variable being correctly set and the compiled
``.o`` and ``.mod`` files be available in the ``${SPS_HOME}/src``
directory.


License
-------

*Copyright 2013 Dan Foreman-Mackey.*

These bindings are available under the `MIT License
<https://raw.github.com/dfm/python-fsps/master/LICENSE.rst>`_. See the `FSPS
homepage <http://www.ucolick.org/~cconroy/FSPS.html>`_ for usage
restrictions and citation requirements.
