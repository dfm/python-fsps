Python FSPS
===========

These are a set of Python binding's to `Charlie Conroy's Flexible Stellar
Population Synthesis (FSPS) Fortran library
<http://www.ucolick.org/~cconroy/FSPS.html>`_. This code is very much inspired
by `Jonathan Sick's version <https://github.com/jonathansick/pySPS>`_ of the
same thing. The difference is that this version is lighter weight: providing
less baked functionality but trading that for fewer dependencies.


Installation
------------

There are only two required dependencies for building this package: FSPS
(obviously) and `NumPy <http://www.numpy.org/>`_. First, m sure that you have
NumPy installed (using your favorite method) and access to the ``f2py``
executable (this should be installed along with NumPy).
Then, you need to follow the directions on `the FSPS page
<http://www.ucolick.org/~cconroy/FSPS.html>`_ to checkout and compile
the FSPS package. In particular, these bindings rely on the value of the
``SPS_HOME`` environment variable being correctly set and the compiled
``.o`` and ``.mod`` files be available in the ``${SPS_HOME}/src``
directory. Finally, clone the `python-fsps repoistory
<https://github.com/dfm/python-fsps>`_: ``git clone
https://github.com/dfm/python-fsps.git`` and run

::

    cd python-fsps
    python setup.py build_fsps
    python setup.py develop

To test the installation, run (this doesn't actually work yet):

::

    python -c 'import fsps;fsps.test()'


API
---

.. autoclass:: fsps.StellarPopulation
   :inherited-members:

.. autofunction:: fsps.find_filter


License
-------

*Copyright 2013 Dan Foreman-Mackey.*

These bindings are available under the `MIT License
<https://raw.github.com/dfm/python-fsps/master/LICENSE.rst>`_. See the `FSPS
homepage <http://www.ucolick.org/~cconroy/FSPS.html>`_ for usage
restrictions and citation requirements.
