Installation
============

Prerequisites
-------------

There are only two required dependencies for building this package: `FSPS
<https://github.com/cconroy20/fsps>`_ (obviously) and `NumPy
<http://www.numpy.org/>`_.
First, make sure that you have NumPy installed (using your favorite method)
and access to the ``f2py`` executable (this should be installed along with
NumPy).  Note that ``f2py`` versions packaged with ``numpy`` v1.10.x do
not work with python-FSPS, please use ``numpy`` v1.9.3 or lower.

Then, you need to follow the directions on `the FSPS page
<https://github.com/cconroy20/fsps>`_ to clone and compile the FSPS
package. Note that Python-FSPS is built against specific versions of
the FSPS Fortran API, so it is important that you have a recent tagged
release version of FSPS.  If you are missing required updates to FSPS
you will be notified of this when you try to build Python-FSPS.
Currently Python-FSPS is built against FSPS v2.6.

Once FSPS is checked out, modify the Makefile as needed and compile
FSPS.  For some compilers (e.g. Intel) it may be necessary to set the
``-fPIC`` flag when compiling FSPS -- please try this if you encounter
Python-FSPS build errors.

These bindings rely on the value of the ``SPS_HOME`` environment
variable being correctly set and the compiled ``.o`` and ``.mod``
files be available in the ``${SPS_HOME}/src`` directory.


Installing development version
------------------------------

Python-FSPS is being actively developed on GitHub so it's usually best
to use the most recent development version.
You can do this by cloning the `python-fsps repository
<https://github.com/dfm/python-fsps>`_ and building::

    git clone https://github.com/dfm/python-fsps.git
    cd python-fsps
    python setup.py install

Installing stable version
-------------------------

After you have FSPS and NumPy installed, you might be able to install the
most recent stable version of python-fsps using pip:

.. code-block:: bash

    pip install fsps

If this does not work, please ``pip uninstall fsps`` and follow the
instructions above for installing the development version.


Testing the installation
------------------------

To test the installation, run (this doesn't actually work yet)::

    python -c 'import fsps;fsps.test()'
