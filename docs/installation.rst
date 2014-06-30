Installation
============

There are only two required dependencies for building this package: `FSPS
<http://people.ucsc.edu/~conroy/FSPS.html>`_ (obviously) and `NumPy
<http://www.numpy.org/>`_.
First, make sure that you have NumPy installed (using your favorite method) and access to the ``f2py`` executable (this should be installed along with NumPy).

Then, you need to follow the directions on `the FSPS page
<http://people.ucsc.edu/~conroy/FSPS.html>`_ to checkout and compile the FSPS
package. Note that Python-FSPS is built against specific versions of the FSPS
Fortran API so it may be necessary to checkout a slightly older revision.
Currently, revision 143 should be used::

   svn update -r 143

Once checked out, modify the Makefile as needed and compile FSPS.
These bindings rely on the value of the ``SPS_HOME`` environment variable
being correctly set and the compiled ``.o`` and ``.mod`` files be available in
the ``${SPS_HOME}/src`` directory.

Finally, clone the `python-fsps repository
<https://github.com/dfm/python-fsps>`_ and build::

    git clone https://github.com/dfm/python-fsps.git
    cd python-fsps
    python setup.py build_fsps
    python setup.py develop

To test the installation, run (this doesn't actually work yet)::

    python -c 'import fsps;fsps.test()'
