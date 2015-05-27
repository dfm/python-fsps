Installation
============

Prerequisites
-------------

There are only two required dependencies for building this package: `FSPS
<https://github.com/cconroy20/fsps>`_ (obviously) and `NumPy
<http://www.numpy.org/>`_.
First, make sure that you have NumPy installed (using your favorite method)
and access to the ``f2py`` executable (this should be installed along with
NumPy).

Then, you need to follow the directions on `the FSPS page
<https://github.com/cconroy20/fsps>`_ to clone and compile the FSPS
package. Note that Python-FSPS is built against specific versions of
the FSPS Fortran API so it may be necessary to checkout a slightly
older commit - you will be notified of this and how to do it when you
try to import the fsps module.

Once checked out, modify the Makefile as needed and compile FSPS.  For
some compilers (e.g. intel) it is necessary to set the ``-fPIC`` flag
when compiling FSPS. These bindings rely on the value of the
``SPS_HOME`` environment variable being correctly set and the compiled
``.o`` and ``.mod`` files be available in the ``${SPS_HOME}/src``
directory.


Installing stable version
-------------------------

After you have FSPS and NumPy installed, you should be able to install the
most recent stable version of python-fsps using pip:

.. code-block:: bash

    pip install fsps


Installing development version
------------------------------

Python-FSPS is being actively developed on GitHub so it's not a bad idea to
use the most recent development version.
You can to this by cloning the `python-fsps repository
<https://github.com/dfm/python-fsps>`_ and building::

    git clone https://github.com/dfm/python-fsps.git
    cd python-fsps
    python setup.py install


Testing the installation
------------------------

To test the installation, run (this doesn't actually work yet)::

    python -c 'import fsps;fsps.test()'
