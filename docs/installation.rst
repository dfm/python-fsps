Installation
============

Prerequisites
-------------

There are only two required dependencies for building this package: `FSPS
<https://github.com/cconroy20/fsps>`_ (obviously) and `NumPy
<https://www.numpy.org/>`_.
First, make sure that you have NumPy installed (using your favorite method)
and access to the ``f2py`` executable (this should be installed along with
NumPy).

Then, you need to follow the directions on `the FSPS page
<https://github.com/cconroy20/fsps>`_ to clone and compile the FSPS
package. Python-FSPS is built against specific versions of the FSPS Fortran
API, so it is important that you have a recent version of FSPS through
git. Currently Python-FSPS is built against FSPS v3.0.

These bindings rely on the value of the ``SPS_HOME`` environment
variable being correctly set and the compiled ``.o`` and ``.mod``
files being available in the ``${SPS_HOME}/src`` directory.

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

After you have FSPS and NumPy installed, you _might_ be able to install the
most recent stable version of python-fsps using pip:

.. code-block:: bash

    pip install fsps

If this does not work, please ``pip uninstall fsps`` and follow the
instructions above for installing the development version.

Troubleshooting
-----------------------
Here are possible fixes for the most common installation issues:

* Version mismatch.  Python-FSPS is built against specific versions of FSPS.
  This can cause problems if you have a different version of FSPS.  Python-FSPS
  attempts to check that you have a compatible version of FSPS during
  installation, and if you do not it will print an error.  This requires that
  you have `git` and that FSPS be under git version control.  If you do not
  have the correct version, try::

      cd $SPS_HOME
      git pull
      cd src
      make clean
      make all

  to get and compile the most recent version.  Then try to install python-FSPS again.

* fPIC.  For some compilers (e.g. Intel) it may be necessary to set the
  ``-fPIC`` flag in the Makefile when compiling FSPS. Please try this if you
  encounter long python-FSPS installation errors that contain the
  text ``can not be used when making a shared object; recompile with -fPIC``.
  For gfortran the flag is not necesary.::

    cd $SPS_HOME/src

  edit the ``F90FLAGS=`` line in the Makefile to include ``-fPIC`` and then::

    make clean
    make all

  Then try to install Python-FSPS again.

* Problems with ``f2py``. The 1.10.x versions of ``numpy`` introduced some
  compatibility problems in
  ``f2py``.  These seem to have been fixed in 1.11.x versions and the fixes
  retroactively added to the 1.10.x versions, but if you encounter long
  complicated error messages that end with something like ``KeyError: 'void'``,
  please consider upgrading your ``numpy`` installation.
