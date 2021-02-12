Installation
============

Prerequisites
-------------

Python-FSPS no longer has any compile time prerequisites (besides a working
Fortran compiler), but (if you're not installing using git, see the "development
version" section below) it does require a clone of the `FSPS
<https://github.com/cconroy20/fsps>`_ project for the data files. To get this
set up, you can clone that repository using:

.. code-block:: bash

    export SPS_HOME="/path/where/you/want/to/download/fsps"
    git clone https://github.com/cconroy20/fsps.git $SPS_HOME

Where, on different systems, you might have to specify the `SPS_HOME`
environment variable. Python-FSPS requires that variable to be set and it will
fail to import if it is set incorrectly.

Python-FSPS is built against specific versions of the FSPS Fortran API, so it is
important that you have a recent version of FSPS through git. Currently
Python-FSPS is built against FSPS v3.0.

Installing stable version
-------------------------

Python-FSPS is available as a package on `PyPI
<https://pypi.org/project/fsps/>`_ that you can install using `pip`:

.. code-block:: bash

    python -m pip install fsps

Installing development version
------------------------------

Python-FSPS is being actively developed on GitHub so you can got there to get
the most recent development version. You can do this by cloning the `python-fsps
repository <https://github.com/dfm/python-fsps>`_ and building:

..code-block:: bash

    git clone --recursive https://github.com/dfm/python-fsps.git
    cd python-fsps
    python -m pip install .

This repository incldes FSPS as a submodule, so if you forget the `--recursive`
flag above, you can get the submodule by running the following commands in the
root directory of the Python-FSPS repository:

.. code-block:: bash

    git submodule init
    git submodule update

If you install Python-FSPS using this method, you don't actually need a separate
FSPS clone and you can just set the `SPS_HOME` variable as:

.. code-block:: bash

    export SPS_HOME=$(pwd)/src/fsps/libfsps

It is recommended that you install using `pip` (even for a local clone), and you
can use `pip install -e .` to install an "editable" version (like you would get
with `setup.py develop`). But if you want to use the `setup.py` script directly,
you'll need to install some prerequisites in advance:

.. code-block:: bash

    python -m pip install numpy "setuptools_scm[toml]"
    python setup.py develop

Troubleshooting
---------------
Here are possible fixes for the most common installation issues:

* Version mismatch.  Python-FSPS is built against specific versions of FSPS.
  This can cause problems if you have a different version of FSPS.  Python-FSPS
  attempts to check that you have a compatible version of FSPS during
  installation, and if you do not it will print an error.  This requires that
  you have `git` and that FSPS be under git version control.  If you do not
  have the correct version, try:

.. code-block:: bash

      cd $SPS_HOME
      git pull

  to get and compile the most recent version. Then try to install Python-FSPS again.
