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

Python-FSPS is built against specific versions of the FSPS Fortran API and data
files, so it is important that you have a recent version of FSPS through git.
Currently Python-FSPS is built against FSPS v3.2.

Installing stable version
-------------------------

Python-FSPS is available as a package on `PyPI
<https://pypi.org/project/fsps/>`_ that you can install using `pip`:

.. code-block:: bash

    python -m pip install fsps

Choosing Stellar Libraries
---------------------------

FSPS can use several different stellar isochrone and spectral libraries, but
switching between these libraries in python-FSPS requires (re-)installing
python-FSPS with specific compiler flags. The available libraries are described
in more detail in the FSPS documentation, but their names are:

* Stellar Isochrone libraries

   - `MIST` (default)
   - `PADOVA`
   - `PARSEC`
   - `BPASS` (in this case the spectral library and SSP parameters cannot be changed)
   - `BASTI`
   - `GENEVA`

* Stellar spectral libraries

   - `MILES` (default)
   - `BASEL`

* Dust emission libraries

   - `DL07` (default)
   - `THEMIS`

Changing any of these libraries requires switching off the relevant default,
which for isochrones is `MIST` and for spectra is `MILES`, and switching on the
desired library. As an example, you can change to Padova isochrones and the
BaSeL low resolution synthetic stellar library by re-installing:

.. code-block:: bash

    pip uninstall fsps
    FFLAGS="-DMIST=0 -DPADOVA=1 -DMILES=0 -DBASEL=1" python -m pip install fsps --no-binary fsps

where the `--no-binary fsps` flag is required to force building from source.

Installing development version
------------------------------

Python-FSPS is being actively developed on GitHub so you can got there to get
the most recent development version. You can do this by cloning the `python-fsps
repository <https://github.com/dfm/python-fsps>`_ and building:

.. code-block:: bash

    git clone --recursive https://github.com/dfm/python-fsps.git
    cd python-fsps
    python -m pip install .

Flags can be prepended to change the stellar libraries as described above. This
repository includes FSPS as a submodule, so if you forget the `--recursive` flag
above, you can get the submodule by running the following commands in the root
directory of the Python-FSPS repository:

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
