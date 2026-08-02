#!/bin/bash

# This script cycles through the different FSPS libraries and runs afe.py to run
# tests and generate the color evolution plots and data.


PYFSPS_DIR=`realpath ../`
TESTDIR=$PYFSPS_DIR/demos
SPECLIBS=( C3K_LR C3K_HR MILES )
ISOCS=( MIST PADOVA PARSEC )

# Do AFE_FLAG=1 for MIST + C3K_LR (the defaults)
cd $PYFSPS_DIR
FFLAGS="-DAFE_FLAG=1"
echo "Recompiling with flags: ${FFLAGS}"
python -m pip uninstall -y fsps
FFLAGS=${FFLAGS} python -m pip install . --no-binary fsps
cd $TESTDIR
python test_lib.py

# Now do the other combinations of isochrones and spectral libraries.  We will
# always use AFE_FLAG=0 (the default), but we will switch the isochrone and
# spectral library.
for iso in "${ISOCS[@]}"; do
    for lib in "${SPECLIBS[@]}"; do
        echo "Running tests for ${lib} and ${iso}"
        cd $PYFSPS_DIR

        FFLAGS=""

        # Switch the isochrone
        if [[ ${iso} != "MIST" ]]; then
            FFLAGS=${FFLAGS}" -DMIST=0"
            FFLAGS=${FFLAGS}" -D${iso}=1"
        fi

        # Switch the spectral library
        if [[ ${lib} != "C3K_LR" ]]; then
            FFLAGS=${FFLAGS}" -DC3K_LR=0"
            FFLAGS=${FFLAGS}" -D${lib}=1"
        fi

        echo $FFLAGS

        # Recompile with flags
        echo "Recompiling with flags: ${FFLAGS}"
        python -m pip uninstall -y fsps
        FFLAGS=${FFLAGS} python -m pip install . --no-binary fsps

        # run the tests
        cd $TESTDIR
        python test_lib.py
    done
done

# Now do BPASS, which has no spectral library attached
cd $PYFSPS_DIR
FFLAGS="-DMIST=0 -DC3K_LR=0 -DBPASS=1"
echo "Recompiling with flags: ${FFLAGS}"
python -m pip uninstall -y fsps
FFLAGS=${FFLAGS} python -m pip install . --no-binary fsps
cd $TESTDIR
python test_lib.py
