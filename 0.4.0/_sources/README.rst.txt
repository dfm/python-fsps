Building the documentation
--------------------------

The docs are now built using GitHub actions so you shouldn't ever need to build
them. But, to test, you can use the following steps:

1. In this directory, run::

    make clean
    make dirhtml

2. Confirm that the results in ``_build/dirhtml`` look sensible.
3. Change to the root directory of the repository and run (changing
   ``VERSION`` to the correct version number)::

    git checkout gh-pages
    mkdir VERSION
    cp -r docs/_build/dirhtml/* VERSION/
    git add VERSION
    git commit -m "Updated docs for version VERSION"
    git push

4. If you want to release this version as the default (stable) version, run::

    rm current
    ln -s VERSION current
    git add current
    git commit -m "Updating stable version to VERSION"
    git push
