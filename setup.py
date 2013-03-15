#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import tempfile

from numpy.distutils.core import setup, Extension


if sys.argv[-1] == "publish":
    os.system("git rev-parse --short HEAD > COMMIT")
    os.system("python setup.py sdist upload")
    sys.exit()

os.system("cp ../fsps/src/*.o fsps")
os.system("cp ../fsps/src/*.mod fsps")

# First, make sure that the f2py interfaces exist.
interface_exists = os.path.exists("fsps/fsps.pyf")
if "interface" in sys.argv or not interface_exists:
    # Generate the Fortran signature/interface.
    cmd = "cd fsps;f2py fsps.f90"
    cmd += " -m _fsps -h fsps.pyf"
    cmd += " --overwrite-signature"
    os.system(cmd)

    if "interface" in sys.argv:
        sys.exit(0)

# Define the Fortran extension.
bart = Extension("fsps._fsps", ["fsps/fsps.pyf"])


setup(
    name="fsps",
    url="https://github.com/dfm/python-fsps",
    version="0.0.1",
    author="Dan Foreman-Mackey",
    author_email="danfm@nyu.edu",
    description="Python bindings for Charlie Conroy's FSPS.",
    long_description=open("README.rst").read(),
    packages=["fsps", ],
    package_data={"": ["README.rst"]},
    include_package_data=True,
    ext_modules=[bart],
    install_requires="numpy",
    classifiers=[
        # "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
)
