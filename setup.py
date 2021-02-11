#!/usr/bin/env python

import os
import glob
import setuptools  # noqa: magic
import numpy.f2py
from numpy.distutils.core import setup, Extension


if not os.path.exists("src/fsps/libfsps/src/sps_vars.f90"):
    raise RuntimeError(
        "It looks like FSPS is not included in your source distribution. "
        "If you're installing using git, don't forget to include the "
        "submodules using a recursive clone or by running: \n"
        "git submodule init\n"
        "git submodule update"
    )


numpy.f2py.run_main(
    [
        "-m",
        "_fsps",
        "-h",
        "src/fsps/_fsps.pyf",
        "--overwrite-signature",
        "src/fsps/fsps.f90",
    ]
)

flags = ["-cpp", "-fPIC"]
ext = Extension(
    "fsps._fsps",
    sources=[
        "src/fsps/libfsps/src/sps_vars.f90",
        "src/fsps/libfsps/src/sps_utils.f90",
    ]
    + [
        f
        for f in glob.glob("src/fsps/libfsps/src/*.f90")
        if os.path.basename(f)
        not in [
            "autosps.f90",
            "simple.f90",
            "lesssimple.f90",
            "sps_vars.f90",
            "sps_utils.f90",
        ]
    ]
    + [
        "src/fsps/fsps.f90",
        "src/fsps/_fsps.pyf",
    ],
    extra_f90_compile_args=flags,
)

# The final setup command. Note: we override the `build_ext` command with our
# custom version from above.
setup(
    name="fsps",
    url="https://github.com/dfm/python-fsps",
    author="Python-FSPS developers",
    author_email="foreman.mackey@gmail.com",
    description="Python bindings for Charlie Conroy's FSPS.",
    long_description=open("README.rst").read(),
    packages=["fsps"],
    package_dir={"": "src"},
    package_data={
        "": ["README.rst", "LICENSE.rst", "AUTHORS.rst"],
        "fsps": [
            "libfsps/data/*",
            "libfsps/dust/*",
            "libfsps/dust/*/*",
            "libfsps/ISOCHRONES/*/*",
            "libfsps/ISOCHRONES/*/*/*",
            "libfsps/nebular/*",
            "libfsps/OUTPUTS/*",
            "libfsps/SPECTRA/*",
            "libfsps/SPECTRA/*/*",
        ],
    },
    include_package_data=True,
    install_requires=["numpy"],
    ext_modules=[ext],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
    zip_safe=False,
)
