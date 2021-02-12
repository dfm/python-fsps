#!/usr/bin/env python

import os
import sys
import glob
import setuptools  # noqa: magic
import numpy.f2py
from numpy.distutils.core import setup, Extension
from numpy.distutils.command.build_ext import build_ext
from distutils.file_util import copy_file


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
    extra_f90_compile_args=["-cpp"],
)


class custom_build_ext(build_ext):
    def run(self):
        build_ext.run(self)

        # HACKS: Copy over any extra DLL files
        # Note: numpy already tries to do this, but it doesn't do a very good job
        # when using the `src/pkgname` layout. We'll fix that here.
        pkg_roots = {
            os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
            for ext in self.extensions
        }
        for pkg_root in pkg_roots:
            shared_lib_dir = pkg_root
            for fn in os.listdir(self.extra_dll_dir):
                if not os.path.isdir(shared_lib_dir):
                    os.makedirs(shared_lib_dir)
                if not fn.lower().endswith(".dll"):
                    continue
                runtime_lib = os.path.join(self.extra_dll_dir, fn)
                copy_file(runtime_lib, shared_lib_dir)


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
    cmdclass={"build_ext": custom_build_ext},
    zip_safe=False,
)
