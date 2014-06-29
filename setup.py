#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess as sp


try:
    from setuptools import setup
    from setuptools.command.build_ext import build_ext
except ImportError:
    from distutils.core import setup
    from distutils.command.build_ext import build_ext


class build_fsps(build_ext):

    def run(self):
        # Generate the Fortran signature/interface.
        cmd = "cd fsps;f2py fsps.f90"
        cmd += " -m _fsps -h fsps.pyf"
        cmd += " --overwrite-signature"
        print("Running: {0}".format(cmd))
        sp.call(cmd, shell=True)

        # Find the FSPS source files.
        fsps_dir = os.path.join(os.environ["SPS_HOME"], "src")
        fns = [f for f in glob.glob(os.path.join(fsps_dir, "*.o"))
               if os.path.basename(f) not in ["autosps.o", "simple.o",
                                              "lesssimple.o"]]

        # Add the interface source files to the file list.
        fns += ["fsps.f90", "fsps.pyf"]

        # Compile the library.
        cmd = "cd fsps;f2py -c -I{0} ".format(fsps_dir)
        cmd += " ".join(fns)
        cmd += " --f90flags=-cpp"
        cmd += " --f90flags=-fPIC"
        print("Running: {0}".format(cmd))
        sp.call(cmd, shell=True)

        # Move the compiled object to the correct path.
        infn = self.get_ext_filename("fsps._fsps")
        outfn = self.get_ext_fullpath("fsps._fsps")
        try:
            os.makedirs(os.path.split(outfn)[0])
        except os.error:
            pass
        cmd = "mv {0} {1}".format(infn, outfn)
        print("Running: {0}".format(cmd))
        sp.call(cmd, shell=True)

    def get_outputs(self):
        return [self.get_ext_fullpath("fsps._fsps")]


if "publish" in sys.argv[-1]:
    os.system("git rev-parse --short HEAD > COMMIT")
    os.system("python setup.py sdist upload")
    sys.exit()


# Hackishly inject a constant into builtins to enable importing of the
# package before the library is built.
if sys.version_info[0] < 3:
    import __builtin__ as builtins
else:
    import builtins
builtins.__FSPS_SETUP__ = True
from fsps import __version__

setup(
    name="fsps",
    url="https://github.com/dfm/python-fsps",
    version=__version__,
    author="Dan Foreman-Mackey",
    author_email="danfm@nyu.edu",
    description="Python bindings for Charlie Conroy's FSPS.",
    long_description=open("README.rst").read(),
    packages=["fsps"],
    package_data={"": ["README.rst", "LICENSE.rst"],
                  "fsps": ["_fsps.so"]},
    include_package_data=True,
    scripts=glob.glob("scripts/*.py"),
    install_requires="numpy",
    cmdclass={
        "build_ext": build_fsps,
    },
    classifiers=[
        # "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
)
