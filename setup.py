#!/usr/bin/env python
# encoding: utf-8

import os
import subprocess
import glob
import shutil
from setuptools import setup, Command

fsps = "fsps"
desc = open("README.md").read()
with open("requirements.txt") as f:
    required = f.readlines()


class BuildFSPS(Command):
    description = "Build FSPS and the fsps extension module"
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        fsps_src_dir = os.path.join(os.environ["SPS_HOME"], "src")

        # Build FSPS.
        cmd = "cd {0};make".format(fsps_src_dir)
        subprocess.call(cmd, shell=True)
        links = glob.glob(os.path.join(fsps_src_dir, "*.o"))
        [links.remove(os.path.join(fsps_src_dir, f)) for f in ["simple.o",
            "lesssimple.o", "autosps.o"]]

        # Copy the modules to the right place.
        modules = glob.glob(os.path.join(fsps_src_dir, "*.mod"))
        for m in modules:
            dest = os.path.join(fsps, os.path.split(m)[1])
            shutil.copy(m, dest)

        # Build f2py signature.
        cmd = "cd {0};f2py {0}.f95 -m _{0} -h {0}.pyf".format(fsps)
        cmd += " --overwrite-signature"
        subprocess.call(cmd, shell=True)

        # Compile the module.
        cmd = "cd {0};f2py -c --fcompiler=gnu95 --f90flags=\"-fPIC\" --f77flags=\"-fPIC\" {0}.pyf {0}.f95 ".format(fsps)
        cmd += " ".join(links)
        subprocess.call(cmd, shell=True)

        # Clean up.
        for p in glob.glob(os.path.join(fsps, "*.mod")):
            os.remove(p)

cmdclass = {"build_fsps": BuildFSPS}

setup(
    name=fsps,
    version="0.0.1",
    author="Daniel Foreman-Mackey",
    author_email="danfm@nyu.edu",
    url="https://github.com/dfm/python-fsps",
    py_modules=[fsps, ],
    cmdclass=cmdclass,
    install_requires=required,
    license="MIT",
    description="Python bindings to Conroy's (FSPS) Fortran code.",
    long_description=desc,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
)
