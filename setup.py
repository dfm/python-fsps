#!/usr/bin/env python

import os
import sys
import glob
import shutil
import subprocess

try:
    from setuptools import setup, Extension
    from setuptools.command.build_ext import build_ext
except ImportError:
    from distutils.core import setup, Extension
    from distutils.command.build_ext import build_ext

def invoke_f2py(files, flags=[], wd=None):
    from numpy.f2py import main

    olddir = os.path.abspath(os.curdir)
    oldargv = list(sys.argv)
    try:
        if wd is not None:
            os.chdir(wd)
        sys.argv = ['f2py']
        sys.argv.extend(files)
        sys.argv.extend(flags)

        main()
    finally:
        sys.argv = oldargv
        os.chdir(olddir)

def is_fortran_program(path):
    with open(path, 'r') as f:
        return any(line.lstrip().upper().startswith('PROGRAM ') for line in f)

class build_fsps(build_ext):

    def run(self):
        # Generate the Fortran signature/interface.
        files = ['fsps.f90']
        flags = "-m _fsps -h fsps.pyf --overwrite-signature".split()
        print("Running f2py on {0} with flags {1}".format(files, flags))
        invoke_f2py(['fsps.f90'], flags, wd='fsps')

        # Compile FSPS
        fsps_dir = os.path.join(os.environ["SPS_HOME"], "src")
        os.environ['F90FLAGS'] = os.environ.get('F90FLAGS', '') + " -fPIC"

        fns = [f.rsplit('.', 1)[0] + '.o' 
                for f in glob.glob(os.path.join(fsps_dir, '*.f90'))
                if not is_fortran_program(f)]

        return_code = subprocess.call(['make', '-C'+fsps_dir] + fns)
        if return_code != 0:
            sys.exit(return_code)

        # Add the interface source files to the file list.
        fns += ["fsps.f90", "fsps.pyf"]

        # Compile the library.
        flags = '-c -I{0} --f90flags=-cpp --f90flags=-fPIC'.format(fsps_dir)
        flags = flags.split()
        print("Running f2py on {0} with flags {1}".format(fns, flags))
        invoke_f2py(fns, flags, wd='fsps')

        # Move the compiled library to the correct directory.
        infn = os.path.abspath(self.get_ext_filename("fsps._fsps"))
        outfn = os.path.abspath(self.get_ext_fullpath("fsps._fsps"))
        if infn != outfn:
            try:
                os.makedirs(os.path.dirname(outfn))
            except os.error:
                pass
            print("Copying {0} to {1}".format(infn, outfn))
            shutil.copyfile(infn, outfn)


if "publish" in sys.argv[-1]:
    os.system("python setup.py sdist upload")
    sys.exit()


# Hackishly inject a constant into builtins to enable importing of the
# package before the library is built.
if sys.version_info[0] < 3:
    import __builtin__ as builtins
else:
    import builtins
builtins.__FSPS_SETUP__ = True
from fsps import __version__  # NOQA

# This is a fake extension that is used to trick distutils into building our
# real library using the `build_fsps` function above even when `install` is
# called.
ext = Extension("fsps._fsps", sources=["fsps/fsps.f90"])

# The final setup command. Note: we override the `build_ext` command with our
# custom version from above.
setup(
    name="fsps",
    url="https://github.com/dfm/python-fsps",
    version=__version__,
    author="Dan Foreman-Mackey",
    author_email="danfm@nyu.edu",
    description="Python bindings for Charlie Conroy's FSPS.",
    long_description=open("README.rst").read(),
    packages=["fsps"],
    package_data={
        "": ["README.rst", "LICENSE.rst", "AUTHORS.rst"],
        "fsps": ["_fsps.so"],
    },
    include_package_data=True,
    ext_modules=[ext],
    scripts=glob.glob("scripts/*.py"),
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
    setup_requires=[
        'numpy==1.16; python_version=="2.7"',
        'numpy; python_version>="3.5"',
        'pygit2',
    ],
)
