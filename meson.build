project(
    'fsps', 
    'c', 'cpp', 'fortran',
    version: run_command(
        ['python3', 'tools/version.py'],
        check : true
    ).stdout().strip(),
    default_options: ['warning_level=0', 'optimization=3'],
    meson_version: '>=0.63.0'
)

# Enable fortran and check arguments
# add_languages('fortran', native: false)
ff = meson.get_compiler('fortran')
f_args = ff.get_supported_arguments('-Wno-conversion', '-cpp', '-fPIC')
add_project_arguments(f_args, language: 'fortran')

# Find and link the quadmath library on Windows
if host_machine.system() == 'windows'
  cc = meson.get_compiler('c')
  quadmath_dep = cc.find_library('quadmath', required: false)
  if quadmath_dep.found()
    add_project_link_arguments('-lquadmath', language: ['c', 'cpp', 'fortran'])
  endif
endif

# Check to see if the submodule is initialized
fs = import('fs')
if not fs.exists('src/fsps/libfsps/src/sps_vars.f90')
  error('It looks like FSPS is not included in your source distribution. If you\'re installing using git, don\'t forget to include the submodules using a recursive clone or by running: \ngit submodule update --init')
endif

# Python setup
py_mod = import('python')
py3 = py_mod.find_installation('python3')

# Handle distribution
meson.add_dist_script('tools/dist.py')

subdir('src/fsps')
