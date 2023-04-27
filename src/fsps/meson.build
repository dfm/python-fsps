# Include directories
incdir_numpy = run_command(py3, ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'], check: true).stdout().strip()
incdir_f2py = run_command(py3, ['../../tools/f2py_include.py'], check: true).stdout().strip()
inc_np = include_directories(incdir_numpy, incdir_f2py)

# Generate the f2py source
f2py_source = custom_target(
    '_fspsmodule.c',
    input: ['fsps.f90'],
    output: ['_fspsmodule.c', '_fsps-f2pywrappers2.f90'],
    command: [py3, '-m', 'numpy.f2py', '@INPUT@', '-m', '_fsps', '--build-dir', '@OUTDIR@'],
)

# These are the source files for libfsps
lib_source = [
    'libfsps/src/sps_vars.f90',
    'libfsps/src/sps_utils.f90',
    'libfsps/src/add_agb_dust.f90',
    'libfsps/src/add_bs.f90',
    'libfsps/src/add_dust.f90',
    'libfsps/src/add_nebular.f90',
    'libfsps/src/add_remnants.f90',
    'libfsps/src/add_xrb.f90',
    'libfsps/src/agn_dust.f90',
    'libfsps/src/attn_curve.f90',
    'libfsps/src/compsp.f90',
    'libfsps/src/csp_gen.f90',
    'libfsps/src/funcint.f90',
    'libfsps/src/get_lumdist.f90',
    'libfsps/src/get_tuniv.f90',
    'libfsps/src/getindx.f90',
    'libfsps/src/getmags.f90',
    'libfsps/src/getspec.f90',
    'libfsps/src/igm_absorb.f90',
    'libfsps/src/imf.f90',
    'libfsps/src/imf_weight.f90',
    'libfsps/src/intsfwght.f90',
    'libfsps/src/linterp.f90',
    'libfsps/src/linterparr.f90',
    'libfsps/src/locate.f90',
    'libfsps/src/mod_gb.f90',
    'libfsps/src/mod_hb.f90',
    'libfsps/src/pz_convol.f90',
    'libfsps/src/sbf.f90',
    'libfsps/src/setup_tabular_sfh.f90',
    'libfsps/src/sfh_weight.f90',
    'libfsps/src/sfhinfo.f90',
    'libfsps/src/sfhlimit.f90',
    'libfsps/src/sfhstat.f90',
    'libfsps/src/smoothspec.f90',
    'libfsps/src/spec_bin.f90',
    'libfsps/src/sps_setup.f90',
    'libfsps/src/ssp_gen.f90',
    'libfsps/src/tsum.f90',
    'libfsps/src/vacairconv.f90',
    'libfsps/src/write_isochrone.f90',
    'libfsps/src/ztinterp.f90',
]

# And finally the definition of the library
py3.extension_module(
    '_fsps',
    'fsps.f90',
    lib_source,
    f2py_source,
    incdir_f2py + '/fortranobject.c',
    include_directories: inc_np,
    install: true,
    subdir: 'fsps',
)

# Now the source files
py3.install_sources(
    [
        '__init__.py',
        'filters.py',
        'fsps_version.py',
        'fsps.py',
        'sps_home.py',
    ],
    pure: false,
    subdir: 'fsps'
)
