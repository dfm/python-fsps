if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/sps_vars.f90")
  message(FATAL_ERROR "The source code in the FSPS submodule was not found."
    "Please run 'git submodule update --init' to initialize the submodule.")
endif()

# Generate the f2py wrappers
set(F2PY_SOURCES
  "${CMAKE_CURRENT_BINARY_DIR}/_fspsmodule.c"
  "${CMAKE_CURRENT_BINARY_DIR}/_fsps-f2pywrappers2.f90"
)
add_custom_command(
  OUTPUT ${F2PY_SOURCES}
  DEPENDS fsps.f90
  VERBATIM
  COMMAND "${Python_EXECUTABLE}" -m numpy.f2py
  "${CMAKE_CURRENT_SOURCE_DIR}/fsps.f90" -m _fsps --lower
  --build-dir "${CMAKE_CURRENT_BINARY_DIR}")

# List out the explicit FSPS sources; we don't use a glob here since
# the 'src' directory includes some executables that we don't want to
# include here.
set(
  FSPS_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/sps_vars.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/sps_utils.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/add_agb_dust.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/add_bs.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/add_dust.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/add_nebular.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/add_remnants.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/add_xrb.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/agn_dust.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/attn_curve.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/compsp.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/csp_gen.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/funcint.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/get_lumdist.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/get_tuniv.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/getindx.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/getmags.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/getspec.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/igm_absorb.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/imf.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/imf_weight.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/intsfwght.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/linterp.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/linterparr.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/locate.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/mod_gb.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/mod_hb.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/pz_convol.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/sbf.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/setup_tabular_sfh.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/sfh_weight.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/sfhinfo.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/sfhlimit.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/sfhstat.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/smoothspec.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/spec_bin.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/sps_setup.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/ssp_gen.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/tsum.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/vacairconv.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/write_isochrone.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/libfsps/src/ztinterp.f90"
)

# Define the Python library
python_add_library(
  _fsps
  MODULE
  "${CMAKE_CURRENT_SOURCE_DIR}/fsps.f90"
  ${FSPS_SOURCES}
  ${F2PY_SOURCES}
  "${F2PY_INCLUDE_DIR}/fortranobject.c"
  WITH_SOABI)
target_link_libraries(_fsps PUBLIC Python::NumPy)
target_include_directories(_fsps PUBLIC "${F2PY_INCLUDE_DIR}")
install(TARGETS _fsps DESTINATION ${SKBUILD_PROJECT_NAME})
