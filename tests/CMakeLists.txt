# Zaku test suite. Labels: unit (no external inputs), python, cling,
# regression (needs the reference production / flux release; skips otherwise).

find_program(ZAKU_PYTHON python3 REQUIRED)
find_program(ZAKU_ROOT_EXE root HINTS ${ROOT_BINDIR} REQUIRED)
get_filename_component(ZAKU_ROOTSYS "${ROOT_BINDIR}" DIRECTORY)

# The tests carry their own ROOT and Zaku environment so they run from a plain
# shell (no setup script sourced).
set(ZAKU_TEST_ENV
  "ROOTSYS=${ZAKU_ROOTSYS}"
  "PATH=${CMAKE_BINARY_DIR}/bin:${ROOT_BINDIR}:$ENV{PATH}"
  "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib:${ROOT_LIBRARY_DIR}:$ENV{LD_LIBRARY_PATH}"
  "ROOT_INCLUDE_PATH=${CMAKE_SOURCE_DIR}/include:${CMAKE_BINARY_DIR}/include:$ENV{ROOT_INCLUDE_PATH}"
  "PYTHONPATH=${CMAKE_SOURCE_DIR}/python:${ROOT_LIBRARY_DIR}:$ENV{PYTHONPATH}"
  "ZAKU_SOURCE_DIR=${CMAKE_SOURCE_DIR}"
  "ZAKU_BINARY_DIR=${CMAKE_BINARY_DIR}"
  "ZAKU_DATA=${CMAKE_SOURCE_DIR}"
)

# zaku_add_cxx_test(<name> <source> <label> [<label>...])
function(zaku_add_cxx_test name source)
  add_executable(${name} ${source})
  target_link_libraries(${name} PRIVATE Zaku)
  target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_compile_options(${name} PRIVATE -Wall -Wextra)
  set_target_properties(${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)
  add_test(NAME ${name} COMMAND ${name})
  set_tests_properties(${name} PROPERTIES ENVIRONMENT "${ZAKU_TEST_ENV}" LABELS "${ARGN}")
endfunction()

# zaku_add_python_test(<name> <script> <label> [<label>...])
function(zaku_add_python_test name script)
  add_test(NAME ${name} COMMAND ${ZAKU_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/${script})
  set_tests_properties(${name} PROPERTIES ENVIRONMENT "${ZAKU_TEST_ENV}" LABELS "${ARGN}")
endfunction()

# zaku_add_cling_test(<name> <macro> <label> [<label>...])
function(zaku_add_cling_test name macro)
  add_test(NAME ${name} COMMAND ${ZAKU_ROOT_EXE} -l -b -q ${CMAKE_CURRENT_SOURCE_DIR}/${macro})
  set_tests_properties(${name} PROPERTIES ENVIRONMENT "${ZAKU_TEST_ENV}" LABELS "${ARGN}")
endfunction()

# --- Stage 0 ----------------------------------------------------------------
zaku_add_cxx_test(test_units    unit/test_units.cc    unit)
zaku_add_cxx_test(test_types    unit/test_types.cc    unit)
zaku_add_cxx_test(test_detector unit/test_detector.cc unit)
zaku_add_python_test(python_import python/test_import.py python unit)
zaku_add_cling_test(cling_autoload cling/test_autoload.C cling unit)
zaku_add_cling_test(cling_explicit_load cling/test_explicit_load.C cling unit)

# The command-line tool must at least report its version.
if(ZAKU_BUILD_APPS)
  add_test(NAME cli_version COMMAND zaku --version)
  set_tests_properties(cli_version PROPERTIES ENVIRONMENT "${ZAKU_TEST_ENV}" LABELS "unit"
                       PASS_REGULAR_EXPRESSION "^${PROJECT_VERSION}")
endif()

# --- Stage 1 ----------------------------------------------------------------
zaku_add_cxx_test(test_flux        unit/test_flux.cc        unit)
zaku_add_cxx_test(test_experiment  unit/test_experiment.cc  unit)
zaku_add_cxx_test(test_manifest    unit/test_manifest.cc    unit)
zaku_add_cxx_test(test_truth       unit/test_truth.cc       unit)
zaku_add_cxx_test(test_flattree_io unit/test_flattree_io.cc unit)
zaku_add_python_test(python_experiment python/test_experiment.py python unit)

# Regression tests need the read-only reference inputs; they skip (exit 77) otherwise.
zaku_add_cxx_test(regression_flux_reference regression/test_flux_reference.cc regression)
zaku_add_cxx_test(regression_join_reference regression/test_join_reference.cc regression)
set_tests_properties(regression_flux_reference regression_join_reference PROPERTIES SKIP_RETURN_CODE 77)

# --- Stage 2 ----------------------------------------------------------------
zaku_add_cxx_test(test_formula_response unit/test_formula_response.cc unit)
zaku_add_cxx_test(test_smearing_reco    unit/test_smearing_reco.cc    unit)
zaku_add_cxx_test(test_binning_samples  unit/test_binning_samples.cc  unit)
zaku_add_cxx_test(test_oscillator       unit/test_oscillator.cc       unit)
zaku_add_python_test(python_detector python/test_detector.py python unit)
zaku_add_cxx_test(regression_recon_reference  regression/test_recon_reference.cc  regression)
zaku_add_cxx_test(regression_select_reference regression/test_select_reference.cc regression)
set_tests_properties(regression_recon_reference regression_select_reference PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 3600)

if(ZAKU_BUILD_APPS)
  add_test(NAME cli_show_response COMMAND zaku ${CMAKE_SOURCE_DIR}/configs/hk_t2k_like.yaml show-response)
  add_test(NAME cli_show_samples COMMAND zaku ${CMAKE_SOURCE_DIR}/configs/legacy_hk_mock.yaml show-samples)
  set_tests_properties(cli_show_response cli_show_samples PROPERTIES ENVIRONMENT "${ZAKU_TEST_ENV}" LABELS "unit")
  set_tests_properties(cli_show_response PROPERTIES PASS_REGULAR_EXPRESSION "9 categories")
  set_tests_properties(cli_show_samples PROPERTIES PASS_REGULAR_EXPRESSION "RHC_1Re0de")
endif()

if(ZAKU_BUILD_APPS)
  add_test(NAME cli_describe COMMAND zaku ${CMAKE_SOURCE_DIR}/configs/hk_t2k_like.yaml describe)
  add_test(NAME cli_list_productions COMMAND zaku ${CMAKE_SOURCE_DIR}/configs/legacy_hk_mock.yaml list-productions)
  set_tests_properties(cli_describe cli_list_productions PROPERTIES ENVIRONMENT "${ZAKU_TEST_ENV}" LABELS "unit")
  set_tests_properties(cli_list_productions PROPERTIES PASS_REGULAR_EXPRESSION "numu_x_numu_200K_FHC")
endif()
