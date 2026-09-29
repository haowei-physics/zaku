# Zaku test suite. Labels: unit (no external inputs), python, cling,
# regression (needs the reference production / flux release; skips otherwise).

find_program(ZAKU_PYTHON python3 REQUIRED)
find_program(ZAKU_ROOT_EXE root HINTS ${ROOT_BINDIR} REQUIRED)

set(ZAKU_TEST_ENV
  "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib:$ENV{LD_LIBRARY_PATH}"
  "ROOT_INCLUDE_PATH=${CMAKE_SOURCE_DIR}/include:${CMAKE_BINARY_DIR}/include:$ENV{ROOT_INCLUDE_PATH}"
  "PYTHONPATH=${CMAKE_SOURCE_DIR}/python:$ENV{PYTHONPATH}"
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
