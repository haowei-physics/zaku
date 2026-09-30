cmake_minimum_required(VERSION 3.20)
project(Zaku VERSION 0.2.0 LANGUAGES CXX
        DESCRIPTION "Zaku Assembles Kinematics and Systematics")

# ---------------------------------------------------------------------------
# Build settings
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

option(ZAKU_BUILD_TESTS "Build the Zaku test suite" ON)
option(ZAKU_BUILD_APPS  "Build the zaku command-line tool" ON)

# The ROOT dictionary (.pcm and .rootmap) must live next to libZaku.so, in the
# build tree as well as in the install tree.
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------
find_package(ROOT 6.30 REQUIRED COMPONENTS Core RIO Tree TreePlayer Hist MathCore Physics Matrix)
find_package(yaml-cpp REQUIRED)
if(TARGET yaml-cpp::yaml-cpp)
  set(ZAKU_YAML_TARGET yaml-cpp::yaml-cpp)
else()
  set(ZAKU_YAML_TARGET yaml-cpp)
endif()
message(STATUS "Zaku: ROOT ${ROOT_VERSION} at ${ROOT_DIR}")
message(STATUS "Zaku: yaml-cpp target ${ZAKU_YAML_TARGET}")

# ---------------------------------------------------------------------------
# Generated version header
# ---------------------------------------------------------------------------
configure_file(cmake/VersionInfo.h.in ${CMAKE_BINARY_DIR}/include/zaku/core/VersionInfo.h @ONLY)

# ---------------------------------------------------------------------------
# Library
# ---------------------------------------------------------------------------
set(ZAKU_SOURCES
  src/core/Version.cc
  src/core/Types.cc
  src/core/Units.cc
  src/core/Detector.cc
  src/core/FluxSpec.cc
  src/core/Production.cc
  src/core/Manifest.cc
  src/core/Provenance.cc
  src/core/Experiment.cc
  src/core/ExperimentYaml.cc
  src/flux/FluxHistograms.cc
  src/flux/FluxPreparer.cc
  src/flux/FluxCovarianceInput.cc
  src/event/FlatEvent.cc
  src/event/FlatTreeCursor.cc
  src/event/WeightCursor.cc
  src/event/EventJoin.cc
  src/osc/Oscillator.cc
  src/detector/Formula.cc
  src/detector/Variables.cc
  src/detector/Response.cc
  src/detector/ResponseYaml.cc
  src/detector/Smearing.cc
  src/detector/Reco.cc
  src/detector/Binning.cc
  src/detector/Sample.cc
  src/stages/RateTable.cc
  src/stages/Reconstruct.cc
  src/stages/Select.cc
  third_party/NuFast/NuFast_LBL.cc
)

# Headers that go into the ROOT dictionary. They are named relative to
# include/ so that the dictionary records them the way user code includes
# them (rootcling resolves them through the target's include directories).
set(ZAKU_DICT_HEADERS
  zaku/core/Version.h
  zaku/core/Types.h
  zaku/core/Units.h
  zaku/core/Detector.h
  zaku/core/FluxSpec.h
  zaku/core/Production.h
  zaku/core/Manifest.h
  zaku/core/Provenance.h
  zaku/core/RunSettings.h
  zaku/core/Experiment.h
  zaku/flux/FluxHistograms.h
  zaku/flux/FluxPreparer.h
  zaku/flux/FluxCovarianceInput.h
  zaku/event/FlatEvent.h
  zaku/event/FlatTreeCursor.h
  zaku/event/WeightCursor.h
  zaku/event/EventJoin.h
  zaku/osc/Oscillator.h
  zaku/detector/Formula.h
  zaku/detector/Variables.h
  zaku/detector/Response.h
  zaku/detector/Smearing.h
  zaku/detector/Reco.h
  zaku/detector/Binning.h
  zaku/detector/Sample.h
  zaku/stages/Stages.h
  zaku/Zaku.h
)

add_library(Zaku SHARED ${ZAKU_SOURCES})
add_library(Zaku::Zaku ALIAS Zaku)
target_include_directories(Zaku
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/NuFast
)
target_link_libraries(Zaku
  PUBLIC  ROOT::Core ROOT::RIO ROOT::Tree ROOT::TreePlayer ROOT::Hist
          ROOT::MathCore ROOT::Physics ROOT::Matrix
  PRIVATE ${ZAKU_YAML_TARGET}
)
target_compile_options(Zaku PRIVATE -Wall -Wextra -Wpedantic)

# ROOT dictionary. NO_CXXMODULE keeps the classic libZaku_rdict.pcm +
# libZaku.rootmap pair (this ROOT is built with runtime C++ modules, which
# would otherwise replace them). NOINSTALL: we install them ourselves below so
# they end up next to the library regardless of CMAKE_INSTALL_LIBDIR.
ROOT_GENERATE_DICTIONARY(G__Zaku ${ZAKU_DICT_HEADERS}
  MODULE Zaku
  LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/src/LinkDef.h
  NO_CXXMODULE
  NOINSTALL
)

# ---------------------------------------------------------------------------
# Command-line tool
# ---------------------------------------------------------------------------
if(ZAKU_BUILD_APPS)
  add_executable(zaku apps/zaku.cc)
  target_link_libraries(zaku PRIVATE Zaku)
  target_compile_options(zaku PRIVATE -Wall -Wextra)
endif()

# ---------------------------------------------------------------------------
# Environment scripts: one for the build tree, one for the install tree
# ---------------------------------------------------------------------------
set(ZAKU_SETUP_LIBDIR   "${CMAKE_BINARY_DIR}/lib")
set(ZAKU_SETUP_BINDIR   "${CMAKE_BINARY_DIR}/bin")
set(ZAKU_SETUP_INCLUDE  "${CMAKE_SOURCE_DIR}/include:${CMAKE_BINARY_DIR}/include")
set(ZAKU_SETUP_PYTHON   "${CMAKE_SOURCE_DIR}/python")
set(ZAKU_SETUP_DATA     "${CMAKE_SOURCE_DIR}")
set(ZAKU_SETUP_ROOT     "${CMAKE_BINARY_DIR}")
configure_file(cmake/setup.sh.in ${CMAKE_BINARY_DIR}/setup.sh @ONLY)

set(ZAKU_SETUP_LIBDIR   "${CMAKE_INSTALL_PREFIX}/lib")
set(ZAKU_SETUP_BINDIR   "${CMAKE_INSTALL_PREFIX}/bin")
set(ZAKU_SETUP_INCLUDE  "${CMAKE_INSTALL_PREFIX}/include")
set(ZAKU_SETUP_PYTHON   "${CMAKE_INSTALL_PREFIX}/python")
set(ZAKU_SETUP_DATA     "${CMAKE_INSTALL_PREFIX}/share/zaku")
set(ZAKU_SETUP_ROOT     "${CMAKE_INSTALL_PREFIX}")
configure_file(cmake/setup.sh.in ${CMAKE_BINARY_DIR}/install-setup.sh @ONLY)

# ---------------------------------------------------------------------------
# Install
# ---------------------------------------------------------------------------
include(CMakePackageConfigHelpers)

install(TARGETS Zaku EXPORT ZakuTargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
)
install(FILES
  ${CMAKE_BINARY_DIR}/lib/libZaku_rdict.pcm
  ${CMAKE_BINARY_DIR}/lib/libZaku.rootmap
  DESTINATION lib
)
install(DIRECTORY include/ DESTINATION include)
install(FILES ${CMAKE_BINARY_DIR}/include/zaku/core/VersionInfo.h DESTINATION include/zaku/core)
install(DIRECTORY python/ DESTINATION python)
install(PROGRAMS ${CMAKE_BINARY_DIR}/install-setup.sh DESTINATION . RENAME setup.sh)
foreach(_dir data configs docs examples)
  if(EXISTS ${CMAKE_SOURCE_DIR}/${_dir})
    install(DIRECTORY ${_dir}/ DESTINATION share/zaku/${_dir})
  endif()
endforeach()
if(ZAKU_BUILD_APPS)
  install(TARGETS zaku RUNTIME DESTINATION bin)
endif()

install(EXPORT ZakuTargets NAMESPACE Zaku:: DESTINATION lib/cmake/Zaku)
configure_package_config_file(cmake/ZakuConfig.cmake.in
  ${CMAKE_BINARY_DIR}/ZakuConfig.cmake
  INSTALL_DESTINATION lib/cmake/Zaku)
write_basic_package_version_file(${CMAKE_BINARY_DIR}/ZakuConfigVersion.cmake
  COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_BINARY_DIR}/ZakuConfig.cmake ${CMAKE_BINARY_DIR}/ZakuConfigVersion.cmake
  DESTINATION lib/cmake/Zaku)

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
if(ZAKU_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
