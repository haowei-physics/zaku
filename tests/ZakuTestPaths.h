// Locations of the read-only reference inputs used by the regression tests.
// Each can be overridden with an environment variable; a missing directory
// makes the test skip (exit code 77).
#pragma once
#include <TSystem.h>

#include <cstdlib>
#include <string>

#include "ZakuTest.h"

namespace zakutest {

inline std::string envOr(const char* var, const char* dflt) {
  const char* v = std::getenv(var);
  return v && *v ? std::string(v) : std::string(dflt);
}

/// The reference production (flat trees, weight files, fakerecon, splines, fakeselect, parameters).
inline std::string referenceDir() {
  return envOr("ZAKU_REFERENCE_DIR", "/home/haowei/shared_storage/neut_output/Mock_FD_Input_Data_2");
}
/// The legacy macro repository (lbl_config, gen cards, src).
inline std::string legacyDir() {
  return envOr("ZAKU_LEGACY_DIR", "/home/haowei/shared_storage/neut_output/mock_lbl_inputs");
}
/// The T2K 2020 public flux release.
inline std::string fluxReleaseDir() {
  return envOr("ZAKU_FLUX_RELEASE_DIR", "/home/haowei/shared_storage/t2kflux_2020_public_release");
}
/// The Zaku source tree (data/, configs/).
inline std::string sourceDir() { return envOr("ZAKU_SOURCE_DIR", "."); }
/// A scratch directory for test outputs.
inline std::string scratchDir() {
  const std::string d = envOr("ZAKU_BINARY_DIR", ".") + "/tests/scratch";
  gSystem->mkdir(d.c_str(), true);
  return d;
}

inline bool exists(const std::string& path) { return !gSystem->AccessPathName(path.c_str()); }

inline void requirePath(const std::string& path, const std::string& what) {
  if (!exists(path)) skip(what + " not available at " + path);
}

}  // namespace zakutest
