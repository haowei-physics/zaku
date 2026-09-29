// Stage 0 Cling check: with libZaku on LD_LIBRARY_PATH the rootmap must make
// ROOT load the library on first use, without an explicit gSystem->Load.
// Run:  root -l -b -q test_autoload.C
#include <cstdio>
#include <string>

#include "zaku/Zaku.h"

int test_autoload() {
  int failures = 0;
  const std::string v = zaku::version();
  if (v.size() < 5) {
    std::printf("FAIL: version '%s'\n", v.c_str());
    ++failures;
  }
  const zaku::Detector det = zaku::Detector::water(188.0);
  const double n = det.nNucleons();
  if (std::abs(n - 1.132e35) / 1.132e35 > 1e-3) {
    std::printf("FAIL: nNucleons %g\n", n);
    ++failures;
  }
  if (zaku::NuSpecies::fromName("nueb").pdg() != -12) {
    std::printf("FAIL: NuSpecies\n");
    ++failures;
  }
  std::printf("cling_autoload: Zaku %s, %d failure(s)\n", v.c_str(), failures);
  if (failures) gSystem->Exit(1);
  return failures;
}
