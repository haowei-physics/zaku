// Stage 0 Cling check: the documented explicit-load path.
// Run:  root -l -b -q test_explicit_load.C
#include <cstdio>

int test_explicit_load() {
  if (gSystem->Load("libZaku") < 0) {
    std::printf("FAIL: gSystem->Load(libZaku)\n");
    gSystem->Exit(1);
  }
  // Everything below is resolved through the dictionary + autoparsed headers.
  gROOT->ProcessLine("#include \"zaku/Zaku.h\"");
  const int rc = static_cast<int>(gROOT->ProcessLine(
      "zaku::Detector::water(188.0).neutTarget().numFreeProtons == 2 ? 0 : 1"));
  std::printf("cling_explicit_load: rc=%d\n", rc);
  if (rc != 0) gSystem->Exit(1);
  return rc;
}
