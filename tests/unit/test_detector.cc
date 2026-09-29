#include "ZakuTest.h"
#include "zaku/core/Detector.h"

int main() {
  using namespace zaku;

  // The legacy chain used 1.132e35 nucleons for 188 kt of water.
  const Detector hk = Detector::water(188.0);
  ZAKU_CHECK_CLOSE(hk.nNucleons(), 1.132e35, 1e-3);
  ZAKU_CHECK_STR_EQ(hk.name(), "water");
  ZAKU_CHECK_CLOSE(hk.fiducialMass_kt(), 188.0, 1e-12);
  ZAKU_CHECK_EQ(static_cast<int>(hk.composition().size()), 2);
  ZAKU_CHECK_STR_EQ(hk.genieTargetString(), "1000080160[0.888889],1000010010[0.111111]");

  const NeutTarget nt = hk.neutTarget();
  ZAKU_CHECK_EQ(nt.numNucleons, 16);
  ZAKU_CHECK_EQ(nt.numBoundProtons, 8);
  ZAKU_CHECK_EQ(nt.numFreeProtons, 2);

  // nuclei: O-16 count = 188e9 g * (16/18) / 16 * N_A
  ZAKU_CHECK_CLOSE(hk.nNuclei(hk.composition()[0]), 188.0e9 * (16.0 / 18.0) / 16.0 * 6.02214076e23, 1e-9);
  ZAKU_CHECK_CLOSE(hk.nNucleiTotal(), hk.nNuclei(hk.composition()[0]) + hk.nNuclei(hk.composition()[1]), 1e-12);

  // argon and scintillator
  const Detector ar = Detector::argon(40.0);
  ZAKU_CHECK_EQ(ar.neutTarget().numNucleons, 40);
  ZAKU_CHECK_EQ(ar.neutTarget().numBoundProtons, 18);
  ZAKU_CHECK_EQ(ar.neutTarget().numFreeProtons, 0);
  ZAKU_CHECK_STR_EQ(ar.genieTargetString(), "1000180400[1.000000]");
  const Detector ch = Detector::scintillator(14.0);
  ZAKU_CHECK_EQ(ch.neutTarget().numNucleons, 12);
  ZAKU_CHECK_EQ(ch.neutTarget().numBoundProtons, 6);
  ZAKU_CHECK_EQ(ch.neutTarget().numFreeProtons, 1);
  ZAKU_CHECK_STR_EQ(Detector::fromMaterial("Water", 1.0).name(), "water");
  ZAKU_CHECK_STR_EQ(Detector::fromMaterial("LAr", 1.0).name(), "argon");
  ZAKU_CHECK_THROWS(Detector::fromMaterial("cheese", 1.0));

  // overrides
  Detector custom("custom", 10.0, {Target(8, 16, 0.5), Target(6, 12, 0.5)});
  ZAKU_CHECK_THROWS(custom.neutTarget());  // two bound nuclei
  custom.setNeutTarget(NeutTarget{16, 8, 0});
  ZAKU_CHECK_EQ(custom.neutTarget().numNucleons, 16);
  custom.setNucleonsOverride(1.0e30);
  ZAKU_CHECK_CLOSE(custom.nNucleons(), 1.0e30, 1e-12);
  custom.setNucleonsOverride(0.0);
  ZAKU_CHECK_CLOSE(custom.nNucleons(), 10.0e9 * 6.02214076e23, 1e-9);

  // validation
  ZAKU_CHECK_THROWS(Detector("bad", -1.0, {Target(8, 16, 1.0)}));
  ZAKU_CHECK_THROWS(Detector("bad", 1.0, {}));
  ZAKU_CHECK_THROWS(Detector("bad", 1.0, {Target(8, 16, 0.5)}));  // fractions do not sum to 1
  ZAKU_CHECK_THROWS(Detector("bad", 1.0, {Target(20, 16, 1.0)})); // Z > A

  ZAKU_CHECK(hk.describe().find("188") != std::string::npos);

  return zakutest::finish();
}
