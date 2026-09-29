// Regression against the reference production: the flat<->weight join and the
// POT weight on numu_x_numu_200K_FHC.
#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <memory>

#include "ZakuTest.h"
#include "ZakuTestPaths.h"
#include "zaku/core/Experiment.h"
#include "zaku/event/EventJoin.h"
#include "zaku/event/FlatTreeCursor.h"
#include "zaku/event/WeightCursor.h"

using namespace zaku;

int run();

int main() {
  try {
    return run();
  } catch (const std::exception& e) {
    std::fprintf(stderr, "EXCEPTION: %s\n", e.what());
    return 1;
  }
}

int run() {
  const std::string ref = zakutest::referenceDir();
  const std::string flatFile = ref + "/numu/numu_x_numu_200K_FHC_flattree.root";
  const std::string weightFile = ref + "/numu/numu_x_numu_200K_FHC_Weight.flat.root";
  const std::string reconFile = ref + "/numu/numu_x_numu_200K_FHC_fakerecon.root";
  zakutest::requirePath(flatFile, "reference flat tree");
  zakutest::requirePath(weightFile, "reference weight file");

  // --- weight metadata ------------------------------------------------------
  ZAKU_CHECK(T2KRWWeightCursor::looksLikeT2KRW(weightFile));
  T2KRWWeightCursor w(weightFile);
  ZAKU_CHECK_EQ(static_cast<int>(w.entries()), 198495);
  ZAKU_CHECK_EQ(static_cast<int>(w.dials().size()), 19);
  ZAKU_CHECK_EQ(w.dialIndex("MaCCQE"), 0);
  ZAKU_CHECK_EQ(w.dialIndex("RadCorr_rate"), 18);
  ZAKU_CHECK_EQ(w.dialIndex("nope"), -1);
  const DialKnots& ma = w.dials()[0];
  ZAKU_CHECK_CLOSE(ma.nominal, 1.21, 1e-6);
  ZAKU_CHECK_EQ(static_cast<int>(ma.knots.size()), 4);
  ZAKU_CHECK_CLOSE(ma.knots[0], 0.818, 1e-6);
  ZAKU_CHECK_CLOSE(ma.knots[3], 1.602, 1e-6);
  ZAKU_CHECK(!ma.nominalIsKnot());
  ZAKU_CHECK_EQ(static_cast<int>(ma.knotsWithNominal().size()), 5);
  ZAKU_CHECK_CLOSE(ma.knotsWithNominal()[2], 1.21, 1e-6);
  const DialKnots& rad = w.dials()[static_cast<std::size_t>(w.dialIndex("RadCorr_rate"))];
  ZAKU_CHECK(rad.nominalIsKnot());
  ZAKU_CHECK_EQ(static_cast<int>(rad.knotsWithNominal().size()), 3);
  const DialKnots& src = w.dials()[static_cast<std::size_t>(w.dialIndex("SRCFrac_O"))];
  ZAKU_CHECK_CLOSE(src.nominal, 0.34, 1e-6);
  ZAKU_CHECK_EQ(static_cast<int>(src.knots.size()), 2);

  ZAKU_CHECK(w.next());
  ZAKU_CHECK_EQ(static_cast<int>(w.entry()), 0);
  ZAKU_CHECK_EQ(static_cast<int>(w.weights(0).size()), 4);
  ZAKU_CHECK_EQ(static_cast<int>(w.weights("MDLSPiEj").size()), 2);
  ZAKU_CHECK(w.Enu() > 0.0 && w.Enu() < 30.0);  // GeV
  ZAKU_CHECK_CLOSE(w.cvWeight(), 1.0, 1e-9);
  ZAKU_CHECK_CLOSE(w.scaleWeight(), 4.456765e-44, 1e-5);
  ZAKU_CHECK(w.goTo(0));

  // --- the join -------------------------------------------------------------
  FlatTreeCursor flat(flatFile);
  ZAKU_CHECK_EQ(static_cast<int>(flat.entries()), 200000);
  EventJoin join(flat, &w);
  long long n = 0, withWeights = 0, dropped = 0;
  double firstScale = 0.0;
  while (join.next()) {
    const FlatEvent& ev = join.event();
    if (n == 0) firstScale = ev.fScaleFactor;
    if (join.hasWeights()) {
      ++withWeights;
      // the matched row carries this event's dial weights
      ZAKU_CHECK_EQ(static_cast<int>(join.weights().weights(0).size()), 4);
    } else {
      ++dropped;
      ZAKU_CHECK(ev.emptyFinalState());
      ZAKU_CHECK(ev.legacyDrop());
    }
    ++n;
  }
  const JoinStats& st = join.stats();
  ZAKU_CHECK_EQ(static_cast<int>(n), 200000);
  ZAKU_CHECK_EQ(static_cast<int>(withWeights), 198495);
  ZAKU_CHECK_EQ(static_cast<int>(dropped), 1505);
  ZAKU_CHECK_EQ(static_cast<int>(st.matched), 198495);
  ZAKU_CHECK_EQ(static_cast<int>(st.unmatchedFlat), 1505);
  ZAKU_CHECK_EQ(static_cast<int>(st.unmatchedEmptyFinalState), 1505);
  ZAKU_CHECK_EQ(static_cast<int>(st.leptonMismatch), 13);
  ZAKU_CHECK(st.maxEnergyRelDiff < 1e-6);
  ZAKU_CHECK_NOTHROW(st.check());
  std::printf("%s\n", st.describe().c_str());

  // --- POT weight with the legacy nucleon count reproduces the stored branch --
  Experiment exp("HK");
  Detector det = Detector::water(188.0);
  det.setNucleonsOverride(1.132e35);
  exp.setDetector(det);
  const std::string legacy = zakutest::legacyDir();
  FluxSpec spec = FluxSpec::t2k(legacy + "/lbl_config/hk/t2kflux_2020_plus250kA_uniform.root", "sk");
  if (zakutest::exists(spec.file)) {
    exp.setFlux(BeamMode::FHC, spec, 6.75e21);
    const double pw = exp.potWeight(BeamMode::FHC, 14, firstScale);
    ZAKU_CHECK_CLOSE(pw, 0.411157, 2e-6);
    if (zakutest::exists(reconFile)) {
      std::unique_ptr<TFile> f(TFile::Open(reconFile.c_str(), "READ"));
      TTree* t = f->Get<TTree>("FlatTree_VARS");
      Double_t stored = 0.0;
      t->SetBranchAddress("POTWeight", &stored);
      t->GetEntry(0);
      ZAKU_CHECK_CLOSE(pw, stored, 1e-12);
      ZAKU_CHECK_EQ(static_cast<int>(t->GetEntries()), 198495);
    }
  }

  // --- a mismatched pair of files is refused --------------------------------
  const std::string otherFlat = ref + "/nue/nue_x_nue_200K_FHC_flattree.root";
  if (zakutest::exists(otherFlat)) {
    FlatTreeCursor flat2(otherFlat);
    T2KRWWeightCursor w2(weightFile);
    EventJoin bad(flat2, &w2);
    ZAKU_CHECK_THROWS(while (bad.next()) {});
  }

  return zakutest::finish();
}
