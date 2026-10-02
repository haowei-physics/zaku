// Regression: the recon stage with the legacy configuration reproduces the
// reference _fakerecon.root files bit for bit (weights, category, smeared
// kinematics, decay-electron tags) for one disappearance and one appearance
// production.
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <cstdlib>
#include <memory>

#include "ZakuTest.h"
#include "ZakuTestPaths.h"
#include "zaku/core/Experiment.h"
#include "zaku/stages/Stages.h"

using namespace zaku;

namespace {

void link(const std::string& target, const std::string& linkPath) {
  const std::string dir = gSystem->GetDirName(linkPath.c_str()).Data();
  gSystem->mkdir(dir.c_str(), true);
  gSystem->Unlink(linkPath.c_str());
  if (gSystem->Symlink(target.c_str(), linkPath.c_str()) != 0) {
    throw std::runtime_error("cannot create symlink " + linkPath);
  }
}

struct ReconCompare {
  long long entries = 0, exactP = 0, exactCos = 0, exactME = 0, exactCat = 0;
  double maxPot = 0, maxOsc = 0, maxEff = 0;
};

ReconCompare compare(const std::string& mine, const std::string& ref) {
  std::unique_ptr<TFile> fm(TFile::Open(mine.c_str(), "READ"));
  std::unique_ptr<TFile> fr(TFile::Open(ref.c_str(), "READ"));
  TTree* tm = fm->Get<TTree>("FlatTree_VARS");
  TTree* tr = fr->Get<TTree>("FlatTree_VARS");
  ReconCompare c;
  if (!tm || !tr) return c;
  c.entries = tr->GetEntries();
  if (tm->GetEntries() != tr->GetEntries()) return c;
  Int_t catM, catR, meM, meR, tmeM, tmeR;
  Double_t potM, potR, oscM, oscR, effM, effR;
  Float_t pM, pR, cM, cR;
  tm->SetBranchAddress("RespCat", &catM);
  tr->SetBranchAddress("RespCat", &catR);
  tm->SetBranchAddress("recon_ME", &meM);
  tr->SetBranchAddress("recon_ME", &meR);
  tm->SetBranchAddress("true_ME", &tmeM);
  tr->SetBranchAddress("true_ME", &tmeR);
  tm->SetBranchAddress("POTWeight", &potM);
  tr->SetBranchAddress("POTWeight", &potR);
  tm->SetBranchAddress("OscWeight", &oscM);
  tr->SetBranchAddress("OscWeight", &oscR);
  tm->SetBranchAddress("SelEffWeight", &effM);
  tr->SetBranchAddress("SelEffWeight", &effR);
  tm->SetBranchAddress("PLep_rec", &pM);
  tr->SetBranchAddress("PLep_rec", &pR);
  tm->SetBranchAddress("CosLep_rec", &cM);
  tr->SetBranchAddress("CosLep_rec", &cR);
  auto rel = [](double a, double b) { return std::fabs(a - b) / (std::fabs(b) > 0 ? std::fabs(b) : 1.0); };
  for (long long i = 0; i < c.entries; ++i) {
    tm->GetEntry(i);
    tr->GetEntry(i);
    if (catM == catR) ++c.exactCat;
    if (meM == meR && tmeM == tmeR) ++c.exactME;
    if (pM == pR) ++c.exactP;
    if (cM == cR) ++c.exactCos;
    c.maxPot = std::max(c.maxPot, rel(potM, potR));
    c.maxOsc = std::max(c.maxOsc, rel(oscM, oscR));
    c.maxEff = std::max(c.maxEff, rel(effM, effR));
  }
  return c;
}

}  // namespace

int run() {
  const std::string ref = zakutest::referenceDir();
  const std::string src = zakutest::sourceDir();
  zakutest::requirePath(ref + "/numu/numu_x_numu_200K_FHC_flattree.root", "reference production");
  setenv("ZAKU_DATA", src.c_str(), 1);

  Experiment exp = Experiment::fromConfig(src + "/configs/legacy_hk_mock.yaml");
  const std::string work = zakutest::scratchDir() + "/legacy_recon";
  exp.setRun(work, exp.run().seed);
  exp.setLegacyRngOrder(true);

  for (const char* name : {"numu_x_numu_200K_FHC", "numu_x_nue_50K_FHC"}) {
    const Production& prod = exp.production(name);
    link(ref + "/" + Layout::flatFile(prod), exp.path(Layout::flatFile(prod)));
    const stages::ReconSummary s = stages::reconstruct(exp, prod, /*force=*/true);
    std::printf("%s\n", s.describe().c_str());
    ZAKU_CHECK_EQ(static_cast<int>(s.seed), 1234);
    const ReconCompare c = compare(exp.path(Layout::reconFile(prod)), ref + "/" + Layout::reconFile(prod));
    const long long expected = std::string(name) == "numu_x_numu_200K_FHC" ? 198495 : 49588;
    ZAKU_CHECK_EQ(static_cast<int>(s.kept), static_cast<int>(expected));
    ZAKU_CHECK_EQ(static_cast<int>(c.entries), static_cast<int>(expected));
    ZAKU_CHECK_MSG(c.exactCat == c.entries, std::string(name) + ": RespCat exact in " + std::to_string(c.exactCat));
    ZAKU_CHECK_MSG(c.exactME == c.entries, std::string(name) + ": decay-electron tags exact in " + std::to_string(c.exactME));
    ZAKU_CHECK_MSG(c.exactP == c.entries, std::string(name) + ": PLep_rec bit-identical in " + std::to_string(c.exactP));
    ZAKU_CHECK_MSG(c.exactCos == c.entries, std::string(name) + ": CosLep_rec bit-identical in " + std::to_string(c.exactCos));
    ZAKU_CHECK_MSG(c.maxPot < 1e-12, std::string(name) + ": POTWeight max rel diff " + std::to_string(c.maxPot));
    ZAKU_CHECK_MSG(c.maxOsc < 1e-12, std::string(name) + ": OscWeight max rel diff " + std::to_string(c.maxOsc));
    ZAKU_CHECK_MSG(c.maxEff < 1e-12, std::string(name) + ": SelEffWeight max rel diff " + std::to_string(c.maxEff));
  }
  // the second run is skipped unless forced
  const stages::ReconSummary again = stages::reconstruct(exp, exp.production("numu_x_nue_50K_FHC"), false);
  ZAKU_CHECK(again.skipped);
  return zakutest::finish();
}

int main() {
  try {
    return run();
  } catch (const std::exception& e) {
    std::fprintf(stderr, "EXCEPTION: %s\n", e.what());
    return 1;
  }
}
