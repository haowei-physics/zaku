// Regression: recon + select with the legacy configuration on all twelve
// reference productions reproduce the reference fakeselect files: entry
// counts per EventType, the per-entry assignment and kinematics, and the
// weighted sample rates. Also checks that split mode gives the same expected
// rates as the legacy draw.
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <cstdlib>
#include <map>
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
  if (gSystem->Symlink(target.c_str(), linkPath.c_str()) != 0) throw std::runtime_error("cannot symlink " + linkPath);
}

struct SelectCompare {
  long long entries = 0, refEntries = 0, exactType = 0;
  std::map<int, long long> counts;
  double maxEnu = 0, maxQ2 = 0, maxTheta = 0;
  std::map<int, double> rates;  // per EventType: sum POT*Osc*SelEff
};

SelectCompare compare(const std::string& mine, const std::string& ref) {
  std::unique_ptr<TFile> fm(TFile::Open(mine.c_str(), "READ"));
  std::unique_ptr<TFile> fr(TFile::Open(ref.c_str(), "READ"));
  TTree* tm = fm->Get<TTree>("FakeSelect");
  TTree* tr = fr->Get<TTree>("FakeSelect");
  SelectCompare c;
  c.entries = tm->GetEntries();
  c.refEntries = tr->GetEntries();
  Int_t typeM, typeR, sidxM, sidxR, pdgM, pdgR;
  Float_t enuM, enuR, q2M, q2R, thM, thR;
  Double_t potR, oscR, effR;
  tm->SetBranchAddress("EventType", &typeM);
  tr->SetBranchAddress("EventType", &typeR);
  tm->SetBranchAddress("SampleIndex", &sidxM);
  tr->SetBranchAddress("SampleIndex", &sidxR);
  tm->SetBranchAddress("PDGLep_rec", &pdgM);
  tr->SetBranchAddress("PDGLep_rec", &pdgR);
  tm->SetBranchAddress("Enu_rec", &enuM);
  tr->SetBranchAddress("Enu_rec", &enuR);
  tm->SetBranchAddress("Q2_rec", &q2M);
  tr->SetBranchAddress("Q2_rec", &q2R);
  tm->SetBranchAddress("Theta_rec", &thM);
  tr->SetBranchAddress("Theta_rec", &thR);
  tr->SetBranchAddress("POTWeight", &potR);
  tr->SetBranchAddress("OscWeight", &oscR);
  tr->SetBranchAddress("SelEffWeight", &effR);
  const long long n = std::min(c.entries, c.refEntries);
  auto rel = [](double a, double b) { return std::fabs(a - b) / (std::fabs(b) > 1e-3 ? std::fabs(b) : 1.0); };
  for (long long i = 0; i < n; ++i) {
    tm->GetEntry(i);
    tr->GetEntry(i);
    ++c.counts[typeM];
    if (typeM == typeR && sidxM == sidxR && pdgM == pdgR) ++c.exactType;
    c.maxEnu = std::max(c.maxEnu, rel(enuM, enuR));
    c.maxQ2 = std::max(c.maxQ2, rel(q2M, q2R));
    c.maxTheta = std::max(c.maxTheta, std::fabs(static_cast<double>(thM) - thR));
    if (typeR >= 0) c.rates[typeR] += potR * oscR * effR;
  }
  return c;
}

}  // namespace

int run() {
  const std::string ref = zakutest::referenceDir();
  const std::string src = zakutest::sourceDir();
  zakutest::requirePath(ref + "/fakeselect_FHC.root", "reference fakeselect files");
  setenv("ZAKU_DATA", src.c_str(), 1);

  Experiment exp = Experiment::fromConfig(src + "/configs/legacy_hk_mock.yaml");
  const std::string work = zakutest::scratchDir() + "/legacy_select";
  exp.setRun(work, exp.run().seed);
  exp.setLegacyRngOrder(true);
  for (const auto& prod : exp.productions()) {
    const std::string flat = ref + "/" + Layout::flatFile(prod);
    zakutest::requirePath(flat, "reference flat tree " + prod.name());
    link(flat, exp.path(Layout::flatFile(prod)));
    const stages::ReconSummary s = stages::reconstruct(exp, prod, false);
    std::printf("%s\n", s.describe().c_str());
  }

  // Reference numbers (legacy FakeSelect printout / file contents).
  const std::map<int, long long> expectedCounts = {{-2, 155409 + 162801}, {0, 202115}, {1, 185165}, {2, 52080}, {3, 200898}, {4, 229539}};
  const std::map<int, double> expectedRates = {{0, 8790.2}, {1, 1882.1}, {2, 299.5}, {3, 12028.8}, {4, 1287.7}};

  std::map<int, long long> counts;
  std::map<int, double> myRates;
  for (BeamMode mode : {BeamMode::FHC, BeamMode::RHC}) {
    const stages::SelectSummary s = stages::select(exp, mode, true);
    std::printf("%s\n", s.describe().c_str());
    ZAKU_CHECK_EQ(static_cast<int>(s.seed), 20250731);
    const SelectCompare c = compare(exp.path(Layout::selectFile(mode)), ref + "/" + Layout::selectFile(mode));
    ZAKU_CHECK_MSG(c.entries == c.refEntries, toString(mode) + " entries " + std::to_string(c.entries) + " vs " + std::to_string(c.refEntries));
    ZAKU_CHECK_MSG(c.exactType == c.refEntries, toString(mode) + ": EventType/SampleIndex/PDGLep_rec exact in " + std::to_string(c.exactType));
    ZAKU_CHECK_MSG(c.maxEnu < 1e-6, toString(mode) + ": Enu_rec max rel diff " + std::to_string(c.maxEnu));
    ZAKU_CHECK_MSG(c.maxQ2 < 1e-6, toString(mode) + ": Q2_rec max rel diff " + std::to_string(c.maxQ2));
    ZAKU_CHECK_MSG(c.maxTheta < 1e-4, toString(mode) + ": Theta_rec max diff " + std::to_string(c.maxTheta));
    for (const auto& kv : c.counts) counts[kv.first] += kv.second;
    // our rate table vs the reference weights
    const std::vector<int> idx = exp.sampleIndices(mode);
    for (std::size_t k = 0; k < idx.size(); ++k) {
      const double mineRate = s.rates.rowTotal(static_cast<int>(k));
      myRates[idx[k]] = mineRate;
      const double refRate = c.rates.count(idx[k]) ? c.rates.at(idx[k]) : 0.0;
      ZAKU_CHECK_MSG(zakutest::close(mineRate, refRate, 1e-9),
                     "rate " + exp.samples()[static_cast<std::size_t>(idx[k])].name + " " + std::to_string(mineRate) +
                         " vs reference entries " + std::to_string(refRate));
    }
    ZAKU_CHECK(zakutest::exists(exp.path(Layout::rateTableFile(mode))));
  }
  for (const auto& kv : expectedCounts) {
    ZAKU_CHECK_MSG(counts[kv.first] == kv.second,
                   "EventType " + std::to_string(kv.first) + " count " + std::to_string(counts[kv.first]) + " vs " + std::to_string(kv.second));
  }
  for (const auto& kv : expectedRates) {
    ZAKU_CHECK_MSG(zakutest::close(myRates[kv.first], kv.second, 1e-3),
                   "EventType " + std::to_string(kv.first) + " rate " + std::to_string(myRates[kv.first]) + " vs " + std::to_string(kv.second));
  }

  // Split mode: same expected rates within MC statistics of the draw.
  exp.setAssignment(Assignment::Split);
  for (BeamMode mode : {BeamMode::FHC, BeamMode::RHC}) {
    const stages::SelectSummary s = stages::select(exp, mode, true);
    const std::vector<int> idx = exp.sampleIndices(mode);
    for (std::size_t k = 0; k < idx.size(); ++k) {
      const double split = s.rates.rowTotal(static_cast<int>(k));
      const double draw = myRates[idx[k]];
      ZAKU_CHECK_MSG(zakutest::close(split, draw, 1e-2),
                     "split vs draw " + exp.samples()[static_cast<std::size_t>(idx[k])].name + " " + std::to_string(split) + " vs " + std::to_string(draw));
    }
    std::printf("split %s:\n%s\n", toString(mode).c_str(), s.rates.text().c_str());
  }
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
