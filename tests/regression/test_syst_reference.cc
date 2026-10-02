// Regression: the splines stage reproduces every response graph of the
// reference _spline.root; the covariances stage reproduces the four reference
// parameter files and the two binning files; select on the reference spline
// files followed by the summary stage gives 77 parameters with finite effects.
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <sstream>

#include "ZakuTest.h"
#include "ZakuTestPaths.h"
#include "zaku/core/Experiment.h"
#include "zaku/stages/Stages.h"
#include "zaku/syst/Covariances.h"

using namespace zaku;

namespace {

void link(const std::string& target, const std::string& linkPath) {
  const std::string dir = gSystem->GetDirName(linkPath.c_str()).Data();
  gSystem->mkdir(dir.c_str(), true);
  gSystem->Unlink(linkPath.c_str());
  if (gSystem->Symlink(target.c_str(), linkPath.c_str()) != 0) throw std::runtime_error("cannot symlink " + linkPath);
}

struct GraphCompare {
  long long entries = 0, refEntries = 0, graphsCompared = 0, sameN = 0, samePoints = 0;
  double maxDiff = 0.0;
};

GraphCompare compareGraphs(const std::string& mine, const std::string& ref, const std::vector<std::string>& dials) {
  std::unique_ptr<TFile> fm(TFile::Open(mine.c_str(), "READ"));
  std::unique_ptr<TFile> fr(TFile::Open(ref.c_str(), "READ"));
  TTree* tm = fm->Get<TTree>("FlatTree_VARS");
  TTree* tr = fr->Get<TTree>("FlatTree_VARS");
  GraphCompare c;
  c.entries = tm->GetEntries();
  c.refEntries = tr->GetEntries();
  std::vector<TGraph*> gm(dials.size(), nullptr), gr(dials.size(), nullptr);
  for (std::size_t d = 0; d < dials.size(); ++d) {
    tm->SetBranchAddress(dials[d].c_str(), &gm[d]);
    tr->SetBranchAddress(dials[d].c_str(), &gr[d]);
  }
  const long long n = std::min(c.entries, c.refEntries);
  for (long long i = 0; i < n; ++i) {
    tm->GetEntry(i);
    tr->GetEntry(i);
    for (std::size_t d = 0; d < dials.size(); ++d) {
      ++c.graphsCompared;
      if (gm[d]->GetN() != gr[d]->GetN()) continue;
      ++c.sameN;
      bool same = true;
      for (int k = 0; k < gm[d]->GetN(); ++k) {
        const double dx = std::fabs(gm[d]->GetX()[k] - gr[d]->GetX()[k]);
        const double dy = std::fabs(gm[d]->GetY()[k] - gr[d]->GetY()[k]);
        c.maxDiff = std::max(c.maxDiff, std::max(dx, dy));
        if (dx > 1e-12 || dy > 1e-12) same = false;
      }
      if (same) ++c.samePoints;
    }
  }
  return c;
}

}  // namespace

int run() {
  const std::string ref = zakutest::referenceDir();
  const std::string src = zakutest::sourceDir();
  zakutest::requirePath(ref + "/numu/numu_x_numu_200K_FHC_spline.root", "reference spline file");
  zakutest::requirePath(ref + "/parameters/xsec/xsecCovHK.root", "reference parameter files");
  setenv("ZAKU_DATA", src.c_str(), 1);

  Experiment exp = Experiment::fromConfig(src + "/configs/legacy_hk_mock.yaml");
  const std::string work = zakutest::scratchDir() + "/legacy_syst";
  exp.setRun(work, exp.run().seed);
  exp.setLegacyRngOrder(true);

  // --- splines on numu_x_numu_200K_FHC (Zaku recon of the reference flat tree) ------------
  const Production& prod = exp.production("numu_x_numu_200K_FHC");
  for (const auto& p : exp.productions()) {
    zakutest::requirePath(ref + "/" + Layout::flatFile(p), "reference flat tree " + p.name());
    zakutest::requirePath(ref + "/" + Layout::weightFile(p), "reference weight file " + p.name());
    link(ref + "/" + Layout::flatFile(p), exp.path(Layout::flatFile(p)));
    link(ref + "/" + Layout::weightFile(p), exp.path(Layout::weightFile(p)));
    // a leftover from an earlier run may be a link to the legacy recon file: start clean
    gSystem->Unlink(exp.path(Layout::reconFile(p)).c_str());
    gSystem->Unlink(exp.path(Layout::splineFile(p)).c_str());
  }
  std::printf("%s\n", stages::reconstruct(exp, prod, false).describe().c_str());
  const stages::SplineSummary ss = stages::makeSplines(exp, prod, /*force=*/true);
  std::printf("%s\n", ss.describe().c_str());
  ZAKU_CHECK_EQ(static_cast<int>(ss.entries), 198495);
  ZAKU_CHECK_EQ(static_cast<int>(ss.dials.size()), 19);
  ZAKU_CHECK_STR_EQ(ss.dials[0], "MaCCQE");
  ZAKU_CHECK_STR_EQ(ss.dials[18], "RadCorr_rate");
  const GraphCompare gc = compareGraphs(ss.outputFile, ref + "/" + Layout::splineFile(prod), ss.dials);
  ZAKU_CHECK_EQ(static_cast<int>(gc.entries), static_cast<int>(gc.refEntries));
  ZAKU_CHECK_MSG(gc.sameN == gc.graphsCompared, "graphs with the same point count: " + std::to_string(gc.sameN) + " of " + std::to_string(gc.graphsCompared));
  ZAKU_CHECK_MSG(gc.samePoints == gc.graphsCompared, "graphs identical: " + std::to_string(gc.samePoints) + " of " + std::to_string(gc.graphsCompared) + ", max diff " + std::to_string(gc.maxDiff));
  // skipped unless forced
  ZAKU_CHECK(stages::makeSplines(exp, prod, false).skipped);

  // --- covariances -------------------------------------------------------------------------
  const stages::CovarianceSummary cs = stages::makeCovariances(exp, /*force=*/true);
  std::printf("%s\n", cs.describe().c_str());
  ZAKU_CHECK_EQ(static_cast<int>(cs.files.size()), 8);
  std::string why;
  const ParameterSet refOsc = ParameterSet::read(ref + "/parameters/oscProb/oscCovHK.root", "osc");
  const ParameterSet myOsc = ParameterSet::read(exp.path(Layout::oscCovarianceFile()), "osc");
  ZAKU_CHECK_MSG(myOsc.sameAs(refOsc, 1e-9, &why), "osc: " + why);
  const ParameterSet refDet = ParameterSet::read(ref + "/parameters/det/detCovHK.root", "det");
  const ParameterSet myDet = ParameterSet::read(exp.path(Layout::detCovarianceFile()), "det");
  ZAKU_CHECK_MSG(myDet.sameAs(refDet, 1e-9, &why), "det: " + why);
  const ParameterSet refXsec = ParameterSet::read(ref + "/parameters/xsec/xsecCovHK.root", "xsec");
  const ParameterSet myXsec = ParameterSet::read(exp.path(Layout::xsecCovarianceFile()), "xsec");
  ZAKU_CHECK_MSG(myXsec.sameAs(refXsec, 1e-9, &why), "xsec: " + why);
  ZAKU_CHECK(myXsec.anyBounds());
  const ParameterSet refFlux = ParameterSet::read(ref + "/parameters/flux/fluxCovHK.root", "flux");
  const ParameterSet myFlux = ParameterSet::read(exp.path(Layout::fluxCovarianceFile(false)), "flux");
  ZAKU_CHECK_MSG(myFlux.sameAs(refFlux, 1e-9, &why), "flux: " + why);
  const ParameterSet refFluxND = ParameterSet::read(ref + "/parameters/flux/fluxCovHK_ND.root", "flux");
  const ParameterSet myFluxND = ParameterSet::read(exp.path(Layout::fluxCovarianceFile(true)), "flux");
  ZAKU_CHECK_MSG(myFluxND.sameAs(refFluxND, 1e-9, &why), "flux ND: " + why);
  double lo = 1.0, hi = 0.0;
  for (std::size_t i = 0; i < myFluxND.size(); ++i) {
    lo = std::min(lo, myFluxND.sigma(static_cast<int>(i)));
    hi = std::max(hi, myFluxND.sigma(static_cast<int>(i)));
  }
  ZAKU_CHECK_CLOSE(lo, 0.0214, 5e-3);
  ZAKU_CHECK_CLOSE(hi, 0.0744, 5e-3);
  ZAKU_CHECK(Binning::fromFile(exp.path(Layout::detBinningFile())).sameAs(Binning::fromFile(ref + "/parameters/det/detBinning.txt"), 1e-9));
  ZAKU_CHECK(Binning::fromFile(exp.path(Layout::fluxBinningFile())).sameAs(Binning::fromFile(ref + "/parameters/flux/fluxBinning.txt"), 1e-9));
  {
    std::unique_ptr<TFile> g(TFile::Open(exp.path(Layout::oscGridFile()).c_str(), "READ"));
    auto* h = g->Get<TH1>("energyBinning");
    ZAKU_CHECK(h != nullptr);
    if (h) {
      ZAKU_CHECK_EQ(h->GetNbinsX(), 2000);
      ZAKU_CHECK_CLOSE(h->GetXaxis()->GetXmax(), 10.0, 1e-12);
    }
  }
  ZAKU_CHECK(stages::makeCovariances(exp, false).skipped);

  // --- recon + splines for every production, select, then the summary ----------------------------
  // (the reference spline files are legacy files without Zaku's recon branches, so the chain is rebuilt)
  for (const auto& p : exp.productions()) {
    std::printf("%s\n", stages::reconstruct(exp, p, false).describe().c_str());
    const stages::SplineSummary s = stages::makeSplines(exp, p, false);
    std::printf("splines %s: %s\n", p.name().c_str(), s.skipped ? "skipped" : std::to_string(s.entries).c_str());
  }
  for (BeamMode mode : exp.beamModes()) {
    const stages::SelectSummary s = stages::select(exp, mode, true);
    std::printf("%s\n", s.describe().c_str());
  }
  const stages::SummaryResult sr = stages::summarise(exp, true);
  std::printf("%s\n", sr.describe().c_str());
  const SummaryTable& t = sr.table;
  ZAKU_CHECK_EQ(static_cast<int>(t.parameters().size()), 77);
  ZAKU_CHECK_EQ(static_cast<int>(t.samples().size()), 5);
  ZAKU_CHECK_EQ(static_cast<int>(t.size()), 77 * 5);
  bool finite = true;
  for (const auto& e : t.entries()) finite = finite && std::isfinite(e.plus) && std::isfinite(e.minus) && e.nominal > 0.0;
  ZAKU_CHECK(finite);
  // detector normalisations move their own sample by exactly +-sigma and no other
  for (std::size_t i = 0; i < myDet.size(); ++i) {
    const std::string& par = myDet.parameter(static_cast<int>(i)).name;
    for (std::size_t s = 0; s < exp.samples().size(); ++s) {
      const SummaryEntry* e = t.find(par, exp.samples()[s].name);
      ZAKU_CHECK(e != nullptr);
      if (!e) continue;
      const double expected = s == i ? myDet.sigma(static_cast<int>(i)) : 0.0;
      ZAKU_CHECK_MSG(std::fabs(e->plusFrac() - expected) < 1e-9, par + " on " + e->sample + ": " + std::to_string(e->plusFrac()));
      ZAKU_CHECK_MSG(std::fabs(e->minusFrac() + expected) < 1e-9, par + " on " + e->sample + " minus");
    }
  }
  // flux normalisations never exceed their sigma; the numu FHC peak bin moves FHC_1Rmu
  for (std::size_t i = 0; i < myFluxND.size(); ++i) {
    for (const auto& smp : t.samples()) {
      const SummaryEntry* e = t.find(myFluxND.parameter(static_cast<int>(i)).name, smp);
      ZAKU_CHECK(e && e->maxFrac() <= myFluxND.sigma(static_cast<int>(i)) + 1e-9);
    }
  }
  // the bin above the oscillation dip (0.7-1 GeV) moves the muon sample; antineutrino-mode bins do not touch it
  const SummaryEntry* above = t.find("FLUX_FHC_NUMU_04", "FHC_1Rmu");
  ZAKU_CHECK(above && above->plusFrac() > 0.003 && std::fabs(above->plusFrac() + above->minusFrac()) < 1e-12);
  const SummaryEntry* other = t.find("FLUX_RHC_NUMUB_00", "FHC_1Rmu");
  ZAKU_CHECK(other && other->plusFrac() == 0.0 && other->minusFrac() == 0.0);
  double fhcNumuTotal = 0.0;
  for (int b = 0; b <= 10; ++b) {
    char name[32];
    std::snprintf(name, sizeof name, "FLUX_FHC_NUMU_%02d", b);
    const SummaryEntry* e = t.find(name, "FHC_1Rmu");
    if (e) fhcNumuTotal += e->plusFrac();
  }
  ZAKU_CHECK_MSG(fhcNumuTotal > 0.02 && fhcNumuTotal < 0.05, "sum of FHC numu flux effects " + std::to_string(fhcNumuTotal));
  // the nominal equals the rate table of the select stage
  const RateTable rates = RateTable::fromCsv([&] {
    std::ifstream in(exp.path(Layout::rateTableFile(BeamMode::FHC)));
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
  }());
  const SummaryEntry* mu = t.find("MaCCQE", "FHC_1Rmu");
  ZAKU_CHECK(mu != nullptr);
  if (mu) {
    // the summary counts events inside the sample's binning range only (no Enu_rec = -1, nothing above 30 GeV)
    ZAKU_CHECK(mu->nominal <= rates.rowTotal(0) + 1e-9);
    ZAKU_CHECK_CLOSE(mu->nominal, rates.rowTotal(0), 1e-3);
    ZAKU_CHECK(mu->plusFrac() > 0.05 && mu->minusFrac() < -0.05);  // MaCCQE up = more CCQE
  }
  const SummaryEntry* dcp = t.find("PMNS_DELTA_CP", "FHC_1Re0de");
  ZAKU_CHECK(dcp && std::fabs(dcp->plusFrac()) > 0.01);   // delta_CP moves the appearance sample
  const SummaryEntry* sign = t.find("PMNS_SIGN_MASS_SQUARED_32", "FHC_1Rmu");
  ZAKU_CHECK(sign && std::fabs(sign->plusFrac()) < 1e-12);  // +-0.1 does not flip the ordering
  ZAKU_CHECK(zakutest::exists(sr.rootFile));
  ZAKU_CHECK(stages::summarise(exp, false).skipped);
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
