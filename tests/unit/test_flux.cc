#include <TAxis.h>
#include <TFile.h>
#include <TH1D.h>
#include <TMatrixD.h>

#include <memory>

#include "ZakuTest.h"
#include "ZakuTestPaths.h"
#include "zaku/core/FluxSpec.h"
#include "zaku/flux/FluxCovarianceInput.h"
#include "zaku/flux/FluxHistograms.h"
#include "zaku/flux/FluxPreparer.h"

using namespace zaku;

namespace {

// A T2K-like variable-bin density: 10 MeV bins to 1 GeV, then 100 MeV bins to 3 GeV.
std::unique_ptr<TH1D> makeVariableFlux(const char* name, double scale) {
  std::vector<double> edges;
  for (int i = 0; i <= 100; ++i) edges.push_back(0.01 * i);
  for (int i = 1; i <= 20; ++i) edges.push_back(1.0 + 0.1 * i);
  auto h = std::make_unique<TH1D>(name, "flux;E_{#nu} [GeV];/cm^{2}/50MeV/1e21POT", static_cast<int>(edges.size()) - 1,
                                  edges.data());
  h->SetDirectory(nullptr);
  for (int i = 1; i <= h->GetNbinsX(); ++i) h->SetBinContent(i, scale * (1.0 + i));  // density, GeV-independent test values
  return h;
}

}  // namespace

int main() {
  const std::string scratch = zakutest::scratchDir();

  // --- FluxSpec ----------------------------------------------------------
  FluxSpec t2k = FluxSpec::t2k("f.root", "sk");
  ZAKU_CHECK_STR_EQ(t2k.histogram(14), "enu_sk_numu");
  ZAKU_CHECK_STR_EQ(t2k.histogram(-14), "enu_sk_numub");
  ZAKU_CHECK_STR_EQ(t2k.histogram(-12), "enu_sk_nueb");
  ZAKU_CHECK_EQ(static_cast<int>(t2k.pdgs().size()), 4);
  ZAKU_CHECK_THROWS(t2k.histogram(16));
  ZAKU_CHECK_STR_EQ(t2k.withFile("g.root").file, "g.root");
  ZAKU_CHECK(t2k.describe().find("enu_sk_numu") != std::string::npos);

  // --- a synthetic GeV flux file ------------------------------------------
  const std::string fluxFile = scratch + "/flux_gev.root";
  {
    TFile f(fluxFile.c_str(), "RECREATE");
    auto hmu = makeVariableFlux("enu_sk_numu", 1.0e5);
    auto he = makeVariableFlux("enu_sk_nue", 1.0e3);
    hmu->Write();
    he->Write();
  }
  FluxSpec spec("", {{14, "enu_sk_numu"}, {12, "enu_sk_nue"}});
  spec.file = fluxFile;
  FluxHistograms fh(spec);
  ZAKU_CHECK(fh.has(14));
  ZAKU_CHECK(!fh.has(-14));
  ZAKU_CHECK_THROWS(fh.histogram(-14));
  const double integral = fh.integral(14);
  ZAKU_CHECK_CLOSE(fh.integralPerPOT(14), integral / 0.05 / 1.0e21, 1e-12);
  // per-bin-content convention
  FluxSpec specCounts = spec;
  specCounts.binNormGeV = 0.0;
  FluxHistograms fhc(specCounts);
  ZAKU_CHECK_CLOSE(fhc.integralPerPOT(14), fh.histogram(14).Integral() / 1.0e21, 1e-12);

  // --- uniform rebinning preserves the density --------------------------
  const std::string uniformFile = scratch + "/flux_uniform.root";
  FluxSpec uspec = FluxPreparer::makeUniform(spec, uniformFile, 0.01, true);
  ZAKU_CHECK_STR_EQ(uspec.file, uniformFile);
  FluxHistograms fu(uspec);
  const TH1D& hu = fu.histogram(14);
  ZAKU_CHECK_EQ(hu.GetNbinsX(), 300);
  ZAKU_CHECK_CLOSE(hu.GetXaxis()->GetXmax(), 3.0, 1e-12);
  ZAKU_CHECK_CLOSE(fu.integral(14), integral, 1e-12);
  ZAKU_CHECK_CLOSE(fu.integralPerPOT(14), fh.integralPerPOT(14), 1e-12);
  // a wide original bin is copied into all the uniform bins it covers
  const TH1D& ho = fh.histogram(14);
  ZAKU_CHECK_CLOSE(hu.GetBinContent(hu.GetXaxis()->FindFixBin(1.55)), ho.GetBinContent(ho.GetXaxis()->FindFixBin(1.55)), 1e-12);
  ZAKU_CHECK_CLOSE(hu.GetBinContent(hu.GetXaxis()->FindFixBin(1.59)), ho.GetBinContent(ho.GetXaxis()->FindFixBin(1.55)), 1e-12);
  // the nue histogram was converted too (allHistograms = true)
  ZAKU_CHECK(fu.has(12));

  std::unique_ptr<TH1D> rb(FluxPreparer::rebinUniform(ho, 0.5, "coarse"));
  ZAKU_CHECK_EQ(rb->GetNbinsX(), 6);
  ZAKU_CHECK_STR_EQ(rb->GetName(), "coarse");
  ZAKU_CHECK_THROWS(FluxPreparer::rebinUniform(ho, -1.0));

  // --- MeV inputs: declared -> converted, undeclared -> rejected ----------
  const std::string mevFile = scratch + "/flux_mev.root";
  {
    TFile f(mevFile.c_str(), "RECREATE");
    TH1D h("enu_sk_numu", "flux;E [MeV]", 300, 0.0, 3000.0);
    for (int i = 1; i <= 300; ++i) h.SetBinContent(i, 7.0);  // density per MeV
    h.Write();
  }
  FluxSpec mev("", {{14, "enu_sk_numu"}});
  mev.file = mevFile;
  ZAKU_CHECK_THROWS(FluxHistograms undeclared(mev));
  mev.units = EnergyUnit::MeV;
  FluxHistograms fm(mev);
  ZAKU_CHECK_CLOSE(fm.histogram(14).GetXaxis()->GetXmax(), 3.0, 1e-12);
  ZAKU_CHECK_CLOSE(fm.histogram(14).GetBinContent(5), 7000.0, 1e-12);        // per GeV
  ZAKU_CHECK_CLOSE(fm.integral(14), 7.0 * 3000.0, 1e-9);                     // density preserved

  // --- missing histogram / file --------------------------------------------
  FluxSpec bad("", {{14, "nope"}});
  bad.file = fluxFile;
  ZAKU_CHECK_THROWS(FluxHistograms b1(bad));
  FluxSpec nofile("", {{14, "enu_sk_numu"}});
  nofile.file = scratch + "/does_not_exist.root";
  ZAKU_CHECK_THROWS(FluxHistograms b2(nofile));

  // --- covariance input on a synthetic file ------------------------------
  const std::string covFile = scratch + "/flux_cov.root";
  {
    TFile f(covFile.c_str(), "RECREATE");
    TMatrixD m(5, 5);
    for (int i = 0; i < 5; ++i) {
      for (int j = 0; j < 5; ++j) m(i, j) = (i == j) ? 0.01 * (i + 1) : 0.001;
    }
    m.Write("total_flux_cov");
    const double nd[] = {0.0, 1.0, 30.0};
    TAxis a1(2, nd);
    a1.Write("nd280_numode_numu_bins");
    const double sk[] = {0.0, 0.5, 1.0, 30.0};
    TAxis a2(3, sk);
    a2.Write("sk_numode_numu_bins");
  }
  FluxCovarianceSpec cs;
  cs.file = covFile;
  FluxCovarianceInput cov(cs);
  ZAKU_CHECK_EQ(static_cast<int>(cov.blocks().size()), 2);
  ZAKU_CHECK_EQ(cov.firstRow("sk"), 2);
  ZAKU_CHECK_EQ(cov.nRows("sk"), 3);
  ZAKU_CHECK_EQ(cov.nRows("nd280"), 2);
  ZAKU_CHECK(cov.blocksFor("sk")[0].beamMode == BeamMode::FHC);
  ZAKU_CHECK_EQ(cov.blocksFor("sk")[0].pdg, 14);
  ZAKU_CHECK_CLOSE(cov.blocksFor("sk")[0].edges[1], 0.5, 1e-12);
  const TMatrixD skm = cov.subMatrix("sk");
  ZAKU_CHECK_EQ(skm.GetNrows(), 3);
  ZAKU_CHECK_CLOSE(skm(0, 0), 0.03, 1e-12);
  const TMatrixD cross = cov.crossMatrix("nd280", "sk");
  ZAKU_CHECK_EQ(cross.GetNrows(), 2);
  ZAKU_CHECK_EQ(cross.GetNcols(), 3);
  ZAKU_CHECK_CLOSE(cross(0, 0), 0.001, 1e-12);
  ZAKU_CHECK_THROWS(cov.blocksFor("dune"));
  ZAKU_CHECK(cov.describe().find("sk FHC numu") != std::string::npos);

  return zakutest::finish();
}
