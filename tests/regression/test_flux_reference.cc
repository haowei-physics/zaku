// Regression against the legacy inputs: the uniform flux files and the flux
// integrals the legacy macros used.
#include <TFile.h>
#include <TH1.h>

#include <memory>

#include "ZakuTest.h"
#include "ZakuTestPaths.h"
#include "zaku/core/FluxSpec.h"
#include "zaku/flux/FluxCovarianceInput.h"
#include "zaku/flux/FluxHistograms.h"
#include "zaku/flux/FluxPreparer.h"

using namespace zaku;

int main() {
  const std::string release = zakutest::fluxReleaseDir();
  const std::string legacy = zakutest::legacyDir();
  zakutest::requirePath(release + "/t2kflux_2020_plus250kA_nominal.root", "T2K flux release");
  zakutest::requirePath(legacy + "/lbl_config/hk/t2kflux_2020_plus250kA_uniform.root", "legacy uniform flux");
  const std::string scratch = zakutest::scratchDir();

  struct Case {
    const char* nominal;
    const char* uniform;
    BeamMode mode;
  };
  const Case cases[] = {{"t2kflux_2020_plus250kA_nominal.root", "t2kflux_2020_plus250kA_uniform.root", BeamMode::FHC},
                        {"t2kflux_2020_minus250kA_nominal.root", "t2kflux_2020_minus250kA_uniform.root", BeamMode::RHC}};

  // Legacy flux_per_pot values (FakeRecon printout, /cm^2/POT).
  const double expectedFHC[] = {1.20736e-14, 8.12493e-16, 1.43387e-16, 2.23263e-17};
  const double expectedRHC[] = {1.09042e-15, 9.76267e-15, 3.17518e-17, 9.81011e-17};
  const int pdgs[] = {14, -14, 12, -12};

  for (const auto& c : cases) {
    FluxSpec spec = FluxSpec::t2k(release + "/" + c.nominal, "sk");
    const std::string out = scratch + "/uniform_" + toString(c.mode) + ".root";
    FluxSpec uspec = FluxPreparer::makeUniform(spec, out, 0.01, /*allHistograms=*/true);

    std::unique_ptr<TFile> ref(TFile::Open((legacy + "/lbl_config/hk/" + c.uniform).c_str(), "READ"));
    std::unique_ptr<TFile> mine(TFile::Open(out.c_str(), "READ"));
    ZAKU_CHECK(ref && !ref->IsZombie());
    ZAKU_CHECK(mine && !mine->IsZombie());
    for (const char* det : {"sk", "nd280"}) {
      for (int pdg : pdgs) {
        const std::string name = "enu_" + std::string(det) + "_" + NuSpecies::fromPdg(pdg).t2kSuffix();
        TH1* hr = ref->Get<TH1>(name.c_str());
        TH1* hm = mine->Get<TH1>(name.c_str());
        ZAKU_CHECK_MSG(hr != nullptr && hm != nullptr, "histogram " + name);
        if (!hr || !hm) continue;
        ZAKU_CHECK_EQ(hm->GetNbinsX(), hr->GetNbinsX());
        double maxDiff = 0.0;
        for (int i = 1; i <= hr->GetNbinsX(); ++i) {
          const double r = hr->GetBinContent(i), m = hm->GetBinContent(i);
          const double d = std::fabs(r - m) / (std::fabs(r) > 0 ? std::fabs(r) : 1.0);
          if (d > maxDiff) maxDiff = d;
        }
        ZAKU_CHECK_MSG(maxDiff < 1e-12, name + " max relative bin difference " + std::to_string(maxDiff));
      }
    }

    FluxHistograms fh(uspec);
    FluxHistograms fnom(spec);
    for (int k = 0; k < 4; ++k) {
      const double expected = c.mode == BeamMode::FHC ? expectedFHC[k] : expectedRHC[k];
      ZAKU_CHECK_MSG(zakutest::close(fh.integralPerPOT(pdgs[k]), expected, 1e-5),
                     toString(c.mode) + " " + NuSpecies::fromPdg(pdgs[k]).name() + " integral " +
                         std::to_string(fh.integralPerPOT(pdgs[k])) + " vs " + std::to_string(expected));
      // the nominal (variable-bin) flux gives the same integral
      ZAKU_CHECK_CLOSE(fnom.integralPerPOT(pdgs[k]), fh.integralPerPOT(pdgs[k]), 1e-9);
    }
  }

  // The T2K covariance file: 100 rows, 16 blocks, SK rows from 50.
  FluxCovarianceSpec cs;
  cs.file = release + "/t2kflux_2020_covariance.root";
  FluxCovarianceInput cov(cs);
  ZAKU_CHECK_EQ(cov.matrix().GetNrows(), 100);
  ZAKU_CHECK_EQ(static_cast<int>(cov.blocks().size()), 16);
  ZAKU_CHECK_EQ(cov.firstRow("sk"), 50);
  ZAKU_CHECK_EQ(cov.nRows("sk"), 50);
  ZAKU_CHECK_EQ(cov.nRows("nd280"), 50);
  const auto sk = cov.blocksFor("sk");
  ZAKU_CHECK_EQ(sk[0].nBins, 11);
  ZAKU_CHECK(sk[0].beamMode == BeamMode::FHC);
  ZAKU_CHECK_EQ(sk[0].pdg, 14);
  ZAKU_CHECK_CLOSE(sk[0].edges[1], 0.4, 1e-12);
  ZAKU_CHECK_CLOSE(sk[0].edges.back(), 30.0, 1e-12);
  ZAKU_CHECK_EQ(sk[1].pdg, -14);
  ZAKU_CHECK_EQ(sk[1].nBins, 5);
  ZAKU_CHECK_EQ(sk[7].pdg, -12);
  ZAKU_CHECK(sk[7].beamMode == BeamMode::RHC);
  ZAKU_CHECK_EQ(sk[7].nBins, 7);
  // fractional sigmas of the unconstrained SK block: 4.69% .. 15.22%
  const TMatrixD skm = cov.subMatrix("sk");
  double smin = 1.0, smax = 0.0;
  for (int i = 0; i < 50; ++i) {
    const double s = std::sqrt(skm(i, i));
    if (s < smin) smin = s;
    if (s > smax) smax = s;
  }
  ZAKU_CHECK_CLOSE(smin, 0.0469, 3e-3);
  ZAKU_CHECK_CLOSE(smax, 0.1522, 3e-3);

  return zakutest::finish();
}
