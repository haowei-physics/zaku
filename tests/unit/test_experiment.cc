#include <TFile.h>
#include <TH1D.h>
#include <TSystem.h>

#include <cstdlib>
#include <fstream>

#include "ZakuTest.h"
#include "ZakuTestPaths.h"
#include "zaku/core/Experiment.h"
#include "zaku/flux/FluxHistograms.h"

using namespace zaku;

namespace {

void writeFlatFlux(const std::string& file, double value) {
  TFile f(file.c_str(), "RECREATE");
  for (const char* name : {"enu_sk_numu", "enu_sk_numub", "enu_sk_nue", "enu_sk_nueb"}) {
    TH1D h(name, name, 300, 0.0, 3.0);
    for (int i = 1; i <= 300; ++i) h.SetBinContent(i, value);
    h.Write();
  }
}

}  // namespace

int main() {
  const std::string scratch = zakutest::scratchDir();
  const std::string fhcFile = scratch + "/exp_flux_fhc.root";
  const std::string rhcFile = scratch + "/exp_flux_rhc.root";
  writeFlatFlux(fhcFile, 2.0e5);
  writeFlatFlux(rhcFile, 1.0e5);

  // --- programmatic construction -------------------------------------------
  Experiment exp("HK");
  exp.setBaseline(295.0).setMatter(2.6, 0.5).setDetector(Detector::water(188.0));
  exp.setFlux(BeamMode::FHC, FluxSpec::t2k(fhcFile, "sk"), 6.75e21);
  exp.setFlux(BeamMode::RHC, FluxSpec::t2k(rhcFile, "sk"), 2.025e22);
  exp.setFluxCovariance(scratch + "/cov.root", "total_flux_cov", 0.8);
  exp.addMcSample(BeamMode::FHC, 14, 14, 200000);
  exp.addMcSample(BeamMode::FHC, 14, 12, 50000, 77);
  exp.addMcSample(BeamMode::RHC, -14, -14, 200000);
  exp.setRun(scratch + "/work", 20250929);

  ZAKU_CHECK_STR_EQ(exp.name(), "HK");
  ZAKU_CHECK_CLOSE(exp.rhoYe(), 1.3, 1e-12);
  ZAKU_CHECK_CLOSE(exp.baseline_km(), 295.0, 1e-12);
  ZAKU_CHECK(exp.hasDetector());
  ZAKU_CHECK_CLOSE(exp.nNucleons(), 1.132e35, 1e-3);
  ZAKU_CHECK(exp.hasBeamMode(BeamMode::FHC));
  ZAKU_CHECK(exp.hasBeamMode(BeamMode::RHC));
  ZAKU_CHECK_EQ(static_cast<int>(exp.beamModes().size()), 2);
  ZAKU_CHECK_CLOSE(exp.pot(BeamMode::RHC), 2.025e22, 1e-12);
  ZAKU_CHECK_EQ(static_cast<int>(exp.mcSamples().size()), 3);
  ZAKU_CHECK_STR_EQ(exp.mcSampleNames()[1], "numu_x_nue_50K_FHC");
  ZAKU_CHECK_EQ(static_cast<int>(exp.mcSamples(BeamMode::RHC).size()), 1);
  ZAKU_CHECK(exp.mcSample("numu_x_nue_50K_FHC").isAppearance());
  ZAKU_CHECK_EQ(static_cast<int>(exp.sampleSeed(exp.mcSample("numu_x_nue_50K_FHC"))), 77);
  ZAKU_CHECK(exp.sampleSeed(exp.mcSample("numu_x_numu_200K_FHC")) != 0u);
  ZAKU_CHECK(exp.sampleSeed(exp.mcSample("numu_x_numu_200K_FHC")) != exp.sampleSeed(exp.mcSample("numubar_x_numubar_200K_RHC")));
  ZAKU_CHECK_THROWS(exp.addMcSample(BeamMode::FHC, 14, 14, 200000));  // duplicate name
  ZAKU_CHECK_THROWS(exp.mcSample("nope"));
  ZAKU_CHECK_STR_EQ(exp.path("a/b.root"), scratch + "/work/a/b.root");
  ZAKU_CHECK_STR_EQ(exp.path("/abs.root"), "/abs.root");

  // flux integrals and the POT weight, legacy expression order
  const double perPot = 2.0e5 * 3.0 / 0.05 / 1.0e21;
  ZAKU_CHECK_CLOSE(exp.fluxIntegralPerPOT(BeamMode::FHC, 14), perPot, 1e-12);
  const double fScale = 4.456765e-44;
  ZAKU_CHECK_CLOSE(exp.potWeight(BeamMode::FHC, 14, fScale), fScale * perPot * 6.75e21 * exp.nNucleons(), 1e-12);
  ZAKU_CHECK_CLOSE(exp.potWeight(exp.mcSample("numu_x_nue_50K_FHC"), fScale), exp.potWeight(BeamMode::FHC, 14, fScale), 1e-12);
  ZAKU_CHECK_THROWS(exp.flux(BeamMode::RHC).histogram(16));

  // setFluxUnits propagates to the loaded fluxes
  exp.setFluxUnits(1.0e20, 0.1);
  ZAKU_CHECK_CLOSE(exp.fluxIntegralPerPOT(BeamMode::FHC, 14), 2.0e5 * 3.0 / 0.1 / 1.0e20, 1e-12);
  exp.setFluxUnits(1.0e21, 0.05);

  // validate: the covariance file does not exist yet (remove any leftover from an earlier run)
  gSystem->Unlink((scratch + "/cov.root").c_str());
  ZAKU_CHECK_THROWS(exp.validate());
  {
    TFile f((scratch + "/cov.root").c_str(), "RECREATE");
  }
  ZAKU_CHECK_NOTHROW(exp.validate());

  // --- YAML round trip ------------------------------------------------------
  exp.setRawSection("response", "tables: []\n");
  const std::string yaml = exp.toYaml();
  ZAKU_CHECK(yaml.find("beam_modes") != std::string::npos);
  Experiment back = Experiment::fromYamlString(yaml);
  ZAKU_CHECK_STR_EQ(back.name(), "HK");
  ZAKU_CHECK_CLOSE(back.baseline_km(), 295.0, 1e-12);
  ZAKU_CHECK_CLOSE(back.nNucleons(), exp.nNucleons(), 1e-12);
  ZAKU_CHECK_STR_EQ(back.detector().name(), "water");
  ZAKU_CHECK_CLOSE(back.pot(BeamMode::FHC), 6.75e21, 1e-12);
  ZAKU_CHECK_STR_EQ(back.flux(BeamMode::RHC).file, rhcFile);
  ZAKU_CHECK_STR_EQ(back.flux(BeamMode::RHC).histogram(-12), "enu_sk_nueb");
  ZAKU_CHECK_CLOSE(back.fluxCovariance().ndConstraint, 0.8, 1e-12);
  ZAKU_CHECK_EQ(static_cast<int>(back.mcSamples().size()), 3);
  ZAKU_CHECK_STR_EQ(back.mcSampleNames()[2], "numubar_x_numubar_200K_RHC");
  ZAKU_CHECK_EQ(static_cast<int>(back.mcSample("numu_x_nue_50K_FHC").seed), 77);
  ZAKU_CHECK_EQ(static_cast<int>(back.run().seed), 20250929);
  ZAKU_CHECK_STR_EQ(back.workDir(), scratch + "/work");
  ZAKU_CHECK(back.rawSection("response").find("tables") != std::string::npos);
  ZAKU_CHECK_STR_EQ(back.toYaml(), yaml);  // idempotent

  // --- YAML from a file with relative paths and an environment variable -----
  setenv("ZAKU_TEST_FLUX_DIR", scratch.c_str(), 1);
  const std::string cfg = scratch + "/exp_config.yaml";
  {
    std::ofstream f(cfg);
    f << "name: FromFile\n"
         "baseline_km: 295\n"
         "matter: {density_g_cm3: 2.6, electron_fraction: 0.5}\n"
         "detector: {material: water, fiducial_mass_kt: 188, n_nucleons: 1.132e35}\n"
         "flux_units: {pot_norm: 1.0e21, bin_norm_GeV: 0.05}\n"
         "beam_modes:\n"
         "  FHC: {pot: 6.75e21, flux: {file: exp_flux_fhc.root, t2k_detector: sk}}\n"
         "  RHC: {pot: 2.025e22, flux: {file: \"${ZAKU_TEST_FLUX_DIR}/exp_flux_rhc.root\", histograms: {numubar: enu_sk_numub}}}\n"
         "mc_samples:\n"
         "  - {beam: FHC, flux: numu, interaction: numu, n: 2000}\n"
         "  - {beam: RHC, flux: numubar, n: 2000}\n"
         "run: {work_dir: work_rel, seed: 5}\n"
         "samples: [{name: a}]\n";
  }
  Experiment fromFile = Experiment::fromConfig(cfg);
  ZAKU_CHECK_STR_EQ(fromFile.configPath(), cfg);
  ZAKU_CHECK_STR_EQ(fromFile.flux(BeamMode::FHC).file, scratch + "/exp_flux_fhc.root");
  ZAKU_CHECK_STR_EQ(fromFile.flux(BeamMode::RHC).file, scratch + "/exp_flux_rhc.root");
  ZAKU_CHECK_STR_EQ(fromFile.flux(BeamMode::FHC).histogram(-12), "enu_sk_nueb");
  ZAKU_CHECK_EQ(fromFile.mcSample("numubar_x_numubar_2K_RHC").interactionPdg, -14);
  ZAKU_CHECK_CLOSE(fromFile.nNucleons(), 1.132e35, 1e-12);
  ZAKU_CHECK_STR_EQ(fromFile.workDir(), scratch + "/work_rel");
  ZAKU_CHECK(!fromFile.rawSection("samples").empty());
  ZAKU_CHECK(fromFile.status().find("numu_x_numu_2K_FHC") != std::string::npos);
  ZAKU_CHECK(fromFile.describe().find("FromFile") != std::string::npos);

  // unknown keys are rejected, missing files reported
  ZAKU_CHECK_THROWS(Experiment::fromYamlString("name: x\nbogus: 1\n"));
  ZAKU_CHECK_THROWS(Experiment::fromYamlString("name: x\ndetector: {material: cheese, fiducial_mass_kt: 1}\n"));
  ZAKU_CHECK_THROWS(Experiment::fromYamlString("name: x\nbeam_modes: {FHC: {pot: 1, flux: {file: a.root}}}\n"));
  ZAKU_CHECK_THROWS(Experiment::fromConfig(scratch + "/missing.yaml"));

  // copy semantics: independent state
  Experiment copy = exp;
  copy.setName("copy").setPOT(BeamMode::FHC, 1.0);
  ZAKU_CHECK_STR_EQ(exp.name(), "HK");
  ZAKU_CHECK_CLOSE(exp.pot(BeamMode::FHC), 6.75e21, 1e-12);
  ZAKU_CHECK_CLOSE(copy.pot(BeamMode::FHC), 1.0, 1e-12);

  return zakutest::finish();
}
