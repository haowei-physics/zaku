#include "ZakuTest.h"
#include "ZakuTestPaths.h"
#include "zaku/core/Experiment.h"
#include "zaku/detector/Binning.h"
#include "zaku/detector/Sample.h"
#include "zaku/osc/Oscillator.h"

using namespace zaku;

int main() {
  const std::string scratch = zakutest::scratchDir();
  const std::string src = zakutest::sourceDir();

  // --- product binning vs the shipped T2K files --------------------------------
  std::vector<double> extra;
  for (int i = 1; i <= 10; ++i) extra.push_back(2.0 + 0.1 * i);
  for (double e : {3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 30.0}) extra.push_back(e);
  BinningAxis enu = BinningAxis::uniform("Enu_rec", 0.0, 2.0, 0.05, extra);
  ZAKU_CHECK_EQ(enu.nBins(), 63);
  BinningAxis theta;
  theta.variable = "Theta_rec";
  theta.edges = {0, 20, 40, 60, 80, 100, 180};
  Binning mu = Binning::product({enu, theta});
  ZAKU_CHECK_EQ(mu.nBins(), 378);
  ZAKU_CHECK_EQ(static_cast<int>(mu.variables().size()), 2);
  ZAKU_CHECK(mu.isRange()[0] && mu.isRange()[1]);
  ZAKU_CHECK_EQ(mu.find({0.07, 25.0}), 1 * 6 + 1);  // Enu bin 1, theta bin 1 (last axis fastest)
  ZAKU_CHECK_EQ(mu.find({20.0, 25.0}), 63 * 6 - 6 + 1);
  ZAKU_CHECK_EQ(mu.find({40.0, 25.0}), -1);
  ZAKU_CHECK_THROWS(mu.find({1.0}));

  const std::string muFile = src + "/data/binning/binning_T2K_FHC_MuLike.txt";
  if (zakutest::exists(muFile)) {
    Binning ref = Binning::fromFile(muFile);
    ZAKU_CHECK_EQ(ref.nBins(), 378);
    ZAKU_CHECK_STR_EQ(ref.variables()[1], "Theta_rec");
    ZAKU_CHECK_MSG(mu.sameAs(ref), "product binning equals the T2K MuLike file");
    // text round trip
    Binning back = Binning::fromText(ref.toText("comment"), "text");
    ZAKU_CHECK(back.sameAs(ref));
    ref.writeFile(scratch + "/mu.txt", "written by the test");
    ZAKU_CHECK(Binning::fromFile(scratch + "/mu.txt").sameAs(ref));
  }
  const std::string eFile = src + "/data/binning/binning_T2K_FHC_ELike.txt";
  if (zakutest::exists(eFile)) {
    BinningAxis p = BinningAxis::uniform("PLep_rec", 0.0, 1.5, 0.1);
    BinningAxis th = BinningAxis::uniform("Theta_rec", 0.0, 140.0, 10.0, {180.0});
    Binning e = Binning::product({p, th});
    ZAKU_CHECK_EQ(e.nBins(), 225);
    ZAKU_CHECK(e.sameAs(Binning::fromFile(eFile)));
  }

  // point variables and mixed binnings
  BinningAxis type;
  type.variable = "EventType";
  type.values = {0, 1, 2};
  Binning pts = Binning::product({type});
  ZAKU_CHECK_EQ(pts.nBins(), 3);
  ZAKU_CHECK(!pts.isRange()[0]);
  ZAKU_CHECK_EQ(pts.find({1.0}), 1);
  ZAKU_CHECK_EQ(pts.find({1.5}), -1);
  const std::string text = pts.toText();
  ZAKU_CHECK(text.find("variables: EventType\n") != std::string::npos);
  Binning mixed = Binning::fromText("variables: isRHC flux_pdg Enu_true Enu_true\n0 14 0.0 0.4\n0 14 0.4 30\n1 -14 0 30\n");
  ZAKU_CHECK_EQ(mixed.nBins(), 3);
  ZAKU_CHECK_EQ(mixed.find({0, 14, 0.2}), 0);
  ZAKU_CHECK_EQ(mixed.find({1, -14, 5.0}), 2);
  ZAKU_CHECK_THROWS(Binning::fromText("0 1 2\n"));                          // no variables line
  ZAKU_CHECK_THROWS(Binning::fromText("variables: a a\n0 1 2\n"));          // wrong column count
  ZAKU_CHECK_THROWS(BinningAxis::uniform("v", 1.0, 0.0, 0.1));
  ZAKU_CHECK_THROWS(Binning::product({}));

  // --- samples, hypotheses, assignment ----------------------------------------------
  SampleDefinition s;
  s.name = "FHC_1Rmu";
  s.title = "FHC 1R#mu";
  s.beamMode = static_cast<int>(BeamMode::FHC);
  s.column = "1Rmu";
  s.hypothesis = Hypothesis::Muon;
  s.binning = mu;
  ZAKU_CHECK(s.appliesTo(BeamMode::FHC));
  ZAKU_CHECK(!s.appliesTo(BeamMode::RHC));
  ZAKU_CHECK_STR_EQ(s.columnName(), "1Rmu");
  SampleDefinition any;
  any.name = "X";
  ZAKU_CHECK(any.appliesTo(BeamMode::RHC));
  ZAKU_CHECK_STR_EQ(any.columnName(), "X");
  ZAKU_CHECK(s.describe().find("378 bins") != std::string::npos);
  ZAKU_CHECK(assignmentFromString("draw") == Assignment::Draw);
  ZAKU_CHECK(assignmentFromString("split") == Assignment::Split);
  ZAKU_CHECK_THROWS(assignmentFromString("maybe"));

  Experiment exp("t");
  exp.addSample(s);
  ZAKU_CHECK_THROWS(exp.addSample(s));  // duplicate
  SampleDefinition r = s;
  r.name = "RHC_1Rmu";
  r.beamMode = static_cast<int>(BeamMode::RHC);
  exp.addSample(r);
  ZAKU_CHECK_EQ(exp.sampleIndex("RHC_1Rmu"), 1);
  ZAKU_CHECK_EQ(exp.sampleIndex("nope"), -1);
  ZAKU_CHECK_EQ(static_cast<int>(exp.sampleIndices(BeamMode::FHC).size()), 1);
  ZAKU_CHECK_EQ(exp.sampleIndices(BeamMode::RHC)[0], 1);
  exp.setAssignment(Assignment::Draw).setEnergyEstimator(EnergyEstimator::Calorimetric).setBindingEnergy(0.03);
  ZAKU_CHECK(exp.assignment() == Assignment::Draw);
  ZAKU_CHECK(exp.energyEstimator() == EnergyEstimator::Calorimetric);
  ZAKU_CHECK_CLOSE(exp.bindingEnergy(), 0.03, 1e-12);
  ZAKU_CHECK_THROWS(exp.response());  // none set
  ZAKU_CHECK_THROWS(exp.oscillator()); // no baseline

  // --- the shipped sample and experiment configs load ---------------------------------
  setenv("ZAKU_DATA", src.c_str(), 1);
  for (const char* cfg : {"configs/hk_t2k_like.yaml", "configs/legacy_hk_mock.yaml", "configs/tiny_test.yaml",
                          "configs/hk_nova_like.yaml"}) {
    const std::string path = src + "/" + cfg;
    if (!zakutest::exists(path)) continue;
    Experiment e = Experiment::fromConfig(path);
    ZAKU_CHECK_MSG(e.hasResponse(), std::string(cfg) + " has a response");
    ZAKU_CHECK_MSG(!e.samples().empty(), std::string(cfg) + " has samples");
    ZAKU_CHECK_MSG(e.samples()[0].binning.nBins() > 0, std::string(cfg) + " sample binning");
    ZAKU_CHECK_MSG(e.sampleIndices(BeamMode::FHC).size() >= 2, std::string(cfg) + " FHC samples");
    // round trip through YAML keeps the detector blocks
    Experiment back = Experiment::fromYamlString(e.toYaml());
    ZAKU_CHECK_MSG(back.samples().size() == e.samples().size(), std::string(cfg) + " samples round trip");
    ZAKU_CHECK_MSG(back.assignment() == e.assignment(), std::string(cfg) + " assignment round trip");
    ZAKU_CHECK_MSG(back.energyEstimator() == e.energyEstimator(), std::string(cfg) + " energy round trip");
    ZAKU_CHECK_MSG(back.response().tables().size() == e.response().tables().size(), std::string(cfg) + " response round trip");
  }
  Experiment legacy = Experiment::fromConfig(src + "/configs/legacy_hk_mock.yaml");
  ZAKU_CHECK(legacy.assignment() == Assignment::Draw);
  ZAKU_CHECK(legacy.run().legacyRngOrder);
  ZAKU_CHECK_EQ(static_cast<int>(legacy.samples().size()), 5);
  ZAKU_CHECK_STR_EQ(legacy.samples()[2].columnName(), "1Re1de");
  ZAKU_CHECK(legacy.samples()[2].hypothesis == Hypothesis::Electron);
  ZAKU_CHECK_EQ(legacy.samples()[3].binning.nBins(), 378);
  ZAKU_CHECK_CLOSE(legacy.oscillation().deltaCP, -1.601, 1e-12);
  ZAKU_CHECK(legacy.response().categories().size() == 9);

  return zakutest::finish();
}
