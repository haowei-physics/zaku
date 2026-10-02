// The four parameter-set builders on synthetic inputs, the dial catalogue,
// and the systematics block of the experiment configuration.
#include <TAxis.h>
#include <TFile.h>
#include <TMatrixD.h>

#include <cmath>
#include <fstream>

#include "ZakuTest.h"
#include "ZakuTestPaths.h"
#include "zaku/core/Experiment.h"
#include "zaku/syst/Covariances.h"

using namespace zaku;

namespace {

SampleDefinition sample(const std::string& name, BeamMode mode, Hypothesis h) {
  SampleDefinition s;
  s.name = name;
  s.title = name;
  s.beamMode = static_cast<int>(mode);
  s.hypothesis = h;
  s.binning = Binning::product({BinningAxis::uniform("Enu_rec", 0.0, 3.0, 0.5)});
  return s;
}

}  // namespace

int main() {
  const std::string scratch = zakutest::scratchDir();
  const std::string src = zakutest::sourceDir();
  setenv("ZAKU_DATA", src.c_str(), 1);

  // --- dial catalogue --------------------------------------------------------
  const DialCatalog cat = DialCatalog::fromFile(src + "/configs/dials/neut_default.yaml");
  ZAKU_CHECK_EQ(static_cast<int>(cat.size()), 19);
  ZAKU_CHECK_EQ(static_cast<int>(cat.fitted().size()), 15);
  ZAKU_CHECK_EQ(static_cast<int>(cat.fitted({"MaCCQE"}).size()), 14);
  ZAKU_CHECK(!cat.dial("MDLSPiEj").fit);
  ZAKU_CHECK_EQ(static_cast<int>(cat.dial("PionFSI_InelProb").tweaks.size()), 4);
  ZAKU_CHECK_CLOSE(cat.dial("PionFSI_InelProb").tweaks[0], -0.9, 1e-12);
  ZAKU_CHECK_THROWS(cat.dial("nope"));
  const DialCatalog code = DialCatalog::neutDefault();
  ZAKU_CHECK(code.names() == cat.names());
  const DialCatalog back = DialCatalog::fromYamlString(cat.toYaml());
  ZAKU_CHECK(back.names() == cat.names());
  ZAKU_CHECK_EQ(static_cast<int>(back.fitted().size()), 15);
  ZAKU_CHECK_THROWS(DialCatalog::fromYamlString("dials: [{name: a}, {name: a}]"));
  ZAKU_CHECK_THROWS(DialCatalog::fromYamlString("dials: [{nam: a}]"));
  ZAKU_CHECK_THROWS(DialCatalog::fromYamlString("bogus: 1"));

  // --- oscillation ----------------------------------------------------------
  OscillationParameters point = OscillationParameters::asimovA();
  OscillationSystematics ocfg;
  const ParameterSet osc = OscCovariance::build(point, ocfg);
  ZAKU_CHECK_EQ(static_cast<int>(osc.size()), 7);
  ZAKU_CHECK_STR_EQ(osc.parameter(0).name, "PMNS_SIN_SQUARED_12");
  ZAKU_CHECK_CLOSE(osc.parameter(0).prior, 0.307, 1e-12);
  ZAKU_CHECK_CLOSE(osc.sigma(0), 0.013, 1e-12);
  ZAKU_CHECK_CLOSE(osc.parameter(4).prior, 2.509e-3, 1e-12);
  ZAKU_CHECK_CLOSE(osc.sigma(5), 1.17, 1e-12);
  ZAKU_CHECK_CLOSE(osc.parameter(6).prior, 1.0, 1e-12);
  ZAKU_CHECK(osc.parameter(6).fixed);
  ZAKU_CHECK(!osc.parameter(5).fixed);
  point.normalOrdering = false;
  point.s23sq = 0.45;
  const ParameterSet oscIO = OscCovariance::build(point, ocfg);
  ZAKU_CHECK_CLOSE(oscIO.parameter(6).prior, -1.0, 1e-12);
  ZAKU_CHECK_CLOSE(oscIO.parameter(2).prior, 0.45, 1e-12);
  const OscillationParameters moved = OscCovariance::moved(point, "PMNS_DELTA_CP", 0.5);
  ZAKU_CHECK_CLOSE(moved.deltaCP, point.deltaCP + 0.5, 1e-12);
  ZAKU_CHECK(OscCovariance::moved(point, "PMNS_SIGN_MASS_SQUARED_32", 2.5).normalOrdering);
  ZAKU_CHECK_THROWS(OscCovariance::moved(point, "nope", 1.0));

  // --- detector -------------------------------------------------------------
  std::vector<SampleDefinition> samples = {sample("FHC_1Rmu", BeamMode::FHC, Hypothesis::Muon),
                                           sample("FHC_1Re", BeamMode::FHC, Hypothesis::Electron),
                                           sample("RHC_1Rmu", BeamMode::RHC, Hypothesis::Muon)};
  DetectorSystematics dcfg;
  dcfg.norms = {{"FHC_1Rmu", 0.02, ""}, {"FHC_1Re", 0.03, "DET_NORM_FHC_1RE_0DE"}, {"RHC_1Rmu", 0.019, ""}};
  dcfg.correlation = 0.5;
  dcfg.pairs = {{"FHC_1Rmu", "RHC_1Rmu"}};
  const ParameterSet det = DetectorCovariance::normalisation(samples, dcfg);
  ZAKU_CHECK_EQ(static_cast<int>(det.size()), 3);
  ZAKU_CHECK_STR_EQ(det.parameter(0).name, "DET_NORM_FHC_1RMU");
  ZAKU_CHECK_STR_EQ(det.parameter(1).name, "DET_NORM_FHC_1RE_0DE");
  ZAKU_CHECK_CLOSE(det.parameter(0).prior, 1.0, 1e-12);
  ZAKU_CHECK_CLOSE(det.covariance()(0, 2), 0.5 * 0.02 * 0.019, 1e-12);
  ZAKU_CHECK_CLOSE(det.covariance()(0, 1), 0.0, 1e-12);
  ZAKU_CHECK(det.hasBinning());
  ZAKU_CHECK_EQ(det.binning().nBins(), 3);
  ZAKU_CHECK_EQ(det.binning().find({1.0}), 1);
  ZAKU_CHECK(!det.binning().isRange()[0]);
  ZAKU_CHECK_STR_EQ(det.binning().variables()[0], "EventType");
  // the text form round-trips, including a multi-line comment
  det.binning().writeFile(scratch + "/detBinning.txt", "line one\nline two");
  ZAKU_CHECK(Binning::fromFile(scratch + "/detBinning.txt").sameAs(det.binning()));
  DetectorSystematics missing = dcfg;
  missing.norms.pop_back();
  ZAKU_CHECK_THROWS(DetectorCovariance::normalisation(samples, missing));
  DetectorSystematics unknown = dcfg;
  unknown.norms.push_back({"nope", 0.1, ""});
  ZAKU_CHECK_THROWS(DetectorCovariance::normalisation(samples, unknown));
  DetectorSystematics badPair = dcfg;
  badPair.pairs = {{"FHC_1Rmu"}};
  ZAKU_CHECK_THROWS(DetectorCovariance::normalisation(samples, badPair));
  DetectorSystematics notPD = dcfg;
  notPD.correlation = 1.0;
  notPD.pairs = {{"FHC_1Rmu", "RHC_1Rmu"}, {"FHC_1Rmu", "FHC_1Re"}, {"FHC_1Re", "RHC_1Rmu"}};
  notPD.correlation = -1.0;
  ZAKU_CHECK_THROWS(DetectorCovariance::normalisation(samples, notPD));

  // --- cross section ----------------------------------------------------------
  DialKnots ma;
  ma.name = "MaCCQE";
  ma.nominal = 1.21;
  ma.knots = {0.818, 1.014, 1.406, 1.602};
  DialKnots inel;
  inel.name = "PionFSI_InelProb";
  inel.nominal = 1.002;
  inel.knots = {0.0111, 0.4515, 1.5525, 2.103};
  DialKnots srcO;
  srcO.name = "SRCFrac_O";
  srcO.nominal = 0.34;
  srcO.knots = {0.4225, 0.505};
  DialKnots rad;
  rad.name = "RadCorr_rate";
  rad.nominal = 0.0;
  rad.knots = {0.0, 0.25, 0.75};
  DialCatalog small;
  small.add("MaCCQE", {-2, -1, 1, 2});
  small.add("PionFSI_InelProb", {-0.9, -0.5, 0.5, 1});
  small.add("SRCFrac_O", {0.25, 0.5});
  small.add("RadCorr_rate", {0, 0.25, 0.75}, false);
  const ParameterSet xsec = XsecCovariance::build(small, {ma, inel, srcO, rad});
  ZAKU_CHECK_EQ(static_cast<int>(xsec.size()), 3);
  ZAKU_CHECK_CLOSE(xsec.parameter(0).prior, 1.21, 1e-12);
  ZAKU_CHECK_CLOSE(xsec.sigma(0), 0.196, 1e-9);
  ZAKU_CHECK_CLOSE(xsec.parameter(0).lower, 0.818, 1e-12);
  ZAKU_CHECK_CLOSE(xsec.parameter(0).upper, 1.602, 1e-12);
  ZAKU_CHECK_CLOSE(xsec.sigma(1), 1.101, 1e-9);
  ZAKU_CHECK_CLOSE(xsec.parameter(1).lower, 0.0111, 1e-12);
  ZAKU_CHECK_CLOSE(xsec.sigma(2), 0.33, 1e-9);
  ZAKU_CHECK_CLOSE(xsec.parameter(2).lower, 0.34, 1e-12);   // the nominal, below every knot
  ZAKU_CHECK_CLOSE(xsec.parameter(2).upper, 0.505, 1e-12);
  ZAKU_CHECK(xsec.anyBounds());
  // exclude on top of the catalogue
  ZAKU_CHECK_EQ(static_cast<int>(XsecCovariance::build(small, {ma, inel, srcO, rad}, {"MaCCQE"}).size()), 2);
  // a fitted dial with no weights, and knots that disagree with the multiples
  ZAKU_CHECK_THROWS(XsecCovariance::build(small, {ma, inel}));
  DialKnots wrong = ma;
  wrong.knots = {0.818, 1.1, 1.406, 1.602};
  ZAKU_CHECK_THROWS(XsecCovariance::build(small, {wrong, inel, srcO}));
  DialKnots fewer = ma;
  fewer.knots = {0.818, 1.602};
  ZAKU_CHECK_THROWS(XsecCovariance::build(small, {fewer, inel, srcO}));
  // RadCorr's tweak at the nominal (multiple 0) is skipped when deriving sigma
  DialCatalog radOnly;
  radOnly.add("RadCorr_rate", {0, 0.25, 0.75});
  ZAKU_CHECK_CLOSE(XsecCovariance::build(radOnly, {rad}).sigma(0), 1.0, 1e-12);

  // --- flux: a synthetic 2+2 covariance file --------------------------------------
  const std::string covFile = scratch + "/syn_flux_cov.root";
  {
    TFile f(covFile.c_str(), "RECREATE");
    TMatrixD full(4, 4);
    // ND block A (rows 0-1), far block C (rows 2-3), cross B
    const double A[2][2] = {{0.04, 0.01}, {0.01, 0.09}};
    const double C[2][2] = {{0.0025, 0.0010}, {0.0010, 0.0049}};
    const double B[2][2] = {{0.006, 0.002}, {0.003, 0.012}};  // rows ND, cols far
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        full(i, j) = A[i][j];
        full(2 + i, 2 + j) = C[i][j];
        full(i, 2 + j) = B[i][j];
        full(2 + j, i) = B[i][j];
      }
    }
    f.WriteTObject(&full, "total_flux_cov");
    const double edges[3] = {0.0, 1.0, 30.0};
    TAxis nd(2, edges);
    f.WriteTObject(&nd, "nd280_numode_numu_bins");
    TAxis sk(2, edges);
    f.WriteTObject(&sk, "sk_anumode_numub_bins");
  }
  FluxCovarianceSpec spec;
  spec.file = covFile;
  FluxCovarianceInput input(spec);
  const FluxCovarianceResult r = FluxCovariance::build(input, "sk", "nd280", 0.8);
  ZAKU_CHECK_EQ(static_cast<int>(r.unconstrained.size()), 2);
  ZAKU_CHECK_STR_EQ(r.unconstrained.parameter(0).name, "FLUX_RHC_NUMUB_00");
  ZAKU_CHECK_STR_EQ(r.constrained.parameter(1).name, "FLUX_RHC_NUMUB_01");
  ZAKU_CHECK_CLOSE(r.unconstrained.parameter(0).prior, 1.0, 1e-12);
  ZAKU_CHECK_CLOSE(r.unconstrained.sigma(0), 0.05, 1e-12);
  ZAKU_CHECK_CLOSE(r.unconstrained.covariance()(0, 1), 0.001, 1e-12);
  // C - s B^T A^-1 B by hand
  {
    TMatrixD A(2, 2), B(2, 2), C(2, 2);
    A(0, 0) = 0.04; A(0, 1) = 0.01; A(1, 0) = 0.01; A(1, 1) = 0.09;
    B(0, 0) = 0.006; B(0, 1) = 0.002; B(1, 0) = 0.003; B(1, 1) = 0.012;
    C(0, 0) = 0.0025; C(0, 1) = 0.0010; C(1, 0) = 0.0010; C(1, 1) = 0.0049;
    TMatrixD Ainv(A);
    Ainv.Invert();
    TMatrixD Bt(TMatrixD::kTransposed, B);
    TMatrixD expected = C - 0.8 * (Bt * Ainv * B);
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) ZAKU_CHECK_CLOSE(r.constrained.covariance()(i, j), expected(i, j), 1e-9);
    }
    ZAKU_CHECK(r.constrained.sigma(0) < r.unconstrained.sigma(0));
  }
  ZAKU_CHECK(r.lambdaMin > 0.0);
  // binning: isRHC and flux_pdg points, Enu_true range
  ZAKU_CHECK_EQ(r.binning.nBins(), 2);
  ZAKU_CHECK_EQ(static_cast<int>(r.binning.variables().size()), 3);
  ZAKU_CHECK_EQ(r.binning.find({1.0, -14.0, 0.5}), 0);
  ZAKU_CHECK_EQ(r.binning.find({1.0, -14.0, 5.0}), 1);
  ZAKU_CHECK_EQ(r.binning.find({0.0, -14.0, 0.5}), -1);
  ZAKU_CHECK(r.constrained.hasBinning());
  // s = 0 reproduces C
  const FluxCovarianceResult r0 = FluxCovariance::build(input, "sk", "nd280", 0.0);
  ZAKU_CHECK(r0.constrained.sameAs(r0.unconstrained, 1e-12));
  ZAKU_CHECK_THROWS(FluxCovariance::build(input, "sk", "nd280", 1.5));
  ZAKU_CHECK_STR_EQ(FluxCovariance::parameterName(BeamMode::FHC, 12, 7), "FLUX_FHC_NUE_07");
  ZAKU_CHECK(!r.describe().empty());

  // --- the systematics block of an experiment --------------------------------------
  Experiment exp = Experiment::fromConfig(src + "/configs/hk_t2k_like.yaml");
  ZAKU_CHECK(exp.hasDialCatalog());
  ZAKU_CHECK_EQ(static_cast<int>(exp.dialCatalog().size()), 19);
  ZAKU_CHECK_EQ(static_cast<int>(exp.systematics().xsec.fitted().size()), 15);
  ZAKU_CHECK_EQ(static_cast<int>(exp.systematics().xsec.exclude.size()), 4);
  ZAKU_CHECK(exp.systematics().xsec.writeUnusedDials);
  ZAKU_CHECK_STR_EQ(exp.systematics().detector.mode, "normalisation");
  ZAKU_CHECK_EQ(static_cast<int>(exp.systematics().detector.norms.size()), 5);
  ZAKU_CHECK_CLOSE(exp.systematics().detector.norm("FHC_1Re1de")->sigma, 0.132, 1e-12);
  ZAKU_CHECK_EQ(static_cast<int>(exp.systematics().detector.pairs.size()), 2);
  ZAKU_CHECK_CLOSE(exp.systematics().oscillation.sigmas.deltaCP, 1.17, 1e-12);
  ZAKU_CHECK_EQ(exp.systematics().grid.bins, 2000);
  ZAKU_CHECK_EQ(static_cast<int>(exp.shifts().size()), 2);
  ZAKU_CHECK(exp.shifts()[1].target == ShiftTarget::HadronicEnergy);
  ZAKU_CHECK_CLOSE(exp.shifts()[0].sigma, 0.024, 1e-12);
  ZAKU_CHECK(exp.describe().find("dials in the catalogue") != std::string::npos);
  ZAKU_CHECK_NOTHROW(exp.validate());
  // round trip
  const std::string yaml = exp.toYaml();
  Experiment again = Experiment::fromYamlString(yaml);
  ZAKU_CHECK_EQ(static_cast<int>(again.systematics().xsec.fitted().size()), 15);
  ZAKU_CHECK_EQ(static_cast<int>(again.shifts().size()), 2);
  ZAKU_CHECK_STR_EQ(again.toYaml(), yaml);
  // legacy names reach the detector set
  Experiment legacy = Experiment::fromConfig(src + "/configs/legacy_hk_mock.yaml");
  const ParameterSet legacyDet = DetectorCovariance::normalisation(legacy.samples(), legacy.systematics().detector);
  ZAKU_CHECK_STR_EQ(legacyDet.parameter(1).name, "DET_NORM_FHC_1RE_0DE");
  ZAKU_CHECK_CLOSE(legacyDet.covariance()(0, 3), 0.0, 1e-15);
  // unknown keys and bad values are rejected
  ZAKU_CHECK_THROWS(Experiment::fromYamlString("name: x\nsystematics: {bogus: 1}\n"));
  ZAKU_CHECK_THROWS(Experiment::fromYamlString("name: x\nsystematics: {detector: {mode: magic}}\n"));
  ZAKU_CHECK_THROWS(Experiment::fromYamlString("name: x\nshifts: [{name: a, target: nope, sigma: 1}]\n"));
  ZAKU_CHECK_THROWS(Experiment::fromYamlString("name: x\nshifts: [{name: a, target: lepton_momentum, sigma: 0}]\n"));
  // validate flags a sample without a detector sigma
  Experiment incomplete = exp;
  SampleDefinition extra = sample("FHC_extra", BeamMode::FHC, Hypothesis::Muon);
  incomplete.addSample(extra);
  ZAKU_CHECK_THROWS(incomplete.validate());

  return zakutest::finish();
}
