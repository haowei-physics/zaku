// The covariances stage: the four GUNDAM parameter-set files, their binning
// files and the oscillation energy grid, from the experiment's systematics
// block. Every matrix is a TMatrixD (see ParameterSet::write).
#include <TFile.h>
#include <TNamed.h>
#include <TSystem.h>

#include <memory>
#include <sstream>
#include <stdexcept>

#include "zaku/core/Experiment.h"
#include "zaku/core/Provenance.h"
#include "zaku/event/WeightCursor.h"
#include "zaku/osc/Oscillator.h"
#include "zaku/stages/Stages.h"
#include "zaku/syst/Covariances.h"
#include "zaku/syst/Splines.h"

namespace zaku {
namespace stages {

namespace {

bool fileExists(const std::string& p) { return !p.empty() && !gSystem->AccessPathName(p.c_str()); }

void ensureDir(const std::string& file) {
  const std::string dir = gSystem->GetDirName(file.c_str()).Data();
  if (!dir.empty()) gSystem->mkdir(dir.c_str(), true);
}

void stampFile(const std::string& file, const std::string& text) {
  std::unique_ptr<TFile> f(TFile::Open(file.c_str(), "UPDATE"));
  if (!f || f->IsZombie()) throw std::runtime_error("covariances: cannot reopen '" + file + "'");
  Provenance::stamp(*f, "covariances", text);
  f->Close();
}

}  // namespace

std::string ParameterSets::describe() const {
  std::ostringstream os;
  if (!osc.empty()) os << osc.describe();
  if (!det.empty()) os << det.describe();
  if (!xsec.empty()) os << xsec.describe() << "  knots from " << knotsSource << "\n";
  if (!fluxUnconstrained.empty()) {
    os << "flux parameter set: " << fluxUnconstrained.size() << " far-detector normalisations\n";
  }
  return os.str();
}

std::string CovarianceSummary::describe() const {
  std::ostringstream os;
  if (skipped) {
    os << "covariances: every parameter file exists, skipped";
    return os.str();
  }
  os << "covariances: wrote " << files.size() << " file(s)\n";
  for (const auto& f : files) os << "  " << f << "\n";
  os << sets.describe();
  return os.str();
}

std::vector<DialKnots> dialKnots(const Experiment& exp, std::string* source) {
  std::vector<std::string> tried;
  for (const auto& p : exp.productions()) {
    const std::string w = exp.path(Layout::weightFile(p));
    if (fileExists(w) && T2KRWWeightCursor::looksLikeT2KRW(w)) {
      T2KRWWeightCursor cur(w);
      if (source) *source = w;
      return cur.dials();
    }
    tried.push_back(w);
  }
  for (const auto& p : exp.productions()) {
    const std::string s = exp.path(Layout::splineFile(p));
    if (!fileExists(s)) continue;
    std::unique_ptr<TFile> f(TFile::Open(s.c_str(), "READ"));
    if (!f || f->IsZombie()) continue;
    TDirectory* z = f->GetDirectory("zaku");
    TNamed* meta = z ? z->Get<TNamed>("dials") : nullptr;
    if (!meta) continue;
    if (source) *source = s;
    return SplineBuilder::fromMetadataText(meta->GetTitle());
  }
  std::string msg = "no weight file (or spline file with dial metadata) exists for any production; run the weight manifest first. Looked for:";
  for (const auto& t : tried) msg += "\n  " + t;
  throw std::runtime_error(msg);
}

ParameterSets buildParameterSets(const Experiment& exp) {
  const SystematicsConfig& cfg = exp.systematics();
  ParameterSets sets;
  if (cfg.oscillation.enabled) sets.osc = OscCovariance::build(exp.oscillation(), cfg.oscillation);
  if (cfg.detector.enabled && cfg.detector.mode == "normalisation") {
    if (exp.samples().empty()) throw std::runtime_error("covariances: the experiment has no samples");
    sets.det = DetectorCovariance::normalisation(exp.samples(), cfg.detector);
  }
  if (cfg.xsec.enabled) {
    if (cfg.xsec.dials.empty()) {
      throw std::runtime_error("covariances: systematics.xsec.dials is empty; give the dial catalogue (e.g. {include: dials/neut_default.yaml})");
    }
    sets.knots = dialKnots(exp, &sets.knotsSource);
    sets.xsec = XsecCovariance::build(cfg.xsec.dials, sets.knots, cfg.xsec.exclude);
  }
  if (cfg.flux.enabled && exp.fluxCovariance().valid()) {
    const FluxCovarianceSpec& spec = exp.fluxCovariance();
    FluxCovarianceInput input(spec);
    const FluxCovarianceResult r = FluxCovariance::build(input, spec.farPrefix, spec.nearPrefix, spec.ndConstraint);
    sets.fluxUnconstrained = r.unconstrained;
    sets.fluxConstrained = r.constrained;
  }
  return sets;
}

CovarianceSummary makeCovariances(const Experiment& exp, bool force) {
  const SystematicsConfig& cfg = exp.systematics();
  CovarianceSummary summary;

  std::vector<std::string> expected;
  if (cfg.oscillation.enabled) {
    expected.push_back(exp.path(Layout::oscCovarianceFile()));
    expected.push_back(exp.path(Layout::oscGridFile()));
  }
  if (cfg.detector.enabled) {
    expected.push_back(exp.path(Layout::detCovarianceFile()));
    expected.push_back(exp.path(Layout::detBinningFile()));
  }
  if (cfg.xsec.enabled) expected.push_back(exp.path(Layout::xsecCovarianceFile()));
  if (cfg.flux.enabled && exp.fluxCovariance().valid()) {
    expected.push_back(exp.path(Layout::fluxCovarianceFile(false)));
    expected.push_back(exp.path(Layout::fluxCovarianceFile(true)));
    expected.push_back(exp.path(Layout::fluxBinningFile()));
  }
  if (expected.empty()) throw std::runtime_error("covariances: every systematics set is disabled; nothing to write");
  bool allExist = true;
  for (const auto& f : expected) allExist = allExist && fileExists(f);
  if (allExist && !force) {
    summary.skipped = true;
    summary.files = expected;
    return summary;
  }

  summary.sets = buildParameterSets(exp);
  if (cfg.detector.enabled && cfg.detector.mode == "shifts") {
    // the detector covariance from the shift universes needs the FakeSelect files
    if (exp.shifts().empty()) throw std::runtime_error("covariances: detector mode 'shifts' needs a `shifts:` list");
    SystematicsSummary sum(exp, summary.sets.xsec, ParameterSet("flux"), ParameterSet("det"), ParameterSet("osc"));
    for (BeamMode mode : exp.beamModes()) {
      const std::string f = exp.path(Layout::selectFile(mode));
      if (!fileExists(f)) throw std::runtime_error("covariances: detector mode 'shifts' needs '" + f + "'; run select first");
      sum.process(mode, f);
    }
    sum.finish();
    summary.sets.det = sum.shiftCovariance();
  }
  const ParameterSets& sets = summary.sets;
  const std::string head = "experiment " + exp.name() + "\n";

  if (!sets.osc.empty()) {
    const std::string f = exp.path(Layout::oscCovarianceFile());
    ensureDir(f);
    sets.osc.write(f);
    stampFile(f, head + sets.osc.describe());
    summary.files.push_back(f);
    const std::string g = exp.path(Layout::oscGridFile());
    Oscillator::writeGrid(g, cfg.grid.bins, cfg.grid.min, cfg.grid.max);
    stampFile(g, head + "oscillation grid: " + std::to_string(cfg.grid.bins) + " bins " + std::to_string(cfg.grid.min) +
                     "-" + std::to_string(cfg.grid.max) + " GeV");
    summary.files.push_back(g);
  }
  if (!sets.det.empty()) {
    const std::string f = exp.path(Layout::detCovarianceFile());
    ensureDir(f);
    sets.det.write(f);
    stampFile(f, head + sets.det.describe());
    summary.files.push_back(f);
    const std::string b = exp.path(Layout::detBinningFile());
    std::ostringstream comment;
    comment << "Detector parameters, one row per parameter, in the order of det_param_cov.\n"
            << "GUNDAM pairs them by position. Do not reorder.";
    for (std::size_t i = 0; i < sets.det.size(); ++i) {
      comment << "\n  row " << i << ": " << sets.det.parameter(static_cast<int>(i)).name << " sigma "
              << 100.0 * sets.det.sigma(static_cast<int>(i)) << "%";
    }
    sets.det.binning().writeFile(b, comment.str());
    summary.files.push_back(b);
  }
  if (!sets.xsec.empty()) {
    const std::string f = exp.path(Layout::xsecCovarianceFile());
    ensureDir(f);
    sets.xsec.write(f);
    stampFile(f, head + "knots from " + sets.knotsSource + "\n" + sets.xsec.describe());
    summary.files.push_back(f);
  }
  if (!sets.fluxUnconstrained.empty()) {
    const std::string f0 = exp.path(Layout::fluxCovarianceFile(false));
    const std::string f1 = exp.path(Layout::fluxCovarianceFile(true));
    ensureDir(f0);
    sets.fluxUnconstrained.write(f0);
    stampFile(f0, head + sets.fluxUnconstrained.title());
    sets.fluxConstrained.write(f1);
    stampFile(f1, head + sets.fluxConstrained.title());
    summary.files.push_back(f0);
    summary.files.push_back(f1);
    const std::string b = exp.path(Layout::fluxBinningFile());
    sets.fluxUnconstrained.binning().writeFile(
        b, "Far-detector flux normalisation bins, in the order of flux_param_cov.\nGUNDAM pairs them by position. Do not reorder.");
    summary.files.push_back(b);
  }
  return summary;
}

}  // namespace stages
}  // namespace zaku
