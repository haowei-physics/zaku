#include "zaku/syst/Covariances.h"

#include <TMatrixD.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace zaku {

// --- SystematicsConfig helpers ---------------------------------------------------

const DetectorNorm* DetectorSystematics::norm(const std::string& sample) const {
  for (const auto& n : norms) {
    if (n.sample == sample) return &n;
  }
  return nullptr;
}

std::string SystematicsConfig::describe() const {
  std::ostringstream os;
  os << "systematics:\n";
  os << "  xsec: " << (xsec.enabled ? "" : "disabled, ") << xsec.dials.size() << " dials in the catalogue, "
     << xsec.fitted().size() << " fitted";
  if (!xsec.exclude.empty()) {
    os << ", excluded";
    for (const auto& e : xsec.exclude) os << " " << e;
  }
  os << (xsec.writeUnusedDials ? ", splines for every dial" : ", splines for fitted dials only") << "\n";
  os << "  flux: " << (flux.enabled ? "enabled" : "disabled") << "\n";
  os << "  detector: " << (detector.enabled ? "" : "disabled, ") << detector.mode << ", " << detector.norms.size()
     << " normalisation(s), correlation " << detector.correlation << " on " << detector.pairs.size() << " pair(s)\n";
  os << "  oscillation: " << (oscillation.enabled ? "" : "disabled, ") << "sigmas s12sq " << oscillation.sigmas.s12sq
     << ", s13sq " << oscillation.sigmas.s13sq << ", s23sq " << oscillation.sigmas.s23sq << ", dm21 "
     << oscillation.sigmas.dmsq21 << ", dm32 " << oscillation.sigmas.dmsq32 << ", dCP " << oscillation.sigmas.deltaCP
     << "; fixed";
  for (const auto& f : oscillation.fixed) os << " " << f;
  os << "\n  grid: " << grid.bins << " bins " << grid.min << "-" << grid.max << " GeV\n";
  return os.str();
}

// --- oscillation ------------------------------------------------------------------

namespace OscCovariance {

const std::vector<std::string>& parameterNames() {
  static const std::vector<std::string> n = {"PMNS_SIN_SQUARED_12",        "PMNS_SIN_SQUARED_13",
                                             "PMNS_SIN_SQUARED_23",        "PMNS_DELTA_MASS_SQUARED_21",
                                             "PMNS_DELTA_MASS_SQUARED_32", "PMNS_DELTA_CP",
                                             "PMNS_SIGN_MASS_SQUARED_32"};
  return n;
}

ParameterSet build(const OscillationParameters& point, const OscillationSystematics& cfg) {
  ParameterSet s("osc");
  s.setTitle("PMNS parameters, priors at the experiment's oscillation point");
  const std::vector<double> priors = {point.s12sq,  point.s13sq,  point.s23sq,  point.dmsq21,
                                      point.dmsq32, point.deltaCP, point.normalOrdering ? 1.0 : -1.0};
  const std::vector<double> sigmas = {cfg.sigmas.s12sq,  cfg.sigmas.s13sq,   cfg.sigmas.s23sq,   cfg.sigmas.dmsq21,
                                      cfg.sigmas.dmsq32, cfg.sigmas.deltaCP, cfg.sigmas.ordering};
  for (std::size_t i = 0; i < parameterNames().size(); ++i) {
    Parameter p;
    p.name = parameterNames()[i];
    p.prior = priors[i];
    p.sigma = sigmas[i];
    p.fixed = std::find(cfg.fixed.begin(), cfg.fixed.end(), p.name) != cfg.fixed.end();
    s.add(p);
  }
  return s;
}

OscillationParameters moved(const OscillationParameters& point, const std::string& name, double delta) {
  OscillationParameters p = point;
  if (name == "PMNS_SIN_SQUARED_12") p.s12sq += delta;
  else if (name == "PMNS_SIN_SQUARED_13") p.s13sq += delta;
  else if (name == "PMNS_SIN_SQUARED_23") p.s23sq += delta;
  else if (name == "PMNS_DELTA_MASS_SQUARED_21") p.dmsq21 += delta;
  else if (name == "PMNS_DELTA_MASS_SQUARED_32") p.dmsq32 += delta;
  else if (name == "PMNS_DELTA_CP") p.deltaCP += delta;
  else if (name == "PMNS_SIGN_MASS_SQUARED_32") p.normalOrdering = ((point.normalOrdering ? 1.0 : -1.0) + delta) >= 0.0;
  else throw std::invalid_argument("OscCovariance::moved: unknown parameter '" + name + "'");
  return p;
}

}  // namespace OscCovariance

// --- detector -----------------------------------------------------------------------

namespace DetectorCovariance {

std::string defaultName(const std::string& sample) {
  std::string s = "DET_NORM_";
  for (char c : sample) s += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return s;
}

ParameterSet normalisation(const std::vector<SampleDefinition>& samples, const DetectorSystematics& cfg) {
  if (cfg.correlation < -1.0 || cfg.correlation > 1.0) {
    throw std::invalid_argument("DetectorCovariance: correlation must be in [-1, 1]");
  }
  ParameterSet s("det");
  s.setTitle("per-sample normalisations");
  std::vector<double> eventTypes;
  for (std::size_t i = 0; i < samples.size(); ++i) {
    const DetectorNorm* n = cfg.norm(samples[i].name);
    if (!n) {
      throw std::runtime_error("DetectorCovariance: systematics.detector.sigmas has no entry for sample '" +
                               samples[i].name + "'");
    }
    Parameter p;
    p.name = n->name.empty() ? defaultName(samples[i].name) : n->name;
    p.prior = 1.0;
    p.sigma = n->sigma;
    s.add(p);
    eventTypes.push_back(static_cast<double>(i));
  }
  for (const auto& n : cfg.norms) {
    bool known = false;
    for (const auto& smp : samples) known = known || smp.name == n.sample;
    if (!known) throw std::runtime_error("DetectorCovariance: '" + n.sample + "' in systematics.detector.sigmas is not a sample");
  }
  for (const auto& pair : cfg.pairs) {
    if (pair.size() != 2) throw std::runtime_error("DetectorCovariance: every entry of `pairs` needs exactly two sample names");
    int a = -1, b = -1;
    for (std::size_t i = 0; i < samples.size(); ++i) {
      if (samples[i].name == pair[0]) a = static_cast<int>(i);
      if (samples[i].name == pair[1]) b = static_cast<int>(i);
    }
    if (a < 0 || b < 0) throw std::runtime_error("DetectorCovariance: pair (" + pair[0] + ", " + pair[1] + ") names an unknown sample");
    s.setCorrelation(a, b, cfg.correlation);
  }
  double lo = 0.0, hi = 0.0;
  if (!s.isPositiveDefinite(&lo, &hi)) {
    throw std::runtime_error("DetectorCovariance: the matrix is not positive definite (lambda_min " + std::to_string(lo) + ")");
  }
  BinningAxis axis;
  axis.variable = "EventType";
  axis.values = eventTypes;
  s.setBinning(Binning::product({axis}));
  return s;
}

}  // namespace DetectorCovariance

// --- cross section --------------------------------------------------------------------

namespace XsecCovariance {

double sigmaFromKnots(const DialSpec& dial, const DialKnots& knots, double* spread) {
  if (dial.tweaks.size() != knots.knots.size()) {
    throw std::runtime_error("XsecCovariance: dial '" + dial.name + "' has " + std::to_string(knots.knots.size()) +
                             " knots in the weight file but " + std::to_string(dial.tweaks.size()) +
                             " tweak multiples in the catalogue");
  }
  double sum = 0.0, lo = 1e300, hi = -1e300;
  for (std::size_t k = 0; k < knots.knots.size(); ++k) {
    if (dial.tweaks[k] == 0.0) continue;  // a tweak at the nominal says nothing about sigma
    const double s = (knots.knots[k] - knots.nominal) / dial.tweaks[k];
    sum += s;
    lo = std::min(lo, s);
    hi = std::max(hi, s);
  }
  std::size_t nUsed = 0;
  for (double t : dial.tweaks) nUsed += t != 0.0 ? 1 : 0;
  if (nUsed == 0) throw std::runtime_error("XsecCovariance: dial '" + dial.name + "' has only a tweak at the nominal");
  const double mean = sum / static_cast<double>(nUsed);
  const double sp = mean != 0.0 ? (hi - lo) / std::fabs(mean) : 0.0;
  if (spread) *spread = sp;
  return mean;
}

ParameterSet build(const DialCatalog& catalogue, const std::vector<DialKnots>& knots,
                   const std::vector<std::string>& exclude) {
  ParameterSet s("xsec");
  s.setTitle(catalogue.generator() + " cross-section dials, sigma from the tweak points");
  for (const auto& name : catalogue.fitted(exclude)) {
    const DialSpec& dial = catalogue.dial(name);
    const DialKnots* k = nullptr;
    for (const auto& kk : knots) {
      if (kk.name == name) k = &kk;
    }
    if (!k) throw std::runtime_error("XsecCovariance: fitted dial '" + name + "' has no weights in the weight file");
    double spread = 0.0;
    const double sigma = sigmaFromKnots(dial, *k, &spread);
    if (spread > 1e-3) {
      throw std::runtime_error("XsecCovariance: the knots of '" + name + "' imply sigmas spread by " +
                               std::to_string(spread) + " relative; the catalogue's tweak multiples do not match the file");
    }
    Parameter p;
    p.name = name;
    p.prior = k->nominal;
    p.sigma = std::fabs(sigma);
    p.hasBounds = true;
    p.lower = std::min(k->nominal, *std::min_element(k->knots.begin(), k->knots.end()));
    p.upper = std::max(k->nominal, *std::max_element(k->knots.begin(), k->knots.end()));
    s.add(p);
  }
  if (s.empty()) throw std::runtime_error("XsecCovariance: no fitted dials");
  return s;
}

}  // namespace XsecCovariance

// --- flux ------------------------------------------------------------------------------

std::string FluxCovarianceResult::describe() const {
  std::ostringstream os;
  os << "flux parameters: " << unconstrained.size() << " far-detector normalisations, ND constraint " << ndConstraint
     << " (constrained eigenvalues " << lambdaMin << " .. " << lambdaMax << ")\n";
  os << std::left << std::setw(20) << "  parameter" << std::right << std::setw(10) << "no ND" << std::setw(10) << "with ND"
     << "\n";
  double sumA = 0.0, sumB = 0.0;
  for (std::size_t i = 0; i < unconstrained.size(); ++i) {
    const double a = unconstrained.sigma(static_cast<int>(i)), b = constrained.sigma(static_cast<int>(i));
    sumA += a;
    sumB += b;
    os << "  " << std::left << std::setw(18) << unconstrained.parameter(static_cast<int>(i)).name << std::right
       << std::fixed << std::setprecision(2) << std::setw(9) << 100.0 * a << "%" << std::setw(9) << 100.0 * b << "%\n";
  }
  if (!unconstrained.empty()) {
    os << "  mean sigma " << 100.0 * sumA / unconstrained.size() << "% -> " << 100.0 * sumB / unconstrained.size() << "%\n";
  }
  os << std::defaultfloat;
  return os.str();
}

namespace FluxCovariance {

std::string parameterName(BeamMode mode, int pdg, int bin) {
  std::string suffix = NuSpecies::fromPdg(pdg).t2kSuffix();
  for (auto& c : suffix) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  std::ostringstream os;
  os << "FLUX_" << toString(mode) << "_" << suffix << "_" << std::setw(2) << std::setfill('0') << bin;
  return os.str();
}

FluxCovarianceResult build(const FluxCovarianceInput& input, const std::string& far, const std::string& near,
                           double s) {
  if (s < 0.0 || s > 1.0) throw std::invalid_argument("FluxCovariance: the ND constraint strength must be in [0, 1]");
  FluxCovarianceResult r;
  r.ndConstraint = s;
  const std::vector<FluxCovarianceBlock> blocks = input.blocksFor(far);
  const TMatrixD C = input.subMatrix(far);
  const int n = C.GetNrows();

  // parameters and binning rows, in block order
  std::vector<Bin> bins;
  for (const auto& b : blocks) {
    for (int k = 0; k < b.nBins; ++k) {
      Parameter p;
      p.name = parameterName(b.beamMode, b.pdg, k);
      p.prior = 1.0;
      p.sigma = 0.0;
      r.unconstrained.add(p);
      r.constrained.add(p);
      Bin bin;
      bin.low = {isRHC(b.beamMode) ? 1.0 : 0.0, static_cast<double>(b.pdg), b.edges[static_cast<std::size_t>(k)]};
      bin.high = {bin.low[0], bin.low[1], b.edges[static_cast<std::size_t>(k) + 1]};
      bins.push_back(bin);
    }
  }
  r.unconstrained.setPrefix("flux");
  r.constrained.setPrefix("flux");
  r.unconstrained.setTitle(far + " flux normalisations, no near-detector constraint");
  r.constrained.setTitle(far + " flux normalisations, ND constraint " + std::to_string(s));
  if (static_cast<int>(r.unconstrained.size()) != n) throw std::logic_error("FluxCovariance: block/matrix size mismatch");
  r.unconstrained.setCovariance(C, 1e-6);

  TMatrixD Cnd(C);
  if (s > 0.0) {
    TMatrixD A = input.subMatrix(near);
    TMatrixD B = input.crossMatrix(near, far);  // rows ND, columns far
    if (A.GetNrows() != B.GetNrows()) throw std::runtime_error("FluxCovariance: ND and cross blocks disagree in size");
    double det = 0.0;
    TMatrixD Ainv(A);
    Ainv.Invert(&det);
    if (det == 0.0) throw std::runtime_error("FluxCovariance: the near-detector block is singular");
    TMatrixD Bt(TMatrixD::kTransposed, B);
    TMatrixD G = Bt * Ainv * B;
    Cnd -= s * G;
  }
  // symmetrise: A^-1 is numerical
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      const double v = 0.5 * (Cnd(i, j) + Cnd(j, i));
      Cnd(i, j) = v;
      Cnd(j, i) = v;
    }
  }
  r.constrained.setCovariance(Cnd, 1e-6);
  if (!r.constrained.isPositiveDefinite(&r.lambdaMin, &r.lambdaMax)) {
    throw std::runtime_error("FluxCovariance: the constrained matrix is not positive definite (lambda_min " +
                             std::to_string(r.lambdaMin) + ")");
  }

  // binning text: isRHC and flux_pdg as points, Enu_true as a range
  std::ostringstream txt;
  txt << "variables: isRHC flux_pdg Enu_true Enu_true\n";
  std::string current;
  for (std::size_t i = 0; i < bins.size(); ++i) {
    const std::string& name = r.unconstrained.parameter(static_cast<int>(i)).name;
    const std::string tag = name.substr(5, name.size() - 8);
    if (tag != current) {
      current = tag;
      txt << "# " << tag << "\n";
    }
    txt << static_cast<int>(bins[i].low[0]) << " " << static_cast<int>(bins[i].low[1]) << " " << std::fixed
        << std::setprecision(4) << std::setw(9) << bins[i].low[2] << " " << std::setw(9) << bins[i].high[2] << "\n";
  }
  r.binning = Binning::fromText(txt.str(), "flux binning");
  r.unconstrained.setBinning(r.binning);
  r.constrained.setBinning(r.binning);
  return r;
}

}  // namespace FluxCovariance

}  // namespace zaku
