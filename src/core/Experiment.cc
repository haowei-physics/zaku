#include "zaku/core/Experiment.h"

#include <TSystem.h>

#include <functional>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "ExperimentImpl.h"
#include "zaku/core/Units.h"

namespace zaku {

namespace {
bool fileExists(const std::string& path) { return !path.empty() && !gSystem->AccessPathName(path.c_str()); }
}  // namespace

Experiment::Experiment(std::string name) : impl_(new Impl) { impl_->name = std::move(name); }
Experiment::~Experiment() = default;
Experiment::Experiment(const Experiment& other) : impl_(new Impl(*other.impl_)) {}
Experiment& Experiment::operator=(const Experiment& other) {
  if (this != &other) impl_.reset(new Impl(*other.impl_));
  return *this;
}
Experiment::Experiment(Experiment&&) noexcept = default;
Experiment& Experiment::operator=(Experiment&&) noexcept = default;

const std::string& Experiment::configPath() const { return impl_->configPath; }

// --- identity, geometry, matter -------------------------------------------

const std::string& Experiment::name() const { return impl_->name; }
Experiment& Experiment::setName(const std::string& n) {
  impl_->name = n;
  return *this;
}
Experiment& Experiment::setBaseline(double km) {
  if (km <= 0.0) throw std::invalid_argument("Experiment::setBaseline: baseline must be positive (km)");
  impl_->baseline_km = km;
  return *this;
}
double Experiment::baseline_km() const { return impl_->baseline_km; }
Experiment& Experiment::setMatter(double density, double ye) {
  if (density < 0.0 || ye < 0.0 || ye > 1.0) {
    throw std::invalid_argument("Experiment::setMatter: density must be >= 0 g/cm3 and 0 <= electron fraction <= 1");
  }
  impl_->density_g_cm3 = density;
  impl_->electronFraction = ye;
  return *this;
}
double Experiment::matterDensity_g_cm3() const { return impl_->density_g_cm3; }
double Experiment::electronFraction() const { return impl_->electronFraction; }
double Experiment::rhoYe() const { return impl_->density_g_cm3 * impl_->electronFraction; }

// --- detector -------------------------------------------------------------

Experiment& Experiment::setDetector(const Detector& det) {
  impl_->detector = det;
  impl_->hasDetector = true;
  return *this;
}
const Detector& Experiment::detector() const {
  if (!impl_->hasDetector) throw std::runtime_error("Experiment '" + impl_->name + "': no detector set");
  return impl_->detector;
}
bool Experiment::hasDetector() const { return impl_->hasDetector; }
double Experiment::nNucleons() const { return detector().nNucleons(); }

// --- beam modes -----------------------------------------------------------

Experiment& Experiment::setFlux(BeamMode mode, const FluxSpec& spec, double pot) {
  Impl::Beam& b = impl_->beams[mode];
  b.flux = spec;
  b.hasFlux = true;
  b.pot = pot;
  impl_->fluxCache.erase(mode);
  return *this;
}
Experiment& Experiment::setPOT(BeamMode mode, double pot) {
  if (pot < 0.0) throw std::invalid_argument("Experiment::setPOT: POT must be >= 0");
  impl_->beams[mode].pot = pot;
  return *this;
}
Experiment& Experiment::setFluxUnits(double potNorm, double binNormGeV) {
  if (potNorm <= 0.0) throw std::invalid_argument("Experiment::setFluxUnits: potNorm must be positive");
  if (binNormGeV < 0.0) throw std::invalid_argument("Experiment::setFluxUnits: binNormGeV must be >= 0");
  impl_->defaultPotNorm = potNorm;
  impl_->defaultBinNormGeV = binNormGeV;
  impl_->fluxUnitsSet = true;
  for (auto& kv : impl_->beams) {
    kv.second.flux.potNorm = potNorm;
    kv.second.flux.binNormGeV = binNormGeV;
  }
  impl_->fluxCache.clear();
  return *this;
}
bool Experiment::hasBeamMode(BeamMode mode) const {
  auto it = impl_->beams.find(mode);
  return it != impl_->beams.end() && it->second.hasFlux;
}
std::vector<BeamMode> Experiment::beamModes() const {
  std::vector<BeamMode> out;
  for (const auto& kv : impl_->beams) {
    if (kv.second.hasFlux) out.push_back(kv.first);
  }
  return out;
}
const FluxSpec& Experiment::flux(BeamMode mode) const {
  if (!hasBeamMode(mode)) throw std::runtime_error("Experiment '" + impl_->name + "': no flux for " + toString(mode));
  return impl_->beams.at(mode).flux;
}
double Experiment::pot(BeamMode mode) const {
  auto it = impl_->beams.find(mode);
  return it == impl_->beams.end() ? 0.0 : it->second.pot;
}
const FluxHistograms& Experiment::fluxHistograms(BeamMode mode) const {
  auto it = impl_->fluxCache.find(mode);
  if (it == impl_->fluxCache.end()) {
    it = impl_->fluxCache.emplace(mode, std::make_shared<FluxHistograms>(flux(mode))).first;
  }
  return *it->second;
}
double Experiment::fluxIntegralPerPOT(BeamMode mode, int pdg) const { return fluxHistograms(mode).integralPerPOT(pdg); }
double Experiment::potWeight(BeamMode mode, int fluxPdg, double fScaleFactor) const {
  // Legacy expression order kept for bit-level reproducibility.
  return fScaleFactor * fluxIntegralPerPOT(mode, fluxPdg) * pot(mode) * nNucleons();
}
double Experiment::potWeight(const McSample& s, double fScaleFactor) const {
  return potWeight(s.beamMode, s.fluxPdg, fScaleFactor);
}

Experiment& Experiment::setFluxCovariance(const FluxCovarianceSpec& spec) {
  impl_->fluxCov = spec;
  return *this;
}
Experiment& Experiment::setFluxCovariance(const std::string& file, const std::string& matrix, double ndConstraint,
                                          const std::string& farPrefix, const std::string& nearPrefix) {
  FluxCovarianceSpec s;
  s.file = file;
  s.matrix = matrix;
  s.ndConstraint = ndConstraint;
  s.farPrefix = farPrefix;
  s.nearPrefix = nearPrefix;
  return setFluxCovariance(s);
}
const FluxCovarianceSpec& Experiment::fluxCovariance() const { return impl_->fluxCov; }

// --- MC samples -----------------------------------------------------------

Experiment& Experiment::addMcSample(BeamMode mode, int fluxPdg, int interactionPdg, long nEvents, unsigned seed) {
  return addMcSample(McSample(mode, fluxPdg, interactionPdg, nEvents, seed));
}
Experiment& Experiment::addMcSample(const McSample& s) {
  if (s.nEvents <= 0) throw std::invalid_argument("Experiment::addMcSample: nEvents must be positive");
  NuSpecies::fromPdg(s.fluxPdg);
  NuSpecies::fromPdg(s.interactionPdg);
  for (const auto& other : impl_->mcSamples) {
    if (other.name() == s.name()) {
      throw std::invalid_argument("Experiment::addMcSample: a sample named '" + s.name() + "' already exists");
    }
  }
  impl_->mcSamples.push_back(s);
  return *this;
}
Experiment& Experiment::clearMcSamples() {
  impl_->mcSamples.clear();
  return *this;
}
const std::vector<McSample>& Experiment::mcSamples() const { return impl_->mcSamples; }
std::vector<std::string> Experiment::mcSampleNames() const {
  std::vector<std::string> out;
  for (const auto& s : impl_->mcSamples) out.push_back(s.name());
  return out;
}
const McSample& Experiment::mcSample(const std::string& name) const {
  for (const auto& s : impl_->mcSamples) {
    if (s.name() == name) return s;
  }
  throw std::out_of_range("Experiment '" + impl_->name + "': no MC sample named '" + name + "'");
}
std::vector<McSample> Experiment::mcSamples(BeamMode mode) const {
  std::vector<McSample> out;
  for (const auto& s : impl_->mcSamples) {
    if (s.beamMode == mode) out.push_back(s);
  }
  return out;
}
unsigned Experiment::sampleSeed(const McSample& s) const {
  if (s.seed != 0) return s.seed;
  // Deterministic per-sample seed from the run seed and the sample name.
  const std::size_t h = std::hash<std::string>{}(s.name() + "#" + std::to_string(impl_->run.seed));
  unsigned seed = static_cast<unsigned>(h & 0x7fffffffu);
  if (seed == 0) seed = 1;
  return seed;
}

// --- run ------------------------------------------------------------------

Experiment& Experiment::setRun(const std::string& workDir, unsigned seed) {
  impl_->run.workDir = workDir;
  impl_->run.seed = seed;
  return *this;
}
Experiment& Experiment::setLegacyRngOrder(bool on) {
  impl_->run.legacyRngOrder = on;
  return *this;
}
const RunSettings& Experiment::run() const { return impl_->run; }
const std::string& Experiment::workDir() const { return impl_->run.workDir; }
std::string Experiment::path(const std::string& relative) const {
  if (!relative.empty() && relative[0] == '/') return relative;
  if (impl_->run.workDir.empty()) return relative;
  return impl_->run.workDir + "/" + relative;
}

// --- inspection -------------------------------------------------------------

std::string Experiment::describe() const {
  std::ostringstream os;
  os << "Experiment '" << impl_->name << "'";
  if (!impl_->configPath.empty()) os << " (from " << impl_->configPath << ")";
  os << "\n  baseline " << impl_->baseline_km << " km, matter " << impl_->density_g_cm3 << " g/cm3 x Ye "
     << impl_->electronFraction << " = " << rhoYe() << "\n";
  if (impl_->hasDetector) os << "  " << impl_->detector.describe() << "\n";
  else os << "  no detector\n";
  for (const auto& kv : impl_->beams) {
    os << "  " << toString(kv.first) << ": POT " << std::scientific << std::setprecision(4) << kv.second.pot
       << std::defaultfloat << ", " << (kv.second.hasFlux ? kv.second.flux.describe() : std::string("no flux")) << "\n";
  }
  if (impl_->fluxCov.valid()) {
    os << "  flux covariance " << impl_->fluxCov.file << ":" << impl_->fluxCov.matrix << " (far '"
       << impl_->fluxCov.farPrefix << "', ND constraint " << impl_->fluxCov.ndConstraint << ")\n";
  }
  os << "  " << impl_->mcSamples.size() << " MC sample(s):\n";
  for (const auto& s : impl_->mcSamples) os << "    " << s.describe() << " [seed " << sampleSeed(s) << "]\n";
  os << "  run: work dir '" << impl_->run.workDir << "', seed " << impl_->run.seed
     << (impl_->run.legacyRngOrder ? ", legacy RNG order" : "") << "\n";
  for (const auto& kv : impl_->raw) os << "  section '" << kv.first << "' kept for a later stage\n";
  return os.str();
}

std::string Experiment::status() const {
  std::ostringstream os;
  os << "status of '" << impl_->name << "' in '" << impl_->run.workDir << "'\n";
  auto mark = [&](const std::string& rel) { return fileExists(path(rel)) ? "ok      " : "missing "; };
  os << "  " << std::left << std::setw(36) << "MC sample" << "generator flat      weight    recon     spline\n";
  for (const auto& s : impl_->mcSamples) {
    os << "  " << std::left << std::setw(36) << s.name() << mark(SampleLayout::generatorFile(s)) << "  "
       << mark(SampleLayout::flatFile(s)) << "  " << mark(SampleLayout::weightFile(s)) << "  "
       << mark(SampleLayout::reconFile(s)) << "  " << mark(SampleLayout::splineFile(s)) << "\n";
  }
  for (const auto& kv : impl_->beams) {
    if (!kv.second.hasFlux) continue;
    os << "  " << std::left << std::setw(36) << ("select " + toString(kv.first)) << mark(SampleLayout::selectFile(kv.first))
       << "\n";
  }
  os << "  " << std::left << std::setw(36) << "parameters" << mark(SampleLayout::parametersDirectory()) << "\n";
  return os.str();
}

void Experiment::validate() const {
  std::vector<std::string> problems;
  if (impl_->name.empty()) problems.push_back("the experiment has no name");
  if (impl_->baseline_km <= 0.0) problems.push_back("baseline_km is not set");
  if (!impl_->hasDetector) problems.push_back("no detector");
  if (beamModes().empty()) problems.push_back("no beam mode with a flux");
  for (const auto& kv : impl_->beams) {
    if (!kv.second.hasFlux) continue;
    if (kv.second.pot <= 0.0) problems.push_back("beam mode " + toString(kv.first) + " has no POT");
    if (!fileExists(kv.second.flux.file)) {
      problems.push_back("flux file for " + toString(kv.first) + " not found: '" + kv.second.flux.file + "'");
    }
  }
  if (impl_->fluxCov.valid() && !fileExists(impl_->fluxCov.file)) {
    problems.push_back("flux covariance file not found: '" + impl_->fluxCov.file + "'");
  }
  for (const auto& s : impl_->mcSamples) {
    if (!hasBeamMode(s.beamMode)) {
      problems.push_back("MC sample " + s.name() + " uses beam mode " + toString(s.beamMode) + " which has no flux");
    } else if (!impl_->beams.at(s.beamMode).flux.has(s.fluxPdg)) {
      problems.push_back("MC sample " + s.name() + " needs a " + s.fluxSpecies().name() + " flux histogram for " +
                         toString(s.beamMode));
    }
  }
  if (impl_->run.workDir.empty()) problems.push_back("run.work_dir is not set");
  if (problems.empty()) return;
  std::string msg = "Experiment '" + impl_->name + "' is incomplete:";
  for (const auto& p : problems) msg += "\n  - " + p;
  throw std::runtime_error(msg);
}

std::string Experiment::rawSection(const std::string& key) const {
  auto it = impl_->raw.find(key);
  return it == impl_->raw.end() ? std::string() : it->second;
}
Experiment& Experiment::setRawSection(const std::string& key, const std::string& yamlText) {
  if (yamlText.empty()) impl_->raw.erase(key);
  else impl_->raw[key] = yamlText;
  return *this;
}

}  // namespace zaku
