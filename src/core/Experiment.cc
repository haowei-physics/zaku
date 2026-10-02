#include "zaku/core/Experiment.h"

#include <TError.h>
#include <TSystem.h>

#include <functional>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "ExperimentImpl.h"
#include "zaku/core/Units.h"
#include "zaku/stages/Stages.h"

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
double Experiment::potWeight(const Production& s, double fScaleFactor) const {
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

// --- productions -----------------------------------------------------------

Experiment& Experiment::addProduction(BeamMode mode, int fluxPdg, int interactionPdg, long nEvents, unsigned seed) {
  return addProduction(Production(mode, fluxPdg, interactionPdg, nEvents, seed));
}
Experiment& Experiment::addProduction(const Production& s) {
  if (s.nEvents <= 0) throw std::invalid_argument("Experiment::addProduction: nEvents must be positive");
  NuSpecies::fromPdg(s.fluxPdg);
  NuSpecies::fromPdg(s.interactionPdg);
  for (const auto& other : impl_->productions) {
    if (other.name() == s.name()) {
      throw std::invalid_argument("Experiment::addProduction: a production named '" + s.name() + "' already exists");
    }
  }
  impl_->productions.push_back(s);
  return *this;
}
Experiment& Experiment::clearProductions() {
  impl_->productions.clear();
  return *this;
}
const std::vector<Production>& Experiment::productions() const { return impl_->productions; }
std::vector<std::string> Experiment::productionNames() const {
  std::vector<std::string> out;
  for (const auto& s : impl_->productions) out.push_back(s.name());
  return out;
}
const Production& Experiment::production(const std::string& name) const {
  for (const auto& s : impl_->productions) {
    if (s.name() == name) return s;
  }
  throw std::out_of_range("Experiment '" + impl_->name + "': no production named '" + name + "'");
}
std::vector<Production> Experiment::productions(BeamMode mode) const {
  std::vector<Production> out;
  for (const auto& s : impl_->productions) {
    if (s.beamMode == mode) out.push_back(s);
  }
  return out;
}
unsigned Experiment::productionSeed(const Production& s) const {
  if (s.seed != 0) return s.seed;
  // Deterministic per-production seed from the run seed and the production name.
  const std::size_t h = std::hash<std::string>{}(s.name() + "#" + std::to_string(impl_->run.seed));
  unsigned seed = static_cast<unsigned>(h & 0x7fffffffu);
  if (seed == 0) seed = 1;
  return seed;
}

// --- detector response and reconstruction -----------------------------------

Experiment& Experiment::setResponse(const Response& r) {
  impl_->response = r;
  impl_->hasResponse = true;
  impl_->responseInclude.clear();
  return *this;
}
const Response& Experiment::response() const {
  if (!impl_->hasResponse) throw std::runtime_error("Experiment '" + impl_->name + "': no response tables set");
  return impl_->response;
}
bool Experiment::hasResponse() const { return impl_->hasResponse; }
Experiment& Experiment::setSmearing(const Smearing& s) {
  impl_->smearing = s;
  impl_->smearingInclude.clear();
  return *this;
}
const Smearing& Experiment::smearing() const { return impl_->smearing; }
Experiment& Experiment::setEnergyEstimator(EnergyEstimator e) {
  impl_->estimator = e;
  return *this;
}
EnergyEstimator Experiment::energyEstimator() const { return impl_->estimator; }
Experiment& Experiment::setAssignment(Assignment a) {
  impl_->assignment = a;
  return *this;
}
Assignment Experiment::assignment() const { return impl_->assignment; }
Experiment& Experiment::setBindingEnergy(double GeV) {
  if (GeV < 0.0) throw std::invalid_argument("Experiment::setBindingEnergy: must be >= 0 GeV");
  impl_->bindingEnergy = GeV;
  return *this;
}
double Experiment::bindingEnergy() const { return impl_->bindingEnergy; }

// --- analysis samples -------------------------------------------------------

Experiment& Experiment::addSample(const SampleDefinition& s) {
  if (s.name.empty()) throw std::invalid_argument("Experiment::addSample: the sample has no name");
  if (sampleIndex(s.name) >= 0) throw std::invalid_argument("Experiment::addSample: a sample named '" + s.name + "' already exists");
  impl_->samples.push_back(s);
  impl_->samplesInclude.clear();
  return *this;
}
Experiment& Experiment::clearSamples() {
  impl_->samples.clear();
  return *this;
}
const std::vector<SampleDefinition>& Experiment::samples() const { return impl_->samples; }
int Experiment::sampleIndex(const std::string& name) const {
  for (std::size_t i = 0; i < impl_->samples.size(); ++i) {
    if (impl_->samples[i].name == name) return static_cast<int>(i);
  }
  return -1;
}
std::vector<int> Experiment::sampleIndices(BeamMode mode) const {
  std::vector<int> out;
  for (std::size_t i = 0; i < impl_->samples.size(); ++i) {
    if (impl_->samples[i].appliesTo(mode)) out.push_back(static_cast<int>(i));
  }
  return out;
}

// --- oscillation ------------------------------------------------------------

Experiment& Experiment::setOscillation(const OscillationParameters& p) {
  impl_->oscillation = p;
  return *this;
}
const OscillationParameters& Experiment::oscillation() const { return impl_->oscillation; }
Oscillator Experiment::oscillator() const {
  if (impl_->baseline_km <= 0.0) throw std::runtime_error("Experiment '" + impl_->name + "': baseline_km is not set");
  return Oscillator(impl_->oscillation, impl_->baseline_km, rhoYe(), 0);
}

// --- systematics -------------------------------------------------------------

Experiment& Experiment::setSystematics(const SystematicsConfig& s) {
  impl_->systematics = s;
  return *this;
}
const SystematicsConfig& Experiment::systematics() const { return impl_->systematics; }
Experiment& Experiment::setDialCatalog(const DialCatalog& c) {
  impl_->systematics.xsec.dials = c;
  impl_->systematics.xsec.dialsInclude.clear();
  return *this;
}
const DialCatalog& Experiment::dialCatalog() const { return impl_->systematics.xsec.dials; }
bool Experiment::hasDialCatalog() const { return !impl_->systematics.xsec.dials.empty(); }
Experiment& Experiment::addShift(const ShiftParameter& s) {
  if (s.name.empty()) throw std::invalid_argument("Experiment::addShift: the shift has no name");
  if (s.sigma <= 0.0) throw std::invalid_argument("Experiment::addShift: '" + s.name + "' needs a positive sigma");
  for (const auto& o : impl_->shifts) {
    if (o.name == s.name) throw std::invalid_argument("Experiment::addShift: a shift named '" + s.name + "' already exists");
  }
  impl_->shifts.push_back(s);
  return *this;
}
Experiment& Experiment::clearShifts() {
  impl_->shifts.clear();
  return *this;
}
const std::vector<ShiftParameter>& Experiment::shifts() const { return impl_->shifts; }

// --- processing stages ------------------------------------------------------

void Experiment::makeSplines(const std::string& productionName, bool force) {
  if (!productionName.empty()) {
    const stages::SplineSummary s = stages::makeSplines(*this, production(productionName), force);
    Info("zaku::Experiment::makeSplines", "%s", s.describe().c_str());
    return;
  }
  for (const auto& p : impl_->productions) {
    const stages::SplineSummary s = stages::makeSplines(*this, p, force);
    Info("zaku::Experiment::makeSplines", "%s", s.describe().c_str());
  }
}

void Experiment::makeCovariances(bool force) {
  const stages::CovarianceSummary s = stages::makeCovariances(*this, force);
  Info("zaku::Experiment::makeCovariances", "%s", s.describe().c_str());
}

void Experiment::summariseSystematics(bool force) {
  const stages::SummaryResult s = stages::summarise(*this, force);
  Info("zaku::Experiment::summariseSystematics", "%s", s.describe().c_str());
}

void Experiment::reconstruct(const std::string& productionName, bool force) {
  if (!productionName.empty()) {
    const stages::ReconSummary s = stages::reconstruct(*this, production(productionName), force);
    Info("zaku::Experiment::reconstruct", "%s", s.describe().c_str());
    return;
  }
  for (const auto& p : impl_->productions) {
    const stages::ReconSummary s = stages::reconstruct(*this, p, force);
    Info("zaku::Experiment::reconstruct", "%s", s.describe().c_str());
  }
}

void Experiment::select(BeamMode mode, bool force) {
  const stages::SelectSummary s = stages::select(*this, mode, force);
  Info("zaku::Experiment::select", "%s", s.describe().c_str());
}

void Experiment::select(bool force) {
  for (BeamMode mode : beamModes()) select(mode, force);
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
  os << "  " << impl_->productions.size() << " production(s):\n";
  for (const auto& s : impl_->productions) os << "    " << s.describe() << " [seed " << productionSeed(s) << "]\n";
  if (impl_->hasResponse) {
    os << "  response: " << impl_->response.categories().size() << " categories, " << impl_->response.tables().size()
       << " tables, columns";
    for (const auto& c : impl_->response.columns()) os << " " << c;
    os << ", assignment " << toString(impl_->assignment) << "\n";
  } else {
    os << "  no response tables\n";
  }
  os << "  " << impl_->smearing.describe() << "\n";
  os << "  energy estimator: " << toString(impl_->estimator) << ", binding energy " << impl_->bindingEnergy << " GeV\n";
  os << "  " << impl_->samples.size() << " sample(s):\n";
  for (std::size_t i = 0; i < impl_->samples.size(); ++i) os << "    [" << i << "] " << impl_->samples[i].describe() << "\n";
  os << "  oscillation: " << impl_->oscillation.describe() << "\n";
  {
    std::istringstream sys(impl_->systematics.describe());
    std::string line;
    while (std::getline(sys, line)) os << "  " << line << "\n";
  }
  if (!impl_->shifts.empty()) {
    os << "  " << impl_->shifts.size() << " shift systematic(s):\n";
    for (const auto& s : impl_->shifts) os << "    " << s.describe() << "\n";
  }
  os << "  run: work dir '" << impl_->run.workDir << "', seed " << impl_->run.seed
     << (impl_->run.legacyRngOrder ? ", legacy RNG order" : "") << "\n";
  for (const auto& kv : impl_->raw) os << "  section '" << kv.first << "' kept for a later stage\n";
  return os.str();
}

std::string Experiment::status() const {
  std::ostringstream os;
  os << "status of '" << impl_->name << "' in '" << impl_->run.workDir << "'\n";
  auto mark = [&](const std::string& rel) { return fileExists(path(rel)) ? "ok      " : "missing "; };
  os << "  " << std::left << std::setw(36) << "production" << "generator flat      weight    recon     spline\n";
  for (const auto& s : impl_->productions) {
    os << "  " << std::left << std::setw(36) << s.name() << mark(Layout::generatorFile(s)) << "  "
       << mark(Layout::flatFile(s)) << "  " << mark(Layout::weightFile(s)) << "  "
       << mark(Layout::reconFile(s)) << "  " << mark(Layout::splineFile(s)) << "\n";
  }
  for (const auto& kv : impl_->beams) {
    if (!kv.second.hasFlux) continue;
    os << "  " << std::left << std::setw(36) << ("select " + toString(kv.first)) << mark(Layout::selectFile(kv.first))
       << "\n";
  }
  os << "  " << std::left << std::setw(36) << "parameters: oscillation" << mark(Layout::oscCovarianceFile()) << "  grid "
     << mark(Layout::oscGridFile()) << "\n";
  os << "  " << std::left << std::setw(36) << "parameters: detector" << mark(Layout::detCovarianceFile()) << "  binning "
     << mark(Layout::detBinningFile()) << "\n";
  os << "  " << std::left << std::setw(36) << "parameters: cross section" << mark(Layout::xsecCovarianceFile()) << "\n";
  os << "  " << std::left << std::setw(36) << "parameters: flux" << mark(Layout::fluxCovarianceFile(false)) << "  ND "
     << mark(Layout::fluxCovarianceFile(true)) << "  binning " << mark(Layout::fluxBinningFile()) << "\n";
  os << "  " << std::left << std::setw(36) << "systematics summary" << mark(Layout::summaryCsvFile()) << "\n";
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
  for (const auto& s : impl_->productions) {
    if (!hasBeamMode(s.beamMode)) {
      problems.push_back("production " + s.name() + " uses beam mode " + toString(s.beamMode) + " which has no flux");
    } else if (!impl_->beams.at(s.beamMode).flux.has(s.fluxPdg)) {
      problems.push_back("production " + s.name() + " needs a " + s.fluxSpecies().name() + " flux histogram for " +
                         toString(s.beamMode));
    }
  }
  if (impl_->run.workDir.empty()) problems.push_back("run.work_dir is not set");
  if (impl_->systematics.detector.enabled && impl_->systematics.detector.mode == "normalisation") {
    for (const auto& s : impl_->samples) {
      if (!impl_->systematics.detector.norm(s.name)) {
        problems.push_back("systematics.detector.sigmas has no entry for sample '" + s.name + "'");
      }
    }
  }
  if (impl_->systematics.detector.mode == "shifts" && impl_->shifts.empty()) {
    problems.push_back("systematics.detector.mode is 'shifts' but no `shifts:` are defined");
  }
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
