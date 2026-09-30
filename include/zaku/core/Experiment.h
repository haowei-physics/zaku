#ifndef ZAKU_CORE_EXPERIMENT_H
#define ZAKU_CORE_EXPERIMENT_H
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "zaku/core/Detector.h"
#include "zaku/core/FluxSpec.h"
#include "zaku/core/Production.h"
#include "zaku/core/RunSettings.h"
#include "zaku/core/Types.h"
#include "zaku/detector/Reco.h"
#include "zaku/detector/Response.h"
#include "zaku/detector/Sample.h"
#include "zaku/detector/Smearing.h"
#include "zaku/osc/Oscillator.h"

namespace zaku {

class FluxHistograms;

/// The one user-facing object: everything the user asserts about a
/// hypothetical experiment, and the pipeline stages as methods.
///
/// Build it in code with the chainable setters, or from a YAML file whose
/// keys map one-to-one onto them (see docs/config_reference.md). Sections the
/// current build does not interpret yet are kept verbatim so a config file
/// stays usable across stages.
class Experiment {
public:
  explicit Experiment(std::string name = "experiment");
  ~Experiment();
  Experiment(const Experiment& other);
  Experiment& operator=(const Experiment& other);
  Experiment(Experiment&&) noexcept;
  Experiment& operator=(Experiment&&) noexcept;

  // --- configuration files ---------------------------------------------
  /// Load a YAML file. Relative paths inside it are resolved against the
  /// file's directory; `${VAR}` is expanded from the environment.
  static Experiment fromConfig(const std::string& yamlPath);
  /// Parse YAML text; `baseDir` resolves relative paths.
  static Experiment fromYamlString(const std::string& yaml, const std::string& baseDir = "");
  std::string toYaml() const;
  void saveConfig(const std::string& path) const;
  /// Path of the file this experiment was loaded from (empty if built in code).
  const std::string& configPath() const;

  // --- identity, geometry, matter ---------------------------------------
  const std::string& name() const;
  Experiment& setName(const std::string& n);
  Experiment& setBaseline(double km);
  double baseline_km() const;
  Experiment& setMatter(double density_g_cm3, double electronFraction = 0.5);
  double matterDensity_g_cm3() const;
  double electronFraction() const;
  /// rho x Y_e, what the oscillation code takes (2.6 x 0.5 = 1.3 for HK).
  double rhoYe() const;

  // --- detector --------------------------------------------------------
  Experiment& setDetector(const Detector& det);
  const Detector& detector() const;
  bool hasDetector() const;
  double nNucleons() const;

  // --- beam modes: flux + exposure -------------------------------------
  Experiment& setFlux(BeamMode mode, const FluxSpec& spec, double pot);
  Experiment& setPOT(BeamMode mode, double pot);
  /// Default normalisation applied to fluxes that do not set their own.
  Experiment& setFluxUnits(double potNorm, double binNormGeV);
  bool hasBeamMode(BeamMode mode) const;
  std::vector<BeamMode> beamModes() const;
  const FluxSpec& flux(BeamMode mode) const;
  double pot(BeamMode mode) const;
  /// Loaded (and cached) histograms of a beam mode.
  const FluxHistograms& fluxHistograms(BeamMode mode) const;
  /// Neutrinos / cm^2 / POT of one species.
  double fluxIntegralPerPOT(BeamMode mode, int pdg) const;
  /// The per-event normalisation weight: fScaleFactor x flux/POT x POT x nucleons.
  double potWeight(BeamMode mode, int fluxPdg, double fScaleFactor) const;
  double potWeight(const Production& production, double fScaleFactor) const;

  Experiment& setFluxCovariance(const FluxCovarianceSpec& spec);
  Experiment& setFluxCovariance(const std::string& file, const std::string& matrix = "total_flux_cov",
                                double ndConstraint = 0.0, const std::string& farPrefix = "sk",
                                const std::string& nearPrefix = "nd280");
  const FluxCovarianceSpec& fluxCovariance() const;

  // --- productions ------------------------------------------------------
  Experiment& addProduction(BeamMode mode, int fluxPdg, int interactionPdg, long nEvents, unsigned seed = 0);
  Experiment& addProduction(const Production& s);
  Experiment& clearProductions();
  const std::vector<Production>& productions() const;
  std::vector<std::string> productionNames() const;
  const Production& production(const std::string& name) const;
  /// productions of one beam mode, in list order (this order defines SampleIndex).
  std::vector<Production> productions(BeamMode mode) const;
  /// The seed a production uses: its own, or one derived from the run seed and its name.
  unsigned productionSeed(const Production& s) const;

  // --- detector response and reconstruction ------------------------------
  /// The response tables: true event -> detection probability and column probabilities.
  Experiment& setResponse(const Response& r);
  const Response& response() const;
  bool hasResponse() const;
  Experiment& setSmearing(const Smearing& s);
  const Smearing& smearing() const;
  /// Which estimator fills Enu_rec (CCQE = T2K-like, Calorimetric = NOvA-like).
  Experiment& setEnergyEstimator(EnergyEstimator e);
  EnergyEstimator energyEstimator() const;
  /// How events are distributed over samples (split weights or one random draw).
  Experiment& setAssignment(Assignment a);
  Assignment assignment() const;
  /// Binding energy used by the CCQE estimator (GeV; T2K oxygen value 0.027).
  Experiment& setBindingEnergy(double GeV);
  double bindingEnergy() const;

  // --- analysis samples ------------------------------------------------
  /// Append a sample; its index in the list is EventType.
  Experiment& addSample(const SampleDefinition& s);
  Experiment& clearSamples();
  const std::vector<SampleDefinition>& samples() const;
  int sampleIndex(const std::string& name) const;
  /// Indices of the samples of one beam mode, in list order.
  std::vector<int> sampleIndices(BeamMode mode) const;

  // --- oscillation ------------------------------------------------------
  Experiment& setOscillation(const OscillationParameters& p);
  const OscillationParameters& oscillation() const;
  /// NuFast for this experiment's baseline and matter.
  Oscillator oscillator() const;

  // --- processing stages ------------------------------------------------
  /// Apply the response, smearing and energy estimators to one production
  /// (or, with an empty name, to all). Reads the flat tree, writes the recon file.
  void reconstruct(const std::string& production = "", bool force = false);
  /// Assign samples and write the FakeSelect tree and rate table of one beam
  /// mode (or all). Reads the spline files when present, else the recon files.
  void select(BeamMode mode, bool force = false);
  void select(bool force = false);

  // --- run -------------------------------------------------------------
  Experiment& setRun(const std::string& workDir, unsigned seed = 20250929);
  Experiment& setLegacyRngOrder(bool on);
  const RunSettings& run() const;
  const std::string& workDir() const;
  /// `workDir/relative` (or `relative` itself when absolute).
  std::string path(const std::string& relative) const;

  // --- inspection -------------------------------------------------------
  std::string describe() const;
  /// Which production files exist for every production and beam mode.
  std::string status() const;
  /// Throw std::runtime_error listing everything a stage needs but is not set.
  void validate() const;

  /// Raw YAML text of a section this build does not interpret (kept for
  /// later stages); empty if absent.
  std::string rawSection(const std::string& key) const;
  Experiment& setRawSection(const std::string& key, const std::string& yamlText);

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace zaku

#endif  // ZAKU_CORE_EXPERIMENT_H
