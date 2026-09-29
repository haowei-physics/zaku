// Private state of zaku::Experiment, shared by Experiment.cc and ExperimentYaml.cc.
#ifndef ZAKU_SRC_CORE_EXPERIMENTIMPL_H
#define ZAKU_SRC_CORE_EXPERIMENTIMPL_H
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "zaku/core/Detector.h"
#include "zaku/core/Experiment.h"
#include "zaku/core/FluxSpec.h"
#include "zaku/core/McSample.h"
#include "zaku/core/RunSettings.h"
#include "zaku/flux/FluxHistograms.h"

namespace zaku {

struct Experiment::Impl {
  struct Beam {
    FluxSpec flux;
    double pot = 0.0;
    bool hasFlux = false;
  };

  std::string name = "experiment";
  std::string configPath;
  double baseline_km = 0.0;
  double density_g_cm3 = 2.6;
  double electronFraction = 0.5;
  bool hasDetector = false;
  Detector detector;
  std::map<BeamMode, Beam> beams;
  double defaultPotNorm = 1.0e21;
  double defaultBinNormGeV = 0.05;
  bool fluxUnitsSet = false;
  FluxCovarianceSpec fluxCov;
  std::vector<McSample> mcSamples;
  RunSettings run;
  std::map<std::string, std::string> raw;  ///< sections kept verbatim for later stages
  /// Loaded flux histograms (immutable once loaded, so sharing between copies is fine).
  mutable std::map<BeamMode, std::shared_ptr<FluxHistograms>> fluxCache;
};

/// The top-level YAML keys this build interprets (in emission order).
const std::vector<std::string>& experimentKnownKeys();
/// Top-level keys that later stages interpret; kept verbatim.
const std::vector<std::string>& experimentRawKeys();

}  // namespace zaku

#endif  // ZAKU_SRC_CORE_EXPERIMENTIMPL_H
