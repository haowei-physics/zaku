#ifndef ZAKU_SYST_COVARIANCES_H
#define ZAKU_SYST_COVARIANCES_H
#include <string>
#include <vector>

#include "zaku/core/DialCatalog.h"
#include "zaku/core/Types.h"
#include "zaku/detector/Sample.h"
#include "zaku/event/WeightCursor.h"
#include "zaku/flux/FluxCovarianceInput.h"
#include "zaku/osc/Oscillator.h"
#include "zaku/syst/ParameterSet.h"
#include "zaku/syst/SystematicsConfig.h"

namespace zaku {

/// The oscillation parameter set: the seven PMNS parameters GUNDAM's
/// tabulated NuOscillator dials take, priors at the experiment's oscillation
/// point, diagonal external constraints.
namespace OscCovariance {
/// PMNS_SIN_SQUARED_12, _13, _23, PMNS_DELTA_MASS_SQUARED_21, _32, PMNS_DELTA_CP, PMNS_SIGN_MASS_SQUARED_32.
const std::vector<std::string>& parameterNames();
ParameterSet build(const OscillationParameters& point, const OscillationSystematics& cfg = OscillationSystematics());
/// The OscillationParameters obtained by moving one named parameter by `delta`.
OscillationParameters moved(const OscillationParameters& point, const std::string& parameter, double delta);
}  // namespace OscCovariance

/// The detector parameter set as per-sample normalisations (legacy), with an
/// optional correlation between paired samples, and the binning that maps
/// parameter i to EventType i.
namespace DetectorCovariance {
std::string defaultName(const std::string& sample);
ParameterSet normalisation(const std::vector<SampleDefinition>& samples, const DetectorSystematics& cfg);
}  // namespace DetectorCovariance

/// The cross-section parameter set from the dial catalogue and the knots of a
/// weight file: sigma = (knot - nominal) / tweak multiple (every knot must
/// imply the same sigma), bounds = the outer knots, diagonal.
namespace XsecCovariance {
ParameterSet build(const DialCatalog& catalogue, const std::vector<DialKnots>& knots,
                   const std::vector<std::string>& exclude = {});
/// sigma implied by the knots of one dial; `spread` is (max - min) / |mean| over the knots.
double sigmaFromKnots(const DialSpec& dial, const DialKnots& knots, double* spread = nullptr);
}  // namespace XsecCovariance

/// The flux parameter sets: the far-detector block of the released
/// covariance, unconstrained (C) and with an emulated near-detector
/// constraint (C - s B^T A^-1 B), and their common binning.
struct FluxCovarianceResult {
  ParameterSet unconstrained;
  ParameterSet constrained;
  Binning binning;
  double ndConstraint = 0.0;
  double lambdaMin = 0.0;   ///< of the constrained matrix
  double lambdaMax = 0.0;
  std::string describe() const;
};

namespace FluxCovariance {
/// FLUX_<FHC|RHC>_<NUMU|NUMUB|NUE|NUEB>_<bin, two digits>.
std::string parameterName(BeamMode mode, int pdg, int bin);
FluxCovarianceResult build(const FluxCovarianceInput& input, const std::string& farPrefix,
                           const std::string& nearPrefix, double ndConstraint);
}  // namespace FluxCovariance

}  // namespace zaku

#endif  // ZAKU_SYST_COVARIANCES_H
