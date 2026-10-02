#ifndef ZAKU_SYST_SYSTEMATICSCONFIG_H
#define ZAKU_SYST_SYSTEMATICSCONFIG_H
#include <string>
#include <vector>

#include "zaku/core/DialCatalog.h"

namespace zaku {

/// Cross-section systematics: the dial catalogue and which dials are fitted.
struct XsecSystematics {
  DialCatalog dials;
  std::string dialsInclude;           ///< file the catalogue came from (emitted as include)
  std::vector<std::string> exclude;   ///< dials taken out of the fit on top of the catalogue's `fit: false`
  bool writeUnusedDials = true;       ///< spline files carry every dial of the weight file, not only the fitted ones
  bool enabled = true;

  std::vector<std::string> fitted() const { return dials.fitted(exclude); }
};

struct FluxSystematics {
  bool enabled = true;
};

/// One per-sample normalisation parameter of the detector set.
struct DetectorNorm {
  std::string sample;   ///< sample name
  double sigma = 0.0;   ///< fractional one-sigma
  std::string name;     ///< parameter name; empty = DET_NORM_<SAMPLE NAME UPPER-CASED>
};

struct DetectorSystematics {
  std::string mode = "normalisation";   ///< normalisation | shifts
  std::vector<DetectorNorm> norms;
  double correlation = 0.0;             ///< rho between the samples of each pair
  std::vector<std::vector<std::string>> pairs;   ///< pairs of sample names sharing `correlation`
  double rangeLow = 0.0;                ///< GUNDAM parametersRange for the normalisations
  double rangeHigh = 2.0;
  bool enabled = true;

  const DetectorNorm* norm(const std::string& sample) const;
};

/// One-sigma widths of the oscillation priors (the legacy MakeCovariances.C table).
struct OscSigmas {
  double s12sq = 0.013;
  double s13sq = 0.0007;
  double s23sq = 0.021;
  double dmsq21 = 1.8e-6;
  double dmsq32 = 3.4e-5;
  double deltaCP = 1.17;
  double ordering = 0.1;   ///< PMNS_SIGN_MASS_SQUARED_32
};

struct OscillationSystematics {
  OscSigmas sigmas;
  std::vector<std::string> fixed = {"PMNS_SIGN_MASS_SQUARED_32"};
  bool enabled = true;
};

/// The true-energy grid of the tabulated oscillation dials (GeV).
struct OscGridSpec {
  int bins = 2000;
  double min = 0.0;
  double max = 10.0;
};

/// The `systematics:` block of an experiment.
struct SystematicsConfig {
  XsecSystematics xsec;
  FluxSystematics flux;
  DetectorSystematics detector;
  OscillationSystematics oscillation;
  OscGridSpec grid;

  std::string describe() const;
};

}  // namespace zaku

#endif  // ZAKU_SYST_SYSTEMATICSCONFIG_H
