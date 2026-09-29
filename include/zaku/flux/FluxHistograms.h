#ifndef ZAKU_FLUX_FLUXHISTOGRAMS_H
#define ZAKU_FLUX_FLUXHISTOGRAMS_H
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "zaku/core/FluxSpec.h"

class TH1;
class TH1D;

namespace zaku {

/// The flux histograms of one beam mode, loaded from a FluxSpec, converted to
/// GeV and owned here (detached from any ROOT directory).
class FluxHistograms {
public:
  explicit FluxHistograms(const FluxSpec& spec);
  ~FluxHistograms();
  FluxHistograms(const FluxHistograms&) = delete;
  FluxHistograms& operator=(const FluxHistograms&) = delete;

  const FluxSpec& spec() const { return spec_; }
  bool has(int pdg) const { return hists_.count(pdg) != 0; }
  std::vector<int> pdgs() const;
  /// The histogram (GeV axis). Throws std::out_of_range if the species is absent.
  const TH1D& histogram(int pdg) const;

  /// Integral of the density over energy: sum(content x width) in the file's units.
  double integral(int pdg) const;
  /// Neutrinos / cm^2 / POT through the detector: integral / binNormGeV / potNorm
  /// (or sum(content) / potNorm when binNormGeV == 0). This is the legacy
  /// `flux_per_pot`: 1.20736e-14 for the T2K FHC numu flux at SK.
  double integralPerPOT(int pdg) const;
  static double integralPerPOT(const TH1& h, double potNorm, double binNormGeV);

  std::string describe() const;

  /// Load one histogram from a file, converting MeV inputs to GeV and
  /// rejecting undeclared MeV axes. Caller owns the result.
  static TH1D* load(const std::string& file, const std::string& name, EnergyUnit units,
                    bool contentIsDensity = true);

private:
  FluxSpec spec_;
  std::map<int, std::unique_ptr<TH1D>> hists_;
};

}  // namespace zaku

#endif  // ZAKU_FLUX_FLUXHISTOGRAMS_H
