#ifndef ZAKU_CORE_FLUXSPEC_H
#define ZAKU_CORE_FLUXSPEC_H
#include <map>
#include <string>
#include <vector>

#include "zaku/core/Types.h"
#include "zaku/core/Units.h"

namespace zaku {

/// Where a beam flux comes from and how its histograms are normalised.
///
/// A flux file holds one TH1 per neutrino species with the energy on the x
/// axis (GeV; MeV inputs must set `units`). The contents are a flux density:
/// neutrinos / cm^2 / `binNormGeV` GeV / `potNorm` protons on target, which is
/// the T2K public-release convention (per 50 MeV per 1e21 POT). Set
/// `binNormGeV = 0` when the contents are per-bin integrals instead.
struct FluxSpec {
  std::string file;                       ///< ROOT file with the histograms
  std::map<int, std::string> histograms;  ///< neutrino PDG code -> histogram name
  EnergyUnit units = EnergyUnit::GeV;     ///< unit of the histogram axes in the file
  double potNorm = 1.0e21;                ///< POT the histograms are normalised to
  double binNormGeV = 0.05;               ///< energy width the density refers to (0 = per-bin content)

  FluxSpec() = default;
  FluxSpec(std::string file_, std::map<int, std::string> histograms_)
      : file(std::move(file_)), histograms(std::move(histograms_)) {}

  bool has(int pdg) const { return histograms.count(pdg) != 0; }
  /// Histogram name for a species; throws std::out_of_range if absent.
  const std::string& histogram(int pdg) const;
  std::vector<int> pdgs() const;
  std::vector<NuSpecies> species() const;

  /// The T2K public-release layout: enu_<detector>_{numu,numub,nue,nueb}.
  static FluxSpec t2k(const std::string& file, const std::string& detector = "sk");

  /// Copy pointing at another file (same histogram names and conventions).
  FluxSpec withFile(const std::string& newFile) const;

  std::string describe() const;
};

/// The flux covariance input: a fractional covariance matrix with one TAxis
/// object per (detector, beam mode, species) block naming the energy bins.
struct FluxCovarianceSpec {
  std::string file;
  std::string matrix = "total_flux_cov";
  std::string farPrefix = "sk";       ///< prefix of the far-detector TAxis names
  std::string nearPrefix = "nd280";   ///< prefix of the near-detector TAxis names
  double ndConstraint = 0.0;          ///< 0 = unconstrained; s in C - s B^T A^-1 B otherwise
  EnergyUnit units = EnergyUnit::GeV; ///< unit of the TAxis edges

  bool valid() const { return !file.empty(); }
};

}  // namespace zaku

#endif  // ZAKU_CORE_FLUXSPEC_H
