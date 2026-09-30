#ifndef ZAKU_OSC_OSCILLATOR_H
#define ZAKU_OSC_OSCILLATOR_H
#include <string>
#include <vector>

#include "zaku/event/FlatEvent.h"

class TH1D;

namespace zaku {

/// Three-flavour oscillation parameters (PDG convention). Δm² in eV², δ in radians.
struct OscillationParameters {
  double s12sq = 0.307;
  double s13sq = 0.0218;
  double s23sq = 0.528;
  double deltaCP = -1.601;
  double dmsq21 = 7.53e-5;
  double dmsq32 = 2.509e-3;   ///< magnitude; the ordering sets the sign of Δm²31
  bool normalOrdering = true;

  /// The T2K "Asimov A" point (normal ordering), the legacy default.
  static OscillationParameters asimovA();
  /// Δm²31 as NuFast wants it: Δm²32 + Δm²21, negative for inverted ordering.
  double dmsq31() const;
  std::string describe() const;
};

/// NuFast matter-effect probabilities for one baseline and matter density.
class Oscillator {
public:
  Oscillator(const OscillationParameters& p, double baseline_km, double rhoYe, int nNewton = 0);

  const OscillationParameters& parameters() const { return params_; }
  double baseline_km() const { return baseline_km_; }
  double rhoYe() const { return rhoYe_; }

  /// P(from -> to) for a neutrino of energy E (GeV). Antineutrinos are
  /// signalled by a negative `from` PDG code (both codes must have the same sign).
  double probability(int fromPdg, int toPdg, double E_GeV) const;
  /// All nine probabilities for one energy; index 0 = e, 1 = mu, 2 = tau.
  std::vector<double> matrix(double E_GeV, bool antineutrino) const;

  /// The legacy OscWeight rule: CC events get P(flux flavour -> interaction
  /// flavour) with the sign of the flux species; NC events get 1 on a
  /// disappearance production and 0 on an appearance production (their rate
  /// is already carried by the disappearance production of the same flux).
  double oscWeight(const FlatEvent& ev, const ProductionContext& ctx) const;

  /// Uniform true-energy grid histogram ("energyBinning") for GUNDAM's
  /// tabulated NuOscillator dials: `nBins` bins on [eMin, eMax] GeV, contents empty.
  static TH1D* gridHistogram(int nBins = 2000, double eMin = 0.0, double eMax = 10.0,
                             const std::string& name = "energyBinning");
  /// Write the grid histogram into a ROOT file (the legacy MakeOscGrid output).
  static void writeGrid(const std::string& file, int nBins = 2000, double eMin = 0.0, double eMax = 10.0);

private:
  OscillationParameters params_;
  double baseline_km_ = 295.0;
  double rhoYe_ = 1.3;
  int nNewton_ = 0;
};

}  // namespace zaku

#endif  // ZAKU_OSC_OSCILLATOR_H
