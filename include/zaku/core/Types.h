#ifndef ZAKU_CORE_TYPES_H
#define ZAKU_CORE_TYPES_H
#include <string>
#include <vector>

// Fundamental value types shared across Zaku. No ROOT, no heavy dependencies:
// these must stay cheap to include everywhere and simple to use from PyROOT.
namespace zaku {

/// Neutrino flavour, valued as the PDG code of the neutrino.
enum class Flavour : int { NuE = 12, NuMu = 14, NuTau = 16 };

/// A neutrino species: flavour plus particle/antiparticle.
struct NuSpecies {
  Flavour flavour = Flavour::NuMu;
  bool isAnti = false;

  NuSpecies() = default;
  NuSpecies(Flavour f, bool anti) : flavour(f), isAnti(anti) {}

  /// PDG code: 12/14/16, negative for antineutrinos.
  int pdg() const { return isAnti ? -static_cast<int>(flavour) : static_cast<int>(flavour); }
  /// Canonical name: numu, numubar, nue, nuebar, nutau, nutaubar.
  std::string name() const;
  /// ROOT-latex label, e.g. "#bar{#nu}_{#mu}".
  std::string latex() const;
  /// Suffix used by the T2K flux release histograms: numu, numub, nue, nueb.
  std::string t2kSuffix() const;

  /// Build from a PDG code. Throws std::invalid_argument for anything but +-12, +-14, +-16.
  static NuSpecies fromPdg(int pdg);
  /// Build from a name. Accepts the canonical names, the T2K suffixes,
  /// "anti-numu"/"anti_numu"/"nu_mu_bar" style spellings, and PDG codes as
  /// text. Case-insensitive. Throws std::invalid_argument otherwise.
  static NuSpecies fromName(const std::string& s);

  bool operator==(const NuSpecies& o) const { return flavour == o.flavour && isAnti == o.isAnti; }
  bool operator!=(const NuSpecies& o) const { return !(*this == o); }
};

/// Horn polarity of the beam: neutrino-enhanced (FHC) or antineutrino-enhanced (RHC).
enum class BeamMode : int { FHC = 0, RHC = 1 };

/// "FHC" or "RHC".
std::string toString(BeamMode m);
/// Accepts FHC/RHC, numode/anumode, nu/anu, plus250kA/minus250kA (case-insensitive).
BeamMode beamModeFromString(const std::string& s);
inline bool isRHC(BeamMode m) { return m == BeamMode::RHC; }
/// {FHC, RHC}
std::vector<BeamMode> allBeamModes();

/// A nuclear target species inside a detector material.
struct Target {
  int Z = 0;
  int A = 0;
  double massFraction = 1.0;  ///< fraction of the detector mass in this nucleus

  Target() = default;
  Target(int z, int a, double f) : Z(z), A(a), massFraction(f) {}

  /// PDG nuclear code 10LZZZAAAI (ground state): O-16 -> 1000080160.
  int nuclearPdg() const { return 1000000000 + Z * 10000 + A * 10; }
  /// Element symbol plus mass number, e.g. "O16", "H1", "Ar40".
  std::string name() const;
};

/// Element symbol for a proton number (1..118), or "Z<n>" if unknown.
std::string elementSymbol(int Z);

}  // namespace zaku

#endif  // ZAKU_CORE_TYPES_H
