#include "zaku/core/Types.h"

#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace zaku {

namespace {

std::string lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return s;
}

/// Strip characters that only differ between spellings: '_', '-', ' '.
std::string squeeze(const std::string& s) {
  std::string out;
  for (char c : s) {
    if (c == '_' || c == '-' || c == ' ') continue;
    out += c;
  }
  return out;
}

const char* const kSymbols[] = {
    "n",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si",
    "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
    "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
    "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
    "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf",
    "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};

}  // namespace

// --- NuSpecies -------------------------------------------------------------

std::string NuSpecies::name() const {
  const char* base = flavour == Flavour::NuE ? "nue" : flavour == Flavour::NuMu ? "numu" : "nutau";
  return isAnti ? std::string(base) + "bar" : std::string(base);
}

std::string NuSpecies::latex() const {
  const char* base = flavour == Flavour::NuE ? "#nu_{e}" : flavour == Flavour::NuMu ? "#nu_{#mu}" : "#nu_{#tau}";
  if (!isAnti) return base;
  const char* abase = flavour == Flavour::NuE     ? "#bar{#nu}_{e}"
                      : flavour == Flavour::NuMu  ? "#bar{#nu}_{#mu}"
                                                  : "#bar{#nu}_{#tau}";
  return abase;
}

std::string NuSpecies::t2kSuffix() const {
  const char* base = flavour == Flavour::NuE ? "nue" : flavour == Flavour::NuMu ? "numu" : "nutau";
  return isAnti ? std::string(base) + "b" : std::string(base);
}

NuSpecies NuSpecies::fromPdg(int pdg) {
  const int a = pdg < 0 ? -pdg : pdg;
  if (a != 12 && a != 14 && a != 16) {
    throw std::invalid_argument("NuSpecies::fromPdg: " + std::to_string(pdg) +
                                " is not a neutrino PDG code (+-12, +-14, +-16)");
  }
  return NuSpecies(static_cast<Flavour>(a), pdg < 0);
}

NuSpecies NuSpecies::fromName(const std::string& input) {
  // Numeric PDG code as text ("14", "-12", " +16 ").
  {
    std::string t = input;
    t.erase(0, t.find_first_not_of(" \t"));
    if (!t.empty() && (t[0] == '-' || t[0] == '+' || std::isdigit(static_cast<unsigned char>(t[0])))) {
      try {
        std::size_t used = 0;
        const int code = std::stoi(t, &used);
        if (t.find_first_not_of(" \t", used) == std::string::npos) return fromPdg(code);
      } catch (const std::invalid_argument&) {
        // not a number: fall through to the name matching below
      } catch (const std::out_of_range&) {
        // not a sensible number either
      }
    }
  }

  std::string s = squeeze(lower(input));
  if (s.empty()) throw std::invalid_argument("NuSpecies::fromName: empty name");

  bool anti = false;
  auto stripPrefix = [&](const std::string& p) {
    if (s.compare(0, p.size(), p) == 0) {
      s.erase(0, p.size());
      return true;
    }
    return false;
  };
  auto stripSuffix = [&](const std::string& p) {
    if (s.size() >= p.size() && s.compare(s.size() - p.size(), p.size(), p) == 0) {
      s.erase(s.size() - p.size());
      return true;
    }
    return false;
  };
  if (stripPrefix("anti")) anti = true;
  if (stripSuffix("bar")) anti = true;
  else if (stripSuffix("b")) anti = true;  // T2K suffix: numub, nueb

  Flavour f;
  if (s == "numu" || s == "nu14" || s == "num" || s == "mu" || s == "muon")
    f = Flavour::NuMu;
  else if (s == "nue" || s == "nu12" || s == "e" || s == "electron")
    f = Flavour::NuE;
  else if (s == "nutau" || s == "nu16" || s == "tau")
    f = Flavour::NuTau;
  else
    throw std::invalid_argument("NuSpecies::fromName: unrecognised neutrino name '" + input + "'");
  return NuSpecies(f, anti);
}

// --- BeamMode --------------------------------------------------------------

std::string toString(BeamMode m) { return m == BeamMode::FHC ? "FHC" : "RHC"; }

BeamMode beamModeFromString(const std::string& input) {
  const std::string s = squeeze(lower(input));
  if (s == "fhc" || s == "numode" || s == "nu" || s == "plus250ka" || s == "neutrino" || s == "forward")
    return BeamMode::FHC;
  if (s == "rhc" || s == "anumode" || s == "anu" || s == "minus250ka" || s == "antineutrino" || s == "reverse")
    return BeamMode::RHC;
  throw std::invalid_argument("beamModeFromString: unrecognised beam mode '" + input + "' (expected FHC or RHC)");
}

std::vector<BeamMode> allBeamModes() { return {BeamMode::FHC, BeamMode::RHC}; }

// --- Target ----------------------------------------------------------------

std::string elementSymbol(int Z) {
  const int n = static_cast<int>(sizeof(kSymbols) / sizeof(kSymbols[0]));
  if (Z >= 0 && Z < n) return kSymbols[Z];
  return "Z" + std::to_string(Z);
}

std::string Target::name() const { return elementSymbol(Z) + std::to_string(A); }

}  // namespace zaku
