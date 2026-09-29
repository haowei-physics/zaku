#ifndef ZAKU_CORE_UNITS_H
#define ZAKU_CORE_UNITS_H
#include <string>

class TH1;

// Zaku uses one unit per quantity everywhere: GeV for energy, momentum and
// mass (GeV^2 for Q^2 and W^2), degrees for reconstructed angles, km for
// baselines, eV^2 for mass splittings, kt for detector mass, POT for exposure.
// Every histogram Zaku reads, writes or references has its kinematic axes in
// GeV. Inputs in MeV must be declared as such and are converted on load.
// See docs/units.md.
namespace zaku {

/// Energy unit of an external input.
enum class EnergyUnit : int { GeV = 0, MeV = 1 };

namespace units {
// Multiply a value expressed in the named unit by the constant to get Zaku's unit.
inline constexpr double GeV = 1.0;
inline constexpr double MeV = 1.0e-3;
inline constexpr double keV = 1.0e-6;
inline constexpr double km = 1.0;
inline constexpr double m = 1.0e-3;
inline constexpr double cm = 1.0e-5;
inline constexpr double kt = 1.0;
inline constexpr double eV2 = 1.0;
inline constexpr double POT = 1.0;

inline constexpr double pi = 3.14159265358979323846;
inline constexpr double deg = pi / 180.0;  ///< degrees -> radians

inline constexpr double kt_to_g = 1.0e9;           ///< grams per kilotonne
inline constexpr double Avogadro = 6.02214076e23;  ///< per mole

// Particle masses in GeV. These are the values the legacy macros used; they
// are kept unchanged so that the reference production is reproducible.
inline constexpr double m_e = 0.000511;
inline constexpr double m_mu = 0.105658;
inline constexpr double m_p = 0.93827;
inline constexpr double m_n = 0.93957;
inline constexpr double m_pi = 0.13957;
inline constexpr double m_pi0 = 0.13498;
}  // namespace units

/// "GeV" or "MeV".
std::string toString(EnergyUnit u);
/// Accepts GeV/MeV in any case. Throws std::invalid_argument otherwise.
EnergyUnit energyUnitFromString(const std::string& s);

/// Convert a value given in `from` to GeV.
double toGeV(double value, EnergyUnit from);
/// Convert a value in GeV to `to`.
double fromGeV(double valueGeV, EnergyUnit to);

/// Return a new histogram (TH1D/TH2D/TH3D, owned by the caller, not attached
/// to any directory) whose selected axes are converted from `from` to GeV.
/// `axes` is a bit mask: 1 = x, 2 = y, 4 = z. If `contentIsDensity` the
/// contents (and errors) are rescaled so that the density per unit energy is
/// preserved, i.e. a per-MeV density becomes a per-GeV density. With
/// `from == EnergyUnit::GeV` this is a plain clone.
TH1* histogramToGeV(const TH1& h, EnergyUnit from, bool contentIsDensity, int axes = 1,
                    const std::string& newName = "");

/// True if the x axis of `h` extends beyond `maxGeV`, which is how an
/// undeclared MeV axis usually shows up (30000 instead of 30).
bool looksLikeMeVAxis(const TH1& h, double maxGeV = 1000.0);

/// Throw std::runtime_error naming `context` if `looksLikeMeVAxis(h)`.
void requireGeVAxis(const TH1& h, const std::string& context, double maxGeV = 1000.0);

}  // namespace zaku

#endif  // ZAKU_CORE_UNITS_H
