// NuFast: fast three-flavour long-baseline oscillation probabilities.
// P. Denton, MIT licence (see LICENSE). Vendored into Zaku; see README.md.
#ifndef ZAKU_THIRD_PARTY_NUFAST_LBL_H
#define ZAKU_THIRD_PARTY_NUFAST_LBL_H
#include <array>
#include <vector>

namespace NuFast {

// Some constants
constexpr double const eVsqkm_to_GeV_over4 = 1e-9 / 1.97327e-7 * 1e3 / 4;
constexpr double const YerhoE2a = 1.52588e-4;

// Probability_Matter_LBL calculates all nine oscillation probabilities across a
// vector of energies, including the matter effect, in an optimized, fast, and
// efficient way. The precision can be controlled with N_Newton. For many
// applications N_Newton=0 will be enough, but many years of DUNE or HK-LBL may
// require N_Newton=1 for a high precision numerical scan. The code is standalone.
//
// Inputs:
//   mixing angles (usual parameterization)
//   phase (usual parameterization, radians)
//   Delta msq's (eV^2)
//   Make Dmsq31 positive/negative for the NO/IO
//   L (km)
//   E (GeV) positive for neutrinos, negative for antineutrinos. Vector
//   rhoYe density (g/cc) times electron fraction, typically around 3 * 0.5 for the Earth's crust
//   N_Newton: number of Newton's method iterations to do. Should be zero, one, two (can be higher);
//             if it is negative, the code will use the exact expression
// Outputs:
//   probs_returned: a vector (over energies) of all nine oscillation probabilities:
//     e.g. probs_returned[i][1][0] is P(nu_mu->nu_e) for energy E[i]; the vector must be allocated first
void Probability_Matter_LBL(double s12sq, double s13sq, double s23sq, double delta, double Dmsq21, double Dmsq31,
                            double L, const std::vector<double>& E, double rhoYe, int N_Newton,
                            std::vector<std::array<std::array<double, 3>, 3>>& probs_returned);

void Probability_Vacuum_LBL(double s12sq, double s13sq, double s23sq, double delta, double Dmsq21, double Dmsq31,
                            double L, const std::vector<double>& E,
                            std::vector<std::array<std::array<double, 3>, 3>>& probs_returned);

}  // namespace NuFast

#endif  // ZAKU_THIRD_PARTY_NUFAST_LBL_H
