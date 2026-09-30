#ifndef ZAKU_DETECTOR_RECO_H
#define ZAKU_DETECTOR_RECO_H
#include <string>

class TRandom;

namespace zaku {

class Smearing;
class VariableSet;
struct FlatEvent;
struct TruthVariables;

/// Which lepton the reconstruction assumes for a sample.
enum class Hypothesis : int { Muon = 0, Electron = 1, True = 2 };
std::string toString(Hypothesis h);
Hypothesis hypothesisFromString(const std::string& s);

/// Which neutrino-energy estimator fills `Enu_rec`.
enum class EnergyEstimator : int { CCQE = 0, Calorimetric = 1 };
std::string toString(EnergyEstimator e);
EnergyEstimator energyEstimatorFromString(const std::string& s);

/// The reconstructed observables of one event, before a sample hypothesis is chosen.
/// Energies GeV, `Theta_rec` degrees.
struct RecoObservables {
  float PLep_rec = 0.f;      ///< smeared outgoing-lepton momentum
  float CosLep_rec = 0.f;    ///< smeared cos(theta) to the beam
  float Theta_rec = 0.f;     ///< acos(CosLep_rec) in degrees (float acos, as the legacy macro)
  double Ehad_rec = 0.0;     ///< smeared visible hadronic energy
  double Elep_rec_mu = 0.0;  ///< calorimetric lepton energy under the muon hypothesis
  double Elep_rec_e = 0.0;   ///< ... under the electron hypothesis
  float Enu_rec_ccqe_mu = -1.f, Q2_rec_ccqe_mu = -1.f;  ///< CCQE formula, muon hypothesis
  float Enu_rec_ccqe_e = -1.f, Q2_rec_ccqe_e = -1.f;    ///< CCQE formula, electron hypothesis
  double Enu_rec_calo_mu = 0.0, Enu_rec_calo_e = 0.0;   ///< Elep_rec + Ehad_rec
  int true_ME = 0;           ///< decay-electron parents in the final state
  int recon_ME = 0;          ///< tagged decay electrons

  /// The values the chosen hypothesis and estimator imply.
  float enuRec(Hypothesis h, EnergyEstimator e, int pdgLepTrue) const;
  float q2Rec(Hypothesis h, int pdgLepTrue) const;
  double elepRec(Hypothesis h, int pdgLepTrue) const;
  /// The PDG code written as PDGLep_rec: +-13 or +-11 by hypothesis and beam mode; the true lepton for True.
  int pdgLepRec(Hypothesis h, bool isRHC, int pdgLepTrue) const;
};

/// The reconstruction formulas, written to reproduce the legacy macros bit for bit.
namespace reco {

/// Outgoing-lepton momentum from its total energy: sqrt(ELep^2 - m^2) with the
/// float ELep squared in float precision, as the legacy macro did.
double leptonMomentum(float ELep, int pdgLep);

/// CCQE neutrino energy and Q2 under a lepton mass hypothesis. The struck
/// nucleon is a neutron in neutrino mode and a proton in antineutrino mode
/// (beam-mode choice, as the legacy macro), with binding energy Eb. Returns
/// false (and -1, -1) when the formula has no positive solution.
bool ccqe(float PLep_rec, float CosLep_rec, double mLep, bool isRHC, double Eb, float& Enu_rec, float& Q2_rec);

/// acos in float precision, converted to degrees.
float thetaDegrees(float cosTheta);

/// Total lepton energy from momentum and mass (float momentum squared, as legacy).
double leptonEnergy(float PLep_rec, double mLep);

/// Run the smearing, the decay-electron tagging and the energy estimators
/// for one event. `legacyStream` draws, in this order, Gaus(p), Gaus(cos) and
/// one Uniform per decay-electron parent (the legacy sequence); `extraStream`
/// draws everything the legacy chain did not have (hadronic and calorimetric
/// lepton energy smearing), so the legacy quantities stay reproducible.
RecoObservables reconstruct(const FlatEvent& ev, const TruthVariables& truth, const VariableSet& vars,
                            const Smearing& smearing, bool isRHC, TRandom& legacyStream, TRandom& extraStream,
                            double bindingEnergy = 0.027);

}  // namespace reco

}  // namespace zaku

#endif  // ZAKU_DETECTOR_RECO_H
