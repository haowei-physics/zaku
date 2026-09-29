#ifndef ZAKU_EVENT_FLATEVENT_H
#define ZAKU_EVENT_FLATEVENT_H
#include <string>
#include <vector>

namespace zaku {

/// One generator event as read from a NUISANCE "GenericVectors" flat tree
/// (`FlatTree_VARS`), the generator-neutral interchange format Zaku works on.
/// Energies and momenta are GeV; `E[]` includes the rest mass. The final-state
/// list contains the outgoing lepton (or, for NC, the outgoing neutrino).
struct FlatEvent {
  long long entry = -1;  ///< entry number in the flat tree

  int Mode = 0;          ///< generator interaction mode (NEUT convention for NEUT files)
  bool cc = false;       ///< charged current
  int PDGnu = 0;
  float Enu_true = 0.f;
  int tgt = 0, tgta = 0, tgtz = 0;
  int PDGLep = 0;        ///< outgoing lepton PDG (NC: the outgoing neutrino); 0 when absent
  float ELep = 0.f;      ///< total lepton energy
  float CosLep = 0.f;
  float Q2 = 0.f, q0 = 0.f, q3 = 0.f;
  float Enu_QE = 0.f, Q2_QE = 0.f;
  float W = 0.f, x = 0.f, y = 0.f;
  double fScaleFactor = 0.0;  ///< <sigma>/N_gen, cm^2 per nucleon per event
  float Weight = 1.f;

  int nfsp = 0;
  std::vector<int> pdg;
  std::vector<float> px, py, pz, E;

  bool flagCCINC = false, flagNCINC = false, flagCCQE = false, flagCC0pi = false, flagCCQELike = false;
  bool flagNCEL = false, flagNC0pi = false, flagCCcoh = false, flagNCcoh = false;
  bool flagCC1pip = false, flagNC1pip = false, flagCC1pim = false, flagNC1pim = false;
  bool flagCC1pi0 = false, flagNC1pi0 = false;

  /// True when the generator produced no final-state particles at all
  /// (NEUT single-pion/gamma modes whose decay failed): 1505 of the 200k
  /// reference events. T2KReWeight writes no weight row for them and the
  /// legacy macros drop them.
  bool emptyFinalState() const { return nfsp == 0; }
  /// The legacy FakeRecon drop rule: Q2 sentinel or a CC event without lepton.
  bool legacyDrop() const { return Q2 < -999.0f || (cc && (PDGLep == 0 || ELep < 0.f)); }

  /// Kinetic energy of final-state particle i (E - m, with m from E^2 - p^2).
  double kineticEnergy(int i) const;
  double momentum(int i) const;
  double mass(int i) const;

  std::string describe() const;
};

/// Per-production context every event inherits: which beam mode and flux
/// species it was generated on, and which species interacted.
struct ProductionContext {
  bool isRHC = false;
  int fluxPdg = 14;
  int interactionPdg = 14;
  bool isAppearance() const { return (fluxPdg < 0 ? -fluxPdg : fluxPdg) != (interactionPdg < 0 ? -interactionPdg : interactionPdg); }
};

/// Derived truth quantities every response table and cut can use. All energies GeV.
struct TruthVariables {
  double PLep_true = 0.0;   ///< outgoing lepton momentum
  double Ehad_true = 0.0;   ///< visible hadronic energy, see docs/physics.md
  double Evis_true = 0.0;   ///< Ehad_true + lepton energy (CC) or Ehad_true (NC)
  double W2 = 0.0;          ///< W^2 (GeV^2)
  int nFS_pip = 0, nFS_pim = 0, nFS_pi0 = 0, nFS_p = 0, nFS_n = 0, nFS_gamma = 0;
  int nFS_mup = 0, nFS_mum = 0, nFS_ep = 0, nFS_em = 0, nFS_kp = 0, nFS_km = 0, nFS_k0 = 0, nFS_other = 0;
  int nFS_pi = 0;           ///< charged pions
  bool isRHC = false;
  int flux_pdg = 14;
  bool isAppearance = false;

  /// Compute from an event. `neutronFraction` is the visible fraction of neutron kinetic energy.
  static TruthVariables compute(const FlatEvent& ev, const ProductionContext& ctx, double neutronFraction = 0.0);

  /// Names of every variable a formula can reference, in a fixed order:
  /// the flat-tree scalars first, then the derived ones above.
  static const std::vector<std::string>& names();
  /// Values in the same order as names().
  std::vector<double> values(const FlatEvent& ev) const;
};

}  // namespace zaku

#endif  // ZAKU_EVENT_FLATEVENT_H
