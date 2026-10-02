#ifndef ZAKU_DETECTOR_SHIFTS_H
#define ZAKU_DETECTOR_SHIFTS_H
#include <string>
#include <vector>

#include "zaku/detector/Reco.h"

namespace zaku {

/// Which reconstructed quantity a shift systematic moves.
enum class ShiftTarget : int {
  LeptonMomentum = 0,   ///< PLep_rec (both estimators follow)
  LeptonAngle = 1,      ///< theta (CosLep_rec, Theta_rec; the CCQE estimators follow)
  LeptonEnergy = 2,     ///< the calorimetric lepton term Elep_rec
  HadronicEnergy = 3,   ///< Ehad_rec (calorimetric estimator follows)
  RecoEnergy = 4        ///< every Enu_rec estimator directly
};
std::string toString(ShiftTarget t);
ShiftTarget shiftTargetFromString(const std::string& s);

enum class ShiftKind : int { Scale = 0, Offset = 1 };
std::string toString(ShiftKind k);
ShiftKind shiftKindFromString(const std::string& s);

/// Which samples a shift applies to, by their lepton hypothesis.
enum class ShiftApplyTo : int { All = 0, Muon = 1, Electron = 2 };
std::string toString(ShiftApplyTo a);
ShiftApplyTo shiftApplyToFromString(const std::string& s);

/// A detector systematic that moves reconstructed quantities rather than
/// scaling event counts: an energy or momentum scale (kind `scale`, sigma as a
/// fraction) or an offset (kind `offset`, sigma in GeV, or degrees for the angle).
struct ShiftParameter {
  std::string name;
  ShiftTarget target = ShiftTarget::LeptonMomentum;
  ShiftKind kind = ShiftKind::Scale;
  double sigma = 0.0;
  ShiftApplyTo applyTo = ShiftApplyTo::All;

  bool appliesTo(Hypothesis h, int pdgLepTrue) const;
  std::string describe() const;
};

/// Applies a shift to the reconstructed observables of one event and
/// re-derives everything downstream (Theta_rec, both CCQE estimators, the
/// calorimetric sums), so a shifted universe is one consistent reconstruction.
namespace ShiftInjector {

/// Apply `nSigma` sigmas of `p` to `ro` in place. `isRHC` and `bindingEnergy`
/// feed the CCQE formula; `cc` keeps the calorimetric lepton term at 0 for NC.
/// Does nothing when the shift does not apply to the sample hypothesis.
void apply(RecoObservables& ro, const ShiftParameter& p, double nSigma, bool isRHC, double bindingEnergy, bool cc,
           Hypothesis sampleHypothesis, int pdgLepTrue);

/// Recompute the derived quantities of `ro` from PLep_rec, CosLep_rec,
/// Ehad_rec and the lepton terms (what apply() calls after moving an input).
void rederive(RecoObservables& ro, bool isRHC, double bindingEnergy, bool cc, bool leptonTermsFromMomentum);

}  // namespace ShiftInjector

}  // namespace zaku

#endif  // ZAKU_DETECTOR_SHIFTS_H
