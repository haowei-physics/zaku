#include "zaku/detector/Shifts.h"

#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "zaku/core/Units.h"

namespace zaku {

std::string toString(ShiftTarget t) {
  switch (t) {
    case ShiftTarget::LeptonMomentum: return "lepton_momentum";
    case ShiftTarget::LeptonAngle: return "lepton_angle";
    case ShiftTarget::LeptonEnergy: return "lepton_energy";
    case ShiftTarget::HadronicEnergy: return "hadronic_energy";
    case ShiftTarget::RecoEnergy: return "reco_energy";
  }
  return "?";
}

ShiftTarget shiftTargetFromString(const std::string& s) {
  if (s == "lepton_momentum" || s == "momentum" || s == "PLep_rec") return ShiftTarget::LeptonMomentum;
  if (s == "lepton_angle" || s == "angle" || s == "theta") return ShiftTarget::LeptonAngle;
  if (s == "lepton_energy" || s == "Elep_rec") return ShiftTarget::LeptonEnergy;
  if (s == "hadronic_energy" || s == "Ehad_rec") return ShiftTarget::HadronicEnergy;
  if (s == "reco_energy" || s == "energy" || s == "Enu_rec") return ShiftTarget::RecoEnergy;
  throw std::invalid_argument("shiftTargetFromString: unknown target '" + s +
                              "' (lepton_momentum, lepton_angle, lepton_energy, hadronic_energy, reco_energy)");
}

std::string toString(ShiftKind k) { return k == ShiftKind::Scale ? "scale" : "offset"; }

ShiftKind shiftKindFromString(const std::string& s) {
  if (s == "scale" || s == "fractional") return ShiftKind::Scale;
  if (s == "offset" || s == "shift" || s == "absolute") return ShiftKind::Offset;
  throw std::invalid_argument("shiftKindFromString: unknown kind '" + s + "' (scale, offset)");
}

std::string toString(ShiftApplyTo a) {
  switch (a) {
    case ShiftApplyTo::All: return "all";
    case ShiftApplyTo::Muon: return "muon";
    case ShiftApplyTo::Electron: return "electron";
  }
  return "all";
}

ShiftApplyTo shiftApplyToFromString(const std::string& s) {
  if (s == "all" || s == "any") return ShiftApplyTo::All;
  if (s == "muon" || s == "mu" || s == "mu-like") return ShiftApplyTo::Muon;
  if (s == "electron" || s == "e" || s == "e-like") return ShiftApplyTo::Electron;
  throw std::invalid_argument("shiftApplyToFromString: unknown value '" + s + "' (all, muon, electron)");
}

bool ShiftParameter::appliesTo(Hypothesis h, int pdgLepTrue) const {
  if (applyTo == ShiftApplyTo::All) return true;
  bool muon;
  if (h == Hypothesis::Muon) muon = true;
  else if (h == Hypothesis::Electron) muon = false;
  else muon = (pdgLepTrue < 0 ? -pdgLepTrue : pdgLepTrue) == 13;
  return applyTo == ShiftApplyTo::Muon ? muon : !muon;
}

std::string ShiftParameter::describe() const {
  std::ostringstream os;
  os << name << ": " << toString(kind) << " of " << toString(target) << ", sigma " << sigma
     << (kind == ShiftKind::Scale ? " (fraction)" : (target == ShiftTarget::LeptonAngle ? " deg" : " GeV"))
     << ", applies to " << toString(applyTo);
  return os.str();
}

namespace ShiftInjector {

void rederive(RecoObservables& ro, bool isRHC, double Eb, bool cc, bool leptonTermsFromMomentum) {
  ro.Theta_rec = reco::thetaDegrees(ro.CosLep_rec);
  reco::ccqe(ro.PLep_rec, ro.CosLep_rec, units::m_mu, isRHC, Eb, ro.Enu_rec_ccqe_mu, ro.Q2_rec_ccqe_mu);
  reco::ccqe(ro.PLep_rec, ro.CosLep_rec, units::m_e, isRHC, Eb, ro.Enu_rec_ccqe_e, ro.Q2_rec_ccqe_e);
  if (leptonTermsFromMomentum) {
    ro.Elep_rec_mu = cc ? reco::leptonEnergy(ro.PLep_rec, units::m_mu) : 0.0;
    ro.Elep_rec_e = cc ? reco::leptonEnergy(ro.PLep_rec, units::m_e) : 0.0;
  }
  ro.Enu_rec_calo_mu = ro.Elep_rec_mu + ro.Ehad_rec;
  ro.Enu_rec_calo_e = ro.Elep_rec_e + ro.Ehad_rec;
}

namespace {
double moved(double value, ShiftKind kind, double delta) { return kind == ShiftKind::Scale ? value * (1.0 + delta) : value + delta; }
}  // namespace

void apply(RecoObservables& ro, const ShiftParameter& p, double nSigma, bool isRHC, double Eb, bool cc,
           Hypothesis hyp, int pdgLepTrue) {
  if (nSigma == 0.0 || !p.appliesTo(hyp, pdgLepTrue)) return;
  const double delta = nSigma * p.sigma;
  switch (p.target) {
    case ShiftTarget::LeptonMomentum: {
      const float pOld = ro.PLep_rec;
      ro.PLep_rec = static_cast<float>(std::max(0.0, moved(pOld, p.kind, delta)));
      // Keep whatever smearing the calorimetric lepton terms carried: scale
      // them by the ratio of the unsmeared energies.
      const double eMuOld = reco::leptonEnergy(pOld, units::m_mu), eEOld = reco::leptonEnergy(pOld, units::m_e);
      if (eMuOld > 0.0) ro.Elep_rec_mu *= reco::leptonEnergy(ro.PLep_rec, units::m_mu) / eMuOld;
      if (eEOld > 0.0) ro.Elep_rec_e *= reco::leptonEnergy(ro.PLep_rec, units::m_e) / eEOld;
      rederive(ro, isRHC, Eb, cc, false);
      return;
    }
    case ShiftTarget::LeptonAngle: {
      const double c = std::min(1.0, std::max(-1.0, static_cast<double>(ro.CosLep_rec)));
      double theta = std::acos(c);  // radians
      theta = p.kind == ShiftKind::Scale ? theta * (1.0 + delta) : theta + delta * TMath::DegToRad();
      theta = std::min(TMath::Pi(), std::max(0.0, theta));
      ro.CosLep_rec = static_cast<float>(std::cos(theta));
      rederive(ro, isRHC, Eb, cc, false);
      return;
    }
    case ShiftTarget::LeptonEnergy:
      ro.Elep_rec_mu = std::max(0.0, moved(ro.Elep_rec_mu, p.kind, delta));
      ro.Elep_rec_e = std::max(0.0, moved(ro.Elep_rec_e, p.kind, delta));
      ro.Enu_rec_calo_mu = ro.Elep_rec_mu + ro.Ehad_rec;
      ro.Enu_rec_calo_e = ro.Elep_rec_e + ro.Ehad_rec;
      return;
    case ShiftTarget::HadronicEnergy:
      ro.Ehad_rec = std::max(0.0, moved(ro.Ehad_rec, p.kind, delta));
      ro.Enu_rec_calo_mu = ro.Elep_rec_mu + ro.Ehad_rec;
      ro.Enu_rec_calo_e = ro.Elep_rec_e + ro.Ehad_rec;
      return;
    case ShiftTarget::RecoEnergy:
      if (ro.Enu_rec_ccqe_mu > 0.f) ro.Enu_rec_ccqe_mu = static_cast<float>(std::max(0.0, moved(ro.Enu_rec_ccqe_mu, p.kind, delta)));
      if (ro.Enu_rec_ccqe_e > 0.f) ro.Enu_rec_ccqe_e = static_cast<float>(std::max(0.0, moved(ro.Enu_rec_ccqe_e, p.kind, delta)));
      ro.Enu_rec_calo_mu = std::max(0.0, moved(ro.Enu_rec_calo_mu, p.kind, delta));
      ro.Enu_rec_calo_e = std::max(0.0, moved(ro.Enu_rec_calo_e, p.kind, delta));
      return;
  }
}

}  // namespace ShiftInjector

}  // namespace zaku
