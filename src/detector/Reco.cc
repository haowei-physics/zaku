#include "zaku/detector/Reco.h"

#include <TMath.h>
#include <TRandom.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "zaku/core/Units.h"
#include "zaku/detector/Smearing.h"
#include "zaku/detector/Variables.h"
#include "zaku/event/FlatEvent.h"

namespace zaku {

std::string toString(Hypothesis h) {
  switch (h) {
    case Hypothesis::Muon: return "muon";
    case Hypothesis::Electron: return "electron";
    case Hypothesis::True: return "true";
  }
  return "muon";
}

Hypothesis hypothesisFromString(const std::string& s) {
  if (s == "muon" || s == "mu" || s == "mu-like") return Hypothesis::Muon;
  if (s == "electron" || s == "e" || s == "e-like") return Hypothesis::Electron;
  if (s == "true" || s == "truth") return Hypothesis::True;
  throw std::invalid_argument("hypothesisFromString: unknown hypothesis '" + s + "' (muon, electron, true)");
}

std::string toString(EnergyEstimator e) { return e == EnergyEstimator::CCQE ? "ccqe" : "calorimetric"; }

EnergyEstimator energyEstimatorFromString(const std::string& s) {
  if (s == "ccqe" || s == "CCQE" || s == "t2k" || s == "t2k_like") return EnergyEstimator::CCQE;
  if (s == "calorimetric" || s == "calo" || s == "nova" || s == "nova_like") return EnergyEstimator::Calorimetric;
  throw std::invalid_argument("energyEstimatorFromString: unknown estimator '" + s + "' (ccqe, calorimetric)");
}

namespace {

double leptonMass(int pdg) {
  const int a = pdg < 0 ? -pdg : pdg;
  if (a == 13) return units::m_mu;
  if (a == 11) return units::m_e;
  return 0.0;
}

bool hypothesisIsMuon(Hypothesis h, int pdgLepTrue) {
  if (h == Hypothesis::Muon) return true;
  if (h == Hypothesis::Electron) return false;
  return (pdgLepTrue < 0 ? -pdgLepTrue : pdgLepTrue) == 13;
}

}  // namespace

float RecoObservables::enuRec(Hypothesis h, EnergyEstimator e, int pdgLepTrue) const {
  const bool mu = hypothesisIsMuon(h, pdgLepTrue);
  if (e == EnergyEstimator::CCQE) return mu ? Enu_rec_ccqe_mu : Enu_rec_ccqe_e;
  return static_cast<float>(mu ? Enu_rec_calo_mu : Enu_rec_calo_e);
}

float RecoObservables::q2Rec(Hypothesis h, int pdgLepTrue) const {
  return hypothesisIsMuon(h, pdgLepTrue) ? Q2_rec_ccqe_mu : Q2_rec_ccqe_e;
}

double RecoObservables::elepRec(Hypothesis h, int pdgLepTrue) const {
  return hypothesisIsMuon(h, pdgLepTrue) ? Elep_rec_mu : Elep_rec_e;
}

int RecoObservables::pdgLepRec(Hypothesis h, bool isRHC, int pdgLepTrue) const {
  if (h == Hypothesis::True) return pdgLepTrue;
  const int base = hypothesisIsMuon(h, pdgLepTrue) ? 13 : 11;
  return isRHC ? -base : base;
}

namespace reco {

double leptonMomentum(float ELep, int pdgLep) {
  const double m = leptonMass(pdgLep);
  if (!(ELep > m)) return 0.0;
  // Legacy: ELep * ELep is a float product, then promoted for the subtraction.
  const float e2 = ELep * ELep;
  return std::sqrt(static_cast<double>(e2) - m * m);
}

double leptonEnergy(float PLep_rec, double mLep) {
  const float p2 = PLep_rec * PLep_rec;
  return std::sqrt(static_cast<double>(p2) + mLep * mLep);
}

bool ccqe(float PLep_rec, float CosLep_rec, double mLep, bool isRHC, double Eb, float& Enu_rec, float& Q2_rec) {
  const double E_lep = leptonEnergy(PLep_rec, mLep);
  const double m_tgt = isRHC ? units::m_p : units::m_n;  // struck nucleon
  const double m_out = isRHC ? units::m_n : units::m_p;  // outgoing nucleon
  const double m_eff = m_tgt - Eb;
  const float pcos = PLep_rec * CosLep_rec;  // float product, as the legacy macro
  const double num = m_out * m_out - m_eff * m_eff - mLep * mLep + 2.0 * m_eff * E_lep;
  const double den = 2.0 * (m_eff - E_lep + pcos);
  if (den != 0.0 && num / den > 0.0) {
    Enu_rec = static_cast<float>(num / den);
    // The legacy Q2 uses the float-rounded Enu_rec.
    Q2_rec = static_cast<float>(-mLep * mLep + 2.0 * Enu_rec * (E_lep - pcos));
    return true;
  }
  Enu_rec = -1.0f;
  Q2_rec = -1.0f;
  return false;
}

float thetaDegrees(float cosTheta) {
  const float c = std::min(1.0f, std::max(-1.0f, cosTheta));
  const float a = std::acos(c);  // float acos, as the legacy macro
  return static_cast<float>(TMath::RadToDeg() * a);
}

RecoObservables reconstruct(const FlatEvent& ev, const TruthVariables& truth, const VariableSet& vars,
                            const Smearing& sm, bool isRHC, TRandom& legacyStream, TRandom& extraStream,
                            double bindingEnergy) {
  RecoObservables r;

  // --- lepton smearing (legacy stream: Gaus(p), Gaus(cos)) ------------------
  const double p_true = leptonMomentum(ev.ELep, ev.PDGLep);
  double p_smear = p_true;
  double cos_smear = ev.CosLep;
  if (sm.leptonMomentum().active()) {
    p_smear = p_true + legacyStream.Gaus(0.0, sm.leptonMomentum().sigma(p_true, vars));
    if (p_smear < 0.0) p_smear = 0.0;
  }
  if (sm.leptonAngle().active()) {
    cos_smear = ev.CosLep + legacyStream.Gaus(0.0, sm.leptonAngle().sigma(ev.CosLep, vars));
    if (cos_smear > 1.0) cos_smear = 1.0;
    if (cos_smear < -1.0) cos_smear = -1.0;
  }
  r.PLep_rec = static_cast<float>(p_smear);
  r.CosLep_rec = static_cast<float>(cos_smear);

  // --- decay electrons (legacy stream: one Uniform per parent, first 100 particles) ---
  const DecayElectronTagging& tag = sm.decayElectronTagging();
  const int nUse = std::min(ev.nfsp, 100);
  for (int ip = 0; ip < nUse; ++ip) {
    double eff;
    if (ev.pdg[static_cast<std::size_t>(ip)] == -13) eff = tag.muPlus;
    else if (ev.pdg[static_cast<std::size_t>(ip)] == 13) eff = tag.muMinus;
    else if (ev.pdg[static_cast<std::size_t>(ip)] == 211) eff = tag.piPlus;
    else continue;
    ++r.true_ME;
    if (legacyStream.Uniform() < eff) ++r.recon_ME;
  }

  // --- angle and CCQE estimators -----------------------------------------
  r.Theta_rec = thetaDegrees(r.CosLep_rec);
  ccqe(r.PLep_rec, r.CosLep_rec, units::m_mu, isRHC, bindingEnergy, r.Enu_rec_ccqe_mu, r.Q2_rec_ccqe_mu);
  ccqe(r.PLep_rec, r.CosLep_rec, units::m_e, isRHC, bindingEnergy, r.Enu_rec_ccqe_e, r.Q2_rec_ccqe_e);

  // --- calorimetric quantities (extra stream) -----------------------------
  double ehad = truth.Ehad_true;
  if (sm.hadronicEnergy().active()) {
    ehad += extraStream.Gaus(0.0, sm.hadronicEnergy().sigma(truth.Ehad_true, vars));
    if (ehad < 0.0) ehad = 0.0;
  }
  r.Ehad_rec = ehad;
  if (ev.cc) {
    const double emu = leptonEnergy(r.PLep_rec, units::m_mu);
    const double ee = leptonEnergy(r.PLep_rec, units::m_e);
    r.Elep_rec_mu = emu;
    r.Elep_rec_e = ee;
    if (sm.muonEnergy().active()) {
      r.Elep_rec_mu = std::max(0.0, emu + extraStream.Gaus(0.0, sm.muonEnergy().sigma(emu, vars)));
    }
    if (sm.electronEnergy().active()) {
      r.Elep_rec_e = std::max(0.0, ee + extraStream.Gaus(0.0, sm.electronEnergy().sigma(ee, vars)));
    }
  } else {
    r.Elep_rec_mu = 0.0;  // the outgoing neutrino is invisible
    r.Elep_rec_e = 0.0;
  }
  r.Enu_rec_calo_mu = r.Elep_rec_mu + r.Ehad_rec;
  r.Enu_rec_calo_e = r.Elep_rec_e + r.Ehad_rec;
  return r;
}

}  // namespace reco

}  // namespace zaku
