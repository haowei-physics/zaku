#include "zaku/event/FlatEvent.h"

#include <cmath>
#include <sstream>
#include <stdexcept>

#include "zaku/core/Units.h"

namespace zaku {

double FlatEvent::momentum(int i) const {
  const double p2 = static_cast<double>(px[i]) * px[i] + static_cast<double>(py[i]) * py[i] +
                    static_cast<double>(pz[i]) * pz[i];
  return std::sqrt(p2 > 0.0 ? p2 : 0.0);
}

double FlatEvent::mass(int i) const {
  const double p = momentum(i);
  const double m2 = static_cast<double>(E[i]) * E[i] - p * p;
  return std::sqrt(m2 > 0.0 ? m2 : 0.0);
}

double FlatEvent::kineticEnergy(int i) const {
  const double ke = static_cast<double>(E[i]) - mass(i);
  return ke > 0.0 ? ke : 0.0;
}

std::string FlatEvent::describe() const {
  std::ostringstream os;
  os << "entry " << entry << ": mode " << Mode << (cc ? " CC " : " NC ") << "nu " << PDGnu << " E=" << Enu_true
     << " GeV, lepton " << PDGLep << " E=" << ELep << " cos=" << CosLep << ", Q2=" << Q2 << ", " << nfsp
     << " final-state particles";
  return os.str();
}

namespace {

double leptonMass(int pdg) {
  switch (pdg < 0 ? -pdg : pdg) {
    case 11: return units::m_e;
    case 13: return units::m_mu;
    case 15: return 1.77686;
    default: return 0.0;
  }
}

bool isNeutrino(int pdg) {
  const int a = pdg < 0 ? -pdg : pdg;
  return a == 12 || a == 14 || a == 16;
}

}  // namespace

TruthVariables TruthVariables::compute(const FlatEvent& ev, const SampleContext& ctx, double neutronFraction) {
  TruthVariables t;
  t.isRHC = ctx.isRHC;
  t.flux_pdg = ctx.fluxPdg;
  t.isAppearance = ctx.isAppearance();

  // Outgoing lepton momentum, as the legacy macro computed it.
  const double ml = leptonMass(ev.PDGLep);
  t.PLep_true = (ev.PDGLep != 0 && ev.ELep > ml) ? std::sqrt(static_cast<double>(ev.ELep) * ev.ELep - ml * ml) : 0.0;

  // The primary lepton (or outgoing NC neutrino) is in the final-state list:
  // exclude the most energetic particle carrying PDGLep exactly once.
  int primary = -1;
  if (ev.PDGLep != 0) {
    for (int i = 0; i < ev.nfsp; ++i) {
      if (ev.pdg[i] == ev.PDGLep && (primary < 0 || ev.E[i] > ev.E[primary])) primary = i;
    }
  }

  double ehad = 0.0;
  for (int i = 0; i < ev.nfsp; ++i) {
    if (i == primary) continue;
    const int pdg = ev.pdg[i];
    if (isNeutrino(pdg)) continue;  // invisible
    const double e = ev.E[i];
    switch (pdg) {
      case 211: ++t.nFS_pip; ehad += e; break;
      case -211: ++t.nFS_pim; ehad += e; break;
      case 111: ++t.nFS_pi0; ehad += e; break;
      case 22: ++t.nFS_gamma; ehad += e; break;
      case 2212: ++t.nFS_p; ehad += ev.kineticEnergy(i); break;
      case 2112: ++t.nFS_n; ehad += neutronFraction * ev.kineticEnergy(i); break;
      case 321: ++t.nFS_kp; ehad += e; break;
      case -321: ++t.nFS_km; ehad += e; break;
      case 311: case -311: case 310: case 130: ++t.nFS_k0; ehad += e; break;
      case -13: ++t.nFS_mup; ehad += ev.kineticEnergy(i); break;
      case 13: ++t.nFS_mum; ehad += ev.kineticEnergy(i); break;
      case -11: ++t.nFS_ep; ehad += e; break;
      case 11: ++t.nFS_em; ehad += e; break;
      default: {
        ++t.nFS_other;
        const int a = pdg < 0 ? -pdg : pdg;
        if (a > 1000000000) break;                       // nuclear fragment: not visible
        if (a > 1000) ehad += ev.kineticEnergy(i);       // baryons deposit their kinetic energy
        else ehad += e;                                   // mesons (eta, ...) their full energy
        break;
      }
    }
  }
  t.nFS_pi = t.nFS_pip + t.nFS_pim;
  t.Ehad_true = ehad;
  t.Evis_true = ehad + (ev.cc ? static_cast<double>(ev.ELep) : 0.0);
  t.W2 = ev.W > 0.f ? static_cast<double>(ev.W) * ev.W : static_cast<double>(ev.W);
  return t;
}

const std::vector<std::string>& TruthVariables::names() {
  static const std::vector<std::string> n = {
      // flat-tree scalars
      "Mode", "cc", "PDGnu", "Enu_true", "tgt", "tgta", "tgtz", "PDGLep", "ELep", "CosLep", "Q2", "q0", "q3",
      "Enu_QE", "Q2_QE", "W", "x", "y", "fScaleFactor", "nfsp",
      "flagCCINC", "flagNCINC", "flagCCQE", "flagCC0pi", "flagCCQELike", "flagNCEL", "flagNC0pi", "flagCCcoh",
      "flagNCcoh", "flagCC1pip", "flagNC1pip", "flagCC1pim", "flagNC1pim", "flagCC1pi0", "flagNC1pi0",
      // derived
      "PLep_true", "Ehad_true", "Evis_true", "W2", "nFS_pip", "nFS_pim", "nFS_pi0", "nFS_p", "nFS_n", "nFS_gamma",
      "nFS_mup", "nFS_mum", "nFS_ep", "nFS_em", "nFS_kp", "nFS_km", "nFS_k0", "nFS_other", "nFS_pi",
      "isRHC", "flux_pdg", "isAppearance"};
  return n;
}

std::vector<double> TruthVariables::values(const FlatEvent& ev) const {
  std::vector<double> v = {
      static_cast<double>(ev.Mode), ev.cc ? 1.0 : 0.0, static_cast<double>(ev.PDGnu), ev.Enu_true,
      static_cast<double>(ev.tgt), static_cast<double>(ev.tgta), static_cast<double>(ev.tgtz),
      static_cast<double>(ev.PDGLep), ev.ELep, ev.CosLep, ev.Q2, ev.q0, ev.q3, ev.Enu_QE, ev.Q2_QE, ev.W, ev.x, ev.y,
      ev.fScaleFactor, static_cast<double>(ev.nfsp),
      ev.flagCCINC ? 1.0 : 0.0, ev.flagNCINC ? 1.0 : 0.0, ev.flagCCQE ? 1.0 : 0.0, ev.flagCC0pi ? 1.0 : 0.0,
      ev.flagCCQELike ? 1.0 : 0.0, ev.flagNCEL ? 1.0 : 0.0, ev.flagNC0pi ? 1.0 : 0.0, ev.flagCCcoh ? 1.0 : 0.0,
      ev.flagNCcoh ? 1.0 : 0.0, ev.flagCC1pip ? 1.0 : 0.0, ev.flagNC1pip ? 1.0 : 0.0, ev.flagCC1pim ? 1.0 : 0.0,
      ev.flagNC1pim ? 1.0 : 0.0, ev.flagCC1pi0 ? 1.0 : 0.0, ev.flagNC1pi0 ? 1.0 : 0.0,
      PLep_true, Ehad_true, Evis_true, W2,
      static_cast<double>(nFS_pip), static_cast<double>(nFS_pim), static_cast<double>(nFS_pi0),
      static_cast<double>(nFS_p), static_cast<double>(nFS_n), static_cast<double>(nFS_gamma),
      static_cast<double>(nFS_mup), static_cast<double>(nFS_mum), static_cast<double>(nFS_ep),
      static_cast<double>(nFS_em), static_cast<double>(nFS_kp), static_cast<double>(nFS_km),
      static_cast<double>(nFS_k0), static_cast<double>(nFS_other), static_cast<double>(nFS_pi),
      isRHC ? 1.0 : 0.0, static_cast<double>(flux_pdg), isAppearance ? 1.0 : 0.0};
  if (v.size() != names().size()) throw std::logic_error("TruthVariables: names/values size mismatch");
  return v;
}

}  // namespace zaku
