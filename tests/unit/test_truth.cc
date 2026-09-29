#include <cmath>

#include "ZakuTest.h"
#include "zaku/core/Units.h"
#include "zaku/event/FlatEvent.h"

using namespace zaku;

namespace {

void addParticle(FlatEvent& ev, int pdg, double px, double py, double pz, double e) {
  ev.pdg.push_back(pdg);
  ev.px.push_back(static_cast<float>(px));
  ev.py.push_back(static_cast<float>(py));
  ev.pz.push_back(static_cast<float>(pz));
  ev.E.push_back(static_cast<float>(e));
  ev.nfsp = static_cast<int>(ev.pdg.size());
}

}  // namespace

int main() {
  // A CC numu event: mu- (primary), a proton, a pi+, a neutron, a photon and a stray neutrino.
  FlatEvent ev;
  ev.entry = 3;
  ev.Mode = 11;
  ev.cc = true;
  ev.PDGnu = 14;
  ev.Enu_true = 1.2f;
  ev.PDGLep = 13;
  ev.ELep = 0.8f;
  ev.CosLep = 0.9f;
  ev.W = 1.3f;
  const double pmu = std::sqrt(0.8 * 0.8 - units::m_mu * units::m_mu);
  addParticle(ev, 13, 0.0, 0.0, pmu, 0.8);                       // primary muon
  addParticle(ev, 2212, 0.0, 0.3, 0.0, std::hypot(0.3, units::m_p));  // proton KE
  addParticle(ev, 211, 0.1, 0.0, 0.0, std::hypot(0.1, units::m_pi));  // pi+ total E
  addParticle(ev, 2112, 0.0, 0.0, 0.2, std::hypot(0.2, units::m_n));  // neutron
  addParticle(ev, 22, 0.05, 0.0, 0.0, 0.05);                      // photon
  addParticle(ev, 14, 0.0, 0.0, 0.5, 0.5);                        // stray neutrino: invisible
  addParticle(ev, 13, 0.0, 0.0, 0.05, std::hypot(0.05, units::m_mu));  // a second, softer muon: hadronic KE

  ZAKU_CHECK_CLOSE(ev.momentum(1), 0.3, 1e-6);
  ZAKU_CHECK_CLOSE(ev.mass(1), units::m_p, 1e-5);
  ZAKU_CHECK_CLOSE(ev.kineticEnergy(1), std::hypot(0.3, units::m_p) - units::m_p, 1e-5);
  ZAKU_CHECK(!ev.emptyFinalState());
  ZAKU_CHECK(!ev.legacyDrop());

  ProductionContext ctx;
  ctx.isRHC = false;
  ctx.fluxPdg = 14;
  ctx.interactionPdg = 14;
  const TruthVariables t = TruthVariables::compute(ev, ctx, 0.0);
  ZAKU_CHECK_CLOSE(t.PLep_true, pmu, 1e-6);
  ZAKU_CHECK_EQ(t.nFS_p, 1);
  ZAKU_CHECK_EQ(t.nFS_pip, 1);
  ZAKU_CHECK_EQ(t.nFS_n, 1);
  ZAKU_CHECK_EQ(t.nFS_gamma, 1);
  ZAKU_CHECK_EQ(t.nFS_mum, 1);  // the soft muon; the primary is excluded
  ZAKU_CHECK_EQ(t.nFS_pi, 1);
  ZAKU_CHECK(!t.isAppearance);
  const double expectedEhad = (std::hypot(0.3, units::m_p) - units::m_p) + std::hypot(0.1, units::m_pi) + 0.05 +
                              (std::hypot(0.05, units::m_mu) - units::m_mu);
  ZAKU_CHECK_CLOSE(t.Ehad_true, expectedEhad, 1e-5);
  ZAKU_CHECK_CLOSE(t.Evis_true, expectedEhad + 0.8, 1e-5);
  ZAKU_CHECK_CLOSE(t.W2, 1.3 * 1.3, 1e-6);
  // neutron fraction
  const TruthVariables tn = TruthVariables::compute(ev, ctx, 0.5);
  ZAKU_CHECK_CLOSE(tn.Ehad_true, expectedEhad + 0.5 * (std::hypot(0.2, units::m_n) - units::m_n), 1e-5);

  // NC event on an appearance sample: the outgoing neutrino is the "lepton" and invisible.
  FlatEvent nc;
  nc.Mode = 31;
  nc.cc = false;
  nc.PDGnu = 12;
  nc.PDGLep = 12;
  nc.ELep = 0.7f;
  addParticle(nc, 12, 0.0, 0.0, 0.7, 0.7);
  addParticle(nc, 111, 0.0, 0.0, 0.3, std::hypot(0.3, units::m_pi0));
  ProductionContext actx;
  actx.fluxPdg = 14;
  actx.interactionPdg = 12;
  const TruthVariables tnc = TruthVariables::compute(nc, actx);
  ZAKU_CHECK(tnc.isAppearance);
  ZAKU_CHECK_CLOSE(tnc.PLep_true, 0.7, 1e-6);
  ZAKU_CHECK_EQ(tnc.nFS_pi0, 1);
  ZAKU_CHECK_CLOSE(tnc.Ehad_true, std::hypot(0.3, units::m_pi0), 1e-6);
  ZAKU_CHECK_CLOSE(tnc.Evis_true, tnc.Ehad_true, 1e-12);  // no lepton term for NC

  // An empty final state (the events T2KReWeight skips).
  FlatEvent empty;
  empty.cc = true;
  empty.PDGLep = 0;
  empty.ELep = -999.9f;
  empty.Q2 = -999.9f;
  ZAKU_CHECK(empty.emptyFinalState());
  ZAKU_CHECK(empty.legacyDrop());

  // names/values are aligned
  const auto& names = TruthVariables::names();
  const auto values = t.values(ev);
  ZAKU_CHECK_EQ(static_cast<int>(names.size()), static_cast<int>(values.size()));
  int idx = -1;
  for (std::size_t i = 0; i < names.size(); ++i) {
    if (names[i] == "Ehad_true") idx = static_cast<int>(i);
  }
  ZAKU_CHECK(idx >= 0);
  ZAKU_CHECK_CLOSE(values[static_cast<std::size_t>(idx)], t.Ehad_true, 1e-12);
  ZAKU_CHECK(ev.describe().find("mode 11") != std::string::npos);

  return zakutest::finish();
}
