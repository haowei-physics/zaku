// ShiftInjector: a shifted universe is one consistent reconstruction.
#include <cmath>

#include "ZakuTest.h"
#include "zaku/core/Units.h"
#include "zaku/detector/Reco.h"
#include "zaku/detector/Shifts.h"

using namespace zaku;

namespace {

RecoObservables nominal() {
  RecoObservables ro;
  ro.PLep_rec = 0.6f;
  ro.CosLep_rec = 0.9f;
  ro.Ehad_rec = 0.3;
  ShiftInjector::rederive(ro, false, 0.027, true, true);
  return ro;
}

}  // namespace

int main() {
  const RecoObservables ro = nominal();
  ZAKU_CHECK(ro.Enu_rec_ccqe_mu > 0.f);
  ZAKU_CHECK_CLOSE(ro.Elep_rec_mu, std::sqrt(0.6 * 0.6 + units::m_mu * units::m_mu), 1e-6);
  ZAKU_CHECK_CLOSE(ro.Enu_rec_calo_mu, ro.Elep_rec_mu + 0.3, 1e-12);
  ZAKU_CHECK_CLOSE(ro.Theta_rec, reco::thetaDegrees(0.9f), 1e-9);

  // momentum scale: p moves, both estimators follow, and -1 sigma undoes +1 sigma to first order
  ShiftParameter pscale;
  pscale.name = "p";
  pscale.target = ShiftTarget::LeptonMomentum;
  pscale.kind = ShiftKind::Scale;
  pscale.sigma = 0.02;
  RecoObservables up = ro;
  ShiftInjector::apply(up, pscale, 1.0, false, 0.027, true, Hypothesis::Muon, 13);
  ZAKU_CHECK_CLOSE(up.PLep_rec, 0.6 * 1.02, 1e-6);
  ZAKU_CHECK(up.Enu_rec_ccqe_mu > ro.Enu_rec_ccqe_mu);
  float e = 0, q = 0;
  reco::ccqe(up.PLep_rec, up.CosLep_rec, units::m_mu, false, 0.027, e, q);
  ZAKU_CHECK_CLOSE(up.Enu_rec_ccqe_mu, e, 1e-9);
  ZAKU_CHECK_CLOSE(up.Elep_rec_mu, reco::leptonEnergy(up.PLep_rec, units::m_mu), 1e-6);
  ZAKU_CHECK_CLOSE(up.Enu_rec_calo_mu, up.Elep_rec_mu + 0.3, 1e-12);
  ZAKU_CHECK_CLOSE(up.Theta_rec, ro.Theta_rec, 1e-12);  // the angle is untouched
  // smearing of the lepton term is preserved as a ratio
  RecoObservables smeared = ro;
  smeared.Elep_rec_mu *= 1.1;
  ShiftInjector::apply(smeared, pscale, 1.0, false, 0.027, true, Hypothesis::Muon, 13);
  ZAKU_CHECK_CLOSE(smeared.Elep_rec_mu, 1.1 * reco::leptonEnergy(smeared.PLep_rec, units::m_mu), 1e-6);
  // zero sigmas or a non-applying hypothesis leave everything alone
  RecoObservables same = ro;
  ShiftInjector::apply(same, pscale, 0.0, false, 0.027, true, Hypothesis::Muon, 13);
  ZAKU_CHECK_CLOSE(same.PLep_rec, ro.PLep_rec, 1e-12);
  pscale.applyTo = ShiftApplyTo::Electron;
  ShiftInjector::apply(same, pscale, 1.0, false, 0.027, true, Hypothesis::Muon, 13);
  ZAKU_CHECK_CLOSE(same.PLep_rec, ro.PLep_rec, 1e-12);
  ShiftInjector::apply(same, pscale, 1.0, false, 0.027, true, Hypothesis::True, 11);  // true electron: applies
  ZAKU_CHECK_CLOSE(same.PLep_rec, 0.6 * 1.02, 1e-6);
  ZAKU_CHECK(pscale.appliesTo(Hypothesis::Electron, 13));
  ZAKU_CHECK(!pscale.appliesTo(Hypothesis::True, 13));

  // angle offset in degrees
  ShiftParameter angle;
  angle.name = "theta";
  angle.target = ShiftTarget::LeptonAngle;
  angle.kind = ShiftKind::Offset;
  angle.sigma = 1.0;  // degrees
  RecoObservables tilted = ro;
  ShiftInjector::apply(tilted, angle, 2.0, false, 0.027, true, Hypothesis::Muon, 13);
  ZAKU_CHECK_CLOSE(tilted.Theta_rec, ro.Theta_rec + 2.0, 1e-4);
  ZAKU_CHECK(tilted.CosLep_rec < ro.CosLep_rec);
  ZAKU_CHECK_CLOSE(tilted.PLep_rec, ro.PLep_rec, 1e-12);
  ZAKU_CHECK(tilted.Enu_rec_ccqe_mu != ro.Enu_rec_ccqe_mu);
  // clamped at pi
  RecoObservables back = ro;
  back.CosLep_rec = -1.0f;
  ShiftInjector::rederive(back, false, 0.027, true, true);
  ShiftInjector::apply(back, angle, 5.0, false, 0.027, true, Hypothesis::Muon, 13);
  ZAKU_CHECK_CLOSE(back.CosLep_rec, -1.0, 1e-9);

  // hadronic energy scale only touches the calorimetric sums
  ShiftParameter had;
  had.name = "had";
  had.target = ShiftTarget::HadronicEnergy;
  had.kind = ShiftKind::Scale;
  had.sigma = 0.05;
  RecoObservables h = ro;
  ShiftInjector::apply(h, had, -1.0, false, 0.027, true, Hypothesis::Muon, 13);
  ZAKU_CHECK_CLOSE(h.Ehad_rec, 0.3 * 0.95, 1e-12);
  ZAKU_CHECK_CLOSE(h.Enu_rec_calo_mu, ro.Elep_rec_mu + 0.3 * 0.95, 1e-12);
  ZAKU_CHECK_CLOSE(h.Enu_rec_ccqe_mu, ro.Enu_rec_ccqe_mu, 1e-12);

  // lepton energy offset (GeV) and reco energy scale
  ShiftParameter lep;
  lep.name = "lep";
  lep.target = ShiftTarget::LeptonEnergy;
  lep.kind = ShiftKind::Offset;
  lep.sigma = 0.01;
  RecoObservables l = ro;
  ShiftInjector::apply(l, lep, 1.0, false, 0.027, true, Hypothesis::Electron, 13);
  ZAKU_CHECK_CLOSE(l.Elep_rec_e, ro.Elep_rec_e + 0.01, 1e-12);
  ZAKU_CHECK_CLOSE(l.Enu_rec_calo_e, ro.Enu_rec_calo_e + 0.01, 1e-12);
  ShiftParameter enu;
  enu.name = "enu";
  enu.target = ShiftTarget::RecoEnergy;
  enu.kind = ShiftKind::Scale;
  enu.sigma = 0.03;
  RecoObservables n = ro;
  ShiftInjector::apply(n, enu, 1.0, false, 0.027, true, Hypothesis::Muon, 13);
  ZAKU_CHECK_CLOSE(n.Enu_rec_ccqe_mu, ro.Enu_rec_ccqe_mu * 1.03, 1e-6);
  ZAKU_CHECK_CLOSE(n.Enu_rec_calo_mu, ro.Enu_rec_calo_mu * 1.03, 1e-12);
  RecoObservables unphysical = ro;
  unphysical.Enu_rec_ccqe_mu = -1.f;
  ShiftInjector::apply(unphysical, enu, 1.0, false, 0.027, true, Hypothesis::Muon, 13);
  ZAKU_CHECK_CLOSE(unphysical.Enu_rec_ccqe_mu, -1.0, 1e-12);  // the sentinel stays

  // NC: the calorimetric lepton term stays 0
  RecoObservables nc = ro;
  ShiftInjector::rederive(nc, false, 0.027, false, true);
  ZAKU_CHECK_CLOSE(nc.Elep_rec_mu, 0.0, 1e-12);
  ZAKU_CHECK_CLOSE(nc.Enu_rec_calo_mu, 0.3, 1e-12);

  // names round trip
  ZAKU_CHECK(shiftTargetFromString(toString(ShiftTarget::HadronicEnergy)) == ShiftTarget::HadronicEnergy);
  ZAKU_CHECK(shiftKindFromString("offset") == ShiftKind::Offset);
  ZAKU_CHECK(shiftApplyToFromString("mu") == ShiftApplyTo::Muon);
  ZAKU_CHECK_THROWS(shiftTargetFromString("nope"));
  ZAKU_CHECK(!pscale.describe().empty());
  return zakutest::finish();
}
