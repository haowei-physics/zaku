#include <TH1D.h>
#include <TMath.h>
#include <TRandom3.h>

#include <cmath>

#include "ZakuTest.h"
#include "ZakuTestPaths.h"
#include "zaku/core/Units.h"
#include "zaku/detector/Reco.h"
#include "zaku/detector/Smearing.h"
#include "zaku/detector/Variables.h"
#include "zaku/event/FlatEvent.h"

using namespace zaku;

int main() {
  VariableSet vars;

  // --- the parametric smearing reproduces the legacy file bit for bit ---------
  Smearing par = Smearing::parametric(0.024, 1.8);
  const std::string legacyFile = zakutest::sourceDir() + "/data/hk/smearing.root";
  if (zakutest::exists(legacyFile)) {
    Smearing leg = Smearing::fromLegacyFile(legacyFile);
    const TH1* a = par.leptonMomentum().histogramPtr();
    const TH1* b = leg.leptonMomentum().histogramPtr();
    ZAKU_CHECK_EQ(a->GetNbinsX(), b->GetNbinsX());
    for (int i = 1; i <= a->GetNbinsX(); ++i) ZAKU_CHECK(a->GetBinContent(i) == b->GetBinContent(i));
    const TH1* c = par.leptonAngle().histogramPtr();
    const TH1* d = leg.leptonAngle().histogramPtr();
    ZAKU_CHECK_EQ(c->GetNbinsX(), d->GetNbinsX());
    int exact = 0;
    for (int i = 1; i <= c->GetNbinsX(); ++i) {
      if (c->GetBinContent(i) == d->GetBinContent(i)) ++exact;
      ZAKU_CHECK_CLOSE(c->GetBinContent(i), d->GetBinContent(i), 1e-12);
    }
    ZAKU_CHECK_MSG(exact == c->GetNbinsX(), "angle histogram bit-identical bins: " + std::to_string(exact));
  }
  par.prepare(vars);
  ZAKU_CHECK_CLOSE(par.leptonMomentum().sigma(0.62, vars), 0.024 * 0.65, 1e-12);   // bin centre 0.65
  ZAKU_CHECK_CLOSE(par.leptonMomentum().sigma(7.0, vars), 0.024 * 4.95, 1e-12);    // clamped to the last bin
  ZAKU_CHECK_CLOSE(par.leptonAngle().sigma(0.0, vars), std::sqrt(1.0 - 0.025 * 0.025) * 1.8 * units::pi / 180.0, 1e-12);

  // --- resolution models ---------------------------------------------------------
  ResolutionModel frac = ResolutionModel::fractional(0.1);
  ResolutionModel abs = ResolutionModel::absolute(0.05);
  ResolutionModel sto = ResolutionModel::stochastic(0.2, 0.15);
  ResolutionModel form = ResolutionModel::formula("0.1 * value + 0.01 * Enu_true");
  form.prepare(vars, "test");
  ZAKU_CHECK_CLOSE(frac.sigma(2.0, vars), 0.2, 1e-12);
  ZAKU_CHECK_CLOSE(abs.sigma(2.0, vars), 0.05, 1e-12);
  ZAKU_CHECK_CLOSE(sto.sigma(1.0, vars), std::sqrt(0.2 * 0.2 + 0.15 * 0.15), 1e-12);
  ZAKU_CHECK_CLOSE(sto.sigma(4.0, vars), std::sqrt(0.16 + 0.36), 1e-12);
  ZAKU_CHECK_CLOSE(sto.sigma(0.0, vars), 0.0, 1e-12);
  vars.set("Enu_true", 2.0);
  ZAKU_CHECK_CLOSE(form.sigma(1.0, vars), 0.1 + 0.02, 1e-12);
  ZAKU_CHECK(!ResolutionModel::none().active());
  ZAKU_CHECK(frac.active());

  // --- the legacy formulas -----------------------------------------------------
  // lepton momentum with the float square
  const float ELep = 0.8f;
  const float e2 = ELep * ELep;
  ZAKU_CHECK(reco::leptonMomentum(ELep, 13) == std::sqrt(static_cast<double>(e2) - units::m_mu * units::m_mu));
  ZAKU_CHECK_CLOSE(reco::leptonMomentum(0.8f, 14), 0.8, 1e-7);   // massless "lepton" (NC neutrino)
  ZAKU_CHECK_CLOSE(reco::leptonMomentum(0.05f, 13), 0.0, 1e-12);  // below the muon mass

  // CCQE estimator: reproduce the legacy arithmetic independently
  float Enu = 0, Q2 = 0;
  const float p = 0.6f, cth = 0.9f;
  ZAKU_CHECK(reco::ccqe(p, cth, units::m_mu, false, 0.027, Enu, Q2));
  {
    const float p2 = p * p;
    const double El = std::sqrt(static_cast<double>(p2) + units::m_mu * units::m_mu);
    const double meff = units::m_n - 0.027;
    const float pc = p * cth;
    const double num = units::m_p * units::m_p - meff * meff - units::m_mu * units::m_mu + 2.0 * meff * El;
    const double den = 2.0 * (meff - El + pc);
    const float expectedE = static_cast<float>(num / den);
    ZAKU_CHECK(Enu == expectedE);
    const float expectedQ2 = static_cast<float>(-units::m_mu * units::m_mu + 2.0 * expectedE * (El - pc));
    ZAKU_CHECK(Q2 == expectedQ2);
  }
  // RHC swaps the nucleons: different answer
  float EnuR = 0, Q2R = 0;
  reco::ccqe(p, cth, units::m_mu, true, 0.027, EnuR, Q2R);
  ZAKU_CHECK(EnuR != Enu);
  // a fast sideways lepton makes the denominator negative: no solution
  float bad = 0, badQ = 0;
  ZAKU_CHECK(!reco::ccqe(5.0f, 0.0f, units::m_mu, false, 0.027, bad, badQ));
  ZAKU_CHECK(bad == -1.0f && badQ == -1.0f);
  // a fast forward one does have a solution (den > 0)
  float fwd = 0, fwdQ = 0;
  ZAKU_CHECK(reco::ccqe(5.0f, 1.0f, units::m_mu, false, 0.027, fwd, fwdQ));
  ZAKU_CHECK(fwd > 5.0f);

  // theta in float precision
  ZAKU_CHECK(reco::thetaDegrees(0.5f) == static_cast<float>(TMath::RadToDeg() * std::acos(0.5f)));
  ZAKU_CHECK_CLOSE(reco::thetaDegrees(1.5f), 0.0, 1e-12);
  ZAKU_CHECK_CLOSE(reco::thetaDegrees(-1.0f), 180.0, 1e-4);

  // --- reconstruct(): legacy stream order and reproducibility -------------------
  FlatEvent ev;
  ev.cc = true;
  ev.PDGnu = 14;
  ev.Enu_true = 1.0f;
  ev.PDGLep = 13;
  ev.ELep = 0.8f;
  ev.CosLep = 0.95f;
  ev.W = 1.2f;
  const double pmu = std::sqrt(0.8 * 0.8 - units::m_mu * units::m_mu);
  auto add = [&](int pdg, double pz, double e) {
    ev.pdg.push_back(pdg);
    ev.px.push_back(0.f);
    ev.py.push_back(0.f);
    ev.pz.push_back(static_cast<float>(pz));
    ev.E.push_back(static_cast<float>(e));
    ev.nfsp = static_cast<int>(ev.pdg.size());
  };
  add(13, pmu, 0.8);
  add(2212, 0.3, std::hypot(0.3, units::m_p));
  add(211, 0.2, std::hypot(0.2, units::m_pi));
  ProductionContext ctx;
  const TruthVariables truth = TruthVariables::compute(ev, ctx);
  vars.loadTruth(ev, truth);

  TRandom3 legacy(1234), extra(99);
  const RecoObservables r1 = reco::reconstruct(ev, truth, vars, par, false, legacy, extra, 0.027);
  // the same result with the same seeds, computed by hand in the legacy order
  TRandom3 check(1234);
  const double ptrue = reco::leptonMomentum(ev.ELep, ev.PDGLep);
  double ps = ptrue + check.Gaus(0.0, par.leptonMomentum().sigma(ptrue, vars));
  if (ps < 0) ps = 0;
  double cs = ev.CosLep + check.Gaus(0.0, par.leptonAngle().sigma(ev.CosLep, vars));
  cs = std::min(1.0, std::max(-1.0, cs));
  int trueME = 0, reconME = 0;
  for (int pdg : {13, 2212, 211}) {
    double eff = pdg == 13 ? 0.83 : (pdg == 211 ? 0.45 : -1.0);
    if (eff < 0) continue;
    ++trueME;
    if (check.Uniform() < eff) ++reconME;
  }
  ZAKU_CHECK(r1.PLep_rec == static_cast<float>(ps));
  ZAKU_CHECK(r1.CosLep_rec == static_cast<float>(cs));
  ZAKU_CHECK_EQ(r1.true_ME, 2);
  ZAKU_CHECK_EQ(r1.recon_ME, reconME);
  ZAKU_CHECK(r1.Theta_rec == reco::thetaDegrees(r1.CosLep_rec));
  float eMu = 0, q2Mu = 0;
  reco::ccqe(r1.PLep_rec, r1.CosLep_rec, units::m_mu, false, 0.027, eMu, q2Mu);
  ZAKU_CHECK(r1.Enu_rec_ccqe_mu == eMu);
  ZAKU_CHECK(r1.Q2_rec_ccqe_mu == q2Mu);
  // no hadronic smearing configured: Ehad_rec equals the truth, lepton energies from the smeared momentum
  ZAKU_CHECK_CLOSE(r1.Ehad_rec, truth.Ehad_true, 1e-12);
  ZAKU_CHECK_CLOSE(r1.Elep_rec_mu, reco::leptonEnergy(r1.PLep_rec, units::m_mu), 1e-12);
  ZAKU_CHECK_CLOSE(r1.Enu_rec_calo_mu, r1.Elep_rec_mu + r1.Ehad_rec, 1e-12);
  ZAKU_CHECK_CLOSE(r1.enuRec(Hypothesis::Muon, EnergyEstimator::CCQE, 13), eMu, 1e-12);
  ZAKU_CHECK_CLOSE(r1.enuRec(Hypothesis::True, EnergyEstimator::Calorimetric, 13), r1.Enu_rec_calo_mu, 1e-6);
  ZAKU_CHECK_EQ(r1.pdgLepRec(Hypothesis::Electron, true, 13), -11);
  ZAKU_CHECK_EQ(r1.pdgLepRec(Hypothesis::Muon, false, 13), 13);
  ZAKU_CHECK_EQ(r1.pdgLepRec(Hypothesis::True, false, 13), 13);

  // the extra stream does not disturb the legacy quantities
  Smearing withHad = par;
  withHad.setHadronicEnergy(ResolutionModel::stochastic(0.2, 0.15));
  withHad.setMuonEnergy(ResolutionModel::fractional(0.035));
  withHad.prepare(vars);
  TRandom3 legacy2(1234), extra2(99);
  const RecoObservables r2 = reco::reconstruct(ev, truth, vars, withHad, false, legacy2, extra2, 0.027);
  ZAKU_CHECK(r2.PLep_rec == r1.PLep_rec);
  ZAKU_CHECK(r2.CosLep_rec == r1.CosLep_rec);
  ZAKU_CHECK_EQ(r2.recon_ME, r1.recon_ME);
  ZAKU_CHECK(r2.Ehad_rec != r1.Ehad_rec);
  ZAKU_CHECK(r2.Elep_rec_mu != r1.Elep_rec_mu);
  ZAKU_CHECK_CLOSE(r2.Elep_rec_e, r1.Elep_rec_e, 1e-12);  // electron energy not smeared

  // NC event: no lepton term in the calorimetric estimator
  FlatEvent nc = ev;
  nc.cc = false;
  nc.PDGLep = 14;
  nc.pdg[0] = 14;
  const TruthVariables tnc = TruthVariables::compute(nc, ctx);
  vars.loadTruth(nc, tnc);
  TRandom3 l3(1), e3(2);
  const RecoObservables r3 = reco::reconstruct(nc, tnc, vars, par, false, l3, e3, 0.027);
  ZAKU_CHECK_CLOSE(r3.Elep_rec_mu, 0.0, 1e-12);
  ZAKU_CHECK_CLOSE(r3.Enu_rec_calo_e, r3.Ehad_rec, 1e-12);

  // --- YAML round trip of the smearing block ------------------------------------
  const std::string yaml = withHad.toYaml();
  Smearing back = Smearing::fromYamlString(yaml);
  back.prepare(vars);
  ZAKU_CHECK_CLOSE(back.leptonMomentum().sigma(0.62, vars), 0.024 * 0.65, 1e-12);
  ZAKU_CHECK_CLOSE(back.leptonAngle().sigma(0.0, vars), par.leptonAngle().sigma(0.0, vars), 1e-9);
  ZAKU_CHECK_CLOSE(back.hadronicEnergy().sigma(1.0, vars), std::sqrt(0.2 * 0.2 + 0.15 * 0.15), 1e-12);
  ZAKU_CHECK_CLOSE(back.muonEnergy().sigma(2.0, vars), 0.07, 1e-12);
  ZAKU_CHECK_THROWS(Smearing::fromYamlString("bogus: 1"));
  Smearing fromText = Smearing::fromYamlString(
      "lepton_momentum: {type: fractional, sigma: 0.05}\nlepton_angle: {type: gaussian_deg, sigma: 3.0}\n"
      "decay_electron_tag: {mu_plus: 0.5}\nneutron_visible_fraction: 0.2\n");
  ZAKU_CHECK_CLOSE(fromText.leptonMomentum().sigma(2.0, vars), 0.1, 1e-12);
  ZAKU_CHECK_CLOSE(fromText.decayElectronTagging().muPlus, 0.5, 1e-12);
  ZAKU_CHECK_CLOSE(fromText.decayElectronTagging().muMinus, 0.83, 1e-12);
  ZAKU_CHECK_CLOSE(fromText.neutronVisibleFraction(), 0.2, 1e-12);

  ZAKU_CHECK(toString(Hypothesis::Electron) == "electron");
  ZAKU_CHECK(hypothesisFromString("mu") == Hypothesis::Muon);
  ZAKU_CHECK(energyEstimatorFromString("calorimetric") == EnergyEstimator::Calorimetric);
  ZAKU_CHECK_THROWS(energyEstimatorFromString("magic"));

  return zakutest::finish();
}
