#include "zaku/osc/Oscillator.h"

#include <TFile.h>
#include <TH1D.h>

#include <array>
#include <memory>
#include <sstream>
#include <stdexcept>

#include "NuFast_LBL.h"

namespace zaku {

OscillationParameters OscillationParameters::asimovA() {
  OscillationParameters p;  // the defaults are the Asimov A point
  return p;
}

double OscillationParameters::dmsq31() const {
  // NuFast wants Dmsq31 = Dmsq32 + Dmsq21, negative for the inverted ordering.
  return normalOrdering ? dmsq32 + dmsq21 : -dmsq32 + dmsq21;
}

std::string OscillationParameters::describe() const {
  std::ostringstream os;
  os << "s12sq " << s12sq << ", s13sq " << s13sq << ", s23sq " << s23sq << ", dCP " << deltaCP << " rad, dm21 "
     << dmsq21 << " eV2, dm32 " << dmsq32 << " eV2, " << (normalOrdering ? "normal" : "inverted") << " ordering";
  return os.str();
}

namespace {
int flavourIndex(int pdg) {
  switch (pdg < 0 ? -pdg : pdg) {
    case 12: return 0;
    case 14: return 1;
    case 16: return 2;
    default: throw std::invalid_argument("Oscillator: " + std::to_string(pdg) + " is not a neutrino PDG code");
  }
}
}  // namespace

Oscillator::Oscillator(const OscillationParameters& p, double baseline_km, double rhoYe, int nNewton)
    : params_(p), baseline_km_(baseline_km), rhoYe_(rhoYe), nNewton_(nNewton) {
  if (baseline_km <= 0.0) throw std::invalid_argument("Oscillator: baseline must be positive");
}

std::vector<double> Oscillator::matrix(double E_GeV, bool antineutrino) const {
  std::vector<double> e(1, antineutrino ? -E_GeV : E_GeV);
  std::vector<std::array<std::array<double, 3>, 3>> probs(1);
  NuFast::Probability_Matter_LBL(params_.s12sq, params_.s13sq, params_.s23sq, params_.deltaCP, params_.dmsq21,
                                 params_.dmsq31(), baseline_km_, e, rhoYe_, nNewton_, probs);
  std::vector<double> out(9);
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) out[static_cast<std::size_t>(a * 3 + b)] = probs[0][static_cast<std::size_t>(a)][static_cast<std::size_t>(b)];
  }
  return out;
}

double Oscillator::probability(int fromPdg, int toPdg, double E_GeV) const {
  const int a = flavourIndex(fromPdg);
  const int b = flavourIndex(toPdg);
  if ((fromPdg < 0) != (toPdg < 0)) return 0.0;  // no neutrino <-> antineutrino transitions
  std::vector<double> e(1, fromPdg < 0 ? -E_GeV : E_GeV);
  std::vector<std::array<std::array<double, 3>, 3>> probs(1);
  NuFast::Probability_Matter_LBL(params_.s12sq, params_.s13sq, params_.s23sq, params_.deltaCP, params_.dmsq21,
                                 params_.dmsq31(), baseline_km_, e, rhoYe_, nNewton_, probs);
  return probs[0][static_cast<std::size_t>(a)][static_cast<std::size_t>(b)];
}

double Oscillator::oscWeight(const FlatEvent& ev, const ProductionContext& ctx) const {
  if (ev.cc) {
    // Legacy: initial flavour from the flux, final from the interaction, the
    // energy sign from the flux species. The flat-tree energy is a float.
    const int a = ((ctx.fluxPdg < 0 ? -ctx.fluxPdg : ctx.fluxPdg) == 12) ? 0 : 1;
    const int b = ((ev.PDGnu < 0 ? -ev.PDGnu : ev.PDGnu) == 12) ? 0 : 1;
    std::vector<double> e(1, ctx.fluxPdg > 0 ? static_cast<double>(ev.Enu_true) : -static_cast<double>(ev.Enu_true));
    std::vector<std::array<std::array<double, 3>, 3>> probs(1);
    NuFast::Probability_Matter_LBL(params_.s12sq, params_.s13sq, params_.s23sq, params_.deltaCP, params_.dmsq21,
                                   params_.dmsq31(), baseline_km_, e, rhoYe_, nNewton_, probs);
    return probs[0][static_cast<std::size_t>(a)][static_cast<std::size_t>(b)];
  }
  // NC: the flavour-blind rate is carried by the disappearance production of this flux.
  return ctx.isAppearance() ? 0.0 : 1.0;
}

TH1D* Oscillator::gridHistogram(int nBins, double eMin, double eMax, const std::string& name) {
  if (nBins <= 0 || eMax <= eMin) throw std::invalid_argument("Oscillator::gridHistogram: invalid grid");
  auto* h = new TH1D(name.c_str(), "oscillation grid;E_{#nu} [GeV];", nBins, eMin, eMax);
  h->SetDirectory(nullptr);
  return h;
}

void Oscillator::writeGrid(const std::string& file, int nBins, double eMin, double eMax) {
  std::unique_ptr<TH1D> h(gridHistogram(nBins, eMin, eMax));
  std::unique_ptr<TFile> f(TFile::Open(file.c_str(), "RECREATE"));
  if (!f || f->IsZombie()) throw std::runtime_error("Oscillator::writeGrid: cannot create '" + file + "'");
  f->WriteTObject(h.get(), h->GetName(), "Overwrite");
  f->Close();
}

}  // namespace zaku
