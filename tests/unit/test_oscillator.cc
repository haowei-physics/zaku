#include <TFile.h>
#include <TH1D.h>

#include <cmath>
#include <memory>

#include "ZakuTest.h"
#include "ZakuTestPaths.h"
#include "zaku/osc/Oscillator.h"

using namespace zaku;

int main() {
  const OscillationParameters p = OscillationParameters::asimovA();
  ZAKU_CHECK_CLOSE(p.s23sq, 0.528, 1e-12);
  ZAKU_CHECK_CLOSE(p.dmsq31(), 2.509e-3 + 7.53e-5, 1e-15);  // NuFast convention, as the legacy macro
  OscillationParameters io = p;
  io.normalOrdering = false;
  ZAKU_CHECK(io.dmsq31() < 0.0);
  ZAKU_CHECK(p.describe().find("normal") != std::string::npos);

  Oscillator osc(p, 295.0, 2.6 * 0.5);
  ZAKU_CHECK_THROWS(Oscillator bad(p, -1.0, 1.3));

  // unitarity and sanity at the T2K peak
  for (double e : {0.3, 0.6, 1.0, 2.5}) {
    const std::vector<double> m = osc.matrix(e, false);
    for (int a = 0; a < 3; ++a) {
      double sum = 0.0;
      for (int b = 0; b < 3; ++b) {
        sum += m[static_cast<std::size_t>(a * 3 + b)];
        ZAKU_CHECK(m[static_cast<std::size_t>(a * 3 + b)] >= -1e-12 && m[static_cast<std::size_t>(a * 3 + b)] <= 1.0 + 1e-12);
      }
      ZAKU_CHECK_CLOSE(sum, 1.0, 1e-9);
    }
  }
  const double pmm = osc.probability(14, 14, 0.6);
  const double pme = osc.probability(14, 12, 0.6);
  ZAKU_CHECK(pmm < 0.2);               // strong disappearance near the oscillation maximum
  ZAKU_CHECK(pme > 0.03 && pme < 0.08);  // appearance of a few per cent
  ZAKU_CHECK_CLOSE(osc.probability(14, 14, 0.6), osc.matrix(0.6, false)[4], 1e-15);
  ZAKU_CHECK_CLOSE(osc.probability(-14, -12, 0.6), osc.matrix(0.6, true)[3], 1e-15);
  ZAKU_CHECK(osc.probability(-14, -12, 0.6) != pme);  // CP and matter asymmetry
  ZAKU_CHECK_CLOSE(osc.probability(14, -12, 0.6), 0.0, 1e-15);
  ZAKU_CHECK_THROWS(osc.probability(13, 12, 0.6));

  // legacy OscWeight rules
  FlatEvent ev;
  ev.cc = true;
  ev.PDGnu = 12;
  ev.Enu_true = 0.6f;
  ProductionContext app;
  app.fluxPdg = 14;
  app.interactionPdg = 12;
  ZAKU_CHECK_CLOSE(osc.oscWeight(ev, app), osc.probability(14, 12, static_cast<double>(0.6f)), 1e-15);
  ProductionContext dis;
  dis.fluxPdg = 14;
  dis.interactionPdg = 14;
  ev.PDGnu = 14;
  ZAKU_CHECK_CLOSE(osc.oscWeight(ev, dis), osc.probability(14, 14, static_cast<double>(0.6f)), 1e-15);
  ev.cc = false;
  ZAKU_CHECK_CLOSE(osc.oscWeight(ev, dis), 1.0, 1e-15);
  ZAKU_CHECK_CLOSE(osc.oscWeight(ev, app), 0.0, 1e-15);
  ProductionContext anti;
  anti.fluxPdg = -14;
  anti.interactionPdg = -14;
  ev.cc = true;
  ev.PDGnu = -14;
  ZAKU_CHECK_CLOSE(osc.oscWeight(ev, anti), osc.probability(-14, -14, static_cast<double>(0.6f)), 1e-15);

  // the grid for NuOscillator
  std::unique_ptr<TH1D> grid(Oscillator::gridHistogram());
  ZAKU_CHECK_EQ(grid->GetNbinsX(), 2000);
  ZAKU_CHECK_CLOSE(grid->GetXaxis()->GetXmax(), 10.0, 1e-12);
  ZAKU_CHECK_STR_EQ(grid->GetName(), "energyBinning");
  const std::string gridFile = zakutest::scratchDir() + "/grid.root";
  Oscillator::writeGrid(gridFile);
  std::unique_ptr<TFile> f(TFile::Open(gridFile.c_str(), "READ"));
  ZAKU_CHECK(f && f->Get<TH1D>("energyBinning") != nullptr);
  ZAKU_CHECK_THROWS(Oscillator::gridHistogram(0));

  return zakutest::finish();
}
