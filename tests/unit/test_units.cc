#include <TH1D.h>
#include <TH2D.h>

#include <memory>

#include "ZakuTest.h"
#include "zaku/core/Units.h"

int main() {
  using namespace zaku;

  // scalar conversions
  ZAKU_CHECK_CLOSE(toGeV(1000.0, EnergyUnit::MeV), 1.0, 1e-12);
  ZAKU_CHECK_CLOSE(toGeV(2.5, EnergyUnit::GeV), 2.5, 1e-12);
  ZAKU_CHECK_CLOSE(fromGeV(1.0, EnergyUnit::MeV), 1000.0, 1e-12);
  ZAKU_CHECK(energyUnitFromString("MeV") == EnergyUnit::MeV);
  ZAKU_CHECK(energyUnitFromString("gev") == EnergyUnit::GeV);
  ZAKU_CHECK_THROWS(energyUnitFromString("keV"));
  ZAKU_CHECK_STR_EQ(toString(EnergyUnit::MeV), "MeV");

  // constants the legacy chain relied on
  ZAKU_CHECK_CLOSE(units::m_mu, 0.105658, 1e-12);
  ZAKU_CHECK_CLOSE(units::m_e, 0.000511, 1e-12);
  ZAKU_CHECK_CLOSE(units::m_p, 0.93827, 1e-12);
  ZAKU_CHECK_CLOSE(units::m_n, 0.93957, 1e-12);
  ZAKU_CHECK_CLOSE(units::kt_to_g, 1.0e9, 1e-12);
  ZAKU_CHECK_CLOSE(30.0 * units::deg, 0.5235987755982988, 1e-12);

  // a density histogram in MeV: content c per MeV on 0..30000 MeV
  TH1D hMeV("hMeV", "flux;E_{#nu} [MeV];#nu / cm^{2} / MeV", 300, 0.0, 30000.0);
  for (int i = 1; i <= 300; ++i) hMeV.SetBinContent(i, 2.0);
  const double integralMeV = hMeV.Integral("width");  // 2 * 30000

  ZAKU_CHECK(looksLikeMeVAxis(hMeV));
  ZAKU_CHECK_THROWS(requireGeVAxis(hMeV, "test"));

  std::unique_ptr<TH1> hGeV(histogramToGeV(hMeV, EnergyUnit::MeV, /*contentIsDensity=*/true, 1, "hGeV"));
  ZAKU_CHECK(hGeV != nullptr);
  ZAKU_CHECK_STR_EQ(hGeV->GetName(), "hGeV");
  ZAKU_CHECK_CLOSE(hGeV->GetXaxis()->GetXmax(), 30.0, 1e-12);
  ZAKU_CHECK_CLOSE(hGeV->GetBinContent(1), 2000.0, 1e-12);          // per GeV
  ZAKU_CHECK_CLOSE(hGeV->Integral("width"), integralMeV, 1e-12);    // density preserved
  ZAKU_CHECK(!looksLikeMeVAxis(*hGeV));
  ZAKU_CHECK_NOTHROW(requireGeVAxis(*hGeV, "test"));
  ZAKU_CHECK_STR_EQ(hGeV->GetXaxis()->GetTitle(), "E_{#nu} [GeV]");

  // a counts histogram: content must not be rescaled
  std::unique_ptr<TH1> hCounts(histogramToGeV(hMeV, EnergyUnit::MeV, /*contentIsDensity=*/false));
  ZAKU_CHECK_CLOSE(hCounts->GetBinContent(7), 2.0, 1e-12);
  ZAKU_CHECK_STR_EQ(hCounts->GetName(), "hMeV");

  // GeV input: plain clone
  std::unique_ptr<TH1> hClone(histogramToGeV(*hGeV, EnergyUnit::GeV, true, 1, "clone"));
  ZAKU_CHECK_CLOSE(hClone->GetBinContent(3), hGeV->GetBinContent(3), 1e-12);
  ZAKU_CHECK_CLOSE(hClone->GetXaxis()->GetXmax(), 30.0, 1e-12);

  // 2-D migration matrix in MeV on both axes, density on both -> factor 1e6
  TH2D h2("h2", "migration;E_{true} [MeV];E_{rec} [MeV]", 10, 0.0, 5000.0, 10, 0.0, 5000.0);
  h2.SetBinContent(3, 4, 1.0e-6);
  std::unique_ptr<TH1> h2GeV(histogramToGeV(h2, EnergyUnit::MeV, true, 1 | 2));
  ZAKU_CHECK_CLOSE(h2GeV->GetXaxis()->GetXmax(), 5.0, 1e-12);
  ZAKU_CHECK_CLOSE(h2GeV->GetYaxis()->GetXmax(), 5.0, 1e-12);
  ZAKU_CHECK_CLOSE(h2GeV->GetBinContent(h2GeV->GetBin(3, 4)), 1.0, 1e-9);
  // only x converted
  std::unique_ptr<TH1> h2x(histogramToGeV(h2, EnergyUnit::MeV, false, 1));
  ZAKU_CHECK_CLOSE(h2x->GetXaxis()->GetXmax(), 5.0, 1e-12);
  ZAKU_CHECK_CLOSE(h2x->GetYaxis()->GetXmax(), 5000.0, 1e-12);

  return zakutest::finish();
}
