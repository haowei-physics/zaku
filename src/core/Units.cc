#include "zaku/core/Units.h"

#include <TAxis.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace zaku {

namespace {

std::string lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return s;
}

std::vector<double> edges(const TAxis& ax, double scale) {
  const int n = ax.GetNbins();
  std::vector<double> e(static_cast<std::size_t>(n) + 1);
  for (int i = 1; i <= n; ++i) e[static_cast<std::size_t>(i - 1)] = ax.GetBinLowEdge(i) * scale;
  e[static_cast<std::size_t>(n)] = ax.GetBinUpEdge(n) * scale;
  return e;
}

std::string relabel(const char* title, EnergyUnit from) {
  std::string t = title ? title : "";
  if (from == EnergyUnit::MeV) {
    const std::string::size_type pos = t.find("MeV");
    if (pos != std::string::npos) t.replace(pos, 3, "GeV");
  }
  return t;
}

}  // namespace

std::string toString(EnergyUnit u) { return u == EnergyUnit::GeV ? "GeV" : "MeV"; }

EnergyUnit energyUnitFromString(const std::string& s) {
  const std::string l = lower(s);
  if (l == "gev") return EnergyUnit::GeV;
  if (l == "mev") return EnergyUnit::MeV;
  throw std::invalid_argument("energyUnitFromString: unrecognised energy unit '" + s + "' (expected GeV or MeV)");
}

double toGeV(double value, EnergyUnit from) { return from == EnergyUnit::MeV ? value * units::MeV : value; }

double fromGeV(double valueGeV, EnergyUnit to) { return to == EnergyUnit::MeV ? valueGeV / units::MeV : valueGeV; }

TH1* histogramToGeV(const TH1& h, EnergyUnit from, bool contentIsDensity, int axes, const std::string& newName) {
  const std::string name = newName.empty() ? std::string(h.GetName()) : newName;
  const double scale = toGeV(1.0, from);
  const int dim = h.GetDimension();

  const bool cx = (axes & 1) != 0;
  const bool cy = (axes & 2) != 0 && dim >= 2;
  const bool cz = (axes & 4) != 0 && dim >= 3;
  const int nConverted = (cx ? 1 : 0) + (cy ? 1 : 0) + (cz ? 1 : 0);
  const double contentFactor =
      (contentIsDensity && from != EnergyUnit::GeV) ? std::pow(1.0 / scale, nConverted) : 1.0;

  const std::vector<double> ex = edges(*h.GetXaxis(), cx ? scale : 1.0);
  const std::vector<double> ey = dim >= 2 ? edges(*h.GetYaxis(), cy ? scale : 1.0) : std::vector<double>{};
  const std::vector<double> ez = dim >= 3 ? edges(*h.GetZaxis(), cz ? scale : 1.0) : std::vector<double>{};

  TH1* out = nullptr;
  if (dim == 1) {
    out = new TH1D(name.c_str(), h.GetTitle(), h.GetNbinsX(), ex.data());
  } else if (dim == 2) {
    out = new TH2D(name.c_str(), h.GetTitle(), h.GetNbinsX(), ex.data(), h.GetNbinsY(), ey.data());
  } else {
    out = new TH3D(name.c_str(), h.GetTitle(), h.GetNbinsX(), ex.data(), h.GetNbinsY(), ey.data(),
                   h.GetNbinsZ(), ez.data());
  }
  out->SetDirectory(nullptr);
  out->Sumw2();

  const int nx = h.GetNbinsX(), ny = dim >= 2 ? h.GetNbinsY() : 1, nz = dim >= 3 ? h.GetNbinsZ() : 1;
  for (int ix = 0; ix <= nx + 1; ++ix) {
    for (int iy = 0; iy <= ny + 1; ++iy) {
      for (int iz = 0; iz <= nz + 1; ++iz) {
        const int bin = h.GetBin(ix, iy, iz);
        out->SetBinContent(bin, h.GetBinContent(bin) * contentFactor);
        out->SetBinError(bin, h.GetBinError(bin) * contentFactor);
      }
    }
  }
  out->SetEntries(h.GetEntries());

  out->GetXaxis()->SetTitle(relabel(h.GetXaxis()->GetTitle(), cx ? from : EnergyUnit::GeV).c_str());
  if (dim >= 2) out->GetYaxis()->SetTitle(relabel(h.GetYaxis()->GetTitle(), cy ? from : EnergyUnit::GeV).c_str());
  if (dim >= 3) out->GetZaxis()->SetTitle(relabel(h.GetZaxis()->GetTitle(), cz ? from : EnergyUnit::GeV).c_str());
  if (dim == 1 && contentIsDensity) out->GetYaxis()->SetTitle(relabel(h.GetYaxis()->GetTitle(), from).c_str());
  return out;
}

bool looksLikeMeVAxis(const TH1& h, double maxGeV) { return h.GetXaxis()->GetXmax() > maxGeV; }

void requireGeVAxis(const TH1& h, const std::string& context, double maxGeV) {
  if (!looksLikeMeVAxis(h, maxGeV)) return;
  throw std::runtime_error(context + ": histogram '" + h.GetName() + "' has an x axis reaching " +
                           std::to_string(h.GetXaxis()->GetXmax()) +
                           ", which cannot be GeV. Declare 'units: MeV' for this input so Zaku converts it.");
}

}  // namespace zaku
