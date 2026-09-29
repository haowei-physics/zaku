#include "zaku/flux/FluxHistograms.h"

#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>

#include <memory>
#include <sstream>
#include <stdexcept>

#include "zaku/core/Units.h"

namespace zaku {

TH1D* FluxHistograms::load(const std::string& file, const std::string& name, EnergyUnit units,
                           bool contentIsDensity) {
  std::unique_ptr<TFile> f(TFile::Open(file.c_str(), "READ"));
  if (!f || f->IsZombie()) throw std::runtime_error("FluxHistograms: cannot open flux file '" + file + "'");
  TH1* h = f->Get<TH1>(name.c_str());
  if (!h) throw std::runtime_error("FluxHistograms: no histogram '" + name + "' in '" + file + "'");
  if (h->GetDimension() != 1) {
    throw std::runtime_error("FluxHistograms: histogram '" + name + "' in '" + file + "' is not one-dimensional");
  }
  // histogramToGeV clones (as TH1D) and converts MeV axes; GeV inputs are copied as they are.
  TH1* out = histogramToGeV(*h, units, contentIsDensity, 1, name);
  requireGeVAxis(*out, "flux file '" + file + "'");
  return static_cast<TH1D*>(out);
}

FluxHistograms::FluxHistograms(const FluxSpec& spec) : spec_(spec) {
  if (spec_.file.empty()) throw std::runtime_error("FluxHistograms: the flux spec has no file");
  if (spec_.histograms.empty()) throw std::runtime_error("FluxHistograms: the flux spec names no histograms");
  if (spec_.potNorm <= 0.0) throw std::runtime_error("FluxHistograms: potNorm must be positive");
  if (spec_.binNormGeV < 0.0) throw std::runtime_error("FluxHistograms: binNormGeV must be >= 0");
  for (const auto& kv : spec_.histograms) {
    hists_[kv.first].reset(load(spec_.file, kv.second, spec_.units, /*contentIsDensity=*/true));
  }
}

FluxHistograms::~FluxHistograms() = default;

std::vector<int> FluxHistograms::pdgs() const {
  std::vector<int> out;
  for (const auto& kv : hists_) out.push_back(kv.first);
  return out;
}

const TH1D& FluxHistograms::histogram(int pdg) const {
  auto it = hists_.find(pdg);
  if (it == hists_.end()) {
    throw std::out_of_range("FluxHistograms: no flux for PDG " + std::to_string(pdg) + " in '" + spec_.file + "'");
  }
  return *it->second;
}

double FluxHistograms::integral(int pdg) const { return histogram(pdg).Integral("width"); }

double FluxHistograms::integralPerPOT(const TH1& h, double potNorm, double binNormGeV) {
  if (binNormGeV > 0.0) return h.Integral("width") / binNormGeV / potNorm;
  return h.Integral() / potNorm;
}

double FluxHistograms::integralPerPOT(int pdg) const {
  return integralPerPOT(histogram(pdg), spec_.potNorm, spec_.binNormGeV);
}

std::string FluxHistograms::describe() const {
  std::ostringstream os;
  os << spec_.describe() << "\n";
  for (const auto& kv : hists_) {
    const TH1D& h = *kv.second;
    os << "  " << NuSpecies::fromPdg(kv.first).name() << ": " << h.GetNbinsX() << " bins ["
       << h.GetXaxis()->GetXmin() << ", " << h.GetXaxis()->GetXmax() << "] GeV, " << integralPerPOT(kv.first)
       << " /cm^2/POT\n";
  }
  return os.str();
}

}  // namespace zaku
