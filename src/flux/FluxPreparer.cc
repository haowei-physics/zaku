#include "zaku/flux/FluxPreparer.h"

#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TKey.h>
#include <TList.h>
#include <TSystem.h>

#include <cmath>
#include <memory>
#include <set>
#include <stdexcept>
#include <vector>

#include "zaku/flux/FluxHistograms.h"

namespace zaku {
namespace FluxPreparer {

TH1D* rebinUniform(const TH1& h, double widthGeV, const std::string& name) {
  if (widthGeV <= 0.0) throw std::invalid_argument("FluxPreparer::rebinUniform: width must be positive");
  const double lo = h.GetXaxis()->GetXmin();
  const double hi = h.GetXaxis()->GetXmax();
  const int nb = static_cast<int>(std::lround((hi - lo) / widthGeV));
  if (nb <= 0) throw std::invalid_argument("FluxPreparer::rebinUniform: histogram range is smaller than the bin width");
  const std::string outName = name.empty() ? std::string(h.GetName()) : name;
  auto* out = new TH1D(outName.c_str(), h.GetTitle(), nb, lo, hi);
  out->SetDirectory(nullptr);
  out->GetXaxis()->SetTitle(h.GetXaxis()->GetTitle());
  out->GetYaxis()->SetTitle(h.GetYaxis()->GetTitle());
  for (int i = 1; i <= nb; ++i) {
    const double centre = lo + (i - 0.5) * widthGeV;
    const int src = h.GetXaxis()->FindFixBin(centre);
    out->SetBinContent(i, h.GetBinContent(src));
  }
  out->SetEntries(h.GetEntries());
  return out;
}

FluxSpec makeUniform(const FluxSpec& spec, const std::string& outFile, double widthGeV, bool allHistograms) {
  std::unique_ptr<TFile> in(TFile::Open(spec.file.c_str(), "READ"));
  if (!in || in->IsZombie()) throw std::runtime_error("FluxPreparer: cannot open flux file '" + spec.file + "'");

  std::vector<std::string> names;
  if (allHistograms) {
    std::set<std::string> seen;
    TIter next(in->GetListOfKeys());
    while (auto* key = static_cast<TKey*>(next())) {
      if (!seen.insert(key->GetName()).second) continue;  // one cycle per name
      TObject* obj = key->ReadObj();
      auto* h = dynamic_cast<TH1*>(obj);
      if (h && h->GetDimension() == 1) names.push_back(key->GetName());
      delete obj;
    }
  } else {
    for (const auto& kv : spec.histograms) names.push_back(kv.second);
  }
  // The spec's own histograms must be there whichever list is used.
  for (const auto& kv : spec.histograms) {
    if (!in->Get<TH1>(kv.second.c_str())) {
      throw std::runtime_error("FluxPreparer: no histogram '" + kv.second + "' in '" + spec.file + "'");
    }
  }
  in->Close();

  const std::string dir = gSystem->GetDirName(outFile.c_str()).Data();
  if (!dir.empty() && gSystem->AccessPathName(dir.c_str())) gSystem->mkdir(dir.c_str(), true);
  std::unique_ptr<TFile> out(TFile::Open(outFile.c_str(), "RECREATE"));
  if (!out || out->IsZombie()) throw std::runtime_error("FluxPreparer: cannot create '" + outFile + "'");
  for (const auto& n : names) {
    std::unique_ptr<TH1D> h(FluxHistograms::load(spec.file, n, spec.units, /*contentIsDensity=*/true));
    std::unique_ptr<TH1D> u(rebinUniform(*h, widthGeV));
    // Write explicitly into the output file: loading the input changed gDirectory.
    out->WriteTObject(u.get(), u->GetName(), "Overwrite");
  }
  out->Close();

  FluxSpec result = spec.withFile(outFile);
  result.units = EnergyUnit::GeV;
  return result;
}

}  // namespace FluxPreparer
}  // namespace zaku
