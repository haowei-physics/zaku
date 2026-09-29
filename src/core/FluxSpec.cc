#include "zaku/core/FluxSpec.h"

#include <sstream>
#include <stdexcept>

namespace zaku {

const std::string& FluxSpec::histogram(int pdg) const {
  auto it = histograms.find(pdg);
  if (it == histograms.end()) {
    throw std::out_of_range("FluxSpec: no histogram for PDG " + std::to_string(pdg) + " in '" + file + "'");
  }
  return it->second;
}

std::vector<int> FluxSpec::pdgs() const {
  std::vector<int> out;
  out.reserve(histograms.size());
  for (const auto& kv : histograms) out.push_back(kv.first);
  return out;
}

std::vector<NuSpecies> FluxSpec::species() const {
  std::vector<NuSpecies> out;
  for (const auto& kv : histograms) out.push_back(NuSpecies::fromPdg(kv.first));
  return out;
}

FluxSpec FluxSpec::t2k(const std::string& file, const std::string& detector) {
  FluxSpec s;
  s.file = file;
  for (int pdg : {14, -14, 12, -12}) {
    s.histograms[pdg] = "enu_" + detector + "_" + NuSpecies::fromPdg(pdg).t2kSuffix();
  }
  return s;
}

FluxSpec FluxSpec::withFile(const std::string& newFile) const {
  FluxSpec s = *this;
  s.file = newFile;
  return s;
}

std::string FluxSpec::describe() const {
  std::ostringstream os;
  os << "flux '" << file << "' [" << toString(units) << ", per " << binNormGeV << " GeV per " << potNorm << " POT]:";
  for (const auto& kv : histograms) os << " " << NuSpecies::fromPdg(kv.first).name() << "=" << kv.second;
  return os.str();
}

}  // namespace zaku
