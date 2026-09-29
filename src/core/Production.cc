#include "zaku/core/Production.h"

#include <sstream>

namespace zaku {

bool Production::isAppearance() const {
  const int a = fluxPdg < 0 ? -fluxPdg : fluxPdg;
  const int b = interactionPdg < 0 ? -interactionPdg : interactionPdg;
  return a != b;
}

NuSpecies Production::fluxSpecies() const { return NuSpecies::fromPdg(fluxPdg); }
NuSpecies Production::interactionSpecies() const { return NuSpecies::fromPdg(interactionPdg); }
std::string Production::name() const { return Layout::productionName(*this); }
std::string Production::directory() const { return Layout::productionDirectory(*this); }

std::string Production::describe() const {
  std::ostringstream os;
  os << name() << ": " << toString(beamMode) << " beam, " << fluxSpecies().name() << " flux -> "
     << interactionSpecies().name() << " interactions, " << nEvents << " events, seed " << seed
     << (isAppearance() ? " (appearance)" : "");
  return os.str();
}

namespace Layout {

namespace {
std::string countTag(long n) {
  if (n > 0 && n % 1000 == 0) return std::to_string(n / 1000) + "K";
  return std::to_string(n);
}
}  // namespace

std::string productionName(const Production& s) {
  return s.fluxSpecies().name() + "_x_" + s.interactionSpecies().name() + "_" + countTag(s.nEvents) + "_" +
         toString(s.beamMode);
}
std::string productionDirectory(const Production& s) { return s.fluxSpecies().name(); }
std::string generatorFile(const Production& s) { return productionDirectory(s) + "/" + productionName(s) + ".root"; }
std::string genieFile(const Production& s) { return productionDirectory(s) + "/" + productionName(s) + ".ghep.root"; }
std::string flatFile(const Production& s) { return productionDirectory(s) + "/" + productionName(s) + "_flattree.root"; }
std::string weightFile(const Production& s) { return productionDirectory(s) + "/" + productionName(s) + "_Weight.flat.root"; }
std::string genieWeightFile(const Production& s, const std::string& dial) {
  return productionDirectory(s) + "/" + productionName(s) + "_rw_" + dial + ".root";
}
std::string reconFile(const Production& s) { return productionDirectory(s) + "/" + productionName(s) + "_fakerecon.root"; }
std::string splineFile(const Production& s) { return productionDirectory(s) + "/" + productionName(s) + "_spline.root"; }
std::string jobDirectory(const Production& s) { return "gen/" + productionName(s); }
std::string selectFile(BeamMode mode) { return "fakeselect_" + toString(mode) + ".root"; }
std::string rateTableFile(BeamMode mode) { return "fakeselect_" + toString(mode) + "_rates.csv"; }
std::string parametersDirectory() { return "parameters"; }
std::string manifestsDirectory() { return "manifests"; }
std::string uniformFluxFile(BeamMode mode) { return "flux/uniform_" + toString(mode) + ".root"; }

}  // namespace Layout

}  // namespace zaku
