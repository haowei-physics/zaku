#include "zaku/core/McSample.h"

#include <sstream>

namespace zaku {

bool McSample::isAppearance() const {
  const int a = fluxPdg < 0 ? -fluxPdg : fluxPdg;
  const int b = interactionPdg < 0 ? -interactionPdg : interactionPdg;
  return a != b;
}

NuSpecies McSample::fluxSpecies() const { return NuSpecies::fromPdg(fluxPdg); }
NuSpecies McSample::interactionSpecies() const { return NuSpecies::fromPdg(interactionPdg); }
std::string McSample::name() const { return SampleLayout::sampleName(*this); }
std::string McSample::directory() const { return SampleLayout::sampleDirectory(*this); }

std::string McSample::describe() const {
  std::ostringstream os;
  os << name() << ": " << toString(beamMode) << " beam, " << fluxSpecies().name() << " flux -> "
     << interactionSpecies().name() << " interactions, " << nEvents << " events, seed " << seed
     << (isAppearance() ? " (appearance)" : "");
  return os.str();
}

namespace SampleLayout {

namespace {
std::string countTag(long n) {
  if (n > 0 && n % 1000 == 0) return std::to_string(n / 1000) + "K";
  return std::to_string(n);
}
}  // namespace

std::string sampleName(const McSample& s) {
  return s.fluxSpecies().name() + "_x_" + s.interactionSpecies().name() + "_" + countTag(s.nEvents) + "_" +
         toString(s.beamMode);
}
std::string sampleDirectory(const McSample& s) { return s.fluxSpecies().name(); }
std::string generatorFile(const McSample& s) { return sampleDirectory(s) + "/" + sampleName(s) + ".root"; }
std::string genieFile(const McSample& s) { return sampleDirectory(s) + "/" + sampleName(s) + ".ghep.root"; }
std::string flatFile(const McSample& s) { return sampleDirectory(s) + "/" + sampleName(s) + "_flattree.root"; }
std::string weightFile(const McSample& s) { return sampleDirectory(s) + "/" + sampleName(s) + "_Weight.flat.root"; }
std::string genieWeightFile(const McSample& s, const std::string& dial) {
  return sampleDirectory(s) + "/" + sampleName(s) + "_rw_" + dial + ".root";
}
std::string reconFile(const McSample& s) { return sampleDirectory(s) + "/" + sampleName(s) + "_fakerecon.root"; }
std::string splineFile(const McSample& s) { return sampleDirectory(s) + "/" + sampleName(s) + "_spline.root"; }
std::string jobDirectory(const McSample& s) { return "gen/" + sampleName(s); }
std::string selectFile(BeamMode mode) { return "fakeselect_" + toString(mode) + ".root"; }
std::string rateTableFile(BeamMode mode) { return "fakeselect_" + toString(mode) + "_rates.csv"; }
std::string parametersDirectory() { return "parameters"; }
std::string manifestsDirectory() { return "manifests"; }
std::string uniformFluxFile(BeamMode mode) { return "flux/uniform_" + toString(mode) + ".root"; }

}  // namespace SampleLayout

}  // namespace zaku
