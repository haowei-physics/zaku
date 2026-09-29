#ifndef ZAKU_CORE_MCSAMPLE_H
#define ZAKU_CORE_MCSAMPLE_H
#include <string>
#include <vector>

#include "zaku/core/Types.h"

namespace zaku {

/// One generator production: a beam mode, the flux species the events are
/// thrown on, the species that interacts (different for appearance samples,
/// which are generated as "100% converted" and weighted by the oscillation
/// probability later), the number of events and the random seed.
struct McSample {
  BeamMode beamMode = BeamMode::FHC;
  int fluxPdg = 14;         ///< species whose flux is sampled
  int interactionPdg = 14;  ///< species that interacts (== fluxPdg for disappearance samples)
  long nEvents = 0;
  unsigned seed = 0;        ///< 0 = derived from the run seed and the sample name

  McSample() = default;
  McSample(BeamMode mode, int flux, int interaction, long n, unsigned s = 0)
      : beamMode(mode), fluxPdg(flux), interactionPdg(interaction), nEvents(n), seed(s) {}

  bool isAppearance() const;
  NuSpecies fluxSpecies() const;
  NuSpecies interactionSpecies() const;
  /// "numu_x_nue_50K_FHC": flux species, interaction species, event count, beam mode.
  std::string name() const;
  /// Subdirectory used for this sample's files: the flux species name.
  std::string directory() const;
  std::string describe() const;
};

/// The one place that knows how the files of a production are named. All
/// paths are relative to the run work directory.
namespace SampleLayout {
/// "<flux>_x_<int>_<N>_<mode>" with N as "200K" for multiples of 1000, else the number.
std::string sampleName(const McSample& s);
std::string sampleDirectory(const McSample& s);
std::string generatorFile(const McSample& s);   ///< <dir>/<name>.root (NEUT) — GENIE: .ghep.root
std::string genieFile(const McSample& s);       ///< <dir>/<name>.ghep.root
std::string flatFile(const McSample& s);        ///< <dir>/<name>_flattree.root
std::string weightFile(const McSample& s);      ///< <dir>/<name>_Weight.flat.root (T2KReWeight)
std::string genieWeightFile(const McSample& s, const std::string& dial);  ///< <dir>/<name>_rw_<dial>.root
std::string reconFile(const McSample& s);       ///< <dir>/<name>_fakerecon.root
std::string splineFile(const McSample& s);      ///< <dir>/<name>_spline.root
std::string jobDirectory(const McSample& s);    ///< gen/<name>: per-job working directory for generator runs
std::string selectFile(BeamMode mode);          ///< fakeselect_<mode>.root
std::string rateTableFile(BeamMode mode);       ///< fakeselect_<mode>_rates.csv
std::string parametersDirectory();              ///< parameters
std::string manifestsDirectory();               ///< manifests
std::string uniformFluxFile(BeamMode mode);     ///< flux/uniform_<mode>.root
}  // namespace SampleLayout

}  // namespace zaku

#endif  // ZAKU_CORE_MCSAMPLE_H
