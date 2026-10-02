#ifndef ZAKU_CORE_PRODUCTION_H
#define ZAKU_CORE_PRODUCTION_H
#include <string>
#include <vector>

#include "zaku/core/Types.h"

namespace zaku {

/// One generator production: a beam mode, the flux species the events are
/// thrown on, the species that interacts (different for appearance
/// productions, which are generated as "100% converted" and weighted by the
/// oscillation probability later), the number of events and the random seed.
///
/// A production is a generator output file, not an analysis sample: samples
/// are the reconstructed and selected event categories defined later.
struct Production {
  BeamMode beamMode = BeamMode::FHC;
  int fluxPdg = 14;         ///< species whose flux is sampled
  int interactionPdg = 14;  ///< species that interacts (== fluxPdg for disappearance productions)
  long nEvents = 0;
  unsigned seed = 0;        ///< 0 = derived from the run seed and the production name

  Production() = default;
  Production(BeamMode mode, int flux, int interaction, long n, unsigned s = 0)
      : beamMode(mode), fluxPdg(flux), interactionPdg(interaction), nEvents(n), seed(s) {}

  bool isAppearance() const;
  NuSpecies fluxSpecies() const;
  NuSpecies interactionSpecies() const;
  /// "numu_x_nue_50K_FHC": flux species, interaction species, event count, beam mode.
  std::string name() const;
  /// Subdirectory used for this production's files: the flux species name.
  std::string directory() const;
  std::string describe() const;
};

/// The one place that knows how the files of a production are named. All
/// paths are relative to the run work directory.
namespace Layout {
/// "<flux>_x_<int>_<N>_<mode>" with N as "200K" for multiples of 1000, else the number.
std::string productionName(const Production& s);
std::string productionDirectory(const Production& s);
std::string generatorFile(const Production& s);   ///< <dir>/<name>.root (NEUT) — GENIE: .ghep.root
std::string genieFile(const Production& s);       ///< <dir>/<name>.ghep.root
std::string flatFile(const Production& s);        ///< <dir>/<name>_flattree.root
std::string weightFile(const Production& s);      ///< <dir>/<name>_Weight.flat.root (T2KReWeight)
std::string genieWeightFile(const Production& s, const std::string& dial);  ///< <dir>/<name>_rw_<dial>.root
std::string reconFile(const Production& s);       ///< <dir>/<name>_fakerecon.root
std::string splineFile(const Production& s);      ///< <dir>/<name>_spline.root
std::string jobDirectory(const Production& s);    ///< gen/<name>: per-job working directory for generator runs
std::string selectFile(BeamMode mode);          ///< fakeselect_<mode>.root
std::string rateTableFile(BeamMode mode);       ///< fakeselect_<mode>_rates.csv
std::string parametersDirectory();              ///< parameters
std::string oscCovarianceFile();                ///< parameters/oscProb/oscCov.root
std::string oscGridFile();                      ///< parameters/oscProb/grid.root
std::string detCovarianceFile();                ///< parameters/det/detCov.root
std::string detBinningFile();                   ///< parameters/det/detBinning.txt
std::string xsecCovarianceFile();               ///< parameters/xsec/xsecCov.root
std::string fluxCovarianceFile(bool ndConstrained);  ///< parameters/flux/fluxCov.root or fluxCov_ND.root
std::string fluxBinningFile();                  ///< parameters/flux/fluxBinning.txt
std::string summaryDirectory();                 ///< summary
std::string summaryCsvFile();                   ///< summary/systematics_summary.csv
std::string summaryRootFile();                  ///< summary/systematics_summary.root
std::string manifestsDirectory();               ///< manifests
std::string uniformFluxFile(BeamMode mode);     ///< flux/uniform_<mode>.root
}  // namespace Layout

}  // namespace zaku

#endif  // ZAKU_CORE_PRODUCTION_H
