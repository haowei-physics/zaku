#ifndef ZAKU_DETECTOR_SAMPLE_H
#define ZAKU_DETECTOR_SAMPLE_H
#include <string>
#include <vector>

#include "zaku/core/Types.h"
#include "zaku/detector/Binning.h"
#include "zaku/detector/Reco.h"

namespace zaku {

/// An analysis sample: a reconstructed and selected event category that
/// GUNDAM fits. Its index in the experiment's sample list is `EventType`.
struct SampleDefinition {
  std::string name;              ///< identifier, e.g. FHC_1Rmu
  std::string title;             ///< display title (ROOT latex), e.g. "FHC 1R#mu"
  int beamMode = -1;             ///< -1 any, else BeamMode
  std::string column;            ///< response-table column feeding it (defaults to name)
  Hypothesis hypothesis = Hypothesis::Muon;
  Binning binning;
  std::string binningFile;       ///< when the binning comes from a GUNDAM txt file
  std::string cut;               ///< optional reco-level cut (Formula) applied after assignment

  const std::string& columnName() const { return column.empty() ? name : column; }
  bool appliesTo(BeamMode mode) const { return beamMode < 0 || beamMode == static_cast<int>(mode); }
  std::string describe() const;
};

/// How events are distributed over samples.
enum class Assignment : int {
  Split = 0,  ///< one entry per column with probability > 0, weighted by the probability (deterministic)
  Draw = 1    ///< one random draw per event against the running sum of the columns (legacy)
};
std::string toString(Assignment a);
Assignment assignmentFromString(const std::string& s);

}  // namespace zaku

#endif  // ZAKU_DETECTOR_SAMPLE_H
