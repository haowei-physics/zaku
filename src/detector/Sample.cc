#include "zaku/detector/Sample.h"

#include <sstream>
#include <stdexcept>

namespace zaku {

std::string SampleDefinition::describe() const {
  std::ostringstream os;
  os << name;
  if (!title.empty()) os << " (" << title << ")";
  os << ": beam " << (beamMode < 0 ? std::string("any") : toString(static_cast<BeamMode>(beamMode))) << ", column "
     << columnName() << ", hypothesis " << toString(hypothesis) << ", " << binning.describe();
  if (!cut.empty()) os << ", cut '" << cut << "'";
  return os.str();
}

std::string toString(Assignment a) { return a == Assignment::Split ? "split" : "draw"; }

Assignment assignmentFromString(const std::string& s) {
  if (s == "split" || s == "weights" || s == "deterministic") return Assignment::Split;
  if (s == "draw" || s == "random" || s == "legacy") return Assignment::Draw;
  throw std::invalid_argument("assignmentFromString: unknown assignment '" + s + "' (split, draw)");
}

}  // namespace zaku
