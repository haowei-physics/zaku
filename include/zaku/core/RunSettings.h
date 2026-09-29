#ifndef ZAKU_CORE_RUNSETTINGS_H
#define ZAKU_CORE_RUNSETTINGS_H
#include <string>

namespace zaku {

/// Where a production lives and how its random numbers are seeded.
struct RunSettings {
  std::string workDir;      ///< every stage reads and writes below this directory
  unsigned seed = 20250929; ///< master seed; per-sample seeds are derived from it
  bool legacyRngOrder = false;  ///< reproduce the legacy macros' random-number sequence exactly

  bool valid() const { return !workDir.empty(); }
};

}  // namespace zaku

#endif  // ZAKU_CORE_RUNSETTINGS_H
