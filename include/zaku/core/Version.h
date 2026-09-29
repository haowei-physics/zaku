#ifndef ZAKU_CORE_VERSION_H
#define ZAKU_CORE_VERSION_H
#include <string>

namespace zaku {

/// Library version, e.g. "0.1.0".
std::string version();
int versionMajor();
int versionMinor();
int versionPatch();

/// Version of the file formats and configuration keys Zaku writes. Bumped
/// whenever an output tree, parameter file or YAML key changes meaning.
/// Every file Zaku writes carries it, so readers can refuse stale inputs.
int schemaVersion();

/// One-line description including the ROOT version Zaku was built against.
std::string buildInfo();

}  // namespace zaku

#endif  // ZAKU_CORE_VERSION_H
