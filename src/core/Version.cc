#include "zaku/core/Version.h"

#include <RVersion.h>

#include "zaku/core/VersionInfo.h"

namespace zaku {

std::string version() { return ZAKU_VERSION_STRING; }
int versionMajor() { return ZAKU_VERSION_MAJOR; }
int versionMinor() { return ZAKU_VERSION_MINOR; }
int versionPatch() { return ZAKU_VERSION_PATCH; }

int schemaVersion() { return 1; }

std::string buildInfo() {
  return std::string("Zaku ") + ZAKU_VERSION_STRING + " (schema " + std::to_string(schemaVersion()) +
         ", ROOT " + ROOT_RELEASE + ")";
}

}  // namespace zaku
