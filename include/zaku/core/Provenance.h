#ifndef ZAKU_CORE_PROVENANCE_H
#define ZAKU_CORE_PROVENANCE_H
#include <string>

class TDirectory;

namespace zaku {

/// What Zaku records in every ROOT file it writes, in a `zaku` directory:
///   zaku/version     TNamed  title = library version
///   zaku/schema      TNamed  title = schema version
///   zaku/stage       TNamed  title = the stage that wrote the file
///   zaku/provenance  TNamed  title = free text (resolved configuration, seeds, inputs)
struct Provenance {
  std::string version;
  int schema = 0;
  std::string stage;
  std::string text;

  /// Write (or overwrite) the stamp into `dir`.
  static void stamp(TDirectory& dir, const std::string& stage, const std::string& text);
  /// Read the stamp back; `schema` is 0 when the file carries none.
  static Provenance read(TDirectory& dir);
  /// Throw std::runtime_error if the file has a stamp with a newer schema than this library.
  static void requireCompatible(TDirectory& dir, const std::string& context);
};

}  // namespace zaku

#endif  // ZAKU_CORE_PROVENANCE_H
