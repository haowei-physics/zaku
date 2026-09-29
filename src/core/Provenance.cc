#include "zaku/core/Provenance.h"

#include <TDirectory.h>
#include <TNamed.h>

#include <stdexcept>

#include "zaku/core/Version.h"

namespace zaku {

namespace {

void writeNamed(TDirectory& d, const char* name, const std::string& value) {
  TNamed n(name, value.c_str());
  d.WriteTObject(&n, name, "Overwrite");
}

std::string readNamed(TDirectory& d, const char* name) {
  TNamed* n = d.Get<TNamed>(name);
  if (!n) return "";
  return n->GetTitle();
}

}  // namespace

void Provenance::stamp(TDirectory& dir, const std::string& stage, const std::string& text) {
  TDirectory* d = dir.GetDirectory("zaku");
  if (!d) d = dir.mkdir("zaku");
  if (!d) throw std::runtime_error("Provenance::stamp: cannot create the 'zaku' directory in " + std::string(dir.GetName()));
  writeNamed(*d, "version", ::zaku::version());
  writeNamed(*d, "schema", std::to_string(::zaku::schemaVersion()));
  writeNamed(*d, "stage", stage);
  writeNamed(*d, "provenance", text);
}

Provenance Provenance::read(TDirectory& dir) {
  Provenance p;
  TDirectory* d = dir.GetDirectory("zaku");
  if (!d) return p;
  p.version = readNamed(*d, "version");
  const std::string schema = readNamed(*d, "schema");
  p.schema = schema.empty() ? 0 : std::stoi(schema);
  p.stage = readNamed(*d, "stage");
  p.text = readNamed(*d, "provenance");
  return p;
}

void Provenance::requireCompatible(TDirectory& dir, const std::string& context) {
  const Provenance p = read(dir);
  if (p.schema > schemaVersion()) {
    throw std::runtime_error(context + ": file was written with Zaku schema " + std::to_string(p.schema) +
                             " (library " + p.version + "), newer than this library's schema " +
                             std::to_string(schemaVersion()));
  }
}

}  // namespace zaku
