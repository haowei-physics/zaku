#include "zaku/core/Manifest.h"

#include <TSystem.h>

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace zaku {

namespace {

bool exists(const std::string& path) { return !gSystem->AccessPathName(path.c_str()); }

std::string jsonEscape(const std::string& s) {
  std::string out;
  for (char c : s) {
    switch (c) {
      case '"': out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\n': out += "\\n"; break;
      case '\t': out += "\\t"; break;
      default: out += c;
    }
  }
  return out;
}

}  // namespace

std::vector<std::string> ManifestJob::missingOutputs() const {
  std::vector<std::string> missing;
  for (const auto& o : outputs) {
    if (!exists(o)) missing.push_back(o);
  }
  return missing;
}

Manifest::Manifest(std::string stage) : stage_(std::move(stage)) {}

void Manifest::add(const ManifestJob& job) { jobs_.push_back(job); }

ManifestJob& Manifest::add(const std::string& name, const std::string& command) {
  ManifestJob j;
  j.name = name;
  j.command = command;
  jobs_.push_back(j);
  return jobs_.back();
}

std::vector<const ManifestJob*> Manifest::pending() const {
  std::vector<const ManifestJob*> out;
  for (const auto& j : jobs_) {
    if (!j.done()) out.push_back(&j);
  }
  return out;
}

int Manifest::nDone() const {
  int n = 0;
  for (const auto& j : jobs_) {
    if (j.done()) ++n;
  }
  return n;
}

std::string Manifest::quote(const std::string& arg) {
  if (arg.empty()) return "''";
  bool safe = true;
  for (char c : arg) {
    if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '/' || c == '.' || c == '_' || c == '-' || c == ':' ||
          c == '=' || c == ',' || c == '+' || c == '@' || c == '%')) {
      safe = false;
      break;
    }
  }
  if (safe) return arg;
  std::string out = "'";
  for (char c : arg) {
    if (c == '\'') out += "'\\''";
    else out += c;
  }
  out += "'";
  return out;
}

std::string Manifest::text() const {
  std::ostringstream os;
  os << "# zaku manifest: stage " << stage_ << "\n";
  os << "# One job per non-comment line; every line is self-contained (it changes into its\n";
  os << "# own working directory). Run them by hand, with e.g.\n";
  os << "#   grep -v '^#' " << stage_ << ".txt | xargs -d '\\n' -P 4 -I{} bash -c '{}'\n";
  os << "# or from a batch array that picks line $TASK_ID of the command list.\n";
  for (const auto& p : preamble_) os << "# env: " << p << "\n";
  for (const auto& j : jobs_) {
    os << "# job: " << j.name;
    if (!j.comment.empty()) os << "  (" << j.comment << ")";
    os << "\n";
    if (!j.outputs.empty()) {
      os << "# outputs:";
      for (const auto& o : j.outputs) os << " " << o;
      os << "\n";
    }
    if (!j.workDir.empty()) os << "cd " << quote(j.workDir) << " && ";
    os << j.command << "\n";
  }
  return os.str();
}

std::string Manifest::json() const {
  std::ostringstream os;
  os << "{\n  \"stage\": \"" << jsonEscape(stage_) << "\",\n  \"preamble\": [";
  for (std::size_t i = 0; i < preamble_.size(); ++i) {
    os << (i ? ", " : "") << "\"" << jsonEscape(preamble_[i]) << "\"";
  }
  os << "],\n  \"jobs\": [\n";
  for (std::size_t i = 0; i < jobs_.size(); ++i) {
    const auto& j = jobs_[i];
    os << "    {\"name\": \"" << jsonEscape(j.name) << "\", \"command\": \"" << jsonEscape(j.command)
       << "\", \"workDir\": \"" << jsonEscape(j.workDir) << "\", \"comment\": \"" << jsonEscape(j.comment)
       << "\", \"outputs\": [";
    for (std::size_t k = 0; k < j.outputs.size(); ++k) {
      os << (k ? ", " : "") << "\"" << jsonEscape(j.outputs[k]) << "\"";
    }
    os << "]}" << (i + 1 < jobs_.size() ? "," : "") << "\n";
  }
  os << "  ]\n}\n";
  return os.str();
}

std::string Manifest::write(const std::string& dir) const {
  if (stage_.empty()) throw std::runtime_error("Manifest::write: the manifest has no stage name");
  if (gSystem->mkdir(dir.c_str(), true) < 0 && !exists(dir)) {
    throw std::runtime_error("Manifest::write: cannot create directory '" + dir + "'");
  }
  const std::string txt = dir + "/" + stage_ + ".txt";
  const std::string js = dir + "/" + stage_ + ".json";
  {
    std::ofstream f(txt);
    if (!f) throw std::runtime_error("Manifest::write: cannot write '" + txt + "'");
    f << text();
  }
  {
    std::ofstream f(js);
    if (!f) throw std::runtime_error("Manifest::write: cannot write '" + js + "'");
    f << json();
  }
  return txt;
}

std::string Manifest::summary() const {
  std::ostringstream os;
  os << stage_ << ": " << jobs_.size() << " job" << (jobs_.size() == 1 ? "" : "s") << ", " << nDone() << " done";
  return os.str();
}

}  // namespace zaku
