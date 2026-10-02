#include "zaku/core/DialCatalog.h"

#include <yaml-cpp/yaml.h>

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace zaku {

std::string DialSpec::describe() const {
  std::ostringstream os;
  os << name << ": tweaks";
  for (double t : tweaks) os << " " << t;
  os << (fit ? "" : " (not fitted)");
  if (!note.empty()) os << " - " << note;
  return os.str();
}

DialCatalog& DialCatalog::add(const DialSpec& d) {
  if (d.name.empty()) throw std::invalid_argument("DialCatalog::add: a dial needs a name");
  if (has(d.name)) throw std::invalid_argument("DialCatalog::add: dial '" + d.name + "' is already listed");
  dials_.push_back(d);
  return *this;
}

DialCatalog& DialCatalog::add(const std::string& name, const std::vector<double>& tweaks, bool fit,
                              const std::string& note) {
  DialSpec d;
  d.name = name;
  d.tweaks = tweaks;
  d.fit = fit;
  d.note = note;
  return add(d);
}

std::vector<std::string> DialCatalog::names() const {
  std::vector<std::string> out;
  for (const auto& d : dials_) out.push_back(d.name);
  return out;
}

std::vector<std::string> DialCatalog::fitted(const std::vector<std::string>& exclude) const {
  std::vector<std::string> out;
  for (const auto& d : dials_) {
    if (!d.fit) continue;
    bool excluded = false;
    for (const auto& e : exclude) {
      if (e == d.name) {
        excluded = true;
        break;
      }
    }
    if (!excluded) out.push_back(d.name);
  }
  return out;
}

int DialCatalog::index(const std::string& name) const {
  for (std::size_t i = 0; i < dials_.size(); ++i) {
    if (dials_[i].name == name) return static_cast<int>(i);
  }
  return -1;
}

const DialSpec& DialCatalog::dial(const std::string& name) const {
  const int i = index(name);
  if (i < 0) throw std::out_of_range("DialCatalog: no dial named '" + name + "'");
  return dials_[static_cast<std::size_t>(i)];
}

DialCatalog DialCatalog::fromFile(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("DialCatalog::fromFile: cannot read '" + path + "'");
  std::stringstream ss;
  ss << in.rdbuf();
  try {
    return fromYamlString(ss.str());
  } catch (const std::exception& e) {
    throw std::runtime_error("DialCatalog::fromFile: '" + path + "': " + e.what());
  }
}

DialCatalog DialCatalog::fromYamlString(const std::string& yaml) {
  YAML::Node root;
  try {
    root = YAML::Load(yaml);
  } catch (const YAML::Exception& e) {
    throw std::runtime_error(std::string("dial catalogue: YAML parse error: ") + e.what());
  }
  DialCatalog c;
  YAML::Node list;
  if (root.IsSequence()) {
    list = root;
  } else if (root.IsMap()) {
    for (const auto& kv : root) {
      const std::string key = kv.first.as<std::string>();
      if (key != "generator" && key != "dials") {
        throw std::runtime_error("dial catalogue: unknown key '" + key + "' (allowed: generator, dials)");
      }
    }
    if (root["generator"]) c.generator_ = root["generator"].as<std::string>();
    list = root["dials"];
  }
  if (!list || !list.IsSequence()) throw std::runtime_error("dial catalogue: 'dials' must be a list");
  for (const auto& n : list) {
    if (!n.IsMap()) throw std::runtime_error("dial catalogue: every dial must be a map {name, tweaks, fit, note}");
    for (const auto& kv : n) {
      const std::string key = kv.first.as<std::string>();
      if (key != "name" && key != "tweaks" && key != "fit" && key != "note") {
        throw std::runtime_error("dial catalogue: unknown key '" + key + "' in a dial (allowed: name, tweaks, fit, note)");
      }
    }
    DialSpec d;
    if (!n["name"]) throw std::runtime_error("dial catalogue: a dial has no 'name'");
    d.name = n["name"].as<std::string>();
    if (n["tweaks"]) {
      for (const auto& t : n["tweaks"]) d.tweaks.push_back(t.as<double>());
    }
    d.fit = n["fit"] ? n["fit"].as<bool>() : true;
    d.note = n["note"] ? n["note"].as<std::string>() : "";
    c.add(d);
  }
  return c;
}

std::string DialCatalog::toYaml() const {
  YAML::Node root;
  root["generator"] = generator_;
  YAML::Node list;
  for (const auto& d : dials_) {
    YAML::Node n;
    n["name"] = d.name;
    YAML::Node tw;
    for (double t : d.tweaks) tw.push_back(t);
    tw.SetStyle(YAML::EmitterStyle::Flow);
    n["tweaks"] = tw;
    if (!d.fit) n["fit"] = false;
    if (!d.note.empty()) n["note"] = d.note;
    n.SetStyle(YAML::EmitterStyle::Flow);
    list.push_back(n);
  }
  root["dials"] = list;
  std::ostringstream os;
  os << YAML::Dump(root) << "\n";
  return os.str();
}

DialCatalog DialCatalog::neutDefault() {
  DialCatalog c;
  c.setGenerator("neut");
  const std::vector<double> two = {-2, -1, 1, 2};
  const std::vector<double> one = {-1, -0.5, 0.5, 1};
  c.add("MaCCQE", two, true, "CCQE axial mass");
  c.add("MaRES", two, true, "resonant axial mass");
  c.add("CA5RES", two, true, "resonant axial form factor at Q2=0");
  c.add("BgSclRES", two, true, "non-resonant I=1/2 background scale");
  c.add("MDLSPiEj", {-3, -2}, false, "single-pion ejection model switch; not a continuous parameter");
  c.add("NucleonFSI_Double", one, true, "nucleon FSI");
  c.add("NucleonFSI_Elastic", one, true, "nucleon FSI");
  c.add("NucleonFSI_Single", one, true, "nucleon FSI");
  c.add("NucleonFSI_Total", one, true, "nucleon FSI");
  c.add("NucleonFSI_PB", one, false, "Pauli blocking in the cascade; excluded by the legacy fit");
  c.add("PionFSI_AbsProb", one, true, "pion FSI");
  c.add("PionFSI_QELowMomProb", one, true, "pion FSI");
  c.add("PionFSI_QEHighMomProb", one, true, "pion FSI");
  c.add("PionFSI_CExLowMomProb", one, true, "pion FSI");
  c.add("PionFSI_CExHighMomProb", one, true, "pion FSI");
  c.add("PionFSI_InelProb", {-0.9, -0.5, 0.5, 1}, true, "pion FSI; -0.9 sigma keeps the probability positive (non-uniform knots)");
  c.add("SRCFrac_C", {0.25, 0.5}, false, "carbon SRC fraction; the detector is water");
  c.add("SRCFrac_O", {0.25, 0.5}, true, "oxygen SRC fraction (one-sided)");
  c.add("RadCorr_rate", {0, 0.25, 0.75}, false, "radiative-correction rate; excluded by the legacy fit");
  return c;
}

std::string DialCatalog::describe() const {
  std::ostringstream os;
  os << generator_ << " dial catalogue: " << dials_.size() << " dials, " << fitted().size() << " fitted\n";
  for (const auto& d : dials_) os << "  " << d.describe() << "\n";
  return os.str();
}

}  // namespace zaku
