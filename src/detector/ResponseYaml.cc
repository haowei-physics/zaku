// YAML representation of zaku::Response (see docs/response.md). yaml-cpp
// stays private to this file.
#include <TAxis.h>
#include <TH1.h>
#include <TSystem.h>
#include <yaml-cpp/yaml.h>

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "zaku/detector/Response.h"

namespace zaku {

namespace {

std::string expandEnv(const std::string& s) {
  std::string out;
  std::size_t i = 0;
  while (i < s.size()) {
    if (s[i] == '$' && i + 1 < s.size() && s[i + 1] == '{') {
      const std::size_t end = s.find('}', i + 2);
      if (end == std::string::npos) throw std::runtime_error("unterminated ${...} in '" + s + "'");
      const std::string var = s.substr(i + 2, end - i - 2);
      const char* val = std::getenv(var.c_str());
      if (!val) throw std::runtime_error("environment variable '" + var + "' used in '" + s + "' is not set");
      out += val;
      i = end + 1;
    } else {
      out += s[i++];
    }
  }
  return out;
}

std::string resolvePath(const std::string& raw, const std::string& baseDir) {
  const std::string s = expandEnv(raw);
  if (s.empty() || s[0] == '/' || baseDir.empty()) return s;
  return baseDir + "/" + s;
}

/// "file:name" -> {file, name}
void splitHistRef(const std::string& ref, std::string& file, std::string& name) {
  const std::string::size_type pos = ref.rfind(':');
  if (pos == std::string::npos || pos == 0 || pos + 1 >= ref.size()) {
    throw std::runtime_error("response: histogram reference '" + ref + "' must be 'file.root:histogram'");
  }
  file = ref.substr(0, pos);
  name = ref.substr(pos + 1);
}

TableSource parseSource(const YAML::Node& n, const std::string& baseDir, const std::string& where) {
  if (!n) throw std::runtime_error("response: '" + where + "' is missing");
  if (n.IsScalar()) {
    // A number is a constant; anything else a formula.
    try {
      const double v = n.as<double>();
      const std::string s = n.as<std::string>();
      // yaml-cpp converts "1e3"-like strings; make sure it really is numeric text
      bool numeric = !s.empty();
      for (char c : s) {
        if (!(std::isdigit(static_cast<unsigned char>(c)) || c == '.' || c == '-' || c == '+' || c == 'e' || c == 'E')) {
          numeric = false;
          break;
        }
      }
      if (numeric) return TableSource::constant(v);
    } catch (const YAML::Exception&) {
      // not a number
    }
    return TableSource::formula(n.as<std::string>());
  }
  if (n.IsSequence()) {
    std::vector<double> v;
    for (const auto& e : n) v.push_back(e.as<double>());
    return TableSource::array(v);
  }
  if (n.IsMap()) {
    if (n["hist"]) {
      std::string file, name;
      splitHistRef(n["hist"].as<std::string>(), file, name);
      EnergyUnit units = EnergyUnit::GeV;
      if (n["units"]) units = energyUnitFromString(n["units"].as<std::string>());
      return TableSource::histogram(resolvePath(file, baseDir), name, units);
    }
    if (n["formula"]) return TableSource::formula(n["formula"].as<std::string>());
    if (n["constant"]) return TableSource::constant(n["constant"].as<double>());
    if (n["array"]) {
      std::vector<double> v;
      for (const auto& e : n["array"]) v.push_back(e.as<double>());
      return TableSource::array(v);
    }
  }
  throw std::runtime_error("response: cannot interpret '" + where + "' (expected a number, a list, a formula string or {hist: ...})");
}

YAML::Node emitSource(const TableSource& s, const std::vector<TableAxis>& axes) {
  YAML::Node n;
  switch (s.kind()) {
    case TableSource::Kind::Constant:
      n = s.constantValue();
      break;
    case TableSource::Kind::Array: {
      for (double v : s.arrayValues()) n.push_back(v);
      n.SetStyle(YAML::EmitterStyle::Flow);
      break;
    }
    case TableSource::Kind::Histogram: {
      const std::string& label = s.histogramLabel();
      if (label.find(':') != std::string::npos) {
        n["hist"] = label;
      } else {
        // An in-memory histogram (e.g. from the legacy file): write its contents
        // as an inline array over the table axes, last axis fastest.
        const TH1* h = s.histogramPtr();
        const int nx = h->GetNbinsX();
        const int ny = h->GetDimension() >= 2 ? h->GetNbinsY() : 1;
        const int nz = h->GetDimension() >= 3 ? h->GetNbinsZ() : 1;
        for (int ix = 1; ix <= nx; ++ix) {
          for (int iy = 1; iy <= ny; ++iy) {
            for (int iz = 1; iz <= nz; ++iz) n.push_back(h->GetBinContent(h->GetBin(ix, iy, iz)));
          }
        }
        n.SetStyle(YAML::EmitterStyle::Flow);
        (void)axes;
      }
      break;
    }
    case TableSource::Kind::Formula:
      n = s.formulaExpression();
      break;
  }
  return n;
}

}  // namespace

Response Response::fromYamlString(const std::string& yaml, const std::string& baseDir) {
  YAML::Node root;
  try {
    root = YAML::Load(yaml);
  } catch (const YAML::Exception& e) {
    throw std::runtime_error(std::string("response: YAML parse error: ") + e.what());
  }
  if (!root.IsMap()) throw std::runtime_error("response: the top level must be a map");
  if (root["legacy_file"]) {
    Response r = fromLegacyFile(resolvePath(root["legacy_file"].as<std::string>(), baseDir));
    if (root["normalise_probabilities"]) r.setNormaliseProbabilities(root["normalise_probabilities"].as<bool>());
    return r;
  }
  Response r;
  if (root["categories"]) {
    for (const auto& c : root["categories"]) {
      if (!c["name"] || !c["select"]) throw std::runtime_error("response: every category needs 'name' and 'select'");
      r.addCategory(c["name"].as<std::string>(), c["select"].as<std::string>());
    }
  }
  if (root["tables"]) {
    int index = 0;
    for (const auto& tn : root["tables"]) {
      const std::string where = "tables[" + std::to_string(index++) + "]";
      for (const auto& kv : tn) {
        const std::string key = kv.first.as<std::string>();
        if (key != "category" && key != "select" && key != "beam" && key != "axes" && key != "efficiency" &&
            key != "probabilities") {
          throw std::runtime_error("response: unknown key '" + key + "' in " + where);
        }
      }
      ResponseTable t;
      if (tn["category"]) t.category = tn["category"].as<std::string>();
      if (tn["select"]) t.select = tn["select"].as<std::string>();
      if (tn["beam"]) {
        const std::string b = tn["beam"].as<std::string>();
        if (b != "any" && b != "all" && !b.empty()) t.beamMode = static_cast<int>(beamModeFromString(b));
      }
      if (tn["axes"]) {
        for (const auto& kv : tn["axes"]) {
          TableAxis a;
          a.variable = kv.first.as<std::string>();
          for (const auto& e : kv.second) a.edges.push_back(e.as<double>());
          if (a.edges.size() < 2) throw std::runtime_error("response: axis '" + a.variable + "' in " + where + " needs at least two edges");
          t.axes.push_back(a);
        }
      }
      t.efficiency = tn["efficiency"] ? parseSource(tn["efficiency"], baseDir, where + ".efficiency") : TableSource::constant(1.0);
      if (tn["probabilities"]) {
        for (const auto& kv : tn["probabilities"]) {
          t.columns.push_back(kv.first.as<std::string>());
          t.probabilities.push_back(parseSource(kv.second, baseDir, where + ".probabilities." + kv.first.as<std::string>()));
        }
      }
      r.addTable(t);
    }
  }
  if (root["normalise_probabilities"]) r.setNormaliseProbabilities(root["normalise_probabilities"].as<bool>());
  return r;
}

Response Response::fromFile(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("Response::fromFile: cannot read '" + path + "'");
  std::stringstream ss;
  ss << in.rdbuf();
  std::string dir = gSystem->GetDirName(path.c_str()).Data();
  if (dir.empty() || dir == ".") dir = gSystem->WorkingDirectory();
  if (dir[0] != '/') dir = std::string(gSystem->WorkingDirectory()) + "/" + dir;
  return fromYamlString(ss.str(), dir);
}

std::string Response::toYaml() const {
  YAML::Node root;
  YAML::Node cats;
  for (const auto& c : categories_) {
    YAML::Node n;
    n["name"] = c.name;
    n["select"] = c.select;
    n.SetStyle(YAML::EmitterStyle::Flow);
    cats.push_back(n);
  }
  root["categories"] = cats;
  YAML::Node tables;
  for (const auto& t : tables_) {
    YAML::Node n;
    if (!t.category.empty()) n["category"] = t.category;
    if (!t.select.empty()) n["select"] = t.select;
    if (t.beamMode >= 0) n["beam"] = toString(static_cast<BeamMode>(t.beamMode));
    if (!t.axes.empty()) {
      YAML::Node axes;
      for (const auto& a : t.axes) {
        YAML::Node e;
        for (double v : a.edges) e.push_back(v);
        e.SetStyle(YAML::EmitterStyle::Flow);
        axes[a.variable] = e;
      }
      n["axes"] = axes;
    }
    n["efficiency"] = emitSource(t.efficiency, t.axes);
    YAML::Node probs;
    for (std::size_t c = 0; c < t.columns.size(); ++c) probs[t.columns[c]] = emitSource(t.probabilities[c], t.axes);
    n["probabilities"] = probs;
    tables.push_back(n);
  }
  root["tables"] = tables;
  root["normalise_probabilities"] = normalise_;
  std::ostringstream os;
  os << "# Zaku response tables (see docs/response.md): for true events in each category and\n"
        "# bin, the detection probability (efficiency) and the probability of each output column.\n"
     << YAML::Dump(root) << "\n";
  return os.str();
}

}  // namespace zaku
