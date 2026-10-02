// YAML builder for zaku::Experiment: Experiment::fromConfig / toYaml.
// yaml-cpp stays private to this file.
#include <TSystem.h>
#include <yaml-cpp/yaml.h>

#include <cstdlib>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>

#include "ExperimentImpl.h"
#include "zaku/core/Experiment.h"

namespace zaku {

const std::vector<std::string>& experimentKnownKeys() {
  static const std::vector<std::string> k = {"name",     "baseline_km",  "matter",   "detector",       "flux_units",
                                             "beam_modes", "flux_covariance", "productions", "response", "smearing",
                                             "energy",   "binding_energy", "samples", "oscillation",  "systematics",
                                             "shifts",   "run"};
  return k;
}

const std::vector<std::string>& experimentRawKeys() {
  static const std::vector<std::string> k = {"generator", "gundam"};
  return k;
}

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
  if (s.size() > 1 && s[0] == '~' && s[1] == '/') return std::string(gSystem->HomeDirectory()) + s.substr(1);
  return baseDir + "/" + s;
}

template <class T>
T get(const YAML::Node& n, const std::string& key, const T& dflt) {
  if (!n[key]) return dflt;
  return n[key].as<T>();
}

template <class T>
T require(const YAML::Node& n, const std::string& key, const std::string& where) {
  if (!n[key]) throw std::runtime_error("config: '" + where + "' needs the key '" + key + "'");
  return n[key].as<T>();
}

void rejectUnknown(const YAML::Node& n, const std::vector<std::string>& allowed, const std::string& where) {
  if (!n.IsMap()) return;
  for (const auto& kv : n) {
    const std::string key = kv.first.as<std::string>();
    bool ok = false;
    for (const auto& a : allowed) {
      if (a == key) {
        ok = true;
        break;
      }
    }
    if (!ok) {
      std::string list;
      for (const auto& a : allowed) list += (list.empty() ? "" : ", ") + a;
      throw std::runtime_error("config: unknown key '" + key + "' in '" + where + "' (allowed: " + list + ")");
    }
  }
}

FluxSpec parseFlux(const YAML::Node& n, const std::string& baseDir, double defaultPotNorm, double defaultBinNormGeV,
                   const std::string& where) {
  rejectUnknown(n, {"file", "histograms", "units", "pot_norm", "bin_norm_GeV", "t2k_detector"}, where);
  FluxSpec f;
  f.file = resolvePath(require<std::string>(n, "file", where), baseDir);
  f.potNorm = defaultPotNorm;
  f.binNormGeV = defaultBinNormGeV;
  if (n["t2k_detector"]) {
    f.histograms = FluxSpec::t2k("", n["t2k_detector"].as<std::string>()).histograms;
  }
  if (n["histograms"]) {
    for (const auto& kv : n["histograms"]) {
      f.histograms[NuSpecies::fromName(kv.first.as<std::string>()).pdg()] = kv.second.as<std::string>();
    }
  }
  if (f.histograms.empty()) throw std::runtime_error("config: '" + where + "' names no histograms");
  if (n["units"]) f.units = energyUnitFromString(n["units"].as<std::string>());
  if (n["pot_norm"]) f.potNorm = n["pot_norm"].as<double>();
  if (n["bin_norm_GeV"]) f.binNormGeV = n["bin_norm_GeV"].as<double>();
  return f;
}

Detector parseDetector(const YAML::Node& n) {
  rejectUnknown(n, {"name", "material", "fiducial_mass_kt", "composition", "n_nucleons", "neut_target"}, "detector");
  const double kt = require<double>(n, "fiducial_mass_kt", "detector");
  Detector det;
  if (n["composition"]) {
    std::vector<Target> comp;
    for (const auto& t : n["composition"]) {
      rejectUnknown(t, {"Z", "A", "mass_fraction"}, "detector.composition");
      comp.emplace_back(require<int>(t, "Z", "composition"), require<int>(t, "A", "composition"),
                        require<double>(t, "mass_fraction", "composition"));
    }
    det = Detector(get<std::string>(n, "name", "custom"), kt, comp);
  } else if (n["material"]) {
    det = Detector::fromMaterial(n["material"].as<std::string>(), kt);
    if (n["name"]) det = Detector(n["name"].as<std::string>(), kt, det.composition());
  } else {
    throw std::runtime_error("config: 'detector' needs either 'material' or 'composition'");
  }
  if (n["n_nucleons"]) det.setNucleonsOverride(n["n_nucleons"].as<double>());
  if (n["neut_target"]) {
    const YAML::Node t = n["neut_target"];
    rejectUnknown(t, {"num_nucleons", "num_bound_protons", "num_free_protons"}, "detector.neut_target");
    NeutTarget nt;
    nt.numNucleons = require<int>(t, "num_nucleons", "neut_target");
    nt.numBoundProtons = require<int>(t, "num_bound_protons", "neut_target");
    nt.numFreeProtons = get<int>(t, "num_free_protons", 0);
    det.setNeutTarget(nt);
  }
  return det;
}

/// A block that is either inline or `{include: file}`. Returns the YAML text
/// of the block and the directory relative paths inside it resolve against.
struct Block {
  std::string text;
  std::string baseDir;
  std::string includePath;  ///< empty when inline
  YAML::Node node;
};

Block loadBlock(const YAML::Node& n, const std::string& baseDir, const std::string& where) {
  Block b;
  if (n.IsMap() && n["include"]) {
    b.includePath = resolvePath(n["include"].as<std::string>(), baseDir);
    std::ifstream in(b.includePath);
    if (!in) throw std::runtime_error("config: '" + where + "' includes '" + b.includePath + "' which cannot be read");
    std::stringstream ss;
    ss << in.rdbuf();
    b.text = ss.str();
    std::string dir = gSystem->GetDirName(b.includePath.c_str()).Data();
    if (dir.empty() || dir == ".") dir = gSystem->WorkingDirectory();
    if (dir[0] != '/') dir = std::string(gSystem->WorkingDirectory()) + "/" + dir;
    b.baseDir = dir;
    try {
      b.node = YAML::Load(b.text);
    } catch (const YAML::Exception& e) {
      throw std::runtime_error("config: '" + b.includePath + "': " + e.what());
    }
  } else {
    b.node = n;
    b.text = YAML::Dump(n);
    b.baseDir = baseDir;
  }
  return b;
}

BinningAxis parseAxis(const std::string& variable, const YAML::Node& n, const std::string& where) {
  BinningAxis a;
  a.variable = variable;
  if (n.IsSequence()) {
    for (const auto& e : n) a.edges.push_back(e.as<double>());
    if (a.edges.size() < 2) throw std::runtime_error("config: '" + where + "': axis '" + variable + "' needs at least two edges");
    return a;
  }
  if (n.IsMap()) {
    rejectUnknown(n, {"min", "max", "step", "extra"}, where + "." + variable);
    std::vector<double> extra;
    if (n["extra"]) {
      for (const auto& e : n["extra"]) extra.push_back(e.as<double>());
    }
    return BinningAxis::uniform(variable, require<double>(n, "min", where), require<double>(n, "max", where),
                                require<double>(n, "step", where), extra);
  }
  throw std::runtime_error("config: '" + where + "': axis '" + variable + "' must be a list of edges or {min, max, step}");
}

SampleDefinition parseSample(const YAML::Node& n, const std::string& baseDir, const std::string& where) {
  rejectUnknown(n, {"name", "title", "beam", "column", "hypothesis", "binning", "cut"}, where);
  SampleDefinition s;
  s.name = require<std::string>(n, "name", where);
  s.title = get<std::string>(n, "title", s.name);
  if (n["beam"]) {
    const std::string b = n["beam"].as<std::string>();
    if (b != "any" && b != "all") s.beamMode = static_cast<int>(beamModeFromString(b));
  }
  s.column = get<std::string>(n, "column", "");
  if (n["hypothesis"]) s.hypothesis = hypothesisFromString(n["hypothesis"].as<std::string>());
  s.cut = get<std::string>(n, "cut", "");
  const YAML::Node bn = n["binning"];
  if (!bn) throw std::runtime_error("config: sample '" + s.name + "' needs a 'binning'");
  rejectUnknown(bn, {"file", "variables", "edges", "values"}, where + ".binning");
  if (bn["file"]) {
    s.binningFile = resolvePath(bn["file"].as<std::string>(), baseDir);
    s.binning = Binning::fromFile(s.binningFile);
  } else {
    if (!bn["variables"]) throw std::runtime_error("config: sample '" + s.name + "' binning needs 'file' or 'variables'");
    std::vector<BinningAxis> axes;
    for (const auto& v : bn["variables"]) {
      const std::string var = v.as<std::string>();
      if (bn["edges"] && bn["edges"][var]) {
        axes.push_back(parseAxis(var, bn["edges"][var], where + ".binning.edges"));
      } else if (bn["values"] && bn["values"][var]) {
        BinningAxis a;
        a.variable = var;
        for (const auto& e : bn["values"][var]) a.values.push_back(e.as<double>());
        axes.push_back(a);
      } else {
        throw std::runtime_error("config: sample '" + s.name + "': no edges or values for binning variable '" + var + "'");
      }
    }
    s.binning = Binning::product(axes);
  }
  return s;
}

YAML::Node emitSample(const SampleDefinition& s) {
  YAML::Node n;
  n["name"] = s.name;
  n["title"] = s.title;
  n["beam"] = s.beamMode < 0 ? std::string("any") : toString(static_cast<BeamMode>(s.beamMode));
  if (!s.column.empty()) n["column"] = s.column;
  n["hypothesis"] = toString(s.hypothesis);
  YAML::Node b;
  if (!s.binningFile.empty()) {
    b["file"] = s.binningFile;
  } else if (!s.binning.axes().empty()) {
    YAML::Node vars, edges, values;
    for (const auto& a : s.binning.axes()) {
      vars.push_back(a.variable);
      YAML::Node list;
      if (a.isRange()) {
        for (double e : a.edges) list.push_back(e);
        list.SetStyle(YAML::EmitterStyle::Flow);
        edges[a.variable] = list;
      } else {
        for (double v : a.values) list.push_back(v);
        list.SetStyle(YAML::EmitterStyle::Flow);
        values[a.variable] = list;
      }
    }
    vars.SetStyle(YAML::EmitterStyle::Flow);
    b["variables"] = vars;
    if (edges.size() > 0) b["edges"] = edges;
    if (values.size() > 0) b["values"] = values;
  }
  n["binning"] = b;
  if (!s.cut.empty()) n["cut"] = s.cut;
  return n;
}

std::vector<std::string> stringList(const YAML::Node& n, const std::string& where) {
  std::vector<std::string> out;
  if (!n) return out;
  if (n.IsScalar()) {
    out.push_back(n.as<std::string>());
    return out;
  }
  if (!n.IsSequence()) throw std::runtime_error("config: '" + where + "' must be a list of names");
  for (const auto& e : n) out.push_back(e.as<std::string>());
  return out;
}

SystematicsConfig parseSystematics(const YAML::Node& n, const std::string& baseDir, const SystematicsConfig& base) {
  rejectUnknown(n, {"xsec", "flux", "detector", "oscillation", "grid"}, "systematics");
  SystematicsConfig c = base;
  if (n["xsec"]) {
    const YAML::Node x = n["xsec"];
    rejectUnknown(x, {"dials", "exclude", "write_unused_dials", "enabled"}, "systematics.xsec");
    if (x["dials"]) {
      const Block b = loadBlock(x["dials"], baseDir, "systematics.xsec.dials");
      c.xsec.dials = DialCatalog::fromYamlString(b.text);
      c.xsec.dialsInclude = b.includePath;
    }
    c.xsec.exclude = stringList(x["exclude"], "systematics.xsec.exclude");
    c.xsec.writeUnusedDials = get<bool>(x, "write_unused_dials", c.xsec.writeUnusedDials);
    c.xsec.enabled = get<bool>(x, "enabled", c.xsec.enabled);
  }
  if (n["flux"]) {
    rejectUnknown(n["flux"], {"enabled"}, "systematics.flux");
    c.flux.enabled = get<bool>(n["flux"], "enabled", c.flux.enabled);
  }
  if (n["detector"]) {
    const YAML::Node d = n["detector"];
    rejectUnknown(d, {"mode", "sigmas", "names", "correlation", "pairs", "range", "enabled"}, "systematics.detector");
    c.detector.mode = get<std::string>(d, "mode", c.detector.mode);
    if (c.detector.mode != "normalisation" && c.detector.mode != "shifts") {
      throw std::runtime_error("config: systematics.detector.mode must be 'normalisation' or 'shifts'");
    }
    if (d["sigmas"]) {
      c.detector.norms.clear();
      for (const auto& kv : d["sigmas"]) {
        DetectorNorm nm;
        nm.sample = kv.first.as<std::string>();
        nm.sigma = kv.second.as<double>();
        c.detector.norms.push_back(nm);
      }
    }
    if (d["names"]) {
      for (const auto& kv : d["names"]) {
        bool found = false;
        for (auto& nm : c.detector.norms) {
          if (nm.sample == kv.first.as<std::string>()) {
            nm.name = kv.second.as<std::string>();
            found = true;
          }
        }
        if (!found) {
          throw std::runtime_error("config: systematics.detector.names: '" + kv.first.as<std::string>() + "' has no sigma");
        }
      }
    }
    c.detector.correlation = get<double>(d, "correlation", c.detector.correlation);
    if (d["pairs"]) {
      c.detector.pairs.clear();
      for (const auto& p : d["pairs"]) c.detector.pairs.push_back(stringList(p, "systematics.detector.pairs[]"));
    }
    if (d["range"]) {
      const std::vector<std::string> r = stringList(d["range"], "systematics.detector.range");
      if (r.size() != 2) throw std::runtime_error("config: systematics.detector.range needs [low, high]");
      c.detector.rangeLow = std::stod(r[0]);
      c.detector.rangeHigh = std::stod(r[1]);
    }
    c.detector.enabled = get<bool>(d, "enabled", c.detector.enabled);
  }
  if (n["oscillation"]) {
    const YAML::Node o = n["oscillation"];
    rejectUnknown(o, {"sigmas", "fixed", "enabled"}, "systematics.oscillation");
    if (o["sigmas"]) {
      const YAML::Node s = o["sigmas"];
      rejectUnknown(s, {"s12sq", "s13sq", "s23sq", "dmsq21", "dmsq32", "delta_cp", "ordering"}, "systematics.oscillation.sigmas");
      c.oscillation.sigmas.s12sq = get<double>(s, "s12sq", c.oscillation.sigmas.s12sq);
      c.oscillation.sigmas.s13sq = get<double>(s, "s13sq", c.oscillation.sigmas.s13sq);
      c.oscillation.sigmas.s23sq = get<double>(s, "s23sq", c.oscillation.sigmas.s23sq);
      c.oscillation.sigmas.dmsq21 = get<double>(s, "dmsq21", c.oscillation.sigmas.dmsq21);
      c.oscillation.sigmas.dmsq32 = get<double>(s, "dmsq32", c.oscillation.sigmas.dmsq32);
      c.oscillation.sigmas.deltaCP = get<double>(s, "delta_cp", c.oscillation.sigmas.deltaCP);
      c.oscillation.sigmas.ordering = get<double>(s, "ordering", c.oscillation.sigmas.ordering);
    }
    if (o["fixed"]) c.oscillation.fixed = stringList(o["fixed"], "systematics.oscillation.fixed");
    c.oscillation.enabled = get<bool>(o, "enabled", c.oscillation.enabled);
  }
  if (n["grid"]) {
    rejectUnknown(n["grid"], {"bins", "min", "max"}, "systematics.grid");
    c.grid.bins = get<int>(n["grid"], "bins", c.grid.bins);
    c.grid.min = get<double>(n["grid"], "min", c.grid.min);
    c.grid.max = get<double>(n["grid"], "max", c.grid.max);
  }
  return c;
}

ShiftParameter parseShift(const YAML::Node& n, const std::string& where) {
  rejectUnknown(n, {"name", "target", "kind", "sigma", "apply_to"}, where);
  ShiftParameter s;
  s.name = require<std::string>(n, "name", where);
  s.target = shiftTargetFromString(require<std::string>(n, "target", where));
  if (n["kind"]) s.kind = shiftKindFromString(n["kind"].as<std::string>());
  s.sigma = require<double>(n, "sigma", where);
  if (n["apply_to"]) s.applyTo = shiftApplyToFromString(n["apply_to"].as<std::string>());
  return s;
}

YAML::Node emitSystematics(const SystematicsConfig& c) {
  YAML::Node n;
  {
    YAML::Node x;
    if (!c.xsec.dialsInclude.empty()) {
      YAML::Node inc;
      inc["include"] = c.xsec.dialsInclude;
      x["dials"] = inc;
    } else if (!c.xsec.dials.empty()) {
      x["dials"] = YAML::Load(c.xsec.dials.toYaml());
    }
    if (!c.xsec.exclude.empty()) {
      YAML::Node ex;
      for (const auto& e : c.xsec.exclude) ex.push_back(e);
      ex.SetStyle(YAML::EmitterStyle::Flow);
      x["exclude"] = ex;
    }
    x["write_unused_dials"] = c.xsec.writeUnusedDials;
    if (!c.xsec.enabled) x["enabled"] = false;
    n["xsec"] = x;
  }
  {
    YAML::Node f;
    f["enabled"] = c.flux.enabled;
    f.SetStyle(YAML::EmitterStyle::Flow);
    n["flux"] = f;
  }
  {
    YAML::Node d;
    d["mode"] = c.detector.mode;
    YAML::Node sig, names;
    for (const auto& nm : c.detector.norms) {
      sig[nm.sample] = nm.sigma;
      if (!nm.name.empty()) names[nm.sample] = nm.name;
    }
    sig.SetStyle(YAML::EmitterStyle::Flow);
    if (!c.detector.norms.empty()) d["sigmas"] = sig;
    if (names.size() > 0) {
      names.SetStyle(YAML::EmitterStyle::Flow);
      d["names"] = names;
    }
    d["correlation"] = c.detector.correlation;
    if (!c.detector.pairs.empty()) {
      YAML::Node pairs;
      for (const auto& p : c.detector.pairs) {
        YAML::Node pair;
        for (const auto& s : p) pair.push_back(s);
        pair.SetStyle(YAML::EmitterStyle::Flow);
        pairs.push_back(pair);
      }
      pairs.SetStyle(YAML::EmitterStyle::Flow);
      d["pairs"] = pairs;
    }
    YAML::Node range;
    range.push_back(c.detector.rangeLow);
    range.push_back(c.detector.rangeHigh);
    range.SetStyle(YAML::EmitterStyle::Flow);
    d["range"] = range;
    if (!c.detector.enabled) d["enabled"] = false;
    n["detector"] = d;
  }
  {
    YAML::Node o, s;
    s["s12sq"] = c.oscillation.sigmas.s12sq;
    s["s13sq"] = c.oscillation.sigmas.s13sq;
    s["s23sq"] = c.oscillation.sigmas.s23sq;
    s["dmsq21"] = c.oscillation.sigmas.dmsq21;
    s["dmsq32"] = c.oscillation.sigmas.dmsq32;
    s["delta_cp"] = c.oscillation.sigmas.deltaCP;
    s["ordering"] = c.oscillation.sigmas.ordering;
    s.SetStyle(YAML::EmitterStyle::Flow);
    o["sigmas"] = s;
    YAML::Node fixed;
    for (const auto& f : c.oscillation.fixed) fixed.push_back(f);
    fixed.SetStyle(YAML::EmitterStyle::Flow);
    o["fixed"] = fixed;
    if (!c.oscillation.enabled) o["enabled"] = false;
    n["oscillation"] = o;
  }
  {
    YAML::Node g;
    g["bins"] = c.grid.bins;
    g["min"] = c.grid.min;
    g["max"] = c.grid.max;
    g.SetStyle(YAML::EmitterStyle::Flow);
    n["grid"] = g;
  }
  return n;
}

YAML::Node emitShift(const ShiftParameter& s) {
  YAML::Node n;
  n["name"] = s.name;
  n["target"] = toString(s.target);
  n["kind"] = toString(s.kind);
  n["sigma"] = s.sigma;
  n["apply_to"] = toString(s.applyTo);
  n.SetStyle(YAML::EmitterStyle::Flow);
  return n;
}

YAML::Node emitFlux(const FluxSpec& f) {
  YAML::Node n;
  n["file"] = f.file;
  YAML::Node h;
  for (const auto& kv : f.histograms) h[NuSpecies::fromPdg(kv.first).name()] = kv.second;
  n["histograms"] = h;
  n["units"] = toString(f.units);
  n["pot_norm"] = f.potNorm;
  n["bin_norm_GeV"] = f.binNormGeV;
  return n;
}

YAML::Node emitDetector(const Detector& d) {
  YAML::Node n;
  n["name"] = d.name();
  n["fiducial_mass_kt"] = d.fiducialMass_kt();
  YAML::Node comp;
  for (const auto& t : d.composition()) {
    YAML::Node e;
    e["Z"] = t.Z;
    e["A"] = t.A;
    e["mass_fraction"] = t.massFraction;
    e.SetStyle(YAML::EmitterStyle::Flow);
    comp.push_back(e);
  }
  n["composition"] = comp;
  if (d.nucleonsOverride() > 0.0) n["n_nucleons"] = d.nucleonsOverride();
  return n;
}

}  // namespace

Experiment Experiment::fromYamlString(const std::string& yaml, const std::string& baseDir) {
  YAML::Node root;
  try {
    root = YAML::Load(yaml);
  } catch (const YAML::Exception& e) {
    throw std::runtime_error(std::string("config: YAML parse error: ") + e.what());
  }
  if (!root.IsMap()) throw std::runtime_error("config: the top level must be a map");

  std::vector<std::string> allowed = experimentKnownKeys();
  allowed.insert(allowed.end(), experimentRawKeys().begin(), experimentRawKeys().end());
  rejectUnknown(root, allowed, "top level");

  Experiment exp(get<std::string>(root, "name", "experiment"));
  Impl& impl = *exp.impl_;

  if (root["baseline_km"]) exp.setBaseline(root["baseline_km"].as<double>());
  if (root["matter"]) {
    rejectUnknown(root["matter"], {"density_g_cm3", "electron_fraction"}, "matter");
    exp.setMatter(get<double>(root["matter"], "density_g_cm3", 2.6), get<double>(root["matter"], "electron_fraction", 0.5));
  }
  if (root["detector"]) exp.setDetector(parseDetector(root["detector"]));
  if (root["flux_units"]) {
    rejectUnknown(root["flux_units"], {"pot_norm", "bin_norm_GeV"}, "flux_units");
    exp.setFluxUnits(get<double>(root["flux_units"], "pot_norm", 1.0e21), get<double>(root["flux_units"], "bin_norm_GeV", 0.05));
  }
  if (root["beam_modes"]) {
    for (const auto& kv : root["beam_modes"]) {
      const BeamMode mode = beamModeFromString(kv.first.as<std::string>());
      const std::string where = "beam_modes." + kv.first.as<std::string>();
      rejectUnknown(kv.second, {"pot", "flux"}, where);
      const YAML::Node fluxNode = kv.second["flux"];
      if (!fluxNode) throw std::runtime_error("config: '" + where + "' needs a 'flux' block");
      exp.setFlux(mode, parseFlux(fluxNode, baseDir, impl.defaultPotNorm, impl.defaultBinNormGeV, where + ".flux"),
                  get<double>(kv.second, "pot", 0.0));
    }
  }
  if (root["flux_covariance"]) {
    const YAML::Node c = root["flux_covariance"];
    rejectUnknown(c, {"file", "matrix", "far_detector_prefix", "near_detector_prefix", "nd_constraint", "units"},
                  "flux_covariance");
    FluxCovarianceSpec s;
    s.file = resolvePath(require<std::string>(c, "file", "flux_covariance"), baseDir);
    s.matrix = get<std::string>(c, "matrix", "total_flux_cov");
    s.farPrefix = get<std::string>(c, "far_detector_prefix", "sk");
    s.nearPrefix = get<std::string>(c, "near_detector_prefix", "nd280");
    s.ndConstraint = get<double>(c, "nd_constraint", 0.0);
    if (c["units"]) s.units = energyUnitFromString(c["units"].as<std::string>());
    exp.setFluxCovariance(s);
  }
  if (root["productions"]) {
    for (const auto& s : root["productions"]) {
      rejectUnknown(s, {"beam", "flux", "interaction", "n", "seed"}, "productions[]");
      Production m;
      m.beamMode = beamModeFromString(require<std::string>(s, "beam", "productions[]"));
      m.fluxPdg = NuSpecies::fromName(require<std::string>(s, "flux", "productions[]")).pdg();
      m.interactionPdg = s["interaction"] ? NuSpecies::fromName(s["interaction"].as<std::string>()).pdg() : m.fluxPdg;
      m.nEvents = require<long>(s, "n", "productions[]");
      m.seed = get<unsigned>(s, "seed", 0u);
      exp.addProduction(m);
    }
  }
  if (root["response"]) {
    const Block b = loadBlock(root["response"], baseDir, "response");
    if (root["response"].IsMap() && !root["response"]["include"]) {
      // inline: assignment may sit next to the tables
    }
    impl.response = Response::fromYamlString(b.text, b.baseDir);
    impl.hasResponse = true;
    impl.responseInclude = b.includePath;
    const YAML::Node outer = root["response"];
    if (outer.IsMap() && outer["assignment"]) exp.setAssignment(assignmentFromString(outer["assignment"].as<std::string>()));
    else if (b.node.IsMap() && b.node["assignment"]) exp.setAssignment(assignmentFromString(b.node["assignment"].as<std::string>()));
  }
  if (root["smearing"]) {
    const Block b = loadBlock(root["smearing"], baseDir, "smearing");
    impl.smearing = Smearing::fromYamlString(b.text, b.baseDir);
    impl.smearingInclude = b.includePath;
  }
  if (root["energy"]) exp.setEnergyEstimator(energyEstimatorFromString(root["energy"].as<std::string>()));
  if (root["binding_energy"]) exp.setBindingEnergy(root["binding_energy"].as<double>());
  if (root["samples"]) {
    const Block b = loadBlock(root["samples"], baseDir, "samples");
    YAML::Node list = b.node;
    if (list.IsMap() && list["samples"]) list = list["samples"];
    if (!list.IsSequence()) throw std::runtime_error("config: 'samples' must be a list (or {include: file} holding one)");
    int i = 0;
    for (const auto& s : list) exp.addSample(parseSample(s, b.baseDir, "samples[" + std::to_string(i++) + "]"));
    impl.samplesInclude = b.includePath;
  }
  if (root["oscillation"]) {
    const YAML::Node o = root["oscillation"];
    rejectUnknown(o, {"preset", "s12sq", "s13sq", "s23sq", "delta_cp", "dmsq21", "dmsq32", "ordering"}, "oscillation");
    OscillationParameters p;
    if (o["preset"]) {
      const std::string preset = o["preset"].as<std::string>();
      if (preset == "asimov_A" || preset == "asimovA" || preset == "t2k_asimov_A") p = OscillationParameters::asimovA();
      else throw std::runtime_error("config: unknown oscillation preset '" + preset + "' (asimov_A)");
    }
    p.s12sq = get<double>(o, "s12sq", p.s12sq);
    p.s13sq = get<double>(o, "s13sq", p.s13sq);
    p.s23sq = get<double>(o, "s23sq", p.s23sq);
    p.deltaCP = get<double>(o, "delta_cp", p.deltaCP);
    p.dmsq21 = get<double>(o, "dmsq21", p.dmsq21);
    p.dmsq32 = get<double>(o, "dmsq32", p.dmsq32);
    if (o["ordering"]) {
      const std::string ord = o["ordering"].as<std::string>();
      if (ord == "normal" || ord == "NO") p.normalOrdering = true;
      else if (ord == "inverted" || ord == "IO") p.normalOrdering = false;
      else throw std::runtime_error("config: oscillation.ordering must be normal or inverted");
    }
    exp.setOscillation(p);
  }
  if (root["systematics"]) impl.systematics = parseSystematics(root["systematics"], baseDir, impl.systematics);
  if (root["shifts"]) {
    if (!root["shifts"].IsSequence()) throw std::runtime_error("config: 'shifts' must be a list");
    int i = 0;
    for (const auto& s : root["shifts"]) exp.addShift(parseShift(s, "shifts[" + std::to_string(i++) + "]"));
  }
  if (root["run"]) {
    const YAML::Node r = root["run"];
    rejectUnknown(r, {"work_dir", "seed", "legacy_rng_order"}, "run");
    exp.setRun(resolvePath(get<std::string>(r, "work_dir", ""), baseDir), get<unsigned>(r, "seed", 20250929u));
    exp.setLegacyRngOrder(get<bool>(r, "legacy_rng_order", false));
  }
  for (const auto& key : experimentRawKeys()) {
    if (root[key]) impl.raw[key] = YAML::Dump(root[key]);
  }
  return exp;
}

Experiment Experiment::fromConfig(const std::string& yamlPath) {
  std::ifstream in(yamlPath);
  if (!in) throw std::runtime_error("Experiment::fromConfig: cannot read '" + yamlPath + "'");
  std::stringstream ss;
  ss << in.rdbuf();
  std::string dir = gSystem->GetDirName(yamlPath.c_str()).Data();
  if (dir.empty() || dir == ".") dir = gSystem->WorkingDirectory();
  if (dir[0] != '/') dir = std::string(gSystem->WorkingDirectory()) + "/" + dir;
  Experiment exp = fromYamlString(ss.str(), dir);
  exp.impl_->configPath = yamlPath;
  return exp;
}

std::string Experiment::toYaml() const {
  const Impl& impl = *impl_;
  YAML::Node root;
  root["name"] = impl.name;
  if (impl.baseline_km > 0.0) root["baseline_km"] = impl.baseline_km;
  {
    YAML::Node m;
    m["density_g_cm3"] = impl.density_g_cm3;
    m["electron_fraction"] = impl.electronFraction;
    m.SetStyle(YAML::EmitterStyle::Flow);
    root["matter"] = m;
  }
  if (impl.hasDetector) root["detector"] = emitDetector(impl.detector);
  {
    YAML::Node u;
    u["pot_norm"] = impl.defaultPotNorm;
    u["bin_norm_GeV"] = impl.defaultBinNormGeV;
    u.SetStyle(YAML::EmitterStyle::Flow);
    root["flux_units"] = u;
  }
  if (!impl.beams.empty()) {
    YAML::Node beams;
    for (const auto& kv : impl.beams) {
      if (!kv.second.hasFlux) continue;
      YAML::Node b;
      b["pot"] = kv.second.pot;
      b["flux"] = emitFlux(kv.second.flux);
      beams[toString(kv.first)] = b;
    }
    root["beam_modes"] = beams;
  }
  if (impl.fluxCov.valid()) {
    YAML::Node c;
    c["file"] = impl.fluxCov.file;
    c["matrix"] = impl.fluxCov.matrix;
    c["far_detector_prefix"] = impl.fluxCov.farPrefix;
    c["near_detector_prefix"] = impl.fluxCov.nearPrefix;
    c["nd_constraint"] = impl.fluxCov.ndConstraint;
    c["units"] = toString(impl.fluxCov.units);
    root["flux_covariance"] = c;
  }
  if (!impl.productions.empty()) {
    YAML::Node list;
    for (const auto& s : impl.productions) {
      YAML::Node m;
      m["beam"] = toString(s.beamMode);
      m["flux"] = s.fluxSpecies().name();
      m["interaction"] = s.interactionSpecies().name();
      m["n"] = static_cast<long long>(s.nEvents);
      if (s.seed != 0) m["seed"] = s.seed;
      m.SetStyle(YAML::EmitterStyle::Flow);
      list.push_back(m);
    }
    root["productions"] = list;
  }
  if (impl.hasResponse) {
    YAML::Node r;
    if (!impl.responseInclude.empty()) {
      r["include"] = impl.responseInclude;
    } else {
      r = YAML::Load(impl.response.toYaml());
    }
    r["assignment"] = toString(impl.assignment);
    root["response"] = r;
  }
  {
    if (!impl.smearingInclude.empty()) {
      YAML::Node s;
      s["include"] = impl.smearingInclude;
      root["smearing"] = s;
    } else {
      root["smearing"] = YAML::Load(impl.smearing.toYaml());
    }
  }
  root["energy"] = toString(impl.estimator);
  root["binding_energy"] = impl.bindingEnergy;
  if (!impl.samples.empty()) {
    if (!impl.samplesInclude.empty()) {
      YAML::Node s;
      s["include"] = impl.samplesInclude;
      root["samples"] = s;
    } else {
      YAML::Node list;
      for (const auto& s : impl.samples) list.push_back(emitSample(s));
      root["samples"] = list;
    }
  }
  {
    YAML::Node o;
    o["s12sq"] = impl.oscillation.s12sq;
    o["s13sq"] = impl.oscillation.s13sq;
    o["s23sq"] = impl.oscillation.s23sq;
    o["delta_cp"] = impl.oscillation.deltaCP;
    o["dmsq21"] = impl.oscillation.dmsq21;
    o["dmsq32"] = impl.oscillation.dmsq32;
    o["ordering"] = impl.oscillation.normalOrdering ? "normal" : "inverted";
    root["oscillation"] = o;
  }
  root["systematics"] = emitSystematics(impl.systematics);
  if (!impl.shifts.empty()) {
    YAML::Node list;
    for (const auto& s : impl.shifts) list.push_back(emitShift(s));
    root["shifts"] = list;
  }
  {
    YAML::Node r;
    r["work_dir"] = impl.run.workDir;
    r["seed"] = impl.run.seed;
    if (impl.run.legacyRngOrder) r["legacy_rng_order"] = true;
    root["run"] = r;
  }
  for (const auto& key : experimentRawKeys()) {
    auto it = impl.raw.find(key);
    if (it != impl.raw.end()) root[key] = YAML::Load(it->second);
  }
  std::ostringstream os;
  os << "# Zaku experiment configuration (see docs/config_reference.md)\n" << YAML::Dump(root) << "\n";
  return os.str();
}

void Experiment::saveConfig(const std::string& path) const {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("Experiment::saveConfig: cannot write '" + path + "'");
  out << toYaml();
}

}  // namespace zaku
