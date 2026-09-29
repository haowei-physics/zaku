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
  static const std::vector<std::string> k = {"name",        "baseline_km",     "matter",     "detector",
                                             "flux_units",  "beam_modes",      "flux_covariance", "mc_samples",
                                             "run"};
  return k;
}

const std::vector<std::string>& experimentRawKeys() {
  static const std::vector<std::string> k = {"samples",  "response",    "smearing",    "energy", "generator",
                                             "shifts",   "oscillation", "systematics", "gundam"};
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
  if (root["mc_samples"]) {
    for (const auto& s : root["mc_samples"]) {
      rejectUnknown(s, {"beam", "flux", "interaction", "n", "seed"}, "mc_samples[]");
      McSample m;
      m.beamMode = beamModeFromString(require<std::string>(s, "beam", "mc_samples[]"));
      m.fluxPdg = NuSpecies::fromName(require<std::string>(s, "flux", "mc_samples[]")).pdg();
      m.interactionPdg = s["interaction"] ? NuSpecies::fromName(s["interaction"].as<std::string>()).pdg() : m.fluxPdg;
      m.nEvents = require<long>(s, "n", "mc_samples[]");
      m.seed = get<unsigned>(s, "seed", 0u);
      exp.addMcSample(m);
    }
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
  if (!impl.mcSamples.empty()) {
    YAML::Node list;
    for (const auto& s : impl.mcSamples) {
      YAML::Node m;
      m["beam"] = toString(s.beamMode);
      m["flux"] = s.fluxSpecies().name();
      m["interaction"] = s.interactionSpecies().name();
      m["n"] = static_cast<long long>(s.nEvents);
      if (s.seed != 0) m["seed"] = s.seed;
      m.SetStyle(YAML::EmitterStyle::Flow);
      list.push_back(m);
    }
    root["mc_samples"] = list;
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
