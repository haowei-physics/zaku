#include "zaku/detector/Smearing.h"

#include <TAxis.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TMath.h>
#include <TSystem.h>
#include <yaml-cpp/yaml.h>

#include <cmath>
#include <cstdlib>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace zaku {

// --- ResolutionModel ---------------------------------------------------------

ResolutionModel ResolutionModel::none() { return ResolutionModel(); }

ResolutionModel ResolutionModel::fractional(double fraction) {
  ResolutionModel m;
  m.kind_ = Kind::Fractional;
  m.a_ = fraction;
  return m;
}

ResolutionModel ResolutionModel::absolute(double sigma) {
  ResolutionModel m;
  m.kind_ = Kind::Absolute;
  m.a_ = sigma;
  return m;
}

ResolutionModel ResolutionModel::histogram(const std::string& file, const std::string& name, EnergyUnit units) {
  std::unique_ptr<TFile> f(TFile::Open(file.c_str(), "READ"));
  if (!f || f->IsZombie()) throw std::runtime_error("ResolutionModel: cannot open '" + file + "'");
  TH1* h = f->Get<TH1>(name.c_str());
  if (!h) throw std::runtime_error("ResolutionModel: no histogram '" + name + "' in '" + file + "'");
  if (h->GetDimension() != 1) throw std::runtime_error("ResolutionModel: '" + name + "' must be one-dimensional");
  // Axis in GeV; the contents are absolute widths, converted along with the axis.
  TH1* c = histogramToGeV(*h, units, /*contentIsDensity=*/false, 1, name);
  if (units == EnergyUnit::MeV) c->Scale(units::MeV);
  return histogram(c, file + ":" + name);
}

ResolutionModel ResolutionModel::histogram(TH1* h, const std::string& label) {
  if (!h) throw std::invalid_argument("ResolutionModel::histogram: null histogram");
  h->SetDirectory(nullptr);
  ResolutionModel m;
  m.kind_ = Kind::Histogram;
  m.hist_.reset(h);
  m.label_ = label;
  return m;
}

ResolutionModel ResolutionModel::stochastic(double a, double b) {
  ResolutionModel m;
  m.kind_ = Kind::Stochastic;
  m.a_ = a;
  m.b_ = b;
  return m;
}

ResolutionModel ResolutionModel::formula(const std::string& expression) {
  ResolutionModel m;
  m.kind_ = Kind::Formula;
  m.expression_ = expression;
  return m;
}

void ResolutionModel::prepare(const VariableSet& vars, const std::string& where) {
  if (kind_ != Kind::Formula) return;
  std::vector<std::string> names = vars.names();
  valueIndex_ = static_cast<int>(names.size());
  names.push_back("value");
  try {
    formula_ = Formula(expression_, names);
  } catch (const std::exception& e) {
    throw std::runtime_error(where + ": " + e.what());
  }
}

double ResolutionModel::sigma(double value, const VariableSet& vars) const {
  switch (kind_) {
    case Kind::None:
      return 0.0;
    case Kind::Fractional:
      return a_ * value;
    case Kind::Absolute:
      return a_;
    case Kind::Histogram: {
      // Legacy lookup: FindBin clamped into [1, nbins].
      const TAxis* ax = hist_->GetXaxis();
      int b = ax->FindFixBin(value);
      if (b < 1) b = 1;
      if (b > ax->GetNbins()) b = ax->GetNbins();
      return hist_->GetBinContent(b);
    }
    case Kind::Stochastic: {
      if (value <= 0.0) return 0.0;
      const double s = a_ * std::sqrt(value);
      const double c = b_ * value;
      return std::sqrt(s * s + c * c);
    }
    case Kind::Formula: {
      std::vector<double> v(vars.values());
      v.push_back(value);
      return formula_.eval(v.data());
    }
  }
  return 0.0;
}

std::string ResolutionModel::describe() const {
  std::ostringstream os;
  switch (kind_) {
    case Kind::None: os << "none"; break;
    case Kind::Fractional: os << "fractional " << a_; break;
    case Kind::Absolute: os << "absolute " << a_; break;
    case Kind::Histogram: os << "hist(" << (label_.empty() ? hist_->GetName() : label_) << ")"; break;
    case Kind::Stochastic: os << "stochastic a=" << a_ << " b=" << b_; break;
    case Kind::Formula: os << "formula(" << expression_ << ")"; break;
  }
  return os.str();
}

// --- Smearing ----------------------------------------------------------------

Smearing::Smearing() = default;

TH1D* Smearing::legacyMomentumHistogram(double fraction, const std::string& name) {
  const bool prev = TH1::AddDirectoryStatus();
  TH1::AddDirectory(false);
  auto* h = new TH1D(name.c_str(), "sigma_p vs p_lep;p_lep (GeV);sigma_p (GeV)", 50, 0.0, 5.0);
  TH1::AddDirectory(prev);
  for (int i = 1; i <= 50; ++i) h->SetBinContent(i, fraction * h->GetBinCenter(i));
  return h;
}

TH1D* Smearing::legacyAngleHistogram(double angle_deg, const std::string& name) {
  const bool prev = TH1::AddDirectoryStatus();
  TH1::AddDirectory(false);
  auto* h = new TH1D(name.c_str(), "sigma_cos vs cos_lep;cos_lep;sigma_cos", 40, -1.0, 1.0);
  TH1::AddDirectory(prev);
  // Same operation order as the legacy MakeSmearing macro: the angle is
  // converted to radians first (TMath::DegToRad()), then scaled by sin(theta).
  const double angleRad = angle_deg * TMath::DegToRad();
  for (int i = 1; i <= 40; ++i) {
    const double c = h->GetBinCenter(i);
    h->SetBinContent(i, std::sqrt(1.0 - c * c) * angleRad);
  }
  return h;
}

Smearing Smearing::parametric(double fraction, double angle_deg) {
  Smearing s;
  s.leptonMomentum_ = ResolutionModel::histogram(legacyMomentumHistogram(fraction));
  s.leptonAngle_ = ResolutionModel::histogram(legacyAngleHistogram(angle_deg));
  return s;
}

Smearing Smearing::fromLegacyFile(const std::string& file) {
  Smearing s;
  s.leptonMomentum_ = ResolutionModel::histogram(file, "h_dp_vs_p");
  s.leptonAngle_ = ResolutionModel::histogram(file, "h_dcos_vs_cos");
  return s;
}

void Smearing::prepare(const VariableSet& vars) {
  leptonMomentum_.prepare(vars, "smearing.lepton_momentum");
  leptonAngle_.prepare(vars, "smearing.lepton_angle");
  hadronicEnergy_.prepare(vars, "smearing.hadronic_energy");
  electronEnergy_.prepare(vars, "smearing.electron_energy");
  muonEnergy_.prepare(vars, "smearing.muon_energy");
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
  return baseDir + "/" + s;
}

ResolutionModel parseModel(const YAML::Node& n, const std::string& baseDir, const std::string& where, bool isAngle) {
  if (!n) return ResolutionModel::none();
  if (n.IsScalar()) {
    // bare number: fractional (or degrees for the angle)
    const double v = n.as<double>();
    if (isAngle) return ResolutionModel::histogram(Smearing::legacyAngleHistogram(v));
    return ResolutionModel::fractional(v);
  }
  if (!n.IsMap()) throw std::runtime_error("smearing: '" + where + "' must be a map");
  if (n["hist"]) {
    const std::string ref = n["hist"].as<std::string>();
    const std::string::size_type pos = ref.rfind(':');
    if (pos == std::string::npos) throw std::runtime_error("smearing: '" + where + ".hist' must be 'file.root:histogram'");
    EnergyUnit units = n["units"] ? energyUnitFromString(n["units"].as<std::string>()) : EnergyUnit::GeV;
    return ResolutionModel::histogram(resolvePath(ref.substr(0, pos), baseDir), ref.substr(pos + 1), units);
  }
  if (n["formula"]) return ResolutionModel::formula(n["formula"].as<std::string>());
  const std::string type = n["type"] ? n["type"].as<std::string>() : std::string("fractional");
  const double sigma = n["sigma"] ? n["sigma"].as<double>() : 0.0;
  if (type == "none") return ResolutionModel::none();
  if (type == "fractional") return ResolutionModel::fractional(sigma);
  if (type == "absolute") return ResolutionModel::absolute(sigma);
  if (type == "legacy_fractional") return ResolutionModel::histogram(Smearing::legacyMomentumHistogram(sigma));
  if (type == "gaussian_deg" || type == "legacy_gaussian_deg") return ResolutionModel::histogram(Smearing::legacyAngleHistogram(sigma));
  if (type == "stochastic") {
    return ResolutionModel::stochastic(n["a"] ? n["a"].as<double>() : 0.0, n["b"] ? n["b"].as<double>() : 0.0);
  }
  throw std::runtime_error("smearing: unknown resolution type '" + type + "' in '" + where +
                           "' (none, fractional, absolute, legacy_fractional, gaussian_deg, stochastic, hist, formula)");
}

YAML::Node emitModel(const ResolutionModel& m) {
  YAML::Node n;
  switch (m.kind()) {
    case ResolutionModel::Kind::None: n["type"] = "none"; break;
    case ResolutionModel::Kind::Fractional: n["type"] = "fractional"; n["sigma"] = m.parameterA(); break;
    case ResolutionModel::Kind::Absolute: n["type"] = "absolute"; n["sigma"] = m.parameterA(); break;
    case ResolutionModel::Kind::Histogram: {
      const std::string& label = m.histogramLabel();
      if (label.find(':') != std::string::npos) {
        n["hist"] = label;
      } else {
        // In-memory legacy-style histogram: recover its parameter from the last bin.
        const TH1* h = m.histogramPtr();
        const int nb = h->GetNbinsX();
        const double c = h->GetXaxis()->GetBinCenter(nb);
        if (h->GetXaxis()->GetXmin() < 0.0) {  // the angle histogram
          n["type"] = "gaussian_deg";
          n["sigma"] = h->GetBinContent(nb) / std::sqrt(1.0 - c * c) * 180.0 / units::pi;
        } else {
          n["type"] = "legacy_fractional";
          n["sigma"] = h->GetBinContent(nb) / c;
        }
      }
      break;
    }
    case ResolutionModel::Kind::Stochastic: n["type"] = "stochastic"; n["a"] = m.parameterA(); n["b"] = m.parameterB(); break;
    case ResolutionModel::Kind::Formula: n["formula"] = m.expression(); break;
  }
  n.SetStyle(YAML::EmitterStyle::Flow);
  return n;
}

}  // namespace

Smearing Smearing::fromYamlString(const std::string& yaml, const std::string& baseDir) {
  YAML::Node root;
  try {
    root = YAML::Load(yaml);
  } catch (const YAML::Exception& e) {
    throw std::runtime_error(std::string("smearing: YAML parse error: ") + e.what());
  }
  if (!root.IsMap()) throw std::runtime_error("smearing: the block must be a map");
  for (const auto& kv : root) {
    const std::string key = kv.first.as<std::string>();
    if (key != "file" && key != "lepton_momentum" && key != "lepton_angle" && key != "hadronic_energy" &&
        key != "electron_energy" && key != "muon_energy" && key != "neutron_visible_fraction" &&
        key != "decay_electron_tag") {
      throw std::runtime_error("smearing: unknown key '" + key + "'");
    }
  }
  Smearing s;
  if (root["file"]) s = fromLegacyFile(resolvePath(root["file"].as<std::string>(), baseDir));
  if (root["lepton_momentum"]) s.leptonMomentum_ = parseModel(root["lepton_momentum"], baseDir, "lepton_momentum", false);
  if (root["lepton_angle"]) s.leptonAngle_ = parseModel(root["lepton_angle"], baseDir, "lepton_angle", true);
  if (root["hadronic_energy"]) s.hadronicEnergy_ = parseModel(root["hadronic_energy"], baseDir, "hadronic_energy", false);
  if (root["electron_energy"]) s.electronEnergy_ = parseModel(root["electron_energy"], baseDir, "electron_energy", false);
  if (root["muon_energy"]) s.muonEnergy_ = parseModel(root["muon_energy"], baseDir, "muon_energy", false);
  if (root["neutron_visible_fraction"]) s.neutronFraction_ = root["neutron_visible_fraction"].as<double>();
  if (root["decay_electron_tag"]) {
    const YAML::Node t = root["decay_electron_tag"];
    if (t["mu_plus"]) s.tagging_.muPlus = t["mu_plus"].as<double>();
    if (t["mu_minus"]) s.tagging_.muMinus = t["mu_minus"].as<double>();
    if (t["pi_plus"]) s.tagging_.piPlus = t["pi_plus"].as<double>();
  }
  return s;
}

std::string Smearing::toYaml() const {
  YAML::Node root;
  root["lepton_momentum"] = emitModel(leptonMomentum_);
  root["lepton_angle"] = emitModel(leptonAngle_);
  root["hadronic_energy"] = emitModel(hadronicEnergy_);
  root["electron_energy"] = emitModel(electronEnergy_);
  root["muon_energy"] = emitModel(muonEnergy_);
  root["neutron_visible_fraction"] = neutronFraction_;
  YAML::Node t;
  t["mu_plus"] = tagging_.muPlus;
  t["mu_minus"] = tagging_.muMinus;
  t["pi_plus"] = tagging_.piPlus;
  t.SetStyle(YAML::EmitterStyle::Flow);
  root["decay_electron_tag"] = t;
  return YAML::Dump(root);
}

std::string Smearing::describe() const {
  std::ostringstream os;
  os << "smearing: lepton momentum " << leptonMomentum_.describe() << ", lepton angle " << leptonAngle_.describe()
     << ", hadronic energy " << hadronicEnergy_.describe() << ", electron energy " << electronEnergy_.describe()
     << ", muon energy " << muonEnergy_.describe() << ", neutron visible fraction " << neutronFraction_
     << ", decay-e tagging mu+ " << tagging_.muPlus << " mu- " << tagging_.muMinus << " pi+ " << tagging_.piPlus;
  return os.str();
}

}  // namespace zaku
