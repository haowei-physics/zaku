#include "zaku/syst/Splines.h"

#include <TGraph.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace zaku {

std::string toString(SplineDecision d) {
  switch (d) {
    case SplineDecision::Full: return "full";
    case SplineDecision::Flat: return "flat";
    case SplineDecision::Rejected: return "rejected";
  }
  return "?";
}

bool SplineDial::nominalInserted() const {
  for (int s : src) {
    if (s < 0) return true;
  }
  return false;
}

bool SplineDial::uniformKnots(double tol) const {
  if (x.size() < 3) return true;
  const double step = x[1] - x[0];
  if (step <= 0.0) return false;
  for (std::size_t k = 2; k < x.size(); ++k) {
    if (std::fabs((x[k] - x[k - 1]) - step) > tol * std::fabs(step)) return false;
  }
  return true;
}

std::string SplineDial::describe() const {
  std::ostringstream os;
  os << std::left << std::setw(24) << name << " nominal " << std::setw(8) << nominal << " " << x.size() << " knots:";
  for (std::size_t k = 0; k < x.size(); ++k) os << " " << x[k] << (src[k] < 0 ? "*" : "");
  os << (uniformKnots() ? "" : "  (non-uniform)");
  return os.str();
}

SplineBuilder::SplineBuilder(const std::vector<DialKnots>& knots) { setup(knots, nullptr); }

SplineBuilder::SplineBuilder(const std::vector<DialKnots>& knots, const std::vector<std::string>& keep) {
  setup(knots, &keep);
}

void SplineBuilder::setup(const std::vector<DialKnots>& knots, const std::vector<std::string>* keep) {
  std::vector<int> chosen;
  if (keep) {
    for (const auto& name : *keep) {
      int found = -1;
      for (std::size_t i = 0; i < knots.size(); ++i) {
        if (knots[i].name == name) {
          found = static_cast<int>(i);
          break;
        }
      }
      if (found < 0) throw std::runtime_error("SplineBuilder: the weight source has no dial '" + name + "'");
      chosen.push_back(found);
    }
  } else {
    for (std::size_t i = 0; i < knots.size(); ++i) chosen.push_back(static_cast<int>(i));
  }
  for (int i : chosen) {
    const DialKnots& k = knots[static_cast<std::size_t>(i)];
    if (k.knots.empty()) throw std::runtime_error("SplineBuilder: dial '" + k.name + "' has no tweak points");
    SplineDial d;
    d.name = k.name;
    d.nominal = k.nominal;
    std::vector<std::pair<double, int>> pt;
    for (std::size_t j = 0; j < k.knots.size(); ++j) pt.emplace_back(k.knots[j], static_cast<int>(j));
    if (!k.nominalIsKnot()) pt.emplace_back(k.nominal, -1);
    std::stable_sort(pt.begin(), pt.end(),
                     [](const std::pair<double, int>& a, const std::pair<double, int>& b) { return a.first < b.first; });
    for (const auto& p : pt) {
      d.x.push_back(p.first);
      d.src.push_back(p.second);
    }
    knots_.push_back(k);
    dials_.push_back(d);
    sourceIndex_.push_back(i);
  }
}

std::vector<std::string> SplineBuilder::dialNames() const {
  std::vector<std::string> out;
  for (const auto& d : dials_) out.push_back(d.name);
  return out;
}

int SplineBuilder::index(const std::string& name) const {
  for (std::size_t i = 0; i < dials_.size(); ++i) {
    if (dials_[i].name == name) return static_cast<int>(i);
  }
  return -1;
}

SplineDecision SplineBuilder::classify(const SplineDial& d, const std::vector<double>& w, double* nominalDeviation) {
  bool allOne = true, allZero = true, bad = false;
  double worst = 0.0;
  for (std::size_t k = 0; k < d.x.size(); ++k) {
    if (d.src[k] < 0) continue;  // the inserted nominal
    const std::size_t s = static_cast<std::size_t>(d.src[k]);
    if (s >= w.size()) throw std::runtime_error("SplineBuilder: dial '" + d.name + "' has fewer weights than knots");
    const double y = w[s];
    if (y != 1.0) allOne = false;
    if (y != 0.0) allZero = false;
    if (!std::isfinite(y) || y < 0.0) bad = true;
    if (std::fabs(d.x[k] - d.nominal) < 1e-9 * (1.0 + std::fabs(d.nominal))) {
      const double dev = std::fabs(y - 1.0);
      worst = std::max(worst, dev);
      if (dev > 1e-6) bad = true;
    }
  }
  if (nominalDeviation) *nominalDeviation = worst;
  if (allZero || bad) return SplineDecision::Rejected;
  if (allOne) return SplineDecision::Flat;
  return SplineDecision::Full;
}

SplineDecision SplineBuilder::build(int i, const std::vector<double>& w, TGraph& g) {
  SplineDial& d = dials_.at(static_cast<std::size_t>(i));
  double dev = 0.0;
  const SplineDecision decision = classify(d, w, &dev);
  d.worstNominalDeviation = std::max(d.worstNominalDeviation, dev);
  if (decision != SplineDecision::Full) {
    g.Set(0);
    if (decision == SplineDecision::Flat) ++d.nFlat;
    else ++d.nRejected;
    return decision;
  }
  g.Set(d.nKnots());
  for (int k = 0; k < d.nKnots(); ++k) {
    const int s = d.src[static_cast<std::size_t>(k)];
    g.SetPoint(k, d.x[static_cast<std::size_t>(k)], s < 0 ? 1.0 : w[static_cast<std::size_t>(s)]);
  }
  ++d.nFull;
  return decision;
}

std::string SplineBuilder::statsText(long long entries) const {
  std::ostringstream os;
  os << std::left << std::setw(24) << "dial" << std::right << std::setw(10) << "full" << std::setw(12) << "no response"
     << std::setw(10) << "rejected" << "\n";
  for (const auto& d : dials_) {
    os << std::left << std::setw(24) << d.name << std::right << std::setw(10) << d.nFull << std::setw(12) << d.nFlat;
    if (entries > 0) os << " (" << std::fixed << std::setprecision(1) << 100.0 * d.nFlat / entries << "%)" << std::defaultfloat;
    os << std::setw(10) << d.nRejected << (d.nRejected ? "  <--" : "") << "\n";
  }
  return os.str();
}

long long SplineBuilder::totalRejected() const {
  long long n = 0;
  for (const auto& d : dials_) n += d.nRejected;
  return n;
}

std::string SplineBuilder::metadataText(const std::vector<DialKnots>& knots) {
  std::ostringstream os;
  os << std::setprecision(17);
  for (const auto& k : knots) {
    os << k.name << "=" << k.nominal << ":";
    for (std::size_t i = 0; i < k.knots.size(); ++i) os << (i ? "," : "") << k.knots[i];
    os << "\n";
  }
  return os.str();
}

std::vector<DialKnots> SplineBuilder::fromMetadataText(const std::string& text) {
  std::vector<DialKnots> out;
  std::istringstream in(text);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    const std::string::size_type eq = line.find('=');
    const std::string::size_type colon = line.find(':', eq == std::string::npos ? 0 : eq);
    if (eq == std::string::npos || colon == std::string::npos) {
      throw std::runtime_error("SplineBuilder::fromMetadataText: cannot parse '" + line + "'");
    }
    DialKnots k;
    k.name = line.substr(0, eq);
    k.nominal = std::stod(line.substr(eq + 1, colon - eq - 1));
    std::istringstream ks(line.substr(colon + 1));
    std::string item;
    while (std::getline(ks, item, ',')) {
      if (!item.empty()) k.knots.push_back(std::stod(item));
    }
    out.push_back(k);
  }
  return out;
}

}  // namespace zaku
