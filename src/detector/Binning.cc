#include "zaku/detector/Binning.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace zaku {

BinningAxis BinningAxis::uniform(const std::string& variable, double min, double max, double step,
                                 const std::vector<double>& extra) {
  if (step <= 0.0 || max <= min) throw std::invalid_argument("BinningAxis::uniform: invalid range or step");
  BinningAxis a;
  a.variable = variable;
  const int n = static_cast<int>(std::lround((max - min) / step));
  for (int i = 0; i <= n; ++i) a.edges.push_back(min + i * step);
  for (double e : extra) {
    if (e <= a.edges.back()) throw std::invalid_argument("BinningAxis::uniform: extra edges must increase");
    a.edges.push_back(e);
  }
  return a;
}

bool Bin::contains(const std::vector<double>& x, const std::vector<bool>& isRange) const {
  for (std::size_t i = 0; i < low.size(); ++i) {
    if (isRange[i]) {
      if (x[i] < low[i] || x[i] >= high[i]) return false;
    } else {
      if (x[i] != low[i]) return false;
    }
  }
  return true;
}

Binning Binning::product(const std::vector<BinningAxis>& axes) {
  if (axes.empty()) throw std::invalid_argument("Binning::product: no axes");
  Binning b;
  b.axes_ = axes;
  for (const auto& a : axes) {
    if (a.nBins() <= 0) throw std::invalid_argument("Binning::product: axis '" + a.variable + "' has no bins");
    b.variables_.push_back(a.variable);
    b.isRange_.push_back(a.isRange());
  }
  std::vector<int> idx(axes.size(), 0);
  while (true) {
    Bin bin;
    for (std::size_t k = 0; k < axes.size(); ++k) {
      const BinningAxis& a = axes[k];
      if (a.isRange()) {
        bin.low.push_back(a.edges[static_cast<std::size_t>(idx[k])]);
        bin.high.push_back(a.edges[static_cast<std::size_t>(idx[k]) + 1]);
      } else {
        bin.low.push_back(a.values[static_cast<std::size_t>(idx[k])]);
        bin.high.push_back(a.values[static_cast<std::size_t>(idx[k])]);
      }
    }
    b.bins_.push_back(bin);
    // advance, last axis fastest
    int k = static_cast<int>(axes.size()) - 1;
    while (k >= 0) {
      if (++idx[static_cast<std::size_t>(k)] < axes[static_cast<std::size_t>(k)].nBins()) break;
      idx[static_cast<std::size_t>(k)] = 0;
      --k;
    }
    if (k < 0) break;
  }
  return b;
}

Binning Binning::fromText(const std::string& text, const std::string& where) {
  Binning b;
  std::istringstream in(text);
  std::string line;
  int lineNo = 0;
  bool haveVariables = false;
  while (std::getline(in, line)) {
    ++lineNo;
    const std::string::size_type hash = line.find('#');
    if (hash != std::string::npos) line.erase(hash);
    std::istringstream ls(line);
    std::string first;
    if (!(ls >> first)) continue;
    if (first == "variables:") {
      std::vector<std::string> vars;
      std::string v;
      while (ls >> v) vars.push_back(v);
      if (vars.empty()) throw std::runtime_error(where + ":" + std::to_string(lineNo) + ": empty variables line");
      if (haveVariables && !b.bins_.empty()) {
        throw std::runtime_error(where + ":" + std::to_string(lineNo) + ": Zaku binnings use one variables line");
      }
      b.variables_.clear();
      b.isRange_.clear();
      for (std::size_t i = 0; i < vars.size(); ++i) {
        if (i + 1 < vars.size() && vars[i + 1] == vars[i]) {
          b.variables_.push_back(vars[i]);
          b.isRange_.push_back(true);
          ++i;
        } else {
          b.variables_.push_back(vars[i]);
          b.isRange_.push_back(false);
        }
      }
      haveVariables = true;
      continue;
    }
    if (!haveVariables) throw std::runtime_error(where + ":" + std::to_string(lineNo) + ": bin line before 'variables:'");
    std::vector<double> numbers;
    numbers.push_back(std::stod(first));
    double v;
    while (ls >> v) numbers.push_back(v);
    std::size_t expected = 0;
    for (bool r : b.isRange_) expected += r ? 2 : 1;
    if (numbers.size() != expected) {
      throw std::runtime_error(where + ":" + std::to_string(lineNo) + ": expected " + std::to_string(expected) +
                               " numbers, found " + std::to_string(numbers.size()));
    }
    Bin bin;
    std::size_t k = 0;
    for (bool r : b.isRange_) {
      if (r) {
        bin.low.push_back(numbers[k]);
        bin.high.push_back(numbers[k + 1]);
        k += 2;
      } else {
        bin.low.push_back(numbers[k]);
        bin.high.push_back(numbers[k]);
        k += 1;
      }
    }
    b.bins_.push_back(bin);
  }
  if (!haveVariables) throw std::runtime_error(where + ": no 'variables:' line");
  return b;
}

Binning Binning::fromFile(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("Binning::fromFile: cannot read '" + path + "'");
  std::stringstream ss;
  ss << in.rdbuf();
  Binning b = fromText(ss.str(), path);
  b.source_ = path;
  return b;
}

int Binning::find(const std::vector<double>& x) const {
  if (x.size() != variables_.size()) throw std::invalid_argument("Binning::find: wrong number of values");
  for (std::size_t i = 0; i < bins_.size(); ++i) {
    if (bins_[i].contains(x, isRange_)) return static_cast<int>(i);
  }
  return -1;
}

std::string Binning::toText(const std::string& comment) const {
  std::ostringstream os;
  if (!comment.empty()) {
    // every line of a multi-line comment gets its own '#'
    std::istringstream cs(comment);
    std::string line;
    while (std::getline(cs, line)) os << "# " << line << "\n";
  }
  os << "variables:";
  for (std::size_t i = 0; i < variables_.size(); ++i) {
    os << " " << variables_[i];
    if (isRange_[i]) os << " " << variables_[i];
  }
  os << "\n";
  os << std::setprecision(10);
  for (const auto& b : bins_) {
    bool first = true;
    for (std::size_t i = 0; i < variables_.size(); ++i) {
      if (!first) os << " ";
      first = false;
      os << b.low[i];
      if (isRange_[i]) os << " " << b.high[i];
    }
    os << "\n";
  }
  return os.str();
}

void Binning::writeFile(const std::string& path, const std::string& comment) const {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("Binning::writeFile: cannot write '" + path + "'");
  out << toText(comment);
}

bool Binning::sameAs(const Binning& other, double tol) const {
  if (variables_ != other.variables_ || isRange_ != other.isRange_ || bins_.size() != other.bins_.size()) return false;
  for (std::size_t i = 0; i < bins_.size(); ++i) {
    for (std::size_t k = 0; k < variables_.size(); ++k) {
      if (std::fabs(bins_[i].low[k] - other.bins_[i].low[k]) > tol) return false;
      if (std::fabs(bins_[i].high[k] - other.bins_[i].high[k]) > tol) return false;
    }
  }
  return true;
}

std::string Binning::describe() const {
  std::ostringstream os;
  os << bins_.size() << " bins in";
  for (std::size_t i = 0; i < variables_.size(); ++i) os << " " << variables_[i] << (isRange_[i] ? "" : "(point)");
  if (!source_.empty()) os << " from " << source_;
  return os.str();
}

}  // namespace zaku
