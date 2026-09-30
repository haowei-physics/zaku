#include "zaku/detector/Response.h"

#include <TAxis.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>

#include <algorithm>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace zaku {

// --- TableAxis ---------------------------------------------------------------

int TableAxis::bin(double x) const {
  const int n = nBins();
  if (n <= 0) return 0;
  auto it = std::upper_bound(edges.begin(), edges.end(), x);
  int b = static_cast<int>(it - edges.begin()) - 1;
  if (b < 0) b = 0;
  if (b >= n) b = n - 1;
  return b;
}

// --- TableSource -------------------------------------------------------------

TableSource::TableSource() = default;

TableSource TableSource::constant(double value) {
  TableSource s;
  s.kind_ = Kind::Constant;
  s.constant_ = value;
  return s;
}

TableSource TableSource::array(const std::vector<double>& values) {
  TableSource s;
  s.kind_ = Kind::Array;
  s.array_ = values;
  return s;
}

TableSource TableSource::histogram(const std::string& file, const std::string& name, EnergyUnit units) {
  std::unique_ptr<TFile> f(TFile::Open(file.c_str(), "READ"));
  if (!f || f->IsZombie()) throw std::runtime_error("TableSource: cannot open '" + file + "'");
  TH1* h = f->Get<TH1>(name.c_str());
  if (!h) throw std::runtime_error("TableSource: no histogram '" + name + "' in '" + file + "'");
  TH1* c = histogramToGeV(*h, units, /*contentIsDensity=*/false, 1 | 2 | 4, name);
  TableSource s = histogram(c, file + ":" + name);
  return s;
}

TableSource TableSource::histogram(TH1* h, const std::string& label) {
  if (!h) throw std::invalid_argument("TableSource::histogram: null histogram");
  h->SetDirectory(nullptr);
  TableSource s;
  s.kind_ = Kind::Histogram;
  s.hist_.reset(h);
  s.label_ = label;
  return s;
}

TableSource TableSource::formula(const std::string& expression) {
  TableSource s;
  s.kind_ = Kind::Formula;
  s.expression_ = expression;
  return s;
}

void TableSource::prepare(const std::vector<TableAxis>& axes, const VariableSet& vars, const std::string& where) {
  axes_ = axes;
  axisIndex_.clear();
  for (const auto& a : axes) {
    const int i = vars.index(a.variable);
    if (i < 0) throw std::runtime_error(where + ": axis variable '" + a.variable + "' is not a known truth variable");
    axisIndex_.push_back(i);
  }
  switch (kind_) {
    case Kind::Constant:
      break;
    case Kind::Array: {
      std::size_t n = 1;
      for (const auto& a : axes) n *= static_cast<std::size_t>(a.nBins());
      if (axes.empty()) n = 1;
      if (array_.size() != n) {
        throw std::runtime_error(where + ": inline array has " + std::to_string(array_.size()) + " values but the axes define " +
                                 std::to_string(n) + " bins");
      }
      break;
    }
    case Kind::Histogram: {
      if (hist_->GetDimension() != static_cast<int>(axes.size())) {
        throw std::runtime_error(where + ": histogram '" + label_ + "' has " + std::to_string(hist_->GetDimension()) +
                                 " dimensions but the table has " + std::to_string(axes.size()) + " axes");
      }
      break;
    }
    case Kind::Formula:
      formula_ = Formula(expression_, vars.names());
      break;
  }
}

double TableSource::evaluate(const VariableSet& vars) const {
  switch (kind_) {
    case Kind::Constant:
      return constant_;
    case Kind::Array: {
      std::size_t flat = 0;
      for (std::size_t k = 0; k < axes_.size(); ++k) {
        const int b = axes_[k].bin(vars.get(axisIndex_[k]));
        flat = flat * static_cast<std::size_t>(axes_[k].nBins()) + static_cast<std::size_t>(b);
      }
      return array_[flat];
    }
    case Kind::Histogram: {
      const TH1& h = *hist_;
      auto clampBin = [](const TAxis* ax, double x) {
        int b = ax->FindFixBin(x);
        if (b < 1) b = 1;
        if (b > ax->GetNbins()) b = ax->GetNbins();
        return b;
      };
      const int bx = clampBin(h.GetXaxis(), vars.get(axisIndex_[0]));
      if (h.GetDimension() == 1) return h.GetBinContent(bx);
      const int by = clampBin(h.GetYaxis(), vars.get(axisIndex_[1]));
      if (h.GetDimension() == 2) return h.GetBinContent(bx, by);
      const int bz = clampBin(h.GetZaxis(), vars.get(axisIndex_[2]));
      return h.GetBinContent(bx, by, bz);
    }
    case Kind::Formula:
      return formula_.eval(vars.data());
  }
  return 0.0;
}

std::string TableSource::describe() const {
  std::ostringstream os;
  switch (kind_) {
    case Kind::Constant: os << constant_; break;
    case Kind::Array: os << "array[" << array_.size() << "]"; break;
    case Kind::Histogram: os << "hist(" << (label_.empty() ? hist_->GetName() : label_) << ")"; break;
    case Kind::Formula: os << "formula(" << expression_ << ")"; break;
  }
  return os.str();
}

// --- ResponseTable / ResponseResult -----------------------------------------

std::string ResponseTable::describe() const {
  std::ostringstream os;
  os << (category.empty() ? "select '" + select + "'" : "category " + category);
  if (beamMode >= 0) os << " [" << toString(static_cast<BeamMode>(beamMode)) << "]";
  os << ": axes";
  for (const auto& a : axes) os << " " << a.variable << "(" << a.nBins() << ")";
  os << ", efficiency " << efficiency.describe();
  for (std::size_t i = 0; i < columns.size(); ++i) os << ", " << columns[i] << "=" << probabilities[i].describe();
  return os.str();
}

double ResponseResult::probabilitySum() const {
  double s = 0.0;
  for (double p : probabilities) s += p;
  return s;
}

// --- Response ----------------------------------------------------------------

Response& Response::addCategory(const std::string& name, const std::string& select) {
  if (categoryIndex(name) >= 0) throw std::invalid_argument("Response: category '" + name + "' already exists");
  categories_.push_back(ResponseCategory{name, select});
  prepared_ = false;
  return *this;
}

Response& Response::addTable(const ResponseTable& table) {
  if (table.columns.size() != table.probabilities.size()) {
    throw std::invalid_argument("Response: a table has " + std::to_string(table.columns.size()) + " columns but " +
                                std::to_string(table.probabilities.size()) + " probability sources");
  }
  if (table.category.empty() && table.select.empty()) {
    throw std::invalid_argument("Response: a table needs a category or a select expression");
  }
  tables_.push_back(table);
  prepared_ = false;
  return *this;
}

int Response::categoryIndex(const std::string& name) const {
  for (std::size_t i = 0; i < categories_.size(); ++i) {
    if (categories_[i].name == name) return static_cast<int>(i);
  }
  return -1;
}

std::vector<std::string> Response::columns() const {
  std::vector<std::string> out;
  for (const auto& t : tables_) {
    for (const auto& c : t.columns) {
      if (std::find(out.begin(), out.end(), c) == out.end()) out.push_back(c);
    }
  }
  return out;
}

void Response::prepare(const VariableSet& vars) {
  selectors_.clear();
  for (const auto& c : categories_) {
    try {
      selectors_.emplace_back(c.select, vars.names());
    } catch (const std::exception& e) {
      throw std::runtime_error("Response category '" + c.name + "': " + e.what());
    }
  }
  tableSelectors_.clear();
  for (std::size_t i = 0; i < tables_.size(); ++i) {
    ResponseTable& t = tables_[i];
    const std::string where = "Response table " + std::to_string(i) + " (" + (t.category.empty() ? t.select : t.category) + ")";
    if (!t.category.empty() && categoryIndex(t.category) < 0) {
      throw std::runtime_error(where + ": unknown category '" + t.category + "'");
    }
    tableSelectors_.emplace_back(t.select, vars.names());
    t.efficiency.prepare(t.axes, vars, where + " efficiency");
    for (std::size_t c = 0; c < t.columns.size(); ++c) {
      t.probabilities[c].prepare(t.axes, vars, where + " column " + t.columns[c]);
    }
  }
  prepared_ = true;
}

int Response::categoryOf(const VariableSet& vars) const {
  for (std::size_t i = 0; i < selectors_.size(); ++i) {
    if (selectors_[i].empty() || selectors_[i].evalBool(vars.data())) return static_cast<int>(i);
  }
  return -1;
}

ResponseResult Response::evaluate(const VariableSet& vars, BeamMode mode) const {
  if (!prepared_) throw std::logic_error("Response::evaluate: call prepare() first");
  ResponseResult r;
  r.category = categoryOf(vars);
  for (std::size_t i = 0; i < tables_.size(); ++i) {
    const ResponseTable& t = tables_[i];
    if (!t.appliesTo(mode)) continue;
    bool match = false;
    if (!t.category.empty()) {
      match = r.category >= 0 && categories_[static_cast<std::size_t>(r.category)].name == t.category;
    } else {
      match = tableSelectors_[i].evalBool(vars.data());
    }
    if (!match) continue;
    r.table = static_cast<int>(i);
    r.efficiency = t.efficiency.evaluate(vars);
    r.probabilities.resize(t.columns.size());
    for (std::size_t c = 0; c < t.columns.size(); ++c) r.probabilities[c] = t.probabilities[c].evaluate(vars);
    const double sum = r.probabilitySum();
    if (sum > 1.0 + 1.0e-9) {
      if (!normalise_) {
        throw std::runtime_error("Response: the column probabilities of table " + std::to_string(i) + " sum to " +
                                 std::to_string(sum) + " > 1 for an event; set normalise_probabilities to rescale");
      }
      for (double& p : r.probabilities) p /= sum;
    }
    break;
  }
  return r;
}

std::vector<std::string> Response::legacyCategoryNames() {
  return {"numuCCQE", "numubarCCQE", "numuCCnonQE", "numubarCCnonQE", "nueCCQE",
          "nuebarCCQE", "nueCCnonQE", "nuebarCCnonQE", "NC"};
}

std::string Response::legacyCategorySelector(const std::string& name) {
  if (name == "NC") return "!cc";
  const bool anti = name.find("bar") != std::string::npos;
  const bool nue = name.compare(0, 3, "nue") == 0;
  const bool qe = name.find("nonQE") == std::string::npos;
  const int pdg = (nue ? 12 : 14) * (anti ? -1 : 1);
  return std::string("cc && PDGnu == ") + std::to_string(pdg) + (qe ? " && flagCCQE" : " && !flagCCQE");
}

Response Response::fromLegacyFile(const std::string& file) {
  std::unique_ptr<TFile> f(TFile::Open(file.c_str(), "READ"));
  if (!f || f->IsZombie()) throw std::runtime_error("Response::fromLegacyFile: cannot open '" + file + "'");
  const std::vector<std::string> cols = {"1Rmu", "1Re0de", "1Re1de"};
  Response r;
  for (const auto& cat : legacyCategoryNames()) r.addCategory(cat, legacyCategorySelector(cat));
  for (const auto& cat : legacyCategoryNames()) {
    TH1* eff = f->Get<TH1>(("h_eff_" + cat).c_str());
    if (!eff) throw std::runtime_error("Response::fromLegacyFile: no h_eff_" + cat + " in '" + file + "'");
    ResponseTable t;
    t.category = cat;
    TableAxis axis;
    axis.variable = "Enu_true";
    for (int i = 1; i <= eff->GetNbinsX() + 1; ++i) axis.edges.push_back(eff->GetXaxis()->GetBinLowEdge(i));
    t.axes.push_back(axis);
    auto clone = [](TH1* h) {
      auto* c = static_cast<TH1*>(h->Clone());
      c->SetDirectory(nullptr);
      return c;
    };
    t.efficiency = TableSource::histogram(clone(eff));
    for (const auto& col : cols) {
      TH1* mix = f->Get<TH1>(("h_mix_" + cat + "_" + col).c_str());
      if (!mix) throw std::runtime_error("Response::fromLegacyFile: no h_mix_" + cat + "_" + col + " in '" + file + "'");
      t.columns.push_back(col);
      t.probabilities.push_back(TableSource::histogram(clone(mix)));
    }
    r.addTable(t);
  }
  r.setNormaliseProbabilities(true);  // the legacy tables were renormalised when built
  return r;
}

std::string Response::describe() const {
  std::ostringstream os;
  os << "response: " << categories_.size() << " categories, " << tables_.size() << " tables, columns";
  for (const auto& c : columns()) os << " " << c;
  os << (normalise_ ? " (rows renormalised)" : "") << "\n";
  for (std::size_t i = 0; i < categories_.size(); ++i) {
    os << "  category " << i << " " << categories_[i].name << ": " << categories_[i].select << "\n";
  }
  for (std::size_t i = 0; i < tables_.size(); ++i) os << "  table " << i << " " << tables_[i].describe() << "\n";
  return os.str();
}

}  // namespace zaku
