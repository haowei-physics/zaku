#include "zaku/detector/Formula.h"

#include <TFormula.h>
#include <TString.h>

#include <algorithm>
#include <atomic>
#include <cctype>
#include <stdexcept>

namespace zaku {

std::string Formula::formulaName(const std::string& variable) {
  if (variable == "x") return "bjorken_x";
  if (variable == "y") return "bjorken_y";
  if (variable == "z") return "var_z";
  if (variable == "t") return "var_t";
  return variable;
}

Formula::Formula(const std::string& expression, const std::vector<std::string>& variables)
    : expression_(expression), variables_(variables) {
  if (expression.empty()) return;
  static std::atomic<int> counter{0};
  const std::string name = "zaku_formula_" + std::to_string(counter++);
  formula_ = std::make_shared<TFormula>(name.c_str(), "", false);
  std::vector<TString> names;
  names.reserve(variables.size());
  for (const auto& v : variables) names.emplace_back(formulaName(v).c_str());
  if (!names.empty()) formula_->AddVariables(names.data(), static_cast<Int_t>(names.size()));
  if (formula_->Compile(expression.c_str()) != 0 || !formula_->IsValid()) {
    throw std::runtime_error("Formula: cannot compile '" + expression + "'");
  }
  // TFormula pre-defines x, y, z, t (and a spare slot) as its first array
  // positions, so user variables do not sit at their insertion index. Map our
  // value order onto TFormula's own positions once, here.
  positions_.assign(variables.size(), -1);
  for (std::size_t i = 0; i < variables.size(); ++i) {
    positions_[i] = formula_->GetVarNumber(formulaName(variables[i]).c_str());
  }
  ndim_ = formula_->GetNdim();
  buffer_.assign(static_cast<std::size_t>(std::max(ndim_, 1)), 0.0);
}

double Formula::eval(const double* values) const {
  if (!formula_) return 0.0;
  for (std::size_t i = 0; i < positions_.size(); ++i) {
    const int p = positions_[i];
    if (p >= 0 && p < ndim_) buffer_[static_cast<std::size_t>(p)] = values[i];
  }
  return formula_->EvalPar(buffer_.data(), nullptr);
}

std::vector<std::string> Formula::referencedVariables() const {
  std::vector<std::string> out;
  for (const auto& v : variables_) {
    const std::string n = formulaName(v);
    std::string::size_type pos = 0;
    bool found = false;
    while ((pos = expression_.find(n, pos)) != std::string::npos) {
      const bool startOk = pos == 0 || !(std::isalnum(static_cast<unsigned char>(expression_[pos - 1])) || expression_[pos - 1] == '_');
      const std::string::size_type end = pos + n.size();
      const bool endOk = end >= expression_.size() || !(std::isalnum(static_cast<unsigned char>(expression_[end])) || expression_[end] == '_');
      if (startOk && endOk) {
        found = true;
        break;
      }
      pos = end;
    }
    if (found) out.push_back(v);
  }
  return out;
}

}  // namespace zaku
