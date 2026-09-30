#ifndef ZAKU_DETECTOR_FORMULA_H
#define ZAKU_DETECTOR_FORMULA_H
#include <memory>
#include <string>
#include <vector>

class TFormula;

namespace zaku {

/// A ROOT TFormula expression over named variables, the same expression
/// syntax GUNDAM uses for its cuts (`cc && PDGnu == 14 && Enu_true < 1.2`,
/// `0.9 * (1 - exp(-Enu_true / 0.4))`, ...). Booleans are 0/1. The variable
/// names are fixed when the formula is built; eval() takes their values in the
/// same order.
class Formula {
public:
  Formula() = default;
  /// Throws std::runtime_error if the expression does not compile.
  Formula(const std::string& expression, const std::vector<std::string>& variables);

  bool empty() const { return expression_.empty(); }
  const std::string& expression() const { return expression_; }
  const std::vector<std::string>& variables() const { return variables_; }

  double eval(const double* values) const;
  double eval(const std::vector<double>& values) const { return eval(values.data()); }
  bool evalBool(const double* values) const { return eval(values) != 0.0; }
  bool evalBool(const std::vector<double>& values) const { return evalBool(values.data()); }

  /// Names among `variables` that the expression mentions (for diagnostics).
  std::vector<std::string> referencedVariables() const;

  /// Sanitise a variable name for TFormula (it reserves x, y, z, t):
  /// "x" -> "bjorken_x", "y" -> "bjorken_y"; everything else unchanged.
  static std::string formulaName(const std::string& variable);

private:
  std::shared_ptr<TFormula> formula_;
  std::string expression_;
  std::vector<std::string> variables_;
  std::vector<int> positions_;          ///< TFormula array position of each variable (-1 if unused)
  int ndim_ = 0;
  mutable std::vector<double> buffer_;  ///< scratch array in TFormula's variable order
};

}  // namespace zaku

#endif  // ZAKU_DETECTOR_FORMULA_H
