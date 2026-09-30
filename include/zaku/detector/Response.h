#ifndef ZAKU_DETECTOR_RESPONSE_H
#define ZAKU_DETECTOR_RESPONSE_H
#include <memory>
#include <string>
#include <vector>

#include "zaku/core/Types.h"
#include "zaku/core/Units.h"
#include "zaku/detector/Formula.h"
#include "zaku/detector/Variables.h"

class TH1;

namespace zaku {

/// A named region of true-event space: the first category whose selector is
/// true owns the event (`RespCat`). Selectors are Formula expressions over the
/// truth variables.
struct ResponseCategory {
  std::string name;
  std::string select;
};

/// One axis of an inline response table: a truth variable and its bin edges (GeV where kinematic).
struct TableAxis {
  std::string variable;
  std::vector<double> edges;
  int nBins() const { return static_cast<int>(edges.size()) - 1; }
  /// Bin index for a value, clamped into [0, nBins) like the legacy FindBin clamp.
  int bin(double x) const;
};

/// Where a table value comes from: a constant, an inline array over the
/// table's axes (flattened, last axis fastest), a ROOT histogram (TH1/TH2/TH3
/// whose axes are the table axes in order), or a formula over the truth
/// variables. Histogram and array lookups clamp to the table range.
class TableSource {
public:
  enum class Kind : int { Constant = 0, Array = 1, Histogram = 2, Formula = 3 };

  TableSource();
  static TableSource constant(double value);
  static TableSource array(const std::vector<double>& values);
  static TableSource histogram(const std::string& file, const std::string& name, EnergyUnit units = EnergyUnit::GeV);
  /// Take ownership of an already loaded histogram (axes in GeV).
  static TableSource histogram(TH1* h, const std::string& label = "");
  static TableSource formula(const std::string& expression);

  Kind kind() const { return kind_; }
  double constantValue() const { return constant_; }
  const std::vector<double>& arrayValues() const { return array_; }
  const std::string& formulaExpression() const { return expression_; }
  const std::string& histogramLabel() const { return label_; }
  const TH1* histogramPtr() const { return hist_.get(); }

  /// Bind to a table's axes and variable set (compiles formulas, checks sizes).
  void prepare(const std::vector<TableAxis>& axes, const VariableSet& vars, const std::string& where);
  /// Evaluate for the current values of `vars`.
  double evaluate(const VariableSet& vars) const;

  std::string describe() const;

private:
  Kind kind_ = Kind::Constant;
  double constant_ = 0.0;
  std::vector<double> array_;
  std::shared_ptr<TH1> hist_;
  std::string label_;
  std::string expression_;
  Formula formula_;
  std::vector<int> axisIndex_;      ///< variable index of each axis
  std::vector<TableAxis> axes_;     ///< copy of the bound axes (for arrays)
};

/// A table: for events of one category (or matching `select`), in one beam
/// mode or any, the detection probability x and the probability y_c of every
/// output column c. Columns are matched to samples by name and beam mode.
struct ResponseTable {
  std::string category;                 ///< name of the category, or empty with `select`
  std::string select;                   ///< inline selector when no category is given
  int beamMode = -1;                    ///< -1 any, else BeamMode
  std::vector<TableAxis> axes;
  TableSource efficiency;
  std::vector<std::string> columns;     ///< output column names, in draw order
  std::vector<TableSource> probabilities;  ///< one per column

  bool appliesTo(BeamMode mode) const { return beamMode < 0 || beamMode == static_cast<int>(mode); }
  std::string describe() const;
};

/// What the response says about one event.
struct ResponseResult {
  int category = -1;                  ///< RespCat: index into categories(), -1 if none matched
  int table = -1;                     ///< index into tables(), -1 if none applies
  double efficiency = 0.0;            ///< x: written as SelEffWeight
  std::vector<double> probabilities;  ///< y_c per column of the table
  double probabilitySum() const;
};

/// The user's detector design as a map from true events to detection and
/// sample probabilities. No geometry, no rings, no PID: only tables.
class Response {
public:
  Response() = default;

  Response& addCategory(const std::string& name, const std::string& select);
  Response& addTable(const ResponseTable& table);
  const std::vector<ResponseCategory>& categories() const { return categories_; }
  const std::vector<ResponseTable>& tables() const { return tables_; }
  int categoryIndex(const std::string& name) const;
  /// Every column name used by any table, in first-seen order.
  std::vector<std::string> columns() const;

  /// True if rows whose probabilities sum above 1 are rescaled (legacy) rather than rejected.
  bool normaliseProbabilities() const { return normalise_; }
  void setNormaliseProbabilities(bool on) { normalise_ = on; }

  /// Compile selectors and bind the tables to a variable set. Must be called
  /// (once) before evaluate(); the same VariableSet layout must be used after.
  void prepare(const VariableSet& vars);
  bool prepared() const { return prepared_; }

  /// Evaluate for an event whose truth variables are loaded in `vars`.
  ResponseResult evaluate(const VariableSet& vars, BeamMode mode) const;
  /// The category of the event (first matching selector), or -1.
  int categoryOf(const VariableSet& vars) const;

  /// The legacy HK response: 9 categories (numu/numubar/nue/nuebar x CCQE/nonQE,
  /// NC) from `h_eff_<cat>` and `h_mix_<cat>_<col>` histograms with the columns
  /// 1Rmu, 1Re0de, 1Re1de, each table binned in Enu_true.
  static Response fromLegacyFile(const std::string& responseFile);
  /// The nine legacy category names in RespCat order.
  static std::vector<std::string> legacyCategoryNames();
  /// Formula selecting one legacy category.
  static std::string legacyCategorySelector(const std::string& name);

  /// YAML representation (see docs/response.md).
  static Response fromYamlString(const std::string& yaml, const std::string& baseDir = "");
  static Response fromFile(const std::string& path);
  std::string toYaml() const;

  std::string describe() const;

private:
  std::vector<ResponseCategory> categories_;
  std::vector<Formula> selectors_;
  std::vector<ResponseTable> tables_;
  std::vector<Formula> tableSelectors_;
  bool normalise_ = false;
  bool prepared_ = false;
};

}  // namespace zaku

#endif  // ZAKU_DETECTOR_RESPONSE_H
