#ifndef ZAKU_DETECTOR_BINNING_H
#define ZAKU_DETECTOR_BINNING_H
#include <string>
#include <vector>

namespace zaku {

/// One dimension of a binning: a variable and either bin edges (a range
/// variable, `[low, high)` bins) or discrete values (a point variable).
struct BinningAxis {
  std::string variable;
  std::vector<double> edges;   ///< range variable when non-empty
  std::vector<double> values;  ///< point variable when non-empty
  bool isRange() const { return !edges.empty(); }
  int nBins() const { return isRange() ? static_cast<int>(edges.size()) - 1 : static_cast<int>(values.size()); }
  /// Uniform edges from min to max in steps of `step`, plus optional extra upper edges.
  static BinningAxis uniform(const std::string& variable, double min, double max, double step,
                             const std::vector<double>& extra = {});
};

/// One bin of a GUNDAM binning: per axis a [low, high) range or a point value.
struct Bin {
  std::vector<double> low;
  std::vector<double> high;   ///< equal to low for point variables
  bool contains(const std::vector<double>& x, const std::vector<bool>& isRange) const;
};

/// A GUNDAM sample binning: the Cartesian product of its axes (last axis
/// fastest), or an explicit list of bins read from a GUNDAM `.txt` file.
///
/// The text format: `#` starts a comment, a `variables:` line lists the
/// columns (a variable named twice is a range and takes two columns, low and
/// high; named once it is a point value), then one bin per line.
class Binning {
public:
  Binning() = default;

  static Binning product(const std::vector<BinningAxis>& axes);
  /// Parse a GUNDAM binning text file (all bins kept in file order).
  static Binning fromFile(const std::string& path);
  static Binning fromText(const std::string& text, const std::string& where = "text");

  const std::vector<std::string>& variables() const { return variables_; }
  const std::vector<bool>& isRange() const { return isRange_; }
  const std::vector<Bin>& bins() const { return bins_; }
  int nBins() const { return static_cast<int>(bins_.size()); }
  bool empty() const { return bins_.empty(); }
  /// The axes when built as a product (empty for file binnings).
  const std::vector<BinningAxis>& axes() const { return axes_; }
  /// Path of the file this binning was read from (empty if built in code).
  const std::string& sourceFile() const { return source_; }

  /// Index of the bin containing `x` (values in variable order), or -1.
  int find(const std::vector<double>& x) const;

  /// The GUNDAM text representation.
  std::string toText(const std::string& comment = "") const;
  void writeFile(const std::string& path, const std::string& comment = "") const;

  /// True if the two binnings have the same variables and bins (to 1e-9).
  bool sameAs(const Binning& other, double tol = 1e-9) const;

  std::string describe() const;

private:
  std::vector<std::string> variables_;
  std::vector<bool> isRange_;
  std::vector<Bin> bins_;
  std::vector<BinningAxis> axes_;
  std::string source_;
};

}  // namespace zaku

#endif  // ZAKU_DETECTOR_BINNING_H
