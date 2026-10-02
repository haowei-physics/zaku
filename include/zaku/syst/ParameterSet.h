#ifndef ZAKU_SYST_PARAMETERSET_H
#define ZAKU_SYST_PARAMETERSET_H
#include <TMatrixD.h>

#include <string>
#include <vector>

#include "zaku/detector/Binning.h"

namespace zaku {

/// One fit parameter: its name, prior, one-sigma width and optional bounds.
struct Parameter {
  std::string name;
  double prior = 0.0;
  double sigma = 0.0;
  double lower = 0.0;
  double upper = 0.0;
  bool hasBounds = false;
  bool fixed = false;   ///< not written to the file; consumed by the GUNDAM writer

  std::string describe() const;
};

/// A GUNDAM parameter set: names, priors, a covariance and optional bounds,
/// plus (for normalisation sets) the binning that pairs dial i with parameter
/// i by position. Written as the legacy MakeCovariances.C did:
///   <prefix>_param_names   TObjArray of TObjString (kSingleKey)
///   <prefix>_param_priors  TVectorD
///   <prefix>_param_cov     TMatrixD (never TMatrixDSym: GUNDAM's reader prefix-matches the class name)
///   <prefix>_param_lower/_upper  TVectorD, only when a parameter has bounds
class ParameterSet {
public:
  ParameterSet() = default;
  explicit ParameterSet(std::string prefix) : prefix_(std::move(prefix)) {}
  // TMatrixD::operator= refuses matrices of different shape instead of
  // resizing, so the copy operations are spelled out.
  ParameterSet(const ParameterSet& o);
  ParameterSet& operator=(const ParameterSet& o);

  const std::string& prefix() const { return prefix_; }
  ParameterSet& setPrefix(const std::string& p) {
    prefix_ = p;
    return *this;
  }
  /// Free text describing the set (written into the provenance of the file).
  const std::string& title() const { return title_; }
  ParameterSet& setTitle(const std::string& t) {
    title_ = t;
    return *this;
  }

  /// Append a parameter; the covariance grows with sigma^2 on the diagonal.
  ParameterSet& add(const Parameter& p);
  std::size_t size() const { return params_.size(); }
  bool empty() const { return params_.empty(); }
  const std::vector<Parameter>& parameters() const { return params_; }
  const Parameter& parameter(int i) const { return params_.at(static_cast<std::size_t>(i)); }
  Parameter& parameter(int i) { return params_.at(static_cast<std::size_t>(i)); }
  int index(const std::string& name) const;
  std::vector<std::string> names() const;
  bool anyBounds() const;

  const TMatrixD& covariance() const { return cov_; }
  /// Replace the covariance (must be size x size and symmetric to `symTol`).
  void setCovariance(const TMatrixD& cov, double symTol = 1e-9);
  /// Set the correlation between two parameters (both directions).
  void setCorrelation(int i, int j, double rho);
  /// sqrt of the diagonal.
  double sigma(int i) const;
  double correlation(int i, int j) const;
  /// Symmetric eigen decomposition; false if the smallest eigenvalue is <= 0.
  bool isPositiveDefinite(double* lambdaMin = nullptr, double* lambdaMax = nullptr) const;

  bool hasBinning() const { return hasBinning_; }
  const Binning& binning() const { return binning_; }
  ParameterSet& setBinning(const Binning& b) {
    binning_ = b;
    hasBinning_ = true;
    return *this;
  }

  /// Write the objects listed above into a (re-created) ROOT file.
  void write(const std::string& file) const;
  /// Read them back. Bounds are read when present; `fixed` is not stored.
  static ParameterSet read(const std::string& file, const std::string& prefix);

  /// True if names match and priors, bounds and covariance agree to `tol`
  /// (absolute on values below 1, relative above). `why` names the first difference.
  bool sameAs(const ParameterSet& other, double tol = 1e-9, std::string* why = nullptr) const;

  std::string describe() const;

private:
  std::string prefix_;
  std::string title_;
  std::vector<Parameter> params_;
  TMatrixD cov_;
  Binning binning_;
  bool hasBinning_ = false;
};

}  // namespace zaku

#endif  // ZAKU_SYST_PARAMETERSET_H
