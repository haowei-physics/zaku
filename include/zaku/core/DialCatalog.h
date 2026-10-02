#ifndef ZAKU_CORE_DIALCATALOG_H
#define ZAKU_CORE_DIALCATALOG_H
#include <string>
#include <vector>

namespace zaku {

/// One reweighting dial of a generator: its name in the reweighting engine,
/// the tweak points (in units of the dial's one-sigma) at which weights are
/// evaluated, and whether it enters the cross-section parameter set of the fit.
struct DialSpec {
  std::string name;
  std::vector<double> tweaks;  ///< sigma multiples, e.g. [-2, -1, 1, 2]
  bool fit = true;             ///< false: weights are produced but the dial is not fitted
  std::string note;            ///< free text (which model the dial belongs to, why it is excluded, ...)

  std::string describe() const;
};

/// The dial catalogue of a generator: which systematics exist and how they
/// are tweaked. The catalogue is data (`configs/dials/*.yaml`); the reweight
/// stage writes it into the engine's card, the systematics stage turns the
/// fitted dials into the cross-section parameter set.
class DialCatalog {
public:
  DialCatalog() = default;

  const std::string& generator() const { return generator_; }
  DialCatalog& setGenerator(const std::string& g) {
    generator_ = g;
    return *this;
  }

  DialCatalog& add(const DialSpec& d);
  DialCatalog& add(const std::string& name, const std::vector<double>& tweaks, bool fit = true,
                   const std::string& note = "");
  const std::vector<DialSpec>& dials() const { return dials_; }
  std::size_t size() const { return dials_.size(); }
  bool empty() const { return dials_.empty(); }
  std::vector<std::string> names() const;
  /// Names of the dials with `fit` set, minus `exclude`, in catalogue order.
  std::vector<std::string> fitted(const std::vector<std::string>& exclude = {}) const;
  bool has(const std::string& name) const { return index(name) >= 0; }
  int index(const std::string& name) const;
  const DialSpec& dial(const std::string& name) const;

  /// YAML: `generator: neut` and `dials: [{name, tweaks, fit, note}]`.
  static DialCatalog fromFile(const std::string& path);
  static DialCatalog fromYamlString(const std::string& yaml);
  std::string toYaml() const;

  /// The 19 NEUT dials of the reference production with the T2KReWeight
  /// `neut_dial_tweaks.toml` tweak points; the four the legacy fit excluded
  /// (MDLSPiEj, NucleonFSI_PB, SRCFrac_C, RadCorr_rate) have `fit: false`.
  static DialCatalog neutDefault();

  std::string describe() const;

private:
  std::string generator_ = "neut";
  std::vector<DialSpec> dials_;
};

}  // namespace zaku

#endif  // ZAKU_CORE_DIALCATALOG_H
