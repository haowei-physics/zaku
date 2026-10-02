#ifndef ZAKU_SYST_SUMMARY_H
#define ZAKU_SYST_SUMMARY_H
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "zaku/core/Experiment.h"
#include "zaku/syst/ParameterSet.h"

class TH1D;

namespace zaku {

/// The +-1 sigma effect of one parameter on one sample's expected event count.
struct SummaryEntry {
  std::string set;        ///< xsec | flux | det | osc | shift
  std::string parameter;
  std::string sample;
  double nominal = 0.0;
  double plus = 0.0;      ///< expected events at prior + 1 sigma
  double minus = 0.0;     ///< ... at prior - 1 sigma

  double plusFrac() const { return nominal != 0.0 ? plus / nominal - 1.0 : 0.0; }
  double minusFrac() const { return nominal != 0.0 ? minus / nominal - 1.0 : 0.0; }
  /// The larger of the two fractional effects, in absolute value.
  double maxFrac() const;
};

/// The table of +-1 sigma effects, one row per (parameter, sample).
class SummaryTable {
public:
  void add(const SummaryEntry& e) { entries_.push_back(e); }
  const std::vector<SummaryEntry>& entries() const { return entries_; }
  std::size_t size() const { return entries_.size(); }
  bool empty() const { return entries_.empty(); }

  std::vector<std::string> sets() const;
  std::vector<std::string> parameters() const;
  std::vector<std::string> samples() const;
  /// The entry of one (parameter, sample), or nullptr.
  const SummaryEntry* find(const std::string& parameter, const std::string& sample) const;
  /// Entries sorted by the largest fractional effect, at most `n`.
  std::vector<SummaryEntry> top(int n) const;

  /// `set,parameter,sample,nominal,plus,minus,plus_frac,minus_frac`
  std::string csv() const;
  static SummaryTable fromCsv(const std::string& text);
  void writeCsv(const std::string& path) const;
  static SummaryTable readCsv(const std::string& path);
  /// Fixed-width table: one row per parameter, one column per sample, "+x.x/-y.y %".
  std::string text() const;

private:
  std::vector<SummaryEntry> entries_;
};

/// A binned spectrum of one sample in one variable (GeV or degrees).
struct Spectrum {
  std::string sample;
  std::string variable;
  std::vector<double> edges;
  std::vector<double> contents;   ///< one per bin; under/overflow are dropped

  void reset(const std::vector<double>& edges_);
  int find(double x) const;   ///< bin index or -1
  void fill(double x, double w);
  double total() const;
  int nBins() const { return static_cast<int>(contents.size()); }
  /// A TH1D with these edges (owned by the caller, detached from any directory).
  TH1D* histogram(const std::string& name, const std::string& title = "") const;
};

/// Computes, from the FakeSelect files of an experiment, the nominal spectrum
/// of every sample and the spectra at +-1 sigma of every parameter of the
/// four parameter sets and of every shift systematic:
///   xsec  - the event's response graph evaluated at prior +- sigma (empty graph = 1);
///   flux  - 1 +- sigma on the events of the parameter's (beam mode, species, E_true) bin;
///   det   - 1 +- sigma on the events of the sample;
///   osc   - the oscillation weight recomputed with the parameter moved;
///   shift - the reconstruction re-derived with the shift applied (events migrate between bins).
/// Spectra are in each sample's first binning variable; the shift universes
/// are also kept in Enu_rec for the detector covariance built from them.
class SystematicsSummary {
public:
  SystematicsSummary(const Experiment& exp, const ParameterSet& xsec, const ParameterSet& flux,
                     const ParameterSet& det, const ParameterSet& osc);
  ~SystematicsSummary();
  SystematicsSummary(const SystematicsSummary&) = delete;
  SystematicsSummary& operator=(const SystematicsSummary&) = delete;

  /// Accumulate one beam mode's FakeSelect file. Needs the response graphs of
  /// the fitted dials when the xsec set is non-empty.
  void process(BeamMode mode, const std::string& fakeSelectFile);
  /// Build the table from what was accumulated.
  SummaryTable finish();
  const SummaryTable& table() const { return table_; }

  const std::vector<Spectrum>& nominal() const { return nominal_; }
  /// The spectrum of a parameter at +1 (plus=true) or -1 sigma for a sample index.
  Spectrum spectrum(const std::string& parameter, int sample, bool plus) const;
  /// Every spectrum as a TH1D: <sample>/nominal, <sample>/<parameter>_plus, _minus.
  void writeHistograms(const std::string& file) const;

  /// The detector covariance implied by the shift systematics: for every
  /// (sample, Enu_rec bin) the fractional half-difference of the +-1 sigma
  /// spectra, summed as delta delta^T over the shifts (fractional, priors 1),
  /// with the binning `EventType Enu_rec Enu_rec`.
  ParameterSet shiftCovariance() const;

  long long entriesProcessed() const { return entries_; }

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
  SummaryTable table_;
  std::vector<Spectrum> nominal_;
  long long entries_ = 0;
};

}  // namespace zaku

#endif  // ZAKU_SYST_SUMMARY_H
