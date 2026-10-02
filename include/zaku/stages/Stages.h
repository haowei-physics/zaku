#ifndef ZAKU_STAGES_STAGES_H
#define ZAKU_STAGES_STAGES_H
#include <string>
#include <vector>

#include "zaku/core/Types.h"
#include "zaku/event/WeightCursor.h"
#include "zaku/syst/ParameterSet.h"
#include "zaku/syst/Summary.h"

class TH2D;

namespace zaku {

class Experiment;
struct Production;

/// The event-rate table a select stage produces: for every sample of a beam
/// mode (plus the rejected remainder), the weighted rate contributed by each
/// production and the number of entries. Weights are POT x oscillation x
/// efficiency x sample probability, i.e. expected events at the exposure.
class RateTable {
public:
  RateTable() = default;
  RateTable(std::vector<std::string> rows, std::vector<std::string> columns);

  const std::vector<std::string>& rows() const { return rows_; }
  const std::vector<std::string>& columns() const { return columns_; }
  void add(int row, int column, double weight);
  double weight(int row, int column) const;
  long long entries(int row, int column) const;
  double rowTotal(int row) const;
  long long rowEntries(int row) const;
  double columnTotal(int column) const;
  double total() const;

  /// Comma-separated table with a total column and an entries column.
  std::string csv() const;
  /// Fixed-width text table (rows = samples, columns = productions).
  std::string text() const;
  /// TH2D (x = production, y = sample) with the weighted rates, owned by the caller.
  TH2D* histogram(const std::string& name = "rates") const;
  /// Read a table back from csv() output.
  static RateTable fromCsv(const std::string& text);
  void writeCsv(const std::string& path) const;

private:
  std::vector<std::string> rows_;
  std::vector<std::string> columns_;
  std::vector<double> weights_;
  std::vector<long long> counts_;
  std::size_t at(int r, int c) const { return static_cast<std::size_t>(r) * columns_.size() + static_cast<std::size_t>(c); }
};

namespace stages {

/// What the recon stage did to one production.
struct ReconSummary {
  std::string production;
  std::string inputFile;
  std::string outputFile;
  long long entries = 0;
  long long kept = 0;
  long long dropped = 0;
  double sumPOTWeight = 0.0;       ///< interactions in the fiducial mass at this exposure
  double sumOscillatedWeight = 0.0;
  unsigned seed = 0;
  bool skipped = false;            ///< output existed and force was off
  std::string describe() const;
};

/// What the select stage did for one beam mode.
struct SelectSummary {
  BeamMode mode = BeamMode::FHC;
  std::string outputFile;
  std::string rateTableFile;
  long long entriesIn = 0;
  long long entriesOut = 0;
  std::vector<std::string> inputFiles;
  RateTable rates;
  unsigned seed = 0;
  bool skipped = false;
  std::string describe() const;
};

/// What the splines stage did to one production.
struct SplineSummary {
  std::string production;
  std::string inputFile;    ///< the recon file
  std::string weightFile;
  std::string outputFile;
  long long entries = 0;
  std::vector<std::string> dials;       ///< graph branches written, in order
  std::vector<long long> nFlat;         ///< per dial: events without response (empty graph)
  std::vector<long long> nRejected;     ///< per dial: unusable responses (empty graph)
  std::string stats;                    ///< the per-dial table
  bool skipped = false;
  std::string describe() const;
};

/// The parameter sets of an experiment, as the covariances stage writes them.
struct ParameterSets {
  ParameterSet osc;
  ParameterSet det;
  ParameterSet xsec;
  ParameterSet fluxUnconstrained;
  ParameterSet fluxConstrained;
  std::vector<DialKnots> knots;   ///< the weight-file knots the xsec set was built from
  std::string knotsSource;        ///< the file they came from
  std::string describe() const;
};

/// What the covariances stage wrote.
struct CovarianceSummary {
  std::vector<std::string> files;
  ParameterSets sets;
  bool skipped = false;
  std::string describe() const;
};

/// What the summary stage produced.
struct SummaryResult {
  std::string csvFile;
  std::string rootFile;
  SummaryTable table;
  long long entries = 0;
  bool skipped = false;
  std::string describe() const;
};

/// Recon stage: flat tree -> recon file for one production.
ReconSummary reconstruct(const Experiment& exp, const Production& production, bool force = false);

/// Splines stage: recon file + weight file -> spline file (one TGraph branch per dial) for one production.
SplineSummary makeSplines(const Experiment& exp, const Production& production, bool force = false);

/// Select stage: recon (or spline) files of one beam mode -> FakeSelect tree + rate table.
SelectSummary select(const Experiment& exp, BeamMode mode, bool force = false);

/// The dial knots of the experiment's productions: from the first weight file
/// that exists, else from the metadata of the first spline file. `source`
/// receives the file used. Throws naming what to run when neither exists.
std::vector<DialKnots> dialKnots(const Experiment& exp, std::string* source = nullptr);

/// Build the parameter sets in memory (what the covariances stage writes and
/// the summary stage evaluates). Sets whose systematics are disabled stay empty.
ParameterSets buildParameterSets(const Experiment& exp);

/// Covariances stage: write the parameter-set files, binnings and the oscillation grid.
CovarianceSummary makeCovariances(const Experiment& exp, bool force = false);

/// Summary stage: the +-1 sigma effect of every parameter on every sample, from the FakeSelect files.
SummaryResult summarise(const Experiment& exp, bool force = false);

/// Names of the branches the recon stage adds to the flat-tree clone, in order.
const std::vector<std::string>& reconBranches();
/// Names of the scalar branches of the FakeSelect tree, in order.
const std::vector<std::string>& selectBranches();

/// The legacy recon seed (1234) or the derived per-production seed.
unsigned reconSeed(const Experiment& exp, const Production& production);
/// The legacy select seed (20250731) or the derived per-beam-mode seed.
unsigned selectSeed(const Experiment& exp, BeamMode mode);

}  // namespace stages

}  // namespace zaku

#endif  // ZAKU_STAGES_STAGES_H
