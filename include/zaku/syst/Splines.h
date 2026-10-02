#ifndef ZAKU_SYST_SPLINES_H
#define ZAKU_SYST_SPLINES_H
#include <string>
#include <vector>

#include "zaku/event/WeightCursor.h"

class TGraph;

namespace zaku {

/// What became of one (event, dial) response.
enum class SplineDecision : int {
  Full = 0,      ///< a graph with every knot
  Flat = 1,      ///< every stored weight is exactly 1: no response, written as an empty graph
  Rejected = 2   ///< unusable (all zero, negative, non-finite, or not 1 at the nominal): empty graph
};
std::string toString(SplineDecision d);

/// One dial as the spline builder sees it: the knot positions sorted
/// ascending with the nominal inserted (weight 1 by construction) unless it
/// coincides with a tweak point, and per knot the index into the weight
/// vector (-1 for the inserted nominal). Counters accumulate over build().
struct SplineDial {
  std::string name;
  double nominal = 0.0;
  std::vector<double> x;
  std::vector<int> src;
  long long nFull = 0;
  long long nFlat = 0;
  long long nRejected = 0;
  double worstNominalDeviation = 0.0;  ///< largest |w - 1| seen at a tweak point equal to the nominal

  int nKnots() const { return static_cast<int>(x.size()); }
  bool nominalInserted() const;
  /// True if the knots (nominal included) are evenly spaced to `tol` relative
  /// to the spacing: GUNDAM's catmull-rom needs this; otherwise akima.
  bool uniformKnots(double tol = 1e-6) const;
  std::string describe() const;
};

/// Builds per-event response graphs from the knot weights of a weight source,
/// with the rules of the legacy MakeSplines.C: knots sorted, nominal inserted
/// with weight 1, flat responses and unusable ones written as empty graphs
/// (GUNDAM treats an empty point list as "no dial"; a graph flat at 0 would
/// silently delete the event through a Shift dial).
class SplineBuilder {
public:
  /// Every dial of `knots`.
  explicit SplineBuilder(const std::vector<DialKnots>& knots);
  /// Only the dials named in `keep` (in `knots` order); unknown names are an error.
  SplineBuilder(const std::vector<DialKnots>& knots, const std::vector<std::string>& keep);

  const std::vector<SplineDial>& dials() const { return dials_; }
  std::size_t size() const { return dials_.size(); }
  std::vector<std::string> dialNames() const;
  int index(const std::string& name) const;
  /// Index of dial `i` in the source's dial list (for WeightSource::weights).
  int sourceIndex(int i) const { return sourceIndex_.at(static_cast<std::size_t>(i)); }
  const DialKnots& knots(int i) const { return knots_.at(static_cast<std::size_t>(i)); }

  /// Fill `graph` for dial `i` from the weights at its tweak points (in the
  /// weight file's order) and update the dial's counters.
  SplineDecision build(int i, const std::vector<double>& weights, TGraph& graph);
  /// The decision alone, without touching the graph or the counters.
  static SplineDecision classify(const SplineDial& dial, const std::vector<double>& weights,
                                 double* nominalDeviation = nullptr);

  /// Per-dial table of flat / rejected counts for `entries` events.
  std::string statsText(long long entries) const;
  long long totalRejected() const;

  /// Knot metadata as one line per dial, `name=nominal:k1,k2,...`, so a spline
  /// file can carry the knots it was built from; and the inverse.
  static std::string metadataText(const std::vector<DialKnots>& knots);
  static std::vector<DialKnots> fromMetadataText(const std::string& text);

private:
  void setup(const std::vector<DialKnots>& knots, const std::vector<std::string>* keep);
  std::vector<DialKnots> knots_;
  std::vector<SplineDial> dials_;
  std::vector<int> sourceIndex_;
};

}  // namespace zaku

#endif  // ZAKU_SYST_SPLINES_H
