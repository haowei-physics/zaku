#ifndef ZAKU_EVENT_EVENTJOIN_H
#define ZAKU_EVENT_EVENTJOIN_H
#include <string>

#include "zaku/event/FlatTreeCursor.h"
#include "zaku/event/WeightCursor.h"

namespace zaku {

struct JoinOptions {
  /// Relative tolerance on the neutrino energy match.
  double energyRelTol = 1.0e-6;
  /// Stop with an error as soon as a weight row cannot be matched to the
  /// current flat-tree event (a desynchronised pair of files).
  bool strict = true;
};

struct JoinStats {
  long long flatEntries = 0;
  long long weightRows = 0;
  long long matched = 0;
  long long unmatchedFlat = 0;           ///< flat entries with no weight row
  long long unmatchedEmptyFinalState = 0;  ///< of which nfsp == 0 (the expected reason)
  long long leptonMismatch = 0;          ///< matched rows whose PDGLep/ELep differ (multi-lepton DIS)
  long long eventNumberMatches = 0;      ///< rows matched through an explicit event number
  double maxEnergyRelDiff = 0.0;

  std::string describe() const;
  /// Throw std::runtime_error unless every weight row was consumed and every
  /// unmatched flat entry has an empty final state.
  void check() const;
};

/// Walks a flat tree and a weight source in lockstep.
///
/// T2KReWeight writes no event index, and skips events with an empty final
/// state, so the rows are matched by a forward cursor on (Mode, PDGnu,
/// neutrino energy). Lepton quantities are deliberately not part of the key:
/// in a handful of multi-lepton DIS events T2KReWeight picks a different
/// "final-state lepton" than NUISANCE. Sources with an explicit event number
/// (GENIE) are matched on it directly.
class EventJoin {
public:
  /// `weights` may be null: the join then just iterates the flat tree.
  EventJoin(FlatTreeCursor& flat, WeightSource* weights, JoinOptions options = JoinOptions());

  /// Advance to the next flat-tree entry; false at the end.
  bool next();
  const FlatEvent& event() const { return flat_.event(); }
  /// True when the current entry has a matching weight row.
  bool hasWeights() const { return hasWeights_; }
  const WeightSource& weights() const;
  WeightSource* weightSource() const { return weights_; }
  long long weightRow() const { return hasWeights_ ? weights_->entry() : -1; }

  const JoinStats& stats() const { return stats_; }
  /// Validate the bookkeeping at the end (also called automatically by next() at the end).
  const JoinStats& finish();

  static bool keyMatches(const FlatEvent& ev, const WeightSource& w, double energyRelTol, double* relDiff = nullptr);

private:
  FlatTreeCursor& flat_;
  WeightSource* weights_ = nullptr;
  JoinOptions options_;
  JoinStats stats_;
  bool hasWeights_ = false;
  bool weightsAhead_ = false;   ///< a weight row is loaded and not yet consumed
  bool weightsDone_ = false;
  bool finished_ = false;
};

}  // namespace zaku

#endif  // ZAKU_EVENT_EVENTJOIN_H
