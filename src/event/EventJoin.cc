#include "zaku/event/EventJoin.h"

#include <cmath>
#include <sstream>
#include <stdexcept>

namespace zaku {

std::string JoinStats::describe() const {
  std::ostringstream os;
  os << "join: " << flatEntries << " flat entries, " << weightRows << " weight rows, " << matched << " matched, "
     << unmatchedFlat << " unmatched (" << unmatchedEmptyFinalState << " with empty final state), " << leptonMismatch
     << " lepton mismatches, max |dE/E| " << maxEnergyRelDiff;
  if (eventNumberMatches) os << ", " << eventNumberMatches << " matched by event number";
  return os.str();
}

void JoinStats::check() const {
  if (weightRows > 0 && matched != weightRows) {
    throw std::runtime_error("EventJoin: " + std::to_string(matched) + " of " + std::to_string(weightRows) +
                             " weight rows were matched to flat-tree entries; the files do not belong together");
  }
  if (unmatchedFlat != unmatchedEmptyFinalState) {
    throw std::runtime_error("EventJoin: " + std::to_string(unmatchedFlat - unmatchedEmptyFinalState) +
                             " flat-tree entries with a non-empty final state have no weight row");
  }
}

EventJoin::EventJoin(FlatTreeCursor& flat, WeightSource* weights, JoinOptions options)
    : flat_(flat), weights_(weights), options_(options) {
  // The join always walks both inputs from their first entry.
  flat_.rewind();
  if (weights_) weights_->rewind();
  stats_.weightRows = weights_ ? weights_->entries() : 0;
}

const WeightSource& EventJoin::weights() const {
  if (!weights_ || !hasWeights_) throw std::logic_error("EventJoin: the current entry has no weights");
  return *weights_;
}

bool EventJoin::keyMatches(const FlatEvent& ev, const WeightSource& w, double energyRelTol, double* relDiff) {
  if (w.Mode() != ev.Mode || w.PDGnu() != ev.PDGnu) return false;
  const double e = ev.Enu_true;
  const double d = std::fabs(w.Enu() - e) / (std::fabs(e) > 1.0e-12 ? std::fabs(e) : 1.0e-12);
  if (relDiff) *relDiff = d;
  return d <= energyRelTol;
}

bool EventJoin::next() {
  if (!flat_.next()) {
    finish();
    return false;
  }
  ++stats_.flatEntries;
  hasWeights_ = false;
  if (!weights_) return true;

  const FlatEvent& ev = flat_.event();
  if (!weightsAhead_ && !weightsDone_) {
    if (weights_->next()) weightsAhead_ = true;
    else weightsDone_ = true;
  }
  if (!weightsAhead_) {
    ++stats_.unmatchedFlat;
    if (ev.emptyFinalState()) ++stats_.unmatchedEmptyFinalState;
    else if (options_.strict) {
      throw std::runtime_error("EventJoin: flat-tree " + ev.describe() + " has no weight row (weight file '" +
                               weights_->fileName() + "' is exhausted)");
    }
    return true;
  }

  bool matched = false;
  double relDiff = 0.0;
  const long long evnum = weights_->eventNumber();
  if (evnum >= 0) {
    if (evnum == ev.entry) {
      matched = true;
      ++stats_.eventNumberMatches;
    } else if (evnum < ev.entry) {
      throw std::runtime_error("EventJoin: weight row for event " + std::to_string(evnum) +
                               " has no flat-tree entry (current entry " + std::to_string(ev.entry) + ")");
    }
  } else {
    matched = keyMatches(ev, *weights_, options_.energyRelTol, &relDiff);
  }

  if (matched) {
    hasWeights_ = true;
    weightsAhead_ = false;
    ++stats_.matched;
    if (relDiff > stats_.maxEnergyRelDiff) stats_.maxEnergyRelDiff = relDiff;
    if (weights_->PDGLep() != ev.PDGLep ||
        std::fabs(weights_->ELep() - ev.ELep) > options_.energyRelTol * std::fabs(ev.ELep) + 1.0e-9) {
      ++stats_.leptonMismatch;
    }
    return true;
  }

  ++stats_.unmatchedFlat;
  if (ev.emptyFinalState()) {
    ++stats_.unmatchedEmptyFinalState;
  } else if (options_.strict) {
    std::ostringstream os;
    os << "EventJoin: flat-tree " << ev.describe() << " does not match weight row " << weights_->entry() << " (mode "
       << weights_->Mode() << ", nu " << weights_->PDGnu() << ", E=" << weights_->Enu() << "); '"
       << weights_->fileName() << "' does not correspond to '" << flat_.fileName() << "'";
    throw std::runtime_error(os.str());
  }
  return true;
}

const JoinStats& EventJoin::finish() {
  if (finished_) return stats_;
  finished_ = true;
  if (weights_) {
    // A loaded-but-unconsumed row, or rows never reached, are leftovers.
    if (weightsAhead_) {
      throw std::runtime_error("EventJoin: weight row " + std::to_string(weights_->entry()) + " of '" +
                               weights_->fileName() + "' matches no flat-tree entry");
    }
    if (options_.strict) stats_.check();
  }
  return stats_;
}

}  // namespace zaku
