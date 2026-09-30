#ifndef ZAKU_DETECTOR_VARIABLES_H
#define ZAKU_DETECTOR_VARIABLES_H
#include <map>
#include <string>
#include <vector>

#include "zaku/event/FlatEvent.h"

namespace zaku {

/// The dictionary of named values a response table or a cut can refer to:
/// the truth variables (see TruthVariables::names()) followed, once the
/// reconstruction has run, by the reconstructed observables. Values are held
/// in one flat vector so formulas can be evaluated without lookups per event.
class VariableSet {
public:
  VariableSet();

  /// Variable names in value order.
  const std::vector<std::string>& names() const { return names_; }
  /// Index of a variable or -1.
  int index(const std::string& name) const;
  bool has(const std::string& name) const { return index(name) >= 0; }
  std::size_t size() const { return names_.size(); }

  /// Add a variable (returns its index; existing names are returned as is).
  int add(const std::string& name);

  double* data() { return values_.data(); }
  const double* data() const { return values_.data(); }
  const std::vector<double>& values() const { return values_; }
  double get(int i) const { return values_[static_cast<std::size_t>(i)]; }
  double get(const std::string& name) const;
  void set(int i, double v) { values_[static_cast<std::size_t>(i)] = v; }
  void set(const std::string& name, double v);

  /// Load the truth variables of an event (the first TruthVariables::names().size() slots).
  void loadTruth(const FlatEvent& ev, const TruthVariables& t);

  /// The reconstructed observable names Zaku adds after the truth variables.
  static const std::vector<std::string>& recoNames();

private:
  std::vector<std::string> names_;
  std::map<std::string, int> index_;
  std::vector<double> values_;
  std::size_t nTruth_ = 0;
};

}  // namespace zaku

#endif  // ZAKU_DETECTOR_VARIABLES_H
