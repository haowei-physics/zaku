#include "zaku/detector/Variables.h"

#include <stdexcept>

#include "zaku/detector/Formula.h"

namespace zaku {

const std::vector<std::string>& VariableSet::recoNames() {
  static const std::vector<std::string> n = {
      "PLep_rec", "CosLep_rec", "Theta_rec", "Ehad_rec", "Elep_rec_mu", "Elep_rec_e",
      "Enu_rec_ccqe_mu", "Q2_rec_ccqe_mu", "Enu_rec_ccqe_e", "Q2_rec_ccqe_e", "Enu_rec_calo_mu", "Enu_rec_calo_e",
      "true_ME", "recon_ME", "RespCat", "SelEffWeight", "POTWeight", "OscWeight",
      // set once a sample hypothesis is chosen
      "Enu_rec", "Q2_rec", "Elep_rec", "EhadFrac_rec", "PDGLep_rec", "EventType", "SampleProb"};
  return n;
}

VariableSet::VariableSet() {
  for (const auto& n : TruthVariables::names()) add(Formula::formulaName(n));
  nTruth_ = names_.size();
  for (const auto& n : recoNames()) add(n);
}

int VariableSet::index(const std::string& name) const {
  auto it = index_.find(Formula::formulaName(name));
  return it == index_.end() ? -1 : it->second;
}

int VariableSet::add(const std::string& name) {
  const std::string n = Formula::formulaName(name);
  auto it = index_.find(n);
  if (it != index_.end()) return it->second;
  const int i = static_cast<int>(names_.size());
  names_.push_back(n);
  index_[n] = i;
  values_.push_back(0.0);
  return i;
}

double VariableSet::get(const std::string& name) const {
  const int i = index(name);
  if (i < 0) throw std::out_of_range("VariableSet: no variable '" + name + "'");
  return values_[static_cast<std::size_t>(i)];
}

void VariableSet::set(const std::string& name, double v) {
  const int i = index(name);
  if (i < 0) throw std::out_of_range("VariableSet: no variable '" + name + "'");
  values_[static_cast<std::size_t>(i)] = v;
}

void VariableSet::loadTruth(const FlatEvent& ev, const TruthVariables& t) {
  const std::vector<double> v = t.values(ev);
  if (v.size() != nTruth_) throw std::logic_error("VariableSet: truth variable count changed");
  for (std::size_t i = 0; i < nTruth_; ++i) values_[i] = v[i];
}

}  // namespace zaku
