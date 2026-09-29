#include "zaku/event/WeightCursor.h"

#include <TBranch.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TObjArray.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace zaku {

// --- DialKnots ---------------------------------------------------------------

bool DialKnots::nominalIsKnot() const {
  const double tol = 1.0e-9 * (1.0 + std::fabs(nominal));
  for (double k : knots) {
    if (std::fabs(k - nominal) <= tol) return true;
  }
  return false;
}

std::vector<double> DialKnots::knotsWithNominal() const {
  std::vector<double> k = knots;
  if (!nominalIsKnot()) k.push_back(nominal);
  std::sort(k.begin(), k.end());
  return k;
}

std::string DialKnots::describe() const {
  std::ostringstream os;
  os << name << ": nominal " << nominal << ", knots";
  for (double k : knots) os << " " << k;
  return os.str();
}

// --- WeightSource ------------------------------------------------------------

std::vector<std::string> WeightSource::dialNames() const {
  std::vector<std::string> out;
  for (const auto& d : dials()) out.push_back(d.name);
  return out;
}

int WeightSource::dialIndex(const std::string& name) const {
  const auto& d = dials();
  for (std::size_t i = 0; i < d.size(); ++i) {
    if (d[i].name == name) return static_cast<int>(i);
  }
  return -1;
}

const std::vector<double>& WeightSource::weights(const std::string& dial) const {
  const int i = dialIndex(dial);
  if (i < 0) throw std::out_of_range("WeightSource: no dial '" + dial + "' in '" + fileName() + "'");
  return weights(i);
}

// --- T2KRWWeightCursor -------------------------------------------------------

struct T2KRWWeightCursor::Buffers {
  std::vector<std::vector<double>*> dialWeights;  // one per dial, owned by ROOT
  Int_t Mode = 0, PDGnu = 0, PDGLep = 0;
  Float_t Q2 = 0;
  Double_t weight = 0, cv_weight = 0;
  TLorentzVector* PISNu = nullptr;
  TLorentzVector* PFSLep = nullptr;
  std::vector<double> empty;
};

bool T2KRWWeightCursor::looksLikeT2KRW(const std::string& file) {
  std::unique_ptr<TFile> f(TFile::Open(file.c_str(), "READ"));
  if (!f || f->IsZombie()) return false;
  return f->Get<TTree>("T2KRWFlatTree") != nullptr && f->Get<TTree>("TweakMetaData") != nullptr;
}

T2KRWWeightCursor::T2KRWWeightCursor(const std::string& file) : fileName_(file), buf_(new Buffers) {
  file_.reset(TFile::Open(file.c_str(), "READ"));
  if (!file_ || file_->IsZombie()) throw std::runtime_error("T2KRWWeightCursor: cannot open '" + file + "'");
  tree_ = file_->Get<TTree>("T2KRWFlatTree");
  if (!tree_) throw std::runtime_error("T2KRWWeightCursor: no tree 'T2KRWFlatTree' in '" + file + "'");
  TTree* meta = file_->Get<TTree>("TweakMetaData");
  if (!meta) throw std::runtime_error("T2KRWWeightCursor: no tree 'TweakMetaData' in '" + file + "'");
  if (meta->GetEntries() < 1) throw std::runtime_error("T2KRWWeightCursor: 'TweakMetaData' is empty in '" + file + "'");

  // Dials: every metadata branch <name> (vector<double>) with a <name>_from partner.
  std::map<std::string, std::vector<double>*> knots;
  std::map<std::string, Double_t> nominal;
  std::vector<std::string> order;
  TObjArray* branches = meta->GetListOfBranches();
  for (int i = 0; i < branches->GetEntries(); ++i) {
    const std::string name = branches->At(i)->GetName();
    const std::string suffix = "_from";
    if (name.size() > suffix.size() && name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0) continue;
    if (!meta->GetBranch((name + suffix).c_str())) continue;
    knots[name] = nullptr;
    nominal[name] = 0.0;
    order.push_back(name);
    meta->SetBranchAddress(name.c_str(), &knots[name]);
    meta->SetBranchAddress((name + suffix).c_str(), &nominal[name]);
  }
  if (order.empty()) throw std::runtime_error("T2KRWWeightCursor: no dials found in 'TweakMetaData' of '" + file + "'");
  meta->GetEntry(0);
  for (const auto& name : order) {
    DialKnots d;
    d.name = name;
    if (knots[name]) d.knots = *knots[name];
    d.nominal = nominal[name];
    dials_.push_back(d);
  }
  meta->ResetBranchAddresses();

  // Weight tree: one vector<double> per dial plus the matching keys.
  Buffers& b = *buf_;
  b.dialWeights.assign(dials_.size(), nullptr);
  for (std::size_t i = 0; i < dials_.size(); ++i) {
    const std::string& name = dials_[i].name;
    if (!tree_->GetBranch(name.c_str())) {
      throw std::runtime_error("T2KRWWeightCursor: dial '" + name + "' has metadata but no weight branch in '" + file + "'");
    }
    tree_->SetBranchAddress(name.c_str(), &b.dialWeights[i]);
  }
  auto bind = [&](const char* name, void* addr) {
    if (!tree_->GetBranch(name)) throw std::runtime_error(std::string("T2KRWWeightCursor: missing branch ") + name);
    tree_->SetBranchAddress(name, addr);
  };
  bind("Mode", &b.Mode);
  bind("PDGnu", &b.PDGnu);
  bind("PDGLep", &b.PDGLep);
  bind("Q2", &b.Q2);
  bind("weight", &b.weight);
  bind("cv_weight", &b.cv_weight);
  bind("PISNu", &b.PISNu);
  bind("PFSLep", &b.PFSLep);
}

T2KRWWeightCursor::~T2KRWWeightCursor() {
  if (tree_) tree_->ResetBranchAddresses();
  if (file_) file_->Close();
}

long long T2KRWWeightCursor::entries() const { return tree_ ? tree_->GetEntries() : 0; }

bool T2KRWWeightCursor::goTo(long long entry) {
  if (entry < 0 || entry >= entries()) return false;
  if (tree_->GetEntry(entry) <= 0) throw std::runtime_error("T2KRWWeightCursor: cannot read row " + std::to_string(entry));
  entry_ = entry;
  return true;
}

bool T2KRWWeightCursor::next() { return goTo(entry_ + 1); }

const std::vector<double>& T2KRWWeightCursor::weights(int i) const {
  if (i < 0 || i >= static_cast<int>(buf_->dialWeights.size())) {
    throw std::out_of_range("T2KRWWeightCursor: dial index " + std::to_string(i) + " out of range");
  }
  const std::vector<double>* w = buf_->dialWeights[static_cast<std::size_t>(i)];
  return w ? *w : buf_->empty;
}

double T2KRWWeightCursor::cvWeight() const { return buf_->cv_weight; }
double T2KRWWeightCursor::scaleWeight() const { return buf_->weight; }
int T2KRWWeightCursor::Mode() const { return buf_->Mode; }
int T2KRWWeightCursor::PDGnu() const { return buf_->PDGnu; }
double T2KRWWeightCursor::Enu() const { return buf_->PISNu ? buf_->PISNu->E() : 0.0; }
int T2KRWWeightCursor::PDGLep() const { return buf_->PDGLep; }
double T2KRWWeightCursor::ELep() const { return buf_->PFSLep ? buf_->PFSLep->E() : 0.0; }

}  // namespace zaku
