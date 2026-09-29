#include "zaku/event/FlatTreeCursor.h"

#include <TBranch.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TObjArray.h>
#include <TTree.h>

#include <stdexcept>

namespace zaku {

namespace {
constexpr int kMaxFSP = 512;
}

struct FlatTreeCursor::Buffers {
  Int_t Mode = 0;
  Char_t cc = 0;
  Int_t PDGnu = 0;
  Float_t Enu_true = 0;
  Int_t tgt = 0, tgta = 0, tgtz = 0;
  Int_t PDGLep = 0;
  Float_t ELep = 0, CosLep = 0, Q2 = 0, q0 = 0, q3 = 0, Enu_QE = 0, Q2_QE = 0, W = 0, x = 0, y = 0;
  Double_t fScaleFactor = 0;
  Float_t Weight = 1;
  Int_t nfsp = 0;
  Int_t pdg[kMaxFSP] = {};
  Float_t px[kMaxFSP] = {}, py[kMaxFSP] = {}, pz[kMaxFSP] = {}, E[kMaxFSP] = {};
  Bool_t flagCCINC = 0, flagNCINC = 0, flagCCQE = 0, flagCC0pi = 0, flagCCQELike = 0, flagNCEL = 0, flagNC0pi = 0,
         flagCCcoh = 0, flagNCcoh = 0, flagCC1pip = 0, flagNC1pip = 0, flagCC1pim = 0, flagNC1pim = 0,
         flagCC1pi0 = 0, flagNC1pi0 = 0;
  bool hasW = false, hasQ0 = false, hasWeight = false;
};

const std::vector<std::string>& FlatTreeCursor::requiredBranches() {
  static const std::vector<std::string> b = {"Mode", "cc",  "PDGnu", "Enu_true", "PDGLep", "ELep", "CosLep",
                                             "Q2",   "fScaleFactor", "nfsp", "pdg", "px", "py", "pz", "E",
                                             "flagCCQE"};
  return b;
}

int FlatTreeCursor::maxFinalStateParticles() { return kMaxFSP; }

FlatTreeCursor::FlatTreeCursor(const std::string& file, const std::string& treeName)
    : fileName_(file), buf_(new Buffers) {
  file_.reset(TFile::Open(file.c_str(), "READ"));
  if (!file_ || file_->IsZombie()) throw std::runtime_error("FlatTreeCursor: cannot open '" + file + "'");
  tree_ = file_->Get<TTree>(treeName.c_str());
  if (!tree_) throw std::runtime_error("FlatTreeCursor: no tree '" + treeName + "' in '" + file + "'");

  std::vector<std::string> missing;
  for (const auto& b : requiredBranches()) {
    if (!tree_->GetBranch(b.c_str())) missing.push_back(b);
  }
  if (!missing.empty()) {
    std::string list;
    for (const auto& m : missing) list += (list.empty() ? "" : ", ") + m;
    throw std::runtime_error("FlatTreeCursor: '" + file + "' is missing required branches: " + list);
  }

  auto bind = [&](const char* name, void* addr, bool required) -> bool {
    if (!tree_->GetBranch(name)) {
      if (required) throw std::runtime_error(std::string("FlatTreeCursor: missing branch ") + name);
      return false;
    }
    tree_->SetBranchAddress(name, addr);
    return true;
  };
  Buffers& b = *buf_;
  bind("Mode", &b.Mode, true);
  bind("cc", &b.cc, true);
  bind("PDGnu", &b.PDGnu, true);
  bind("Enu_true", &b.Enu_true, true);
  bind("tgt", &b.tgt, false);
  bind("tgta", &b.tgta, false);
  bind("tgtz", &b.tgtz, false);
  bind("PDGLep", &b.PDGLep, true);
  bind("ELep", &b.ELep, true);
  bind("CosLep", &b.CosLep, true);
  bind("Q2", &b.Q2, true);
  b.hasQ0 = bind("q0", &b.q0, false);
  bind("q3", &b.q3, false);
  bind("Enu_QE", &b.Enu_QE, false);
  bind("Q2_QE", &b.Q2_QE, false);
  b.hasW = bind("W", &b.W, false);
  bind("x", &b.x, false);
  bind("y", &b.y, false);
  bind("fScaleFactor", &b.fScaleFactor, true);
  b.hasWeight = bind("Weight", &b.Weight, false);
  bind("nfsp", &b.nfsp, true);
  bind("pdg", b.pdg, true);
  bind("px", b.px, true);
  bind("py", b.py, true);
  bind("pz", b.pz, true);
  bind("E", b.E, true);
  bind("flagCCINC", &b.flagCCINC, false);
  bind("flagNCINC", &b.flagNCINC, false);
  bind("flagCCQE", &b.flagCCQE, true);
  bind("flagCC0pi", &b.flagCC0pi, false);
  bind("flagCCQELike", &b.flagCCQELike, false);
  bind("flagNCEL", &b.flagNCEL, false);
  bind("flagNC0pi", &b.flagNC0pi, false);
  bind("flagCCcoh", &b.flagCCcoh, false);
  bind("flagNCcoh", &b.flagNCcoh, false);
  bind("flagCC1pip", &b.flagCC1pip, false);
  bind("flagNC1pip", &b.flagNC1pip, false);
  bind("flagCC1pim", &b.flagCC1pim, false);
  bind("flagNC1pim", &b.flagNC1pim, false);
  bind("flagCC1pi0", &b.flagCC1pi0, false);
  bind("flagNC1pi0", &b.flagNC1pi0, false);
}

FlatTreeCursor::~FlatTreeCursor() {
  if (tree_) tree_->ResetBranchAddresses();
  if (file_) file_->Close();
}

long long FlatTreeCursor::entries() const { return tree_ ? tree_->GetEntries() : 0; }

void FlatTreeCursor::copyBuffers() {
  const Buffers& b = *buf_;
  FlatEvent& e = event_;
  e.Mode = b.Mode;
  e.cc = b.cc != 0;
  e.PDGnu = b.PDGnu;
  e.Enu_true = b.Enu_true;
  e.tgt = b.tgt;
  e.tgta = b.tgta;
  e.tgtz = b.tgtz;
  e.PDGLep = b.PDGLep;
  e.ELep = b.ELep;
  e.CosLep = b.CosLep;
  e.Q2 = b.Q2;
  e.q0 = b.q0;
  e.q3 = b.q3;
  e.Enu_QE = b.Enu_QE;
  e.Q2_QE = b.Q2_QE;
  e.W = b.W;
  e.x = b.x;
  e.y = b.y;
  e.fScaleFactor = b.fScaleFactor;
  e.Weight = b.Weight;
  if (b.nfsp < 0 || b.nfsp > kMaxFSP) {
    throw std::runtime_error("FlatTreeCursor: entry " + std::to_string(event_.entry) + " has nfsp=" +
                             std::to_string(b.nfsp) + ", beyond the supported " + std::to_string(kMaxFSP));
  }
  e.nfsp = b.nfsp;
  e.pdg.assign(b.pdg, b.pdg + b.nfsp);
  e.px.assign(b.px, b.px + b.nfsp);
  e.py.assign(b.py, b.py + b.nfsp);
  e.pz.assign(b.pz, b.pz + b.nfsp);
  e.E.assign(b.E, b.E + b.nfsp);
  e.flagCCINC = b.flagCCINC;
  e.flagNCINC = b.flagNCINC;
  e.flagCCQE = b.flagCCQE;
  e.flagCC0pi = b.flagCC0pi;
  e.flagCCQELike = b.flagCCQELike;
  e.flagNCEL = b.flagNCEL;
  e.flagNC0pi = b.flagNC0pi;
  e.flagCCcoh = b.flagCCcoh;
  e.flagNCcoh = b.flagNCcoh;
  e.flagCC1pip = b.flagCC1pip;
  e.flagNC1pip = b.flagNC1pip;
  e.flagCC1pim = b.flagCC1pim;
  e.flagNC1pim = b.flagNC1pim;
  e.flagCC1pi0 = b.flagCC1pi0;
  e.flagNC1pi0 = b.flagNC1pi0;
}

bool FlatTreeCursor::goTo(long long entry) {
  if (entry < 0 || entry >= entries()) return false;
  if (tree_->GetEntry(entry) <= 0) throw std::runtime_error("FlatTreeCursor: cannot read entry " + std::to_string(entry));
  event_.entry = entry;
  copyBuffers();
  next_ = entry + 1;
  return true;
}

bool FlatTreeCursor::next() {
  if (next_ >= entries()) return false;
  return goTo(next_);
}

std::vector<std::string> FlatTreeCursor::branchNames() const {
  std::vector<std::string> out;
  TObjArray* branches = tree_->GetListOfBranches();
  for (int i = 0; i < branches->GetEntries(); ++i) out.push_back(branches->At(i)->GetName());
  return out;
}

bool FlatTreeCursor::hasBranch(const std::string& name) const { return tree_->GetBranch(name.c_str()) != nullptr; }

TH1* FlatTreeCursor::fluxHistogram() const {
  TH1* h = file_->Get<TH1>("FlatTree_FLUX");
  if (!h) return nullptr;
  auto* c = static_cast<TH1*>(h->Clone());
  c->SetDirectory(nullptr);
  return c;
}

TH1* FlatTreeCursor::eventRateHistogram() const {
  TH1* h = file_->Get<TH1>("FlatTree_EVT");
  if (!h) return nullptr;
  auto* c = static_cast<TH1*>(h->Clone());
  c->SetDirectory(nullptr);
  return c;
}

TTree* FlatTreeCursor::cloneTo(TDirectory& dir, const std::string& name) {
  dir.cd();
  TTree* clone = tree_->CloneTree(0);
  if (!clone) throw std::runtime_error("FlatTreeCursor: CloneTree failed for '" + fileName_ + "'");
  if (!name.empty()) clone->SetName(name.c_str());
  clone->SetDirectory(&dir);
  return clone;
}

}  // namespace zaku
