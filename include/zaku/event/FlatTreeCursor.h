#ifndef ZAKU_EVENT_FLATTREECURSOR_H
#define ZAKU_EVENT_FLATTREECURSOR_H
#include <memory>
#include <string>
#include <vector>

#include "zaku/event/FlatEvent.h"

class TDirectory;
class TFile;
class TH1;
class TTree;

namespace zaku {

/// Sequential reader of a NUISANCE flat tree with classic branch addresses,
/// so that an output tree cloned with `cloneTo()` copies every input branch on
/// `Fill()` (TTreeReader cannot do that). Missing required branches are fatal.
class FlatTreeCursor {
public:
  explicit FlatTreeCursor(const std::string& file, const std::string& treeName = "FlatTree_VARS");
  ~FlatTreeCursor();
  FlatTreeCursor(const FlatTreeCursor&) = delete;
  FlatTreeCursor& operator=(const FlatTreeCursor&) = delete;

  const std::string& fileName() const { return fileName_; }
  long long entries() const;
  long long entry() const { return event_.entry; }
  /// Load the next entry; false at the end.
  bool next();
  /// Load a specific entry.
  bool goTo(long long entry);
  /// Position before the first entry, so that next() loads entry 0.
  void rewind() {
    next_ = 0;
    event_.entry = -1;
  }
  const FlatEvent& event() const { return event_; }

  TTree* tree() const { return tree_; }
  TFile* file() const { return file_.get(); }
  std::vector<std::string> branchNames() const;
  bool hasBranch(const std::string& name) const;

  /// Clone of `FlatTree_FLUX` / `FlatTree_EVT` (owned by the caller, detached), or nullptr.
  TH1* fluxHistogram() const;
  TH1* eventRateHistogram() const;

  /// `tree->CloneTree(0)` written into `dir`, sharing this cursor's branch
  /// addresses: after each next(), `clone->Fill()` copies the input event.
  /// New branches can be added to the clone before filling.
  TTree* cloneTo(TDirectory& dir, const std::string& name = "");

  static const std::vector<std::string>& requiredBranches();
  /// Maximum final-state multiplicity the reader accepts.
  static int maxFinalStateParticles();

private:
  struct Buffers;
  void copyBuffers();

  std::string fileName_;
  std::unique_ptr<TFile> file_;
  TTree* tree_ = nullptr;
  std::unique_ptr<Buffers> buf_;
  FlatEvent event_;
  long long next_ = 0;
};

}  // namespace zaku

#endif  // ZAKU_EVENT_FLATTREECURSOR_H
