// Round trip through the flat-tree cursor: a synthetic NUISANCE-style tree is
// written, read back, cloned with an extra branch, and the clone must carry
// every input branch (including ones the cursor never bound) plus the new one.
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

#include <memory>

#include "ZakuTest.h"
#include "ZakuTestPaths.h"
#include "zaku/core/Provenance.h"
#include "zaku/core/Version.h"
#include "zaku/event/FlatTreeCursor.h"

using namespace zaku;

namespace {

void writeSyntheticFlatTree(const std::string& file, int nEntries) {
  TFile f(file.c_str(), "RECREATE");
  TTree t("FlatTree_VARS", "synthetic");
  Int_t Mode, PDGnu, PDGLep, nfsp, tgt = 1000080160;
  Char_t cc;
  Float_t Enu_true, ELep, CosLep, Q2, q0, W, x, y, Enu_QE, Q2_QE, q3;
  Double_t fScaleFactor = 3.0e-43;
  Bool_t flagCCQE;
  const int kMax = 8;
  Int_t pdg[kMax];
  Float_t px[kMax], py[kMax], pz[kMax], E[kMax];
  t.Branch("Mode", &Mode, "Mode/I");
  t.Branch("cc", &cc, "cc/B");
  t.Branch("PDGnu", &PDGnu, "PDGnu/I");
  t.Branch("Enu_true", &Enu_true, "Enu_true/F");
  t.Branch("tgt", &tgt, "tgt/I");
  t.Branch("PDGLep", &PDGLep, "PDGLep/I");
  t.Branch("ELep", &ELep, "ELep/F");
  t.Branch("CosLep", &CosLep, "CosLep/F");
  t.Branch("Q2", &Q2, "Q2/F");
  t.Branch("q0", &q0, "q0/F");   // never bound by the cursor: must still be cloned
  t.Branch("q3", &q3, "q3/F");
  t.Branch("W", &W, "W/F");
  t.Branch("x", &x, "x/F");
  t.Branch("y", &y, "y/F");
  t.Branch("Enu_QE", &Enu_QE, "Enu_QE/F");
  t.Branch("Q2_QE", &Q2_QE, "Q2_QE/F");
  t.Branch("fScaleFactor", &fScaleFactor, "fScaleFactor/D");
  t.Branch("nfsp", &nfsp, "nfsp/I");
  t.Branch("px", px, "px[nfsp]/F");
  t.Branch("py", py, "py[nfsp]/F");
  t.Branch("pz", pz, "pz[nfsp]/F");
  t.Branch("E", E, "E[nfsp]/F");
  t.Branch("pdg", pdg, "pdg[nfsp]/I");
  t.Branch("flagCCQE", &flagCCQE, "flagCCQE/O");
  for (int i = 0; i < nEntries; ++i) {
    Mode = (i % 3 == 0) ? 1 : 11;
    cc = 1;
    PDGnu = 14;
    Enu_true = 0.5f + 0.01f * i;
    PDGLep = 13;
    ELep = 0.3f + 0.005f * i;
    CosLep = 0.9f;
    Q2 = 0.1f;
    q0 = 100.0f + i;
    q3 = 0.2f;
    W = 1.1f;
    x = 0.3f;
    y = 0.4f;
    Enu_QE = Enu_true;
    Q2_QE = Q2;
    flagCCQE = (Mode == 1);
    nfsp = (i == 4) ? 0 : 2;  // one empty final state
    pdg[0] = 13;
    pdg[1] = 2212;
    for (int k = 0; k < 2; ++k) {
      px[k] = 0.0f;
      py[k] = 0.0f;
      pz[k] = 0.3f;
      E[k] = (k == 0) ? ELep : 1.0f;
    }
    t.Fill();
  }
  t.Write();
  TH1D flux("FlatTree_FLUX", "flux", 10, 0.0, 10.0);
  flux.SetBinContent(3, 42.0);
  flux.Write();
  TH1D evt("FlatTree_EVT", "evt", 10, 0.0, 10.0);
  evt.Write();
}

}  // namespace

int main() {
  const std::string scratch = zakutest::scratchDir();
  const std::string in = scratch + "/synthetic_flat.root";
  const std::string out = scratch + "/synthetic_clone.root";
  writeSyntheticFlatTree(in, 10);

  {
    FlatTreeCursor cur(in);
    ZAKU_CHECK_EQ(static_cast<int>(cur.entries()), 10);
    ZAKU_CHECK(cur.hasBranch("q0"));
    ZAKU_CHECK(!cur.hasBranch("nope"));
    ZAKU_CHECK_EQ(static_cast<int>(cur.branchNames().size()), 24);

    std::unique_ptr<TH1> flux(cur.fluxHistogram());
    ZAKU_CHECK(flux != nullptr);
    ZAKU_CHECK_CLOSE(flux->GetBinContent(3), 42.0, 1e-12);

    std::unique_ptr<TFile> fout(TFile::Open(out.c_str(), "RECREATE"));
    TTree* clone = cur.cloneTo(*fout);
    Double_t extra = 0.0;
    clone->Branch("POTWeight", &extra, "POTWeight/D");

    int n = 0, empty = 0;
    while (cur.next()) {
      const FlatEvent& ev = cur.event();
      ZAKU_CHECK_EQ(static_cast<int>(ev.entry), n);
      if (ev.emptyFinalState()) ++empty;
      ZAKU_CHECK_CLOSE(ev.Enu_true, 0.5 + 0.01 * n, 1e-6);
      ZAKU_CHECK(ev.cc);
      ZAKU_CHECK_EQ(ev.PDGLep, 13);
      if (!ev.emptyFinalState()) ZAKU_CHECK_EQ(ev.pdg[1], 2212);
      ZAKU_CHECK_CLOSE(ev.q0, 100.0 + n, 1e-6);  // optional branch bound
      extra = 0.5 * n;
      clone->Fill();
      ++n;
    }
    ZAKU_CHECK_EQ(n, 10);
    ZAKU_CHECK_EQ(empty, 1);
    ZAKU_CHECK(cur.goTo(7));
    ZAKU_CHECK_EQ(static_cast<int>(cur.entry()), 7);
    ZAKU_CHECK(!cur.goTo(10));
    fout->cd();
    clone->Write();
    Provenance::stamp(*fout, "test", "synthetic");
    fout->Close();
  }

  // The clone must have every input branch with the input values, plus the new one.
  {
    std::unique_ptr<TFile> f(TFile::Open(out.c_str(), "READ"));
    TTree* t = f->Get<TTree>("FlatTree_VARS");
    ZAKU_CHECK(t != nullptr);
    ZAKU_CHECK_EQ(static_cast<int>(t->GetEntries()), 10);
    ZAKU_CHECK(t->GetBranch("q0") != nullptr);
    ZAKU_CHECK(t->GetBranch("POTWeight") != nullptr);
    Float_t q0 = 0, Enu = 0;
    Double_t pw = 0;
    Int_t nfsp = 0;
    t->SetBranchAddress("q0", &q0);
    t->SetBranchAddress("Enu_true", &Enu);
    t->SetBranchAddress("POTWeight", &pw);
    t->SetBranchAddress("nfsp", &nfsp);
    t->GetEntry(6);
    ZAKU_CHECK_CLOSE(q0, 106.0, 1e-6);
    ZAKU_CHECK_CLOSE(Enu, 0.56, 1e-6);
    ZAKU_CHECK_CLOSE(pw, 3.0, 1e-12);
    ZAKU_CHECK_EQ(nfsp, 2);
    t->GetEntry(4);
    ZAKU_CHECK_EQ(nfsp, 0);
    ZAKU_CHECK(f->Get<TH1>("FlatTree_FLUX") == nullptr);  // not copied by cloneTo itself

    const Provenance p = Provenance::read(*f);
    ZAKU_CHECK_STR_EQ(p.version, version());
    ZAKU_CHECK_EQ(p.schema, schemaVersion());
    ZAKU_CHECK_STR_EQ(p.stage, "test");
    ZAKU_CHECK_STR_EQ(p.text, "synthetic");
    ZAKU_CHECK_NOTHROW(Provenance::requireCompatible(*f, "test"));
  }

  // A file without the required branches is rejected.
  {
    TFile f((scratch + "/bad_flat.root").c_str(), "RECREATE");
    TTree t("FlatTree_VARS", "bad");
    Int_t mode = 0;
    t.Branch("Mode", &mode, "Mode/I");
    t.Fill();
    t.Write();
  }
  ZAKU_CHECK_THROWS(FlatTreeCursor bad(scratch + "/bad_flat.root"));
  ZAKU_CHECK_THROWS(FlatTreeCursor missing(scratch + "/no_such_file.root"));

  return zakutest::finish();
}
