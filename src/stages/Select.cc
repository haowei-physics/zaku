// The select stage: recon (or spline) files of one beam mode -> the FakeSelect
// tree GUNDAM reads, plus the rate table.
//
// Events are distributed over the experiment's samples from the column
// probabilities the recon stage stored: either one entry per column carrying
// the probability as a weight (split), or one random draw per event against
// the running sum of the columns (draw, the legacy behaviour). The energy
// estimator is evaluated under the assigned sample's lepton hypothesis.
#include <TBranch.h>
#include <TClass.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TTree.h>

#include <memory>
#include <sstream>
#include <stdexcept>

#include "zaku/core/Experiment.h"
#include "zaku/core/Provenance.h"
#include "zaku/detector/Formula.h"
#include "zaku/detector/Reco.h"
#include "zaku/detector/Variables.h"
#include "zaku/event/FlatTreeCursor.h"
#include "zaku/stages/Stages.h"

namespace zaku {
namespace stages {

namespace {

bool fileExists(const std::string& p) { return !p.empty() && !gSystem->AccessPathName(p.c_str()); }
constexpr int kMaxColumns = 64;
constexpr int kNotSelected = -2;

/// "0=a,b,c;1=d,e" -> [[a,b,c],[d,e]]
std::vector<std::vector<std::string>> parseColumns(const std::string& text) {
  std::vector<std::vector<std::string>> out;
  std::istringstream in(text);
  std::string table;
  while (std::getline(in, table, ';')) {
    const std::string::size_type eq = table.find('=');
    if (eq == std::string::npos) continue;
    std::vector<std::string> cols;
    std::istringstream cs(table.substr(eq + 1));
    std::string c;
    while (std::getline(cs, c, ',')) {
      if (!c.empty()) cols.push_back(c);
    }
    out.push_back(cols);
  }
  return out;
}

/// The recon-stage branches the select stage reads from each input file.
struct ReconBuffers {
  Int_t RespCat = -1, true_ME = 0, recon_ME = 0, ProductionIndex = 0, RespTable = -1, nSampleProb = 0;
  Double_t POTWeight = 0, OscWeight = 0, SelEffWeight = 0;
  Float_t PLep_rec = 0, CosLep_rec = 0, Theta_rec = 0;
  Long64_t EventIndex = -1;
  Float_t PLep_true = 0, Ehad_true = 0, Evis_true = 0, W2 = 0;
  Float_t Ehad_rec = 0, Elep_rec_mu = 0, Elep_rec_e = 0, Enu_rec_ccqe_mu = 0, Q2_rec_ccqe_mu = 0, Enu_rec_ccqe_e = 0,
          Q2_rec_ccqe_e = 0, Enu_rec_calo_mu = 0, Enu_rec_calo_e = 0;
  Double_t SampleProb[kMaxColumns] = {};

  void bind(TTree* t) {
    auto b = [&](const char* name, void* addr) {
      if (!t->GetBranch(name)) {
        throw std::runtime_error(std::string("select: input is missing the recon branch '") + name +
                                 "'; rerun the recon stage");
      }
      t->SetBranchAddress(name, addr);
    };
    b("RespCat", &RespCat);
    b("true_ME", &true_ME);
    b("recon_ME", &recon_ME);
    b("ProductionIndex", &ProductionIndex);
    b("RespTable", &RespTable);
    b("nSampleProb", &nSampleProb);
    b("POTWeight", &POTWeight);
    b("OscWeight", &OscWeight);
    b("SelEffWeight", &SelEffWeight);
    b("PLep_rec", &PLep_rec);
    b("CosLep_rec", &CosLep_rec);
    b("Theta_rec", &Theta_rec);
    b("EventIndex", &EventIndex);
    b("PLep_true", &PLep_true);
    b("Ehad_true", &Ehad_true);
    b("Evis_true", &Evis_true);
    b("W2", &W2);
    b("Ehad_rec", &Ehad_rec);
    b("Elep_rec_mu", &Elep_rec_mu);
    b("Elep_rec_e", &Elep_rec_e);
    b("Enu_rec_ccqe_mu", &Enu_rec_ccqe_mu);
    b("Q2_rec_ccqe_mu", &Q2_rec_ccqe_mu);
    b("Enu_rec_ccqe_e", &Enu_rec_ccqe_e);
    b("Q2_rec_ccqe_e", &Q2_rec_ccqe_e);
    b("Enu_rec_calo_mu", &Enu_rec_calo_mu);
    b("Enu_rec_calo_e", &Enu_rec_calo_e);
    b("SampleProb", SampleProb);
  }
};

}  // namespace

const std::vector<std::string>& selectBranches() {
  static const std::vector<std::string> b = {
      // legacy, in the legacy order
      "isRHC", "EventType", "SampleIndex", "flux_pdg", "RespCat", "PDGLep_rec", "recon_ME", "Enu_rec", "Q2_rec",
      "POTWeight", "OscWeight", "SelEffWeight", "PLep_rec", "CosLep_rec", "Theta_rec", "Enu_true", "PDGnu", "cc",
      // Zaku additions
      "Mode", "EventIndex", "ProductionIndex", "isAppearance", "SampleProb", "RespTable", "true_ME", "PLep_true",
      "Ehad_true", "Evis_true", "W2", "Ehad_rec", "Elep_rec", "EhadFrac_rec", "Enu_rec_ccqe_mu", "Enu_rec_ccqe_e",
      "Enu_rec_calo_mu", "Enu_rec_calo_e"};
  return b;
}

unsigned selectSeed(const Experiment& exp, BeamMode mode) {
  if (exp.run().legacyRngOrder) return 20250731u;  // the legacy FakeSelect seed
  const std::size_t h = std::hash<std::string>{}("select#" + toString(mode) + "#" + std::to_string(exp.run().seed));
  unsigned s = static_cast<unsigned>(h & 0x7fffffffu);
  return s == 0 ? 1u : s;
}

std::string SelectSummary::describe() const {
  std::ostringstream os;
  if (skipped) {
    os << "select " << toString(mode) << ": output exists, skipped (" << outputFile << ")";
    return os.str();
  }
  os << "select " << toString(mode) << ": " << entriesIn << " entries in, " << entriesOut << " out, seed " << seed
     << " -> " << outputFile << "\n" << rates.text();
  return os.str();
}

SelectSummary select(const Experiment& exp, BeamMode mode, bool force) {
  SelectSummary summary;
  summary.mode = mode;
  summary.outputFile = exp.path(Layout::selectFile(mode));
  summary.rateTableFile = exp.path(Layout::rateTableFile(mode));
  summary.seed = selectSeed(exp, mode);
  if (fileExists(summary.outputFile) && !force) {
    summary.skipped = true;
    return summary;
  }

  const std::vector<Production> productions = exp.productions(mode);
  if (productions.empty()) throw std::runtime_error("select: no productions for beam mode " + toString(mode));
  const std::vector<int> sampleIdx = exp.sampleIndices(mode);
  if (sampleIdx.empty()) throw std::runtime_error("select: no samples defined for beam mode " + toString(mode));
  const std::vector<SampleDefinition>& samples = exp.samples();

  for (const auto& p : productions) {
    const std::string spline = exp.path(Layout::splineFile(p));
    const std::string recon = exp.path(Layout::reconFile(p));
    if (fileExists(spline)) summary.inputFiles.push_back(spline);
    else if (fileExists(recon)) summary.inputFiles.push_back(recon);
    else throw std::runtime_error("select: no recon file for " + p.name() + " ('" + recon + "'); run the recon stage first");
  }

  // --- formulas ---------------------------------------------------------------
  VariableSet vars;
  std::vector<Formula> cuts(samples.size());
  for (std::size_t i = 0; i < samples.size(); ++i) {
    if (!samples[i].cut.empty()) cuts[i] = Formula(samples[i].cut, vars.names());
  }
  const int iEnu = vars.index("Enu_rec"), iQ2 = vars.index("Q2_rec"), iElep = vars.index("Elep_rec"),
            iEhadFrac = vars.index("EhadFrac_rec"), iPdgLep = vars.index("PDGLep_rec"), iType = vars.index("EventType"),
            iProb = vars.index("SampleProb");

  TRandom3 rng(summary.seed);

  // --- output -----------------------------------------------------------------
  std::unique_ptr<TFile> out(TFile::Open(summary.outputFile.c_str(), "RECREATE"));
  if (!out || out->IsZombie()) throw std::runtime_error("select: cannot create '" + summary.outputFile + "'");
  auto* tree = new TTree("FakeSelect", "sample-assigned fake data");

  Bool_t isRHCb = isRHC(mode), isAppearance = false;
  Int_t EventType = kNotSelected, SampleIndex = 0, flux_pdg = 0, RespCat = -1, PDGLep_rec = 0, recon_ME = 0, PDGnu = 0,
        Mode = 0, ProductionIndex = 0, RespTable = -1, true_ME = 0;
  Char_t cc = 0;
  Float_t Enu_rec = -1, Q2_rec = -1, PLep_rec = 0, CosLep_rec = 0, Theta_rec = 0, Enu_true = 0;
  Double_t POTWeight = 0, OscWeight = 0, SelEffWeight = 0, SampleProb = 1.0;
  Long64_t EventIndex = -1;
  Float_t PLep_true = 0, Ehad_true = 0, Evis_true = 0, W2 = 0, Ehad_rec = 0, Elep_rec = 0, EhadFrac_rec = 0,
          Enu_rec_ccqe_mu = 0, Enu_rec_ccqe_e = 0, Enu_rec_calo_mu = 0, Enu_rec_calo_e = 0;

  tree->Branch("isRHC", &isRHCb, "isRHC/O");
  tree->Branch("EventType", &EventType, "EventType/I");
  tree->Branch("SampleIndex", &SampleIndex, "SampleIndex/I");
  tree->Branch("flux_pdg", &flux_pdg, "flux_pdg/I");
  tree->Branch("RespCat", &RespCat, "RespCat/I");
  tree->Branch("PDGLep_rec", &PDGLep_rec, "PDGLep_rec/I");
  tree->Branch("recon_ME", &recon_ME, "recon_ME/I");
  tree->Branch("Enu_rec", &Enu_rec, "Enu_rec/F");
  tree->Branch("Q2_rec", &Q2_rec, "Q2_rec/F");
  tree->Branch("POTWeight", &POTWeight, "POTWeight/D");
  tree->Branch("OscWeight", &OscWeight, "OscWeight/D");
  tree->Branch("SelEffWeight", &SelEffWeight, "SelEffWeight/D");
  tree->Branch("PLep_rec", &PLep_rec, "PLep_rec/F");
  tree->Branch("CosLep_rec", &CosLep_rec, "CosLep_rec/F");
  tree->Branch("Theta_rec", &Theta_rec, "Theta_rec/F");
  tree->Branch("Enu_true", &Enu_true, "Enu_true/F");
  tree->Branch("PDGnu", &PDGnu, "PDGnu/I");
  tree->Branch("cc", &cc, "cc/B");
  tree->Branch("Mode", &Mode, "Mode/I");
  tree->Branch("EventIndex", &EventIndex, "EventIndex/L");
  tree->Branch("ProductionIndex", &ProductionIndex, "ProductionIndex/I");
  tree->Branch("isAppearance", &isAppearance, "isAppearance/O");
  tree->Branch("SampleProb", &SampleProb, "SampleProb/D");
  tree->Branch("RespTable", &RespTable, "RespTable/I");
  tree->Branch("true_ME", &true_ME, "true_ME/I");
  tree->Branch("PLep_true", &PLep_true, "PLep_true/F");
  tree->Branch("Ehad_true", &Ehad_true, "Ehad_true/F");
  tree->Branch("Evis_true", &Evis_true, "Evis_true/F");
  tree->Branch("W2", &W2, "W2/F");
  tree->Branch("Ehad_rec", &Ehad_rec, "Ehad_rec/F");
  tree->Branch("Elep_rec", &Elep_rec, "Elep_rec/F");
  tree->Branch("EhadFrac_rec", &EhadFrac_rec, "EhadFrac_rec/F");
  tree->Branch("Enu_rec_ccqe_mu", &Enu_rec_ccqe_mu, "Enu_rec_ccqe_mu/F");
  tree->Branch("Enu_rec_ccqe_e", &Enu_rec_ccqe_e, "Enu_rec_ccqe_e/F");
  tree->Branch("Enu_rec_calo_mu", &Enu_rec_calo_mu, "Enu_rec_calo_mu/F");
  tree->Branch("Enu_rec_calo_e", &Enu_rec_calo_e, "Enu_rec_calo_e/F");

  // Response graphs (one TGraph branch per dial, unsplit) are carried through
  // by class; the branch set must be identical in every input file.
  std::vector<std::string> graphNames;
  std::vector<TGraph*> graphBuf;

  // --- rate table --------------------------------------------------------------
  std::vector<std::string> rowNames;
  for (int i : sampleIdx) rowNames.push_back(samples[static_cast<std::size_t>(i)].name);
  rowNames.push_back("rejected");
  std::vector<std::string> colNames;
  for (const auto& p : productions) colNames.push_back(p.name());
  summary.rates = RateTable(rowNames, colNames);
  auto rowOf = [&](int sample) {
    for (std::size_t k = 0; k < sampleIdx.size(); ++k) {
      if (sampleIdx[k] == sample) return static_cast<int>(k);
    }
    return static_cast<int>(sampleIdx.size());  // rejected
  };

  const EnergyEstimator estimator = exp.energyEstimator();
  const Assignment assignment = exp.assignment();

  // --- loop over productions -------------------------------------------------------
  for (std::size_t ip = 0; ip < productions.size(); ++ip) {
    const Production& prod = productions[ip];
    FlatTreeCursor cur(summary.inputFiles[ip]);
    TTree* in = cur.tree();
    Provenance::requireCompatible(*cur.file(), "select: " + summary.inputFiles[ip]);
    ReconBuffers rb;
    rb.bind(in);

    // response columns of this file
    std::vector<std::vector<std::string>> tableColumns;
    if (TDirectory* z = cur.file()->GetDirectory("zaku")) {
      if (auto* n = z->Get<TNamed>("response_columns")) tableColumns = parseColumns(n->GetTitle());
    }
    if (tableColumns.empty()) {
      throw std::runtime_error("select: '" + summary.inputFiles[ip] + "' carries no response column names; rerun recon");
    }
    // column -> sample index (for this beam mode), per table
    std::vector<std::vector<int>> columnSample(tableColumns.size());
    for (std::size_t t = 0; t < tableColumns.size(); ++t) {
      for (const auto& col : tableColumns[t]) {
        int s = -1;
        for (int i : sampleIdx) {
          if (samples[static_cast<std::size_t>(i)].columnName() == col) {
            s = i;
            break;
          }
        }
        columnSample[t].push_back(s);
      }
    }

    // graph branches
    std::vector<std::string> names;
    TObjArray* branches = in->GetListOfBranches();
    for (int ib = 0; ib < branches->GetEntries(); ++ib) {
      auto* br = static_cast<TBranch*>(branches->At(ib));
      TClass* cl = TClass::GetClass(br->GetClassName());
      if (cl && cl->InheritsFrom(TGraph::Class())) names.push_back(br->GetName());
    }
    if (ip == 0) {
      graphNames = names;
      graphBuf.assign(names.size(), nullptr);
      for (std::size_t g = 0; g < graphNames.size(); ++g) {
        tree->Branch(graphNames[g].c_str(), "TGraph", &graphBuf[g], 32000, 0);
      }
    } else if (names != graphNames) {
      throw std::runtime_error("select: the response-graph branches of '" + summary.inputFiles[ip] +
                               "' differ from those of the first input; every production must carry the same dials");
    }
    for (std::size_t g = 0; g < graphNames.size(); ++g) in->SetBranchAddress(graphNames[g].c_str(), &graphBuf[g]);

    ProductionContext ctx;
    ctx.isRHC = isRHC(mode);
    ctx.fluxPdg = prod.fluxPdg;
    ctx.interactionPdg = prod.interactionPdg;
    const int col = static_cast<int>(ip);

    while (cur.next()) {
      ++summary.entriesIn;
      const FlatEvent& ev = cur.event();
      const TruthVariables truth = TruthVariables::compute(ev, ctx, exp.smearing().neutronVisibleFraction());
      vars.loadTruth(ev, truth);
      vars.set("PLep_rec", rb.PLep_rec);
      vars.set("CosLep_rec", rb.CosLep_rec);
      vars.set("Theta_rec", rb.Theta_rec);
      vars.set("Ehad_rec", rb.Ehad_rec);
      vars.set("Elep_rec_mu", rb.Elep_rec_mu);
      vars.set("Elep_rec_e", rb.Elep_rec_e);
      vars.set("Enu_rec_ccqe_mu", rb.Enu_rec_ccqe_mu);
      vars.set("Q2_rec_ccqe_mu", rb.Q2_rec_ccqe_mu);
      vars.set("Enu_rec_ccqe_e", rb.Enu_rec_ccqe_e);
      vars.set("Q2_rec_ccqe_e", rb.Q2_rec_ccqe_e);
      vars.set("Enu_rec_calo_mu", rb.Enu_rec_calo_mu);
      vars.set("Enu_rec_calo_e", rb.Enu_rec_calo_e);
      vars.set("true_ME", rb.true_ME);
      vars.set("recon_ME", rb.recon_ME);
      vars.set("RespCat", rb.RespCat);
      vars.set("SelEffWeight", rb.SelEffWeight);
      vars.set("POTWeight", rb.POTWeight);
      vars.set("OscWeight", rb.OscWeight);

      const std::vector<int>* colSample = (rb.RespTable >= 0 && rb.RespTable < static_cast<int>(columnSample.size()))
                                              ? &columnSample[static_cast<std::size_t>(rb.RespTable)]
                                              : nullptr;
      const int nCols = colSample ? std::min<int>(rb.nSampleProb, static_cast<int>(colSample->size())) : 0;

      // --- assignment: (sample index or -1, probability) pairs -----------------
      std::vector<std::pair<int, double>> entries;
      if (assignment == Assignment::Draw) {
        const double r = rng.Uniform();  // one draw per event, always (legacy)
        double acc = 0.0;
        int chosen = -1;
        for (int c = 0; c < nCols; ++c) {
          acc += rb.SampleProb[c];
          if (r < acc) {
            chosen = (*colSample)[static_cast<std::size_t>(c)];
            break;
          }
        }
        entries.emplace_back(chosen, 1.0);
      } else {
        double assigned = 0.0;
        for (int c = 0; c < nCols; ++c) {
          const int s = (*colSample)[static_cast<std::size_t>(c)];
          const double p = rb.SampleProb[c];
          if (s >= 0 && p > 0.0) {
            entries.emplace_back(s, p);
            assigned += p;
          }
        }
        if (1.0 - assigned > 1.0e-12 || entries.empty()) entries.emplace_back(-1, std::max(0.0, 1.0 - assigned));
      }

      for (auto& e : entries) {
        int sample = e.first;
        const double prob = e.second;
        // The rejected remainder is reconstructed under the electron hypothesis (legacy).
        const Hypothesis hyp = sample >= 0 ? samples[static_cast<std::size_t>(sample)].hypothesis : Hypothesis::Electron;
        RecoObservables ro;
        ro.Enu_rec_ccqe_mu = rb.Enu_rec_ccqe_mu;
        ro.Q2_rec_ccqe_mu = rb.Q2_rec_ccqe_mu;
        ro.Enu_rec_ccqe_e = rb.Enu_rec_ccqe_e;
        ro.Q2_rec_ccqe_e = rb.Q2_rec_ccqe_e;
        ro.Enu_rec_calo_mu = rb.Enu_rec_calo_mu;
        ro.Enu_rec_calo_e = rb.Enu_rec_calo_e;
        ro.Elep_rec_mu = rb.Elep_rec_mu;
        ro.Elep_rec_e = rb.Elep_rec_e;
        Enu_rec = ro.enuRec(hyp, estimator, ev.PDGLep);
        Q2_rec = ro.q2Rec(hyp, ev.PDGLep);
        Elep_rec = static_cast<Float_t>(ro.elepRec(hyp, ev.PDGLep));
        EhadFrac_rec = (Elep_rec + rb.Ehad_rec) > 0.f ? rb.Ehad_rec / (Elep_rec + rb.Ehad_rec) : 0.f;
        PDGLep_rec = ro.pdgLepRec(hyp, ctx.isRHC, ev.PDGLep);
        vars.set(iEnu, Enu_rec);
        vars.set(iQ2, Q2_rec);
        vars.set(iElep, Elep_rec);
        vars.set(iEhadFrac, EhadFrac_rec);
        vars.set(iPdgLep, PDGLep_rec);
        vars.set(iType, sample >= 0 ? sample : kNotSelected);
        vars.set(iProb, prob);
        if (sample >= 0 && !cuts[static_cast<std::size_t>(sample)].empty() &&
            !cuts[static_cast<std::size_t>(sample)].evalBool(vars.data())) {
          sample = -1;  // failed the sample cut: rejected
          vars.set(iType, kNotSelected);
        }

        EventType = sample >= 0 ? sample : kNotSelected;
        SampleIndex = rb.ProductionIndex;
        ProductionIndex = rb.ProductionIndex;
        flux_pdg = prod.fluxPdg;
        isAppearance = ctx.isAppearance();
        RespCat = rb.RespCat;
        RespTable = rb.RespTable;
        recon_ME = rb.recon_ME;
        true_ME = rb.true_ME;
        POTWeight = rb.POTWeight;
        OscWeight = rb.OscWeight;
        SelEffWeight = rb.SelEffWeight;
        SampleProb = prob;
        PLep_rec = rb.PLep_rec;
        CosLep_rec = rb.CosLep_rec;
        Theta_rec = rb.Theta_rec;
        Enu_true = ev.Enu_true;
        PDGnu = ev.PDGnu;
        cc = ev.cc ? 1 : 0;
        Mode = ev.Mode;
        EventIndex = rb.EventIndex;
        PLep_true = rb.PLep_true;
        Ehad_true = rb.Ehad_true;
        Evis_true = rb.Evis_true;
        W2 = rb.W2;
        Ehad_rec = rb.Ehad_rec;
        Enu_rec_ccqe_mu = rb.Enu_rec_ccqe_mu;
        Enu_rec_ccqe_e = rb.Enu_rec_ccqe_e;
        Enu_rec_calo_mu = rb.Enu_rec_calo_mu;
        Enu_rec_calo_e = rb.Enu_rec_calo_e;

        tree->Fill();
        ++summary.entriesOut;
        summary.rates.add(rowOf(sample), col, POTWeight * OscWeight * SelEffWeight * prob);
      }
    }
    for (std::size_t g = 0; g < graphNames.size(); ++g) in->ResetBranchAddress(in->GetBranch(graphNames[g].c_str()));
  }

  // --- trailer -----------------------------------------------------------------
  out->cd();
  tree->Write("", TObject::kOverwrite);
  std::unique_ptr<TH2D> h(summary.rates.histogram("rates"));
  out->WriteTObject(h.get(), "rates", "Overwrite");
  std::ostringstream prov;
  prov << "experiment " << exp.name() << "\nbeam mode " << toString(mode) << "\nassignment " << toString(assignment)
       << "\nenergy estimator " << toString(estimator) << "\nseed " << summary.seed
       << (exp.run().legacyRngOrder ? " (legacy RNG order)" : "") << "\ninputs";
  for (const auto& f : summary.inputFiles) prov << "\n  " << f;
  prov << "\nsamples";
  for (int i : sampleIdx) prov << "\n  [" << i << "] " << samples[static_cast<std::size_t>(i)].describe();
  Provenance::stamp(*out, "select", prov.str());
  out->Close();
  summary.rates.writeCsv(summary.rateTableFile);
  return summary;
}

}  // namespace stages
}  // namespace zaku
