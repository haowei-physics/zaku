// The recon stage: flat tree -> recon file for one production.
//
// Per event it computes the normalisation weight, the oscillation weight, the
// truth variables, the response (category, efficiency, column probabilities),
// the smeared lepton kinematics, the decay-electron tags and both energy
// estimators, and writes everything next to a clone of the input event.
//
// Random numbers: the "legacy stream" draws exactly what the legacy FakeRecon
// macro drew, in its order (Gaus(p), Gaus(cos), one Uniform per decay-electron
// parent), so the legacy production is reproducible; everything new draws
// from a second stream.
#include <TFile.h>
#include <TH1.h>
#include <TNamed.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TTree.h>

#include <memory>
#include <sstream>
#include <stdexcept>

#include "zaku/core/Experiment.h"
#include "zaku/core/Provenance.h"
#include "zaku/detector/Reco.h"
#include "zaku/detector/Variables.h"
#include "zaku/event/FlatTreeCursor.h"
#include "zaku/osc/Oscillator.h"
#include "zaku/stages/Stages.h"

namespace zaku {
namespace stages {

namespace {
bool fileExists(const std::string& p) { return !p.empty() && !gSystem->AccessPathName(p.c_str()); }
constexpr int kMaxColumns = 64;
}  // namespace

const std::vector<std::string>& reconBranches() {
  static const std::vector<std::string> b = {
      // legacy, in the legacy order
      "RespCat", "true_ME", "recon_ME", "POTWeight", "OscWeight", "SelEffWeight", "PLep_rec", "CosLep_rec",
      // Zaku additions
      "EventIndex", "ProductionIndex", "isRHC", "flux_pdg", "isAppearance", "RespTable", "Theta_rec",
      "PLep_true", "Ehad_true", "Evis_true", "W2", "nFS_pip", "nFS_pim", "nFS_pi0", "nFS_p", "nFS_n", "nFS_gamma",
      "nFS_mup", "nFS_mum", "nFS_other", "Ehad_rec", "Elep_rec_mu", "Elep_rec_e", "Enu_rec_ccqe_mu", "Q2_rec_ccqe_mu",
      "Enu_rec_ccqe_e", "Q2_rec_ccqe_e", "Enu_rec_calo_mu", "Enu_rec_calo_e", "nSampleProb", "SampleProb"};
  return b;
}

unsigned reconSeed(const Experiment& exp, const Production& production) {
  // The legacy macro used one constant seed for every production.
  if (exp.run().legacyRngOrder) return exp.run().seed;
  return exp.productionSeed(production);
}

std::string ReconSummary::describe() const {
  std::ostringstream os;
  if (skipped) {
    os << "recon " << production << ": output exists, skipped (" << outputFile << ")";
    return os.str();
  }
  os << "recon " << production << ": " << kept << " of " << entries << " entries kept (" << dropped
     << " without final state), " << sumPOTWeight << " interactions in the fiducial mass, " << sumOscillatedWeight
     << " after oscillation, seed " << seed << " -> " << outputFile;
  return os.str();
}

ReconSummary reconstruct(const Experiment& exp, const Production& production, bool force) {
  ReconSummary summary;
  summary.production = production.name();
  summary.inputFile = exp.path(Layout::flatFile(production));
  summary.outputFile = exp.path(Layout::reconFile(production));
  summary.seed = reconSeed(exp, production);

  if (fileExists(summary.outputFile) && !force) {
    summary.skipped = true;
    return summary;
  }
  if (!fileExists(summary.inputFile)) {
    throw std::runtime_error("recon " + production.name() + ": flat tree '" + summary.inputFile +
                             "' does not exist; run the flatten manifest for this production first");
  }
  if (!exp.hasResponse()) throw std::runtime_error("recon: the experiment has no response tables (setResponse)");
  if (!exp.hasBeamMode(production.beamMode)) {
    throw std::runtime_error("recon: no flux for beam mode " + toString(production.beamMode));
  }
  exp.detector();  // throws if unset

  // --- constants of this production ----------------------------------------
  const BeamMode mode = production.beamMode;
  ProductionContext ctx;
  ctx.isRHC = isRHC(mode);
  ctx.fluxPdg = production.fluxPdg;
  ctx.interactionPdg = production.interactionPdg;
  int productionIndex = 0;
  for (const auto& p : exp.productions(mode)) {
    if (p.name() == production.name()) break;
    ++productionIndex;
  }
  const double fluxPerPot = exp.fluxIntegralPerPOT(mode, production.fluxPdg);
  const double pot = exp.pot(mode);
  const double nNucleons = exp.nNucleons();
  const Oscillator osc = exp.oscillator();

  VariableSet vars;
  Response response = exp.response();
  response.prepare(vars);
  Smearing smearing = exp.smearing();
  smearing.prepare(vars);

  TRandom3 legacyStream(summary.seed);
  TRandom3 extraStream(summary.seed + 1000003u);

  // --- input / output --------------------------------------------------------
  FlatTreeCursor cur(summary.inputFile);
  summary.entries = cur.entries();
  const std::string dir = gSystem->GetDirName(summary.outputFile.c_str()).Data();
  if (!dir.empty()) gSystem->mkdir(dir.c_str(), true);
  std::unique_ptr<TFile> out(TFile::Open(summary.outputFile.c_str(), "RECREATE"));
  if (!out || out->IsZombie()) throw std::runtime_error("recon: cannot create '" + summary.outputFile + "'");
  TTree* tree = cur.cloneTo(*out);

  Int_t RespCat = -1, true_ME = 0, recon_ME = 0, ProductionIndex = productionIndex, flux_pdg = production.fluxPdg;
  Int_t RespTable = -1, nSampleProb = 0;
  Double_t POTWeight = 0, OscWeight = 0, SelEffWeight = 0;
  Float_t PLep_rec = 0, CosLep_rec = 0, Theta_rec = 0;
  Long64_t EventIndex = -1;
  Bool_t isRHCb = ctx.isRHC, isAppearance = ctx.isAppearance();
  Float_t PLep_true = 0, Ehad_true = 0, Evis_true = 0, W2 = 0;
  Int_t nFS_pip = 0, nFS_pim = 0, nFS_pi0 = 0, nFS_p = 0, nFS_n = 0, nFS_gamma = 0, nFS_mup = 0, nFS_mum = 0, nFS_other = 0;
  Float_t Ehad_rec = 0, Elep_rec_mu = 0, Elep_rec_e = 0, Enu_rec_ccqe_mu = 0, Q2_rec_ccqe_mu = 0, Enu_rec_ccqe_e = 0,
          Q2_rec_ccqe_e = 0, Enu_rec_calo_mu = 0, Enu_rec_calo_e = 0;
  Double_t SampleProb[kMaxColumns] = {};

  tree->Branch("RespCat", &RespCat, "RespCat/I");
  tree->Branch("true_ME", &true_ME, "true_ME/I");
  tree->Branch("recon_ME", &recon_ME, "recon_ME/I");
  tree->Branch("POTWeight", &POTWeight, "POTWeight/D");
  tree->Branch("OscWeight", &OscWeight, "OscWeight/D");
  tree->Branch("SelEffWeight", &SelEffWeight, "SelEffWeight/D");
  tree->Branch("PLep_rec", &PLep_rec, "PLep_rec/F");
  tree->Branch("CosLep_rec", &CosLep_rec, "CosLep_rec/F");
  tree->Branch("EventIndex", &EventIndex, "EventIndex/L");
  tree->Branch("ProductionIndex", &ProductionIndex, "ProductionIndex/I");
  tree->Branch("isRHC", &isRHCb, "isRHC/O");
  tree->Branch("flux_pdg", &flux_pdg, "flux_pdg/I");
  tree->Branch("isAppearance", &isAppearance, "isAppearance/O");
  tree->Branch("RespTable", &RespTable, "RespTable/I");
  tree->Branch("Theta_rec", &Theta_rec, "Theta_rec/F");
  tree->Branch("PLep_true", &PLep_true, "PLep_true/F");
  tree->Branch("Ehad_true", &Ehad_true, "Ehad_true/F");
  tree->Branch("Evis_true", &Evis_true, "Evis_true/F");
  tree->Branch("W2", &W2, "W2/F");
  tree->Branch("nFS_pip", &nFS_pip, "nFS_pip/I");
  tree->Branch("nFS_pim", &nFS_pim, "nFS_pim/I");
  tree->Branch("nFS_pi0", &nFS_pi0, "nFS_pi0/I");
  tree->Branch("nFS_p", &nFS_p, "nFS_p/I");
  tree->Branch("nFS_n", &nFS_n, "nFS_n/I");
  tree->Branch("nFS_gamma", &nFS_gamma, "nFS_gamma/I");
  tree->Branch("nFS_mup", &nFS_mup, "nFS_mup/I");
  tree->Branch("nFS_mum", &nFS_mum, "nFS_mum/I");
  tree->Branch("nFS_other", &nFS_other, "nFS_other/I");
  tree->Branch("Ehad_rec", &Ehad_rec, "Ehad_rec/F");
  tree->Branch("Elep_rec_mu", &Elep_rec_mu, "Elep_rec_mu/F");
  tree->Branch("Elep_rec_e", &Elep_rec_e, "Elep_rec_e/F");
  tree->Branch("Enu_rec_ccqe_mu", &Enu_rec_ccqe_mu, "Enu_rec_ccqe_mu/F");
  tree->Branch("Q2_rec_ccqe_mu", &Q2_rec_ccqe_mu, "Q2_rec_ccqe_mu/F");
  tree->Branch("Enu_rec_ccqe_e", &Enu_rec_ccqe_e, "Enu_rec_ccqe_e/F");
  tree->Branch("Q2_rec_ccqe_e", &Q2_rec_ccqe_e, "Q2_rec_ccqe_e/F");
  tree->Branch("Enu_rec_calo_mu", &Enu_rec_calo_mu, "Enu_rec_calo_mu/F");
  tree->Branch("Enu_rec_calo_e", &Enu_rec_calo_e, "Enu_rec_calo_e/F");
  tree->Branch("nSampleProb", &nSampleProb, "nSampleProb/I");
  tree->Branch("SampleProb", SampleProb, "SampleProb[nSampleProb]/D");

  // --- event loop ------------------------------------------------------------
  while (cur.next()) {
    const FlatEvent& ev = cur.event();
    if (ev.legacyDrop()) {
      ++summary.dropped;
      continue;
    }
    ++summary.kept;
    EventIndex = ev.entry;

    POTWeight = exp.potWeight(mode, production.fluxPdg, ev.fScaleFactor);
    // (identical to fScaleFactor * fluxPerPot * pot * nNucleons in the legacy order)
    (void)fluxPerPot;
    (void)pot;
    (void)nNucleons;
    OscWeight = osc.oscWeight(ev, ctx);

    const TruthVariables truth = TruthVariables::compute(ev, ctx, smearing.neutronVisibleFraction());
    vars.loadTruth(ev, truth);

    const ResponseResult rr = response.evaluate(vars, mode);
    RespCat = rr.category;
    RespTable = rr.table;
    SelEffWeight = rr.efficiency;
    if (rr.probabilities.size() > static_cast<std::size_t>(kMaxColumns)) {
      throw std::runtime_error("recon: a response table has more than 64 columns");
    }
    nSampleProb = static_cast<Int_t>(rr.probabilities.size());
    for (std::size_t c = 0; c < rr.probabilities.size(); ++c) SampleProb[c] = rr.probabilities[c];

    const RecoObservables ro =
        reco::reconstruct(ev, truth, vars, smearing, ctx.isRHC, legacyStream, extraStream, exp.bindingEnergy());
    PLep_rec = ro.PLep_rec;
    CosLep_rec = ro.CosLep_rec;
    Theta_rec = ro.Theta_rec;
    true_ME = ro.true_ME;
    recon_ME = ro.recon_ME;
    Ehad_rec = static_cast<Float_t>(ro.Ehad_rec);
    Elep_rec_mu = static_cast<Float_t>(ro.Elep_rec_mu);
    Elep_rec_e = static_cast<Float_t>(ro.Elep_rec_e);
    Enu_rec_ccqe_mu = ro.Enu_rec_ccqe_mu;
    Q2_rec_ccqe_mu = ro.Q2_rec_ccqe_mu;
    Enu_rec_ccqe_e = ro.Enu_rec_ccqe_e;
    Q2_rec_ccqe_e = ro.Q2_rec_ccqe_e;
    Enu_rec_calo_mu = static_cast<Float_t>(ro.Enu_rec_calo_mu);
    Enu_rec_calo_e = static_cast<Float_t>(ro.Enu_rec_calo_e);

    PLep_true = static_cast<Float_t>(truth.PLep_true);
    Ehad_true = static_cast<Float_t>(truth.Ehad_true);
    Evis_true = static_cast<Float_t>(truth.Evis_true);
    W2 = static_cast<Float_t>(truth.W2);
    nFS_pip = truth.nFS_pip;
    nFS_pim = truth.nFS_pim;
    nFS_pi0 = truth.nFS_pi0;
    nFS_p = truth.nFS_p;
    nFS_n = truth.nFS_n;
    nFS_gamma = truth.nFS_gamma;
    nFS_mup = truth.nFS_mup;
    nFS_mum = truth.nFS_mum;
    nFS_other = truth.nFS_other;

    summary.sumPOTWeight += POTWeight;
    summary.sumOscillatedWeight += POTWeight * OscWeight;
    tree->Fill();
  }

  if (summary.kept > 0 && (summary.sumPOTWeight <= 0.0 || summary.sumOscillatedWeight <= 0.0)) {
    throw std::runtime_error("recon " + production.name() + ": the weights sum to zero; a weight was never assigned");
  }

  // --- trailer -----------------------------------------------------------------
  out->cd();
  std::unique_ptr<TH1> flux(cur.fluxHistogram());
  std::unique_ptr<TH1> evt(cur.eventRateHistogram());
  if (flux) out->WriteTObject(flux.get(), "FlatTree_FLUX", "Overwrite");
  if (evt) out->WriteTObject(evt.get(), "FlatTree_EVT", "Overwrite");
  tree->Write("", TObject::kOverwrite);

  // Column names of every response table, so the select stage can map them to samples.
  std::ostringstream cols;
  for (std::size_t i = 0; i < response.tables().size(); ++i) {
    if (i) cols << ";";
    cols << i << "=";
    const auto& t = response.tables()[i];
    for (std::size_t c = 0; c < t.columns.size(); ++c) cols << (c ? "," : "") << t.columns[c];
  }
  std::ostringstream prov;
  prov << "experiment " << exp.name() << "\nproduction " << production.describe() << "\ninput " << summary.inputFile
       << "\nseed " << summary.seed << (exp.run().legacyRngOrder ? " (legacy RNG order)" : "") << "\nresponse columns "
       << cols.str() << "\nsmearing " << smearing.describe() << "\nenergy estimator " << toString(exp.energyEstimator())
       << "\noscillation " << osc.parameters().describe() << ", baseline " << osc.baseline_km() << " km, rhoYe "
       << osc.rhoYe();
  Provenance::stamp(*out, "recon", prov.str());
  TDirectory* z = out->GetDirectory("zaku");
  TNamed columns("response_columns", cols.str().c_str());
  z->WriteTObject(&columns, "response_columns", "Overwrite");
  out->Close();
  return summary;
}

}  // namespace stages
}  // namespace zaku
