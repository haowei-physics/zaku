// The splines stage: recon file + weight file -> spline file for one
// production. A copy of the recon tree with one unsplit TGraph branch per
// dial, each graph the event's weight as a function of the dial value (see
// SplineBuilder for the rules). GUNDAM reads the graph as a TObject from the
// tree it takes events from, so the graphs must live in this tree.
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TNamed.h>
#include <TSystem.h>
#include <TTree.h>

#include <memory>
#include <sstream>
#include <stdexcept>

#include "zaku/core/Experiment.h"
#include "zaku/core/Provenance.h"
#include "zaku/event/EventJoin.h"
#include "zaku/event/FlatTreeCursor.h"
#include "zaku/event/WeightCursor.h"
#include "zaku/stages/Stages.h"
#include "zaku/syst/Splines.h"

namespace zaku {
namespace stages {

namespace {
bool fileExists(const std::string& p) { return !p.empty() && !gSystem->AccessPathName(p.c_str()); }
}  // namespace

std::string SplineSummary::describe() const {
  std::ostringstream os;
  if (skipped) {
    os << "splines " << production << ": output exists, skipped (" << outputFile << ")";
    return os.str();
  }
  long long rejected = 0;
  for (long long r : nRejected) rejected += r;
  os << "splines " << production << ": " << entries << " entries, " << dials.size() << " dials -> " << outputFile;
  if (rejected) os << "\n  " << rejected << " dial-event responses were unusable and written as empty graphs (weight 1)";
  os << "\n" << stats;
  return os.str();
}

SplineSummary makeSplines(const Experiment& exp, const Production& production, bool force) {
  SplineSummary summary;
  summary.production = production.name();
  summary.inputFile = exp.path(Layout::reconFile(production));
  summary.weightFile = exp.path(Layout::weightFile(production));
  summary.outputFile = exp.path(Layout::splineFile(production));
  if (fileExists(summary.outputFile) && !force) {
    summary.skipped = true;
    return summary;
  }
  if (!fileExists(summary.inputFile)) {
    throw std::runtime_error("splines " + production.name() + ": recon file '" + summary.inputFile +
                             "' does not exist; run the recon stage first");
  }
  if (!fileExists(summary.weightFile)) {
    throw std::runtime_error("splines " + production.name() + ": weight file '" + summary.weightFile +
                             "' does not exist; run the weight manifest for this production first");
  }

  // --- inputs -------------------------------------------------------------------
  FlatTreeCursor cur(summary.inputFile);
  Provenance::requireCompatible(*cur.file(), "splines: " + summary.inputFile);
  std::unique_ptr<WeightSource> weights;
  if (T2KRWWeightCursor::looksLikeT2KRW(summary.weightFile)) {
    weights.reset(new T2KRWWeightCursor(summary.weightFile));
  } else {
    throw std::runtime_error("splines: '" + summary.weightFile + "' is not a T2KReWeight weight file (T2KRWFlatTree + TweakMetaData)");
  }

  const SystematicsConfig& cfg = exp.systematics();
  std::unique_ptr<SplineBuilder> builder;
  if (cfg.xsec.writeUnusedDials || cfg.xsec.dials.empty()) {
    builder.reset(new SplineBuilder(weights->dials()));
  } else {
    builder.reset(new SplineBuilder(weights->dials(), cfg.xsec.fitted()));
  }
  for (const auto& name : builder->dialNames()) {
    if (cur.hasBranch(name)) {
      throw std::runtime_error("splines: dial '" + name + "' collides with a branch of '" + summary.inputFile +
                               "'; GUNDAM resolves leaves by name, so the output would read back wrong");
    }
  }
  // every fitted dial must be available in the weight file
  for (const auto& name : cfg.xsec.fitted()) {
    if (weights->dialIndex(name) < 0) {
      throw std::runtime_error("splines: fitted dial '" + name + "' is not in '" + summary.weightFile +
                               "'; the dial catalogue and the reweight card disagree");
    }
  }

  // --- output -------------------------------------------------------------------
  const std::string dir = gSystem->GetDirName(summary.outputFile.c_str()).Data();
  if (!dir.empty()) gSystem->mkdir(dir.c_str(), true);
  std::unique_ptr<TFile> out(TFile::Open(summary.outputFile.c_str(), "RECREATE"));
  if (!out || out->IsZombie()) throw std::runtime_error("splines: cannot create '" + summary.outputFile + "'");
  TTree* tree = cur.cloneTo(*out);
  std::vector<TGraph*> graphs(builder->size(), nullptr);
  for (std::size_t d = 0; d < builder->size(); ++d) {
    graphs[d] = new TGraph();
    // Split level 0: the graph must reach GUNDAM as one object.
    tree->Branch(builder->dials()[d].name.c_str(), "TGraph", &graphs[d], 32000, 0);
  }

  // --- event loop -----------------------------------------------------------------
  EventJoin join(cur, weights.get());
  while (join.next()) {
    if (!join.hasWeights()) {
      throw std::runtime_error("splines " + production.name() + ": recon entry " + std::to_string(cur.entry()) +
                               " has no weight row in '" + summary.weightFile + "'");
    }
    for (std::size_t d = 0; d < builder->size(); ++d) {
      builder->build(static_cast<int>(d), weights->weights(builder->sourceIndex(static_cast<int>(d))), *graphs[d]);
    }
    tree->Fill();
    ++summary.entries;
  }
  join.finish();

  // --- trailer -----------------------------------------------------------------------
  out->cd();
  std::unique_ptr<TH1> flux(cur.fluxHistogram());
  std::unique_ptr<TH1> evt(cur.eventRateHistogram());
  if (flux) out->WriteTObject(flux.get(), "FlatTree_FLUX", "Overwrite");
  if (evt) out->WriteTObject(evt.get(), "FlatTree_EVT", "Overwrite");
  tree->Write("", TObject::kOverwrite);
  std::string columns;
  if (TDirectory* z = cur.file()->GetDirectory("zaku")) {
    if (auto* n = z->Get<TNamed>("response_columns")) columns = n->GetTitle();
  }
  std::ostringstream prov;
  prov << "experiment " << exp.name() << "\nproduction " << production.describe() << "\nrecon " << summary.inputFile
       << "\nweights " << summary.weightFile << "\n" << join.stats().describe() << "\ndials";
  for (const auto& d : builder->dials()) prov << "\n  " << d.describe();
  Provenance::stamp(*out, "splines", prov.str());
  TDirectory* z = out->GetDirectory("zaku");
  TNamed dialsMeta("dials", SplineBuilder::metadataText(weights->dials()).c_str());
  z->WriteTObject(&dialsMeta, "dials", "Overwrite");
  if (!columns.empty()) {
    TNamed cols("response_columns", columns.c_str());
    z->WriteTObject(&cols, "response_columns", "Overwrite");
  }
  out->Close();
  for (auto* g : graphs) delete g;

  summary.dials = builder->dialNames();
  for (const auto& d : builder->dials()) {
    summary.nFlat.push_back(d.nFlat);
    summary.nRejected.push_back(d.nRejected);
  }
  summary.stats = builder->statsText(summary.entries);
  return summary;
}

}  // namespace stages
}  // namespace zaku
