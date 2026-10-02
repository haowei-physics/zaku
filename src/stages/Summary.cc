// The summary stage: the +-1 sigma effect of every parameter on every
// sample's expected event count, from the FakeSelect files.
#include <TSystem.h>

#include <sstream>
#include <stdexcept>

#include "zaku/core/Experiment.h"
#include "zaku/stages/Stages.h"
#include "zaku/syst/Summary.h"

namespace zaku {
namespace stages {

namespace {
bool fileExists(const std::string& p) { return !p.empty() && !gSystem->AccessPathName(p.c_str()); }
}  // namespace

std::string SummaryResult::describe() const {
  std::ostringstream os;
  if (skipped) {
    os << "summary: output exists, skipped (" << csvFile << ")";
    return os.str();
  }
  os << "summary: " << table.parameters().size() << " parameters x " << table.samples().size() << " samples from "
     << entries << " FakeSelect entries -> " << csvFile << ", " << rootFile << "\n" << table.text();
  return os.str();
}

SummaryResult summarise(const Experiment& exp, bool force) {
  SummaryResult result;
  result.csvFile = exp.path(Layout::summaryCsvFile());
  result.rootFile = exp.path(Layout::summaryRootFile());
  if (fileExists(result.csvFile) && fileExists(result.rootFile) && !force) {
    result.skipped = true;
    result.table = SummaryTable::readCsv(result.csvFile);
    return result;
  }
  std::vector<std::pair<BeamMode, std::string>> inputs;
  for (BeamMode mode : exp.beamModes()) {
    const std::string f = exp.path(Layout::selectFile(mode));
    if (!fileExists(f)) throw std::runtime_error("summary: '" + f + "' does not exist; run the select stage first");
    inputs.emplace_back(mode, f);
  }
  if (inputs.empty()) throw std::runtime_error("summary: the experiment has no beam modes");

  const ParameterSets sets = buildParameterSets(exp);
  SystematicsSummary sum(exp, sets.xsec, sets.fluxConstrained, sets.det, sets.osc);
  for (const auto& in : inputs) sum.process(in.first, in.second);
  result.table = sum.finish();
  result.entries = sum.entriesProcessed();

  gSystem->mkdir(exp.path(Layout::summaryDirectory()).c_str(), true);
  result.table.writeCsv(result.csvFile);
  sum.writeHistograms(result.rootFile);
  return result;
}

}  // namespace stages
}  // namespace zaku
