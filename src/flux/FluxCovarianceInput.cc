#include "zaku/flux/FluxCovarianceInput.h"

#include <TAxis.h>
#include <TFile.h>
#include <TKey.h>
#include <TList.h>

#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>

#include "zaku/core/Units.h"

namespace zaku {

namespace {

std::vector<std::string> split(const std::string& s, char sep) {
  std::vector<std::string> out;
  std::string cur;
  for (char c : s) {
    if (c == sep) {
      out.push_back(cur);
      cur.clear();
    } else {
      cur += c;
    }
  }
  out.push_back(cur);
  return out;
}

}  // namespace

std::string FluxCovarianceBlock::describe() const {
  std::ostringstream os;
  os << detector << " " << toString(beamMode) << " " << NuSpecies::fromPdg(pdg).name() << ": rows [" << firstRow
     << ", " << firstRow + nBins << "), edges";
  for (double e : edges) os << " " << e;
  return os.str();
}

FluxCovarianceInput::FluxCovarianceInput(const FluxCovarianceSpec& spec) : spec_(spec) {
  if (!spec_.valid()) throw std::runtime_error("FluxCovarianceInput: no covariance file given");
  std::unique_ptr<TFile> f(TFile::Open(spec_.file.c_str(), "READ"));
  if (!f || f->IsZombie()) throw std::runtime_error("FluxCovarianceInput: cannot open '" + spec_.file + "'");
  auto* m = f->Get<TMatrixD>(spec_.matrix.c_str());
  if (!m) {
    throw std::runtime_error("FluxCovarianceInput: no TMatrixD '" + spec_.matrix + "' in '" + spec_.file + "'");
  }
  matrix_.ResizeTo(*m);
  matrix_ = *m;

  const double scale = toGeV(1.0, spec_.units);
  std::set<std::string> seen;
  int row = 0;
  TIter next(f->GetListOfKeys());
  while (auto* key = static_cast<TKey*>(next())) {
    if (std::string(key->GetClassName()) != "TAxis") continue;
    if (!seen.insert(key->GetName()).second) continue;
    std::unique_ptr<TAxis> ax(static_cast<TAxis*>(key->ReadObj()));
    const std::string name = key->GetName();
    const std::vector<std::string> tok = split(name, '_');
    if (tok.size() != 4 || tok[3] != "bins") {
      throw std::runtime_error("FluxCovarianceInput: cannot parse TAxis name '" + name +
                               "' (expected <detector>_<numode|anumode>_<species>_bins)");
    }
    FluxCovarianceBlock b;
    b.axisName = name;
    b.detector = tok[0];
    b.beamMode = beamModeFromString(tok[1]);
    b.pdg = NuSpecies::fromName(tok[2]).pdg();
    b.nBins = ax->GetNbins();
    for (int i = 1; i <= b.nBins; ++i) b.edges.push_back(ax->GetBinLowEdge(i) * scale);
    b.edges.push_back(ax->GetBinUpEdge(b.nBins) * scale);
    b.firstRow = row;
    row += b.nBins;
    blocks_.push_back(b);
  }
  if (blocks_.empty()) throw std::runtime_error("FluxCovarianceInput: no TAxis objects in '" + spec_.file + "'");
  if (row != matrix_.GetNrows() || matrix_.GetNrows() != matrix_.GetNcols()) {
    throw std::runtime_error("FluxCovarianceInput: the TAxis blocks cover " + std::to_string(row) +
                             " rows but the matrix is " + std::to_string(matrix_.GetNrows()) + "x" +
                             std::to_string(matrix_.GetNcols()));
  }
}

std::vector<FluxCovarianceBlock> FluxCovarianceInput::blocksFor(const std::string& det) const {
  std::vector<FluxCovarianceBlock> out;
  for (const auto& b : blocks_) {
    if (b.detector == det) out.push_back(b);
  }
  if (out.empty()) throw std::runtime_error("FluxCovarianceInput: no blocks for detector '" + det + "'");
  return out;
}

int FluxCovarianceInput::firstRow(const std::string& det) const { return blocksFor(det).front().firstRow; }

int FluxCovarianceInput::nRows(const std::string& det) const {
  int n = 0;
  for (const auto& b : blocksFor(det)) n += b.nBins;
  return n;
}

TMatrixD FluxCovarianceInput::subMatrix(const std::string& det) const {
  const int r0 = firstRow(det);
  const int n = nRows(det);
  return matrix_.GetSub(r0, r0 + n - 1, r0, r0 + n - 1);
}

TMatrixD FluxCovarianceInput::crossMatrix(const std::string& a, const std::string& b) const {
  const int ra = firstRow(a), na = nRows(a);
  const int rb = firstRow(b), nb = nRows(b);
  return matrix_.GetSub(ra, ra + na - 1, rb, rb + nb - 1);
}

std::string FluxCovarianceInput::describe() const {
  std::ostringstream os;
  os << "flux covariance '" << spec_.file << ":" << spec_.matrix << "' " << matrix_.GetNrows() << "x"
     << matrix_.GetNcols() << "\n";
  for (const auto& b : blocks_) os << "  " << b.describe() << "\n";
  return os.str();
}

}  // namespace zaku
