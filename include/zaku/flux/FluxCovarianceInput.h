#ifndef ZAKU_FLUX_FLUXCOVARIANCEINPUT_H
#define ZAKU_FLUX_FLUXCOVARIANCEINPUT_H
#include <TMatrixD.h>

#include <string>
#include <vector>

#include "zaku/core/FluxSpec.h"
#include "zaku/core/Types.h"

namespace zaku {

/// One block of the flux covariance: a detector, a beam mode and a species,
/// with its energy bin edges (GeV) and its row range in the full matrix.
struct FluxCovarianceBlock {
  std::string detector;   ///< "sk", "nd280", ...
  BeamMode beamMode = BeamMode::FHC;
  int pdg = 14;
  int firstRow = 0;
  int nBins = 0;
  std::vector<double> edges;
  std::string axisName;   ///< the TAxis object the block came from

  std::string describe() const;
};

/// Reader for a T2K-style flux covariance file: a fractional TMatrixT<double>
/// plus TAxis objects named `<detector>_<numode|anumode>_<species>_bins`
/// whose order in the file defines the row order of the matrix.
class FluxCovarianceInput {
public:
  explicit FluxCovarianceInput(const FluxCovarianceSpec& spec);

  const FluxCovarianceSpec& spec() const { return spec_; }
  const TMatrixD& matrix() const { return matrix_; }
  const std::vector<FluxCovarianceBlock>& blocks() const { return blocks_; }

  /// Blocks of one detector, in row order.
  std::vector<FluxCovarianceBlock> blocksFor(const std::string& detectorPrefix) const;
  int firstRow(const std::string& detectorPrefix) const;
  int nRows(const std::string& detectorPrefix) const;
  /// The square sub-matrix of one detector.
  TMatrixD subMatrix(const std::string& detectorPrefix) const;
  /// The rectangular block rows=detA, columns=detB.
  TMatrixD crossMatrix(const std::string& detectorA, const std::string& detectorB) const;

  std::string describe() const;

private:
  FluxCovarianceSpec spec_;
  TMatrixD matrix_;
  std::vector<FluxCovarianceBlock> blocks_;
};

}  // namespace zaku

#endif  // ZAKU_FLUX_FLUXCOVARIANCEINPUT_H
