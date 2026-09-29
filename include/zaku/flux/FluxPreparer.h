#ifndef ZAKU_FLUX_FLUXPREPARER_H
#define ZAKU_FLUX_FLUXPREPARER_H
#include <string>

#include "zaku/core/FluxSpec.h"

class TH1;
class TH1D;

namespace zaku {

/// Generator-ready copies of the flux histograms.
///
/// NEUT (and GENIE with NOWIDTH) sample the energy from the raw bin contents,
/// so a variable-width density histogram would under-weight the wide bins.
/// `makeUniform` rebins every histogram to uniform bins of `widthGeV`, taking
/// for each new bin the density of the original bin that contains its centre.
/// The density, and therefore Integral("width"), is preserved exactly. This
/// is what the legacy `*_uniform.root` files are.
namespace FluxPreparer {

/// Rebin `h` to uniform bins of `widthGeV` over its full range (GeV axis).
/// Caller owns the result. `name` defaults to the input name.
TH1D* rebinUniform(const TH1& h, double widthGeV, const std::string& name = "");

/// Write uniform-bin copies of the spec's histograms (or of every TH1 in the
/// file when `allHistograms`) into `outFile` and return a FluxSpec pointing
/// at it. Existing output files are overwritten.
FluxSpec makeUniform(const FluxSpec& spec, const std::string& outFile, double widthGeV = 0.01,
                     bool allHistograms = true);

}  // namespace FluxPreparer

}  // namespace zaku

#endif  // ZAKU_FLUX_FLUXPREPARER_H
