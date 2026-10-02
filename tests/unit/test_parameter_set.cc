// ParameterSet: construction, correlations, positive-definiteness, the ROOT
// file layout GUNDAM reads (TMatrixD, TObjArray of names) and the round trip.
#include <TFile.h>
#include <TKey.h>
#include <TMatrixD.h>
#include <TObjArray.h>

#include <memory>

#include "ZakuTest.h"
#include "ZakuTestPaths.h"
#include "zaku/syst/ParameterSet.h"

using namespace zaku;

int main() {
  const std::string scratch = zakutest::scratchDir();

  ParameterSet s("det");
  Parameter a;
  a.name = "A";
  a.prior = 1.0;
  a.sigma = 0.02;
  Parameter b = a;
  b.name = "B";
  b.sigma = 0.03;
  Parameter c = a;
  c.name = "C";
  c.sigma = 0.1;
  c.hasBounds = true;
  c.lower = 0.0;
  c.upper = 2.0;
  s.add(a).add(b).add(c);
  ZAKU_CHECK_EQ(static_cast<int>(s.size()), 3);
  ZAKU_CHECK_EQ(s.index("B"), 1);
  ZAKU_CHECK_EQ(s.index("nope"), -1);
  ZAKU_CHECK_CLOSE(s.sigma(1), 0.03, 1e-12);
  ZAKU_CHECK_CLOSE(s.covariance()(2, 2), 0.01, 1e-12);
  ZAKU_CHECK(s.anyBounds());
  ZAKU_CHECK_THROWS(s.add(a));  // duplicate name
  Parameter neg = a;
  neg.name = "N";
  neg.sigma = -1;
  ZAKU_CHECK_THROWS(s.add(neg));

  // correlation and positive definiteness
  s.setCorrelation(0, 1, 0.5);
  ZAKU_CHECK_CLOSE(s.covariance()(0, 1), 0.5 * 0.02 * 0.03, 1e-12);
  ZAKU_CHECK_CLOSE(s.covariance()(1, 0), 0.5 * 0.02 * 0.03, 1e-12);
  ZAKU_CHECK_CLOSE(s.correlation(0, 1), 0.5, 1e-12);
  double lo = 0, hi = 0;
  ZAKU_CHECK(s.isPositiveDefinite(&lo, &hi));
  ZAKU_CHECK(lo > 0 && hi >= lo);
  ParameterSet bad = s;
  bad.setCorrelation(0, 1, 1.0);
  bad.setCorrelation(0, 2, 1.0);
  bad.setCorrelation(1, 2, -1.0);
  ZAKU_CHECK(!bad.isPositiveDefinite());
  ZAKU_CHECK_THROWS(s.setCorrelation(0, 1, 1.5));

  // a non-symmetric matrix is refused; a wrong size too
  TMatrixD m(3, 3);
  m(0, 0) = 1;
  m(1, 1) = 1;
  m(2, 2) = 1;
  m(0, 1) = 0.1;
  ZAKU_CHECK_THROWS(s.setCovariance(m));
  m(1, 0) = 0.1;
  ZAKU_CHECK_NOTHROW(s.setCovariance(m));
  ZAKU_CHECK_CLOSE(s.sigma(0), 1.0, 1e-12);  // sigma follows the matrix
  TMatrixD wrong(2, 2);
  ZAKU_CHECK_THROWS(s.setCovariance(wrong));
  s.setCovariance(TMatrixD(3, 3));
  s.parameter(0).sigma = 0.02;  // restore a sensible set
  s.setCovariance([&] {
    TMatrixD r(3, 3);
    r(0, 0) = 0.02 * 0.02;
    r(1, 1) = 0.03 * 0.03;
    r(2, 2) = 0.01;
    r(0, 1) = r(1, 0) = 0.5 * 0.02 * 0.03;
    return r;
  }());

  // binning
  BinningAxis axis;
  axis.variable = "EventType";
  axis.values = {0, 1, 2};
  s.setBinning(Binning::product({axis}));
  ZAKU_CHECK(s.hasBinning());
  ZAKU_CHECK_EQ(s.binning().nBins(), 3);

  // --- file layout ---------------------------------------------------------
  const std::string file = scratch + "/paramset.root";
  s.write(file);
  {
    std::unique_ptr<TFile> f(TFile::Open(file.c_str(), "READ"));
    ZAKU_CHECK(f && !f->IsZombie());
    auto* names = f->Get<TObjArray>("det_param_names");
    ZAKU_CHECK(names != nullptr);
    if (names) ZAKU_CHECK_EQ(names->GetEntries(), 3);
    TKey* covKey = f->GetKey("det_param_cov");
    ZAKU_CHECK(covKey != nullptr);
    if (covKey) ZAKU_CHECK_STR_EQ(covKey->GetClassName(), "TMatrixT<double>");  // never TMatrixTSym
    ZAKU_CHECK(f->GetKey("det_param_lower") != nullptr);
    ZAKU_CHECK(f->GetKey("det_param_upper") != nullptr);
    ZAKU_CHECK(f->GetKey("det_param_priors") != nullptr);
  }
  const ParameterSet back = ParameterSet::read(file, "det");
  std::string why;
  ZAKU_CHECK_MSG(back.sameAs(s, 1e-12, &why), why);
  ZAKU_CHECK(back.parameter(2).hasBounds);
  ZAKU_CHECK_CLOSE(back.parameter(2).upper, 2.0, 1e-12);
  ZAKU_CHECK(back.parameter(0).hasBounds);  // bounds vectors cover every parameter once any has bounds
  ZAKU_CHECK_CLOSE(back.parameter(0).lower, 1.0, 1e-12);  // prior when unbounded
  ZAKU_CHECK_THROWS(ParameterSet::read(file, "flux"));

  // a set without bounds writes no bounds vectors
  ParameterSet nb("osc");
  nb.add(a);
  nb.write(scratch + "/paramset_nb.root");
  {
    std::unique_ptr<TFile> f(TFile::Open((scratch + "/paramset_nb.root").c_str(), "READ"));
    ZAKU_CHECK(f->GetKey("osc_param_lower") == nullptr);
  }
  const ParameterSet nbBack = ParameterSet::read(scratch + "/paramset_nb.root", "osc");
  ZAKU_CHECK(!nbBack.anyBounds());

  // assignment into a default-constructed set must carry the matrix (TMatrixD::operator= does not resize)
  ParameterSet assigned;
  assigned = s;
  ZAKU_CHECK_EQ(assigned.covariance().GetNrows(), 3);
  ZAKU_CHECK(assigned.sameAs(s, 1e-12, &why));
  assigned = ParameterSet("empty");
  ZAKU_CHECK_EQ(assigned.covariance().GetNrows(), 0);

  // sameAs detects differences
  ParameterSet other = s;
  other.parameter(1).prior = 1.1;
  TMatrixD cov = other.covariance();
  other.setCovariance(cov);
  ZAKU_CHECK(!other.sameAs(s, 1e-12, &why));
  ZAKU_CHECK(why.find("prior") != std::string::npos);

  ZAKU_CHECK(!s.describe().empty());
  ZAKU_CHECK_THROWS(ParameterSet("x").write(scratch + "/empty.root"));
  return zakutest::finish();
}
