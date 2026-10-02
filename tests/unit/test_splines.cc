// SplineBuilder: knot ordering with the inserted nominal, the legacy
// flat/rejected rules, the uniformity flag and the metadata round trip.
#include <TGraph.h>

#include <cmath>
#include <limits>

#include "ZakuTest.h"
#include "zaku/syst/Splines.h"

using namespace zaku;

int main() {
  DialKnots ma;
  ma.name = "MaCCQE";
  ma.nominal = 1.21;
  ma.knots = {0.818, 1.014, 1.406, 1.602};
  DialKnots rad;
  rad.name = "RadCorr_rate";
  rad.nominal = 0.0;
  rad.knots = {0.0, 0.25, 0.75};
  DialKnots inel;
  inel.name = "PionFSI_InelProb";
  inel.nominal = 1.002;
  inel.knots = {0.0111, 0.4515, 1.5525, 2.103};
  DialKnots src;
  src.name = "SRCFrac_O";
  src.nominal = 0.34;
  src.knots = {0.4225, 0.505};

  SplineBuilder all({ma, rad, inel, src});
  ZAKU_CHECK_EQ(static_cast<int>(all.size()), 4);
  ZAKU_CHECK_EQ(all.index("RadCorr_rate"), 1);
  ZAKU_CHECK_EQ(all.sourceIndex(1), 1);

  // nominal inserted in the middle, sorted
  const SplineDial& dm = all.dials()[0];
  ZAKU_CHECK_EQ(dm.nKnots(), 5);
  ZAKU_CHECK(dm.nominalInserted());
  ZAKU_CHECK_CLOSE(dm.x[2], 1.21, 1e-12);
  ZAKU_CHECK_EQ(dm.src[2], -1);
  ZAKU_CHECK_EQ(dm.src[0], 0);
  ZAKU_CHECK_EQ(dm.src[4], 3);
  ZAKU_CHECK(dm.uniformKnots());
  // nominal already a tweak: nothing inserted
  const SplineDial& dr = all.dials()[1];
  ZAKU_CHECK_EQ(dr.nKnots(), 3);
  ZAKU_CHECK(!dr.nominalInserted());
  ZAKU_CHECK(!dr.uniformKnots());  // 0, 0.25, 0.75
  // the -0.9 sigma knot makes the spacing non-uniform
  ZAKU_CHECK(!all.dials()[2].uniformKnots());
  // one-sided: nominal below every knot
  const SplineDial& ds = all.dials()[3];
  ZAKU_CHECK_EQ(ds.nKnots(), 3);
  ZAKU_CHECK_CLOSE(ds.x[0], 0.34, 1e-12);
  ZAKU_CHECK_EQ(ds.src[0], -1);
  ZAKU_CHECK(ds.uniformKnots());  // 0.34, 0.4225, 0.505

  // --- decisions -----------------------------------------------------------
  TGraph g;
  ZAKU_CHECK(all.build(0, {0.7, 0.87, 1.10, 1.17}, g) == SplineDecision::Full);
  ZAKU_CHECK_EQ(g.GetN(), 5);
  ZAKU_CHECK_CLOSE(g.GetX()[2], 1.21, 1e-12);
  ZAKU_CHECK_CLOSE(g.GetY()[2], 1.0, 1e-12);
  ZAKU_CHECK_CLOSE(g.GetY()[0], 0.7, 1e-12);
  ZAKU_CHECK_CLOSE(g.GetY()[4], 1.17, 1e-12);
  ZAKU_CHECK(all.build(0, {1.0, 1.0, 1.0, 1.0}, g) == SplineDecision::Flat);
  ZAKU_CHECK_EQ(g.GetN(), 0);
  ZAKU_CHECK(all.build(0, {0.0, 0.0, 0.0, 0.0}, g) == SplineDecision::Rejected);  // cv_weight == 0 events
  ZAKU_CHECK(all.build(0, {0.7, -0.1, 1.1, 1.2}, g) == SplineDecision::Rejected);  // negative
  ZAKU_CHECK(all.build(0, {0.7, std::numeric_limits<double>::quiet_NaN(), 1.1, 1.2}, g) == SplineDecision::Rejected);
  ZAKU_CHECK(all.build(0, {1.0, 1.0, 1.0, 1.2}, g) == SplineDecision::Full);  // one knot differs: real response
  ZAKU_CHECK(all.build(0, {1.3, 1.3, 1.3, 1.3}, g) == SplineDecision::Full);  // flat at != 1 is a real rescaling
  ZAKU_CHECK_EQ(g.GetN(), 5);
  ZAKU_CHECK_EQ(static_cast<int>(all.dials()[0].nFull), 3);
  ZAKU_CHECK_EQ(static_cast<int>(all.dials()[0].nFlat), 1);
  ZAKU_CHECK_EQ(static_cast<int>(all.dials()[0].nRejected), 3);
  ZAKU_CHECK_EQ(static_cast<int>(all.totalRejected()), 3);

  // nominal among the tweaks: the stored weight there must be 1
  ZAKU_CHECK(all.build(1, {1.0, 1.1, 1.3}, g) == SplineDecision::Full);
  ZAKU_CHECK_EQ(g.GetN(), 3);
  ZAKU_CHECK(all.build(1, {1.0 + 1e-7, 1.1, 1.3}, g) == SplineDecision::Full);  // within 1e-6
  ZAKU_CHECK(all.build(1, {1.01, 1.1, 1.3}, g) == SplineDecision::Rejected);
  ZAKU_CHECK(all.dials()[1].worstNominalDeviation > 0.009);
  // too few weights is an error
  ZAKU_CHECK_THROWS(all.build(0, {1.0, 1.0}, g));

  // --- subset --------------------------------------------------------------
  SplineBuilder some({ma, rad, inel, src}, {"SRCFrac_O", "MaCCQE"});
  ZAKU_CHECK_EQ(static_cast<int>(some.size()), 2);
  ZAKU_CHECK_STR_EQ(some.dialNames()[0], "SRCFrac_O");
  ZAKU_CHECK_EQ(some.sourceIndex(0), 3);
  ZAKU_CHECK_EQ(some.sourceIndex(1), 0);
  ZAKU_CHECK_THROWS(SplineBuilder({ma}, {"nope"}));
  DialKnots empty;
  empty.name = "empty";
  ZAKU_CHECK_THROWS(SplineBuilder({empty}));

  // --- metadata round trip ---------------------------------------------------
  const std::string meta = SplineBuilder::metadataText({ma, rad, inel, src});
  const std::vector<DialKnots> back = SplineBuilder::fromMetadataText(meta);
  ZAKU_CHECK_EQ(static_cast<int>(back.size()), 4);
  ZAKU_CHECK_STR_EQ(back[2].name, "PionFSI_InelProb");
  ZAKU_CHECK_CLOSE(back[2].nominal, 1.002, 1e-15);
  ZAKU_CHECK_EQ(static_cast<int>(back[2].knots.size()), 4);
  ZAKU_CHECK_CLOSE(back[2].knots[0], 0.0111, 1e-15);
  ZAKU_CHECK_CLOSE(back[1].nominal, 0.0, 1e-15);
  ZAKU_CHECK_THROWS(SplineBuilder::fromMetadataText("garbage"));

  ZAKU_CHECK(!all.statsText(10).empty());
  ZAKU_CHECK_STR_EQ(toString(SplineDecision::Flat), "flat");
  return zakutest::finish();
}
