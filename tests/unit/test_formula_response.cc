#include <TFile.h>
#include <TH1D.h>

#include "ZakuTest.h"
#include "ZakuTestPaths.h"
#include "zaku/detector/Formula.h"
#include "zaku/detector/Response.h"
#include "zaku/detector/Variables.h"

using namespace zaku;

namespace {

FlatEvent ccqeEvent(int pdgnu, double enu, bool qe) {
  FlatEvent ev;
  ev.cc = true;
  ev.PDGnu = pdgnu;
  ev.Enu_true = static_cast<float>(enu);
  ev.flagCCQE = qe;
  ev.PDGLep = pdgnu > 0 ? (std::abs(pdgnu) == 14 ? 13 : 11) : (std::abs(pdgnu) == 14 ? -13 : -11);
  ev.ELep = 0.5f;
  ev.x = 0.3f;
  ev.y = 0.4f;
  ev.W = 1.2f;
  return ev;
}

}  // namespace

int main() {
  const std::string scratch = zakutest::scratchDir();

  // --- Formula over the truth dictionary ------------------------------------
  VariableSet vars;
  ZAKU_CHECK(vars.has("Enu_true"));
  ZAKU_CHECK(vars.has("x"));  // sanitised to bjorken_x internally
  ZAKU_CHECK(vars.has("Enu_rec"));
  ZAKU_CHECK_EQ(vars.index("x"), vars.index("bjorken_x"));
  ZAKU_CHECK_STR_EQ(Formula::formulaName("x"), "bjorken_x");
  ZAKU_CHECK_STR_EQ(Formula::formulaName("Enu_true"), "Enu_true");

  Formula f("cc && PDGnu == 14 && Enu_true < 1.2 && bjorken_x > 0.2", vars.names());
  FlatEvent ev = ccqeEvent(14, 0.7, true);
  ProductionContext ctx;
  vars.loadTruth(ev, TruthVariables::compute(ev, ctx));
  ZAKU_CHECK(f.evalBool(vars.data()));
  vars.set("Enu_true", 2.0);
  ZAKU_CHECK(!f.evalBool(vars.data()));
  const auto refs = f.referencedVariables();
  ZAKU_CHECK_EQ(static_cast<int>(refs.size()), 4);
  ZAKU_CHECK_THROWS(Formula("nonsense(((", vars.names()));
  ZAKU_CHECK(Formula().empty());
  Formula g("0.9 * (1 - exp(-Enu_true / 0.4))", vars.names());
  vars.set("Enu_true", 0.4);
  ZAKU_CHECK_CLOSE(g.eval(vars.data()), 0.9 * (1 - std::exp(-1.0)), 1e-12);

  // --- table sources agree ---------------------------------------------------
  TableAxis axis;
  axis.variable = "Enu_true";
  axis.edges = {0.0, 0.5, 1.0, 2.0, 30.0};
  ZAKU_CHECK_EQ(axis.bin(-1.0), 0);   // clamped
  ZAKU_CHECK_EQ(axis.bin(0.7), 1);
  ZAKU_CHECK_EQ(axis.bin(50.0), 3);   // clamped
  const std::vector<double> values = {0.1, 0.2, 0.3, 0.4};
  TableSource arr = TableSource::array(values);
  auto* h = new TH1D("h", "h", 4, axis.edges.data());
  for (int i = 1; i <= 4; ++i) h->SetBinContent(i, values[static_cast<std::size_t>(i - 1)]);
  TableSource hist = TableSource::histogram(h, "");
  TableSource form = TableSource::formula("Enu_true < 0.5 ? 0.1 : (Enu_true < 1.0 ? 0.2 : (Enu_true < 2.0 ? 0.3 : 0.4))");
  TableSource cst = TableSource::constant(0.25);
  const std::vector<TableAxis> axes = {axis};
  arr.prepare(axes, vars, "arr");
  hist.prepare(axes, vars, "hist");
  form.prepare(axes, vars, "form");
  cst.prepare(axes, vars, "cst");
  for (double e : {0.1, 0.7, 1.5, 5.0, 40.0}) {
    vars.set("Enu_true", e);
    ZAKU_CHECK_CLOSE(arr.evaluate(vars), hist.evaluate(vars), 1e-12);
    ZAKU_CHECK_CLOSE(arr.evaluate(vars), form.evaluate(vars), 1e-12);
    ZAKU_CHECK_CLOSE(cst.evaluate(vars), 0.25, 1e-12);
  }
  TableSource wrong = TableSource::array({1.0, 2.0});
  ZAKU_CHECK_THROWS(wrong.prepare(axes, vars, "wrong"));

  // --- a small response -------------------------------------------------------
  Response r;
  r.addCategory("numuCC", "cc && abs(PDGnu) == 14");
  r.addCategory("NC", "!cc");
  ZAKU_CHECK_THROWS(r.addCategory("NC", "1"));
  ResponseTable t;
  t.category = "numuCC";
  t.axes = axes;
  t.efficiency = TableSource::array({1.0, 0.9, 0.8, 0.7});
  t.columns = {"A", "B"};
  t.probabilities = {TableSource::array({0.5, 0.6, 0.7, 0.8}), TableSource::constant(0.1)};
  r.addTable(t);
  ResponseTable nc;
  nc.category = "NC";
  nc.efficiency = TableSource::constant(1.0);
  nc.columns = {"A"};
  nc.probabilities = {TableSource::constant(0.05)};
  r.addTable(nc);
  ResponseTable bad;
  bad.category = "numuCC";
  bad.columns = {"A"};
  ZAKU_CHECK_THROWS(r.addTable(bad));  // no probability sources
  ZAKU_CHECK_THROWS(r.evaluate(vars, BeamMode::FHC));  // not prepared
  r.prepare(vars);
  ZAKU_CHECK_EQ(static_cast<int>(r.columns().size()), 2);

  ev = ccqeEvent(14, 0.7, true);
  vars.loadTruth(ev, TruthVariables::compute(ev, ctx));
  ResponseResult res = r.evaluate(vars, BeamMode::FHC);
  ZAKU_CHECK_EQ(res.category, 0);
  ZAKU_CHECK_EQ(res.table, 0);
  ZAKU_CHECK_CLOSE(res.efficiency, 0.9, 1e-12);
  ZAKU_CHECK_EQ(static_cast<int>(res.probabilities.size()), 2);
  ZAKU_CHECK_CLOSE(res.probabilities[0], 0.6, 1e-12);
  ZAKU_CHECK_CLOSE(res.probabilities[1], 0.1, 1e-12);
  ZAKU_CHECK_CLOSE(res.probabilitySum(), 0.7, 1e-12);
  FlatEvent ncEv = ccqeEvent(14, 0.7, false);
  ncEv.cc = false;
  vars.loadTruth(ncEv, TruthVariables::compute(ncEv, ctx));
  res = r.evaluate(vars, BeamMode::RHC);
  ZAKU_CHECK_EQ(res.category, 1);
  ZAKU_CHECK_EQ(res.table, 1);
  FlatEvent nueEv = ccqeEvent(12, 0.7, true);  // no category
  vars.loadTruth(nueEv, TruthVariables::compute(nueEv, ctx));
  res = r.evaluate(vars, BeamMode::FHC);
  ZAKU_CHECK_EQ(res.category, -1);
  ZAKU_CHECK_EQ(res.table, -1);

  // probabilities above 1 are an error unless normalised
  Response over;
  over.addCategory("all", "1");
  ResponseTable ot;
  ot.category = "all";
  ot.efficiency = TableSource::constant(1.0);
  ot.columns = {"A", "B"};
  ot.probabilities = {TableSource::constant(0.7), TableSource::constant(0.6)};
  over.addTable(ot);
  over.prepare(vars);
  ZAKU_CHECK_THROWS(over.evaluate(vars, BeamMode::FHC));
  over.setNormaliseProbabilities(true);
  res = over.evaluate(vars, BeamMode::FHC);
  ZAKU_CHECK_CLOSE(res.probabilitySum(), 1.0, 1e-12);

  // --- YAML round trip --------------------------------------------------------
  const std::string yaml = r.toYaml();
  Response back = Response::fromYamlString(yaml);
  back.prepare(vars);
  vars.loadTruth(ev, TruthVariables::compute(ev, ctx));
  ResponseResult res2 = back.evaluate(vars, BeamMode::FHC);
  ZAKU_CHECK_CLOSE(res2.efficiency, 0.9, 1e-12);
  ZAKU_CHECK_CLOSE(res2.probabilities[0], 0.6, 1e-12);
  ZAKU_CHECK_STR_EQ(back.toYaml(), yaml);
  ZAKU_CHECK_THROWS(Response::fromYamlString("tables: [{category: a, bogus: 1}]"));

  // --- the legacy HK file -----------------------------------------------------
  const std::string legacyFile = zakutest::sourceDir() + "/data/hk/response.root";
  if (zakutest::exists(legacyFile)) {
    Response hk = Response::fromLegacyFile(legacyFile);
    ZAKU_CHECK_EQ(static_cast<int>(hk.categories().size()), 9);
    ZAKU_CHECK_EQ(static_cast<int>(hk.tables().size()), 9);
    ZAKU_CHECK_STR_EQ(hk.categories()[8].name, "NC");
    ZAKU_CHECK_STR_EQ(hk.columns()[2], "1Re1de");
    hk.prepare(vars);
    // numu CCQE at 0.7 GeV: h_eff bin 3 (0.6-1.2) = 1.12235, 1Rmu = 0.93853
    vars.loadTruth(ev, TruthVariables::compute(ev, ctx));
    ResponseResult hr = hk.evaluate(vars, BeamMode::FHC);
    ZAKU_CHECK_EQ(hr.category, 0);
    ZAKU_CHECK_CLOSE(hr.efficiency, 1.12235, 1e-9);
    ZAKU_CHECK_CLOSE(hr.probabilities[0], 0.93853, 1e-9);
    // numubar CCnonQE -> category 3
    FlatEvent nb = ccqeEvent(-14, 0.2, false);
    vars.loadTruth(nb, TruthVariables::compute(nb, ctx));
    hr = hk.evaluate(vars, BeamMode::FHC);
    ZAKU_CHECK_EQ(hr.category, 3);
    ZAKU_CHECK_CLOSE(hr.efficiency, 0.16641, 1e-9);
    // nuebar CCQE -> 5, nue CCnonQE -> 6, NC -> 8
    FlatEvent eb = ccqeEvent(-12, 5.0, true);
    vars.loadTruth(eb, TruthVariables::compute(eb, ctx));
    ZAKU_CHECK_EQ(hk.evaluate(vars, BeamMode::FHC).category, 5);
    FlatEvent en = ccqeEvent(12, 1.0, false);
    vars.loadTruth(en, TruthVariables::compute(en, ctx));
    ZAKU_CHECK_EQ(hk.evaluate(vars, BeamMode::FHC).category, 6);
    vars.loadTruth(ncEv, TruthVariables::compute(ncEv, ctx));
    hr = hk.evaluate(vars, BeamMode::FHC);
    ZAKU_CHECK_EQ(hr.category, 8);
    ZAKU_CHECK_CLOSE(hr.efficiency, 0.0063, 1e-9);
    // the YAML export reproduces the histogram lookups
    Response hy = Response::fromYamlString(hk.toYaml());
    hy.prepare(vars);
    for (double e : {0.1, 0.5, 0.9, 2.0, 10.0}) {
      vars.set("Enu_true", e);
      const ResponseResult a = hk.evaluate(vars, BeamMode::RHC);
      const ResponseResult b = hy.evaluate(vars, BeamMode::RHC);
      ZAKU_CHECK_CLOSE(a.efficiency, b.efficiency, 1e-12);
      for (std::size_t c = 0; c < a.probabilities.size(); ++c) ZAKU_CHECK_CLOSE(a.probabilities[c], b.probabilities[c], 1e-12);
    }
    const std::string shipped = zakutest::sourceDir() + "/configs/response/hk_2025.yaml";
    if (zakutest::exists(shipped)) {
      Response hs = Response::fromFile(shipped);
      hs.prepare(vars);
      vars.loadTruth(ev, TruthVariables::compute(ev, ctx));
      ZAKU_CHECK_CLOSE(hs.evaluate(vars, BeamMode::FHC).efficiency, 1.12235, 1e-9);
    }
  }

  return zakutest::finish();
}
