"""Stage 2 python check: response tables, smearing, estimators and samples from PyROOT."""
import os
import sys


def main() -> int:
    import zaku

    failures = 0

    def check(cond, what):
        nonlocal failures
        if not cond:
            failures += 1
            print("FAIL:", what)

    src = os.environ.get("ZAKU_SOURCE_DIR", ".")
    os.environ["ZAKU_DATA"] = src

    resp = zaku.Response.fromFile(os.path.join(src, "configs/response/hk_2025.yaml"))
    check(resp.categories().size() == 9, "9 categories")
    check(zaku.strings(resp.columns()) == ["1Rmu", "1Re0de", "1Re1de"], "columns")
    vars_ = zaku.VariableSet()
    resp.prepare(vars_)
    ev = zaku.FlatEvent()
    ev.cc = True
    ev.PDGnu = 14
    ev.Enu_true = 0.7
    ev.flagCCQE = True
    ev.PDGLep = 13
    ev.ELep = 0.5
    ctx = zaku.ProductionContext()
    vars_.loadTruth(ev, zaku.TruthVariables.compute(ev, ctx, 0.0))
    r = resp.evaluate(vars_, zaku.BeamMode.FHC)
    check(r.category == 0 and abs(r.efficiency - 1.12235) < 1e-9, "numu CCQE efficiency %g" % r.efficiency)
    check(abs(r.probabilities[0] - 0.93853) < 1e-9, "1Rmu probability")

    sm = zaku.Smearing.parametric(0.024, 1.8)
    sm.prepare(vars_)
    check(abs(sm.leptonMomentum().sigma(0.62, vars_) - 0.024 * 0.65) < 1e-12, "sigma_p")

    enu = zaku.reco.leptonMomentum(0.8, 13)
    check(0.79 < enu < 0.8, "lepton momentum")
    check(abs(zaku.reco.thetaDegrees(0.0) - 90.0) < 1e-3, "theta")

    b = zaku.Binning.fromFile(os.path.join(src, "data/binning/binning_T2K_FHC_ELike.txt"))
    check(b.nBins() == 225, "ELike bins")

    exp = zaku.Experiment.fromConfig(os.path.join(src, "configs/hk_t2k_like.yaml"))
    check(exp.samples().size() == 5, "5 samples")
    check(str(exp.samples()[0].name) == "FHC_1Rmu", "first sample")
    check(exp.assignment() == zaku.Assignment.Split, "split assignment")
    check(exp.energyEstimator() == zaku.EnergyEstimator.CCQE, "ccqe estimator")
    osc = exp.oscillator()
    p = osc.probability(14, 12, 0.6)
    check(0.03 < p < 0.08, "appearance probability %g" % p)
    check("FHC_1Rmu" in str(exp.describe()), "describe mentions samples")

    print("test_detector: %d failure(s)" % failures)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
