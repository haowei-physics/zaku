"""Stage 3 python check: dial catalogue, parameter sets and the systematics block from PyROOT."""
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
    scratch = os.path.join(os.environ.get("ZAKU_BINARY_DIR", "."), "tests", "scratch")
    os.makedirs(scratch, exist_ok=True)
    os.environ["ZAKU_DATA"] = src

    cat = zaku.DialCatalog.fromFile(os.path.join(src, "configs/dials/neut_default.yaml"))
    check(cat.size() == 19, "19 dials")
    check(cat.fitted().size() == 15, "15 fitted dials")
    check(str(cat.dial("MaCCQE").name) == "MaCCQE", "dial lookup")

    osc = zaku.OscCovariance.build(zaku.OscillationParameters.asimovA())
    check(osc.size() == 7, "7 oscillation parameters")
    check(abs(osc.sigma(0) - 0.013) < 1e-12, "s12 sigma")
    path = os.path.join(scratch, "py_osc.root")
    osc.write(path)
    back = zaku.ParameterSet.read(path, "osc")
    check(back.sameAs(osc, 1e-12), "parameter set round trip")

    exp = zaku.Experiment.fromConfig(os.path.join(src, "configs/hk_t2k_like.yaml"))
    check(exp.hasDialCatalog(), "experiment has a dial catalogue")
    check(exp.systematics().xsec.fitted().size() == 15, "fitted dials from the config")
    check(exp.shifts().size() == 2, "two shift systematics")
    check(exp.shifts()[0].target == zaku.ShiftTarget.LeptonMomentum, "shift target")
    det = zaku.DetectorCovariance.normalisation(exp.samples(), exp.systematics().detector)
    check(det.size() == 5, "5 detector normalisations")
    check(det.binning().nBins() == 5, "detector binning")

    ro = zaku.RecoObservables()
    ro.PLep_rec = 0.6
    ro.CosLep_rec = 0.9
    ro.Ehad_rec = 0.3
    zaku.ShiftInjector.rederive(ro, False, 0.027, True, True)
    zaku.ShiftInjector.apply(ro, exp.shifts()[0], 1.0, False, 0.027, True, zaku.Hypothesis.Muon, 13)
    check(abs(ro.PLep_rec - 0.6 * 1.024) < 1e-6, "momentum scale applied")

    knots = zaku.DialKnots()
    knots.name = "MaCCQE"
    knots.nominal = 1.21
    for k in (0.818, 1.014, 1.406, 1.602):
        knots.knots.push_back(k)
    builder = zaku.SplineBuilder([knots])
    check(builder.dials()[0].nKnots() == 5, "nominal inserted")

    print("test_systematics: %d failure(s)" % failures)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
