"""Stage 1 python check: Experiment, FluxSpec and the YAML round trip from PyROOT."""
import os
import sys
import tempfile


def main() -> int:
    import ROOT
    import zaku

    failures = 0

    def check(cond, what):
        nonlocal failures
        if not cond:
            failures += 1
            print("FAIL:", what)

    tmp = tempfile.mkdtemp(prefix="zaku_py_")
    flux_file = os.path.join(tmp, "flux.root")
    f = ROOT.TFile(flux_file, "RECREATE")
    for name in ("enu_sk_numu", "enu_sk_numub", "enu_sk_nue", "enu_sk_nueb"):
        h = ROOT.TH1D(name, name, 300, 0.0, 3.0)
        for i in range(1, 301):
            h.SetBinContent(i, 1.0e5)
        h.Write()
    f.Close()

    exp = zaku.Experiment("PyHK")
    exp.setBaseline(295.0).setMatter(2.6, 0.5).setDetector(zaku.Detector.water(188.0))
    exp.setFlux(zaku.BeamMode.FHC, zaku.FluxSpec.t2k(flux_file, "sk"), 6.75e21)
    exp.addMcSample(zaku.BeamMode.FHC, 14, 14, 2000)
    exp.addMcSample(zaku.BeamMode.FHC, 14, 12, 1000)
    exp.setRun(os.path.join(tmp, "work"), 42)

    check(str(exp.name()) == "PyHK", "name")
    check(abs(exp.rhoYe() - 1.3) < 1e-12, "rhoYe")
    names = zaku.strings(exp.mcSampleNames())
    check(names == ["numu_x_numu_2K_FHC", "numu_x_nue_1K_FHC"], "sample names %s" % names)
    per_pot = 1.0e5 * 3.0 / 0.05 / 1.0e21
    check(abs(exp.fluxIntegralPerPOT(zaku.BeamMode.FHC, 14) - per_pot) / per_pot < 1e-12, "flux integral")
    pw = exp.potWeight(zaku.BeamMode.FHC, 14, 3.0e-43)
    check(abs(pw - 3.0e-43 * per_pot * 6.75e21 * exp.nNucleons()) / pw < 1e-12, "potWeight")

    yaml_text = str(exp.toYaml())
    check("beam_modes" in yaml_text and "mc_samples" in yaml_text, "toYaml content")
    cfg = os.path.join(tmp, "exp.yaml")
    exp.saveConfig(cfg)
    back = zaku.Experiment.fromConfig(cfg)
    check(str(back.name()) == "PyHK", "fromConfig name")
    check(abs(back.pot(zaku.BeamMode.FHC) - 6.75e21) < 1.0, "fromConfig POT")
    check(zaku.strings(back.mcSampleNames()) == names, "fromConfig samples")
    check(str(back.toYaml()) == yaml_text, "round trip idempotent")
    check("numu_x_numu_2K_FHC" in str(back.status()), "status")

    spec = zaku.FluxSpec.t2k("x.root", "nd280")
    check(str(spec.histogram(-12)) == "enu_nd280_nueb", "FluxSpec.t2k")
    s = zaku.McSample(zaku.BeamMode.RHC, -14, -12, 50000)
    check(str(s.name()) == "numubar_x_nuebar_50K_RHC" and s.isAppearance(), "McSample")
    check(str(zaku.SampleLayout.flatFile(s)) == "numubar/numubar_x_nuebar_50K_RHC_flattree.root", "SampleLayout")

    m = zaku.Manifest("generate")
    job = m.add("j", "echo hi")
    check("echo hi" in str(m.text()), "Manifest text")

    print("test_experiment: %d failure(s)" % failures)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
