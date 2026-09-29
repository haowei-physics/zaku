"""Stage 0 python check: the dictionary loads and the core API is usable from PyROOT."""
import sys


def main() -> int:
    import zaku

    failures = 0

    def check(cond, what):
        nonlocal failures
        if not cond:
            failures += 1
            print("FAIL:", what)

    ver = zaku.version()
    check(isinstance(str(ver), str) and len(str(ver)) >= 5, "version string")
    check(zaku.versionMajor() >= 0, "versionMajor")
    check("ROOT" in str(zaku.buildInfo()), "buildInfo mentions ROOT")

    det = zaku.Detector.water(188.0)
    n = det.nNucleons()
    check(abs(n - 1.132e35) / 1.132e35 < 1e-3, "water(188).nNucleons() = %g" % n)
    check(str(det.genieTargetString()) == "1000080160[0.888889],1000010010[0.111111]", "genie target string")
    nt = det.neutTarget()
    check((nt.numNucleons, nt.numBoundProtons, nt.numFreeProtons) == (16, 8, 2), "neut target")

    sp = zaku.NuSpecies.fromName("numubar")
    check(sp.pdg() == -14, "NuSpecies.fromName('numubar').pdg()")
    check(str(sp.name()) == "numubar", "NuSpecies.name()")
    check(zaku.beamModeFromString("rhc") == zaku.BeamMode.RHC, "beamModeFromString")
    check(str(zaku.toString(zaku.BeamMode.FHC)) == "FHC", "toString(BeamMode)")

    check(abs(zaku.toGeV(1000.0, zaku.EnergyUnit.MeV) - 1.0) < 1e-12, "toGeV")
    check(abs(zaku.units.m_mu - 0.105658) < 1e-12, "units.m_mu")

    import ROOT
    h = ROOT.TH1D("h", "h", 10, 0.0, 30000.0)
    check(zaku.looksLikeMeVAxis(h), "looksLikeMeVAxis")
    hg = zaku.histogramToGeV(h, zaku.EnergyUnit.MeV, False, 1, "hg")
    check(abs(hg.GetXaxis().GetXmax() - 30.0) < 1e-12, "histogramToGeV axis")

    v = ROOT.std.vector("std::string")()
    v.push_back("a")
    v.push_back("b")
    check(zaku.strings(v) == ["a", "b"], "strings() conversion")

    # Fallback attribute access to the C++ namespace.
    check(zaku.elementSymbol(8) == "O" or str(zaku.elementSymbol(8)) == "O", "namespace fallback")

    print("python_import: %d failure(s)" % failures)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
