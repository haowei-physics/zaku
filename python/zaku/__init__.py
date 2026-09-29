"""Zaku - Zaku Assembles Kinematics and Systematics.

Thin Python layer over the ROOT dictionary of libZaku. Importing this module
loads the library and exposes the C++ namespace ``zaku`` as module attributes,
so that::

    import zaku
    det = zaku.Detector.water(188.0)
    print(zaku.version(), det.nNucleons())

is the same object model as in C++ or on the ROOT prompt. Nothing here
re-implements physics; the helpers only make PyROOT objects nicer to handle.

libZaku.so must be findable: source the ``setup.sh`` generated in the build or
install tree (it sets LD_LIBRARY_PATH, ROOT_INCLUDE_PATH and PYTHONPATH).
"""
from __future__ import annotations

import os as _os

import ROOT as _ROOT

_LIBRARY = "libZaku"

_status = _ROOT.gSystem.Load(_LIBRARY)
if _status < 0:
    raise ImportError(
        "zaku: could not load %s (gSystem.Load returned %d). Source the Zaku setup.sh so that "
        "libZaku.so is on LD_LIBRARY_PATH." % (_LIBRARY, _status)
    )

# Classes are autoparsed from the dictionary on demand, but free functions in
# the namespace are not: parse the umbrella header once (ROOT_INCLUDE_PATH,
# set by setup.sh, makes it findable).
if not _ROOT.gInterpreter.Declare('#include "zaku/Zaku.h"'):
    raise ImportError(
        "zaku: could not parse zaku/Zaku.h. Source the Zaku setup.sh so that the headers are on "
        "ROOT_INCLUDE_PATH."
    )

_ns = _ROOT.zaku

# Explicit re-exports of the stable public API (the module-level __getattr__
# below covers everything else in the namespace).
version = _ns.version
versionMajor = _ns.versionMajor
versionMinor = _ns.versionMinor
versionPatch = _ns.versionPatch
schemaVersion = _ns.schemaVersion
buildInfo = _ns.buildInfo

Flavour = _ns.Flavour
NuSpecies = _ns.NuSpecies
BeamMode = _ns.BeamMode
Target = _ns.Target
Detector = _ns.Detector
NeutTarget = _ns.NeutTarget
EnergyUnit = _ns.EnergyUnit
units = _ns.units

toGeV = _ns.toGeV
fromGeV = _ns.fromGeV
histogramToGeV = _ns.histogramToGeV
requireGeVAxis = _ns.requireGeVAxis
looksLikeMeVAxis = _ns.looksLikeMeVAxis

# Stage 1: experiment, flux, events
Experiment = _ns.Experiment
FluxSpec = _ns.FluxSpec
FluxCovarianceSpec = _ns.FluxCovarianceSpec
Production = _ns.Production
Layout = _ns.Layout
Manifest = _ns.Manifest
ManifestJob = _ns.ManifestJob
Provenance = _ns.Provenance
RunSettings = _ns.RunSettings
FluxHistograms = _ns.FluxHistograms
FluxPreparer = _ns.FluxPreparer
FluxCovarianceInput = _ns.FluxCovarianceInput
FlatEvent = _ns.FlatEvent
ProductionContext = _ns.ProductionContext
TruthVariables = _ns.TruthVariables
FlatTreeCursor = _ns.FlatTreeCursor
DialKnots = _ns.DialKnots
T2KRWWeightCursor = _ns.T2KRWWeightCursor
EventJoin = _ns.EventJoin
JoinOptions = _ns.JoinOptions
JoinStats = _ns.JoinStats

DATA_DIR = _os.environ.get("ZAKU_DATA", "")


def __getattr__(name: str):
    """Fall back to the C++ namespace for anything not re-exported above."""
    try:
        return getattr(_ns, name)
    except AttributeError as err:
        raise AttributeError("module 'zaku' has no attribute '%s'" % name) from err


def strings(sequence) -> list:
    """Convert a C++ ``std::vector<std::string>`` (or any iterable of C++
    strings / bytes) into a list of Python ``str``."""
    out = []
    for item in sequence:
        if isinstance(item, bytes):
            out.append(item.decode())
        else:
            out.append(str(item))
    return out


def to_numpy(hist):
    """Return ``(edges, contents, errors)`` numpy arrays for a 1-D ROOT histogram."""
    import numpy as np

    n = hist.GetNbinsX()
    edges = np.array([hist.GetXaxis().GetBinLowEdge(i) for i in range(1, n + 2)])
    contents = np.array([hist.GetBinContent(i) for i in range(1, n + 1)])
    errors = np.array([hist.GetBinError(i) for i in range(1, n + 1)])
    return edges, contents, errors


__all__ = [
    "version", "versionMajor", "versionMinor", "versionPatch", "schemaVersion", "buildInfo",
    "Flavour", "NuSpecies", "BeamMode", "Target", "Detector", "NeutTarget", "EnergyUnit", "units",
    "toGeV", "fromGeV", "histogramToGeV", "requireGeVAxis", "looksLikeMeVAxis",
    "Experiment", "FluxSpec", "FluxCovarianceSpec", "Production", "Layout", "Manifest", "ManifestJob",
    "Provenance", "RunSettings", "FluxHistograms", "FluxPreparer", "FluxCovarianceInput",
    "FlatEvent", "ProductionContext", "TruthVariables", "FlatTreeCursor", "DialKnots", "T2KRWWeightCursor",
    "EventJoin", "JoinOptions", "JoinStats",
    "DATA_DIR", "strings", "to_numpy",
]
