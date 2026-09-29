# Zaku

**Zaku Assembles Kinematics and Systematics** — a C++17 / ROOT library that
builds fake far-detector inputs for [GUNDAM](https://github.com/gundam-organization/gundam)
sensitivity studies: generator configuration (NEUT or GENIE), a simulation-free
detector response, systematic splines and covariances, and the GUNDAM
configuration itself.

One object, `zaku::Experiment`, holds everything the user asserts about a
hypothetical experiment (flux and POT, baseline, detector mass, response
tables, smearing, samples, generator choices) and exposes the pipeline stages
as methods. The same classes work from C++, the ROOT prompt, ROOT macros and
PyROOT/Jupyter through a ROOT dictionary; the `zaku` command is a thin wrapper.

Zaku never runs NEUT, GENIE, T2KReWeight, NUISANCE or GUNDAM itself. It writes
their inputs and a manifest of exact command lines; you run them (by hand or
on a batch system). Zaku's own stages are plain in-process ROOT file processing.

## Status

Stages 0 and 1 of the implementation plan are in place: build system and ROOT
dictionary, unit conventions, the core value types (`NuSpecies`, `BeamMode`,
`Target`, `Detector`), the `Experiment` object with its YAML builder
([docs/config_reference.md](docs/config_reference.md)), flux loading and
generator-ready rebinning, the flux covariance reader, and the event layer
(flat-tree cursor, T2KReWeight weight cursor, and the join between them,
[docs/pipeline.md](docs/pipeline.md)). The processing stages are added one at
a time; `zaku --help` lists which are available in a given build.

```bash
zaku configs/hk_t2k_like.yaml describe          # what the experiment asserts
zaku configs/hk_t2k_like.yaml list-productions  # generator productions, in index order
zaku configs/hk_t2k_like.yaml status            # which production files exist
```

## Build

Requires ROOT ≥ 6.30 (with PyROOT), yaml-cpp, CMake ≥ 3.20 and a C++17 compiler.

```bash
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$PWD/install
cmake --build build -j
ctest --test-dir build --output-on-failure
cmake --install build
source install/setup.sh          # or: source build/setup.sh to use the build tree
```

`setup.sh` puts `libZaku.so` (with its dictionary and rootmap) on the dynamic
path, makes the headers visible to Cling/PyROOT, and adds the python shim and
the `zaku` command to your paths.

## Quick start

C++ / ROOT prompt. Classes autoload from the rootmap on first use; free
functions such as `zaku::version()` need the umbrella header, which
`setup.sh` makes findable through `ROOT_INCLUDE_PATH`:

```cpp
root [0] #include "zaku/Zaku.h"
root [1] auto det = zaku::Detector::water(188.0);
root [2] det.nNucleons()
(double) 1.1321625e+35
```

Python:

```python
import zaku
det = zaku.Detector.water(188.0)
print(zaku.version(), det.nNucleons(), det.genieTargetString())
```

## Layout

```
include/zaku/   public headers (umbrella zaku/Zaku.h)
src/            implementation and the dictionary LinkDef
apps/           the zaku command-line tool
python/zaku/    PyROOT shim (import zaku)
tests/          CTest suite: unit, python, cling, regression
docs/           user documentation (start with docs/units.md)
```

## Units

GeV everywhere, on every histogram axis. See [docs/units.md](docs/units.md).

## Licence

NuFast (vendored later under `data/third_party`) is MIT-licensed by
P. Denton. Zaku's own code: see `LICENSE` once added.
