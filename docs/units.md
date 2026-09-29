# Units

Zaku uses exactly one unit per physical quantity, everywhere: in the C++ API,
in every tree branch it writes, in every YAML key it reads, and on every axis
of every histogram it reads, writes or references. There is no per-call unit
argument except at the single point where an external input enters.

| Quantity | Unit | Notes |
|---|---|---|
| Energy, momentum, mass | GeV | all histogram axes, branches, YAML edges, response-table axes, smearing tables, migration matrices, flux histograms, oscillation grid |
| Q², W² | GeV² | branch names say which: `Q2`, `W2` |
| Reconstructed angle `Theta_rec` | degrees | as in the T2K binning files; cosines (`CosLep_*`) are dimensionless |
| Δm² | eV² | |
| Baseline | km | |
| Matter | g/cm³ × electron fraction (ρY_e) | HK default 2.6 × 0.5 |
| Detector mass | kt | the nucleon count is derived: 188 kt of water = 1.132e35 nucleons |
| Exposure | POT | fluxes are normalised per `pot_norm` POT (T2K release: 1e21) |
| Flux histogram content | ν / cm² / `bin_norm_GeV` / `pot_norm` POT | T2K: per 50 MeV per 1e21 POT; `bin_norm_GeV: 0` means per-bin integral |
| Cross-section scale (`fScaleFactor`) | cm² / nucleon / event | NUISANCE convention |
| Covariances | fractional | flux and detector sets; cross-section sets in parameter units |

## Rules

1. **GeV on every histogram axis.** A histogram whose x axis reaches beyond
   1000 is treated as an undeclared MeV input and rejected
   (`zaku::requireGeVAxis`). Inputs in MeV are legal but must say so with
   `units: MeV` next to the file reference; Zaku converts them on load with
   `zaku::histogramToGeV`, dividing the axis edges by 1000 and, for density
   histograms, multiplying the contents by 1000 so the density per GeV is
   preserved.
2. **Generator cards always say GeV.** NEUT cards carry the `"gev"` flag on the
   flux histogram; GENIE fluxes are written in GeV.
3. **T2KReWeight four-vectors are GeV** (`PISNu`, `PFSLep`), verified on the
   reference production; no conversion is applied.
4. **Particle masses** are the values the legacy macros used
   (`m_mu = 0.105658`, `m_e = 0.000511`, `m_p = 0.93827`, `m_n = 0.93957` GeV),
   kept unchanged so the reference production is reproducible.
5. **Angles** are converted to degrees only for `Theta_rec`; everything
   internal uses cosines or radians (`zaku::units::deg` converts degrees to
   radians).

## API

```cpp
#include "zaku/core/Units.h"
zaku::toGeV(1000.0, zaku::EnergyUnit::MeV);            // 1.0
zaku::units::m_mu;                                       // 0.105658 GeV
TH1* h = zaku::histogramToGeV(hMeV, zaku::EnergyUnit::MeV, /*contentIsDensity=*/true);
zaku::requireGeVAxis(*h, "flux file");                   // throws if the axis still looks like MeV
```

```python
import zaku
zaku.toGeV(1000.0, zaku.EnergyUnit.MeV)
zaku.units.m_mu
```
