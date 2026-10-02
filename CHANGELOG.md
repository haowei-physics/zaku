# Changelog

Versioning: the minor version advances with each completed implementation
stage, the patch version with fixes or renames in between. The library version
is `zaku::version()`; the file-format/config schema is `zaku::schemaVersion()`.

## 0.3.0 — 2026-09-30

Stage 3: systematics.

- `DialCatalog` (`configs/dials/neut_default.yaml`): the generator's dials with
  their tweak points and `fit` flag; `systematics.xsec.dials` in the config.
- `ParameterSet`: the GUNDAM parameter-set file layout (`TMatrixD`, `TObjArray`
  names, priors, bounds, binning) with read/write/compare.
- `SplineBuilder` and the `splines` stage: per-event `TGraph` response
  branches with the legacy rules (nominal inserted, flat and unusable responses
  as empty graphs); spline files carry `zaku/dials`.
- The four parameter sets: `OscCovariance`, `DetectorCovariance`
  (normalisations with pair correlations, or from shift universes),
  `XsecCovariance` (sigma from the knots, bounds from the outer knots),
  `FluxCovariance` (far block, `C − s·BᵀA⁻¹B`, binning); the `covariances`
  stage writes them plus the oscillation grid under `parameters/`.
- `ShiftParameter`/`ShiftInjector` (`shifts:`): scales and offsets of the
  reconstructed quantities with everything downstream re-derived.
- `SystematicsSummary`/`SummaryTable` and the `summary` stage: ±1σ effect of
  every parameter on every sample, csv + ROOT histograms.
- `Experiment`: `setSystematics/systematics`, `setDialCatalog`, `addShift`,
  `makeSplines`, `makeCovariances`, `summariseSystematics`; `status` lists the
  parameter files; YAML keys `systematics` and `shifts` are now interpreted.
- CLI: `splines`, `covariances`, `summary`, `list-dials`, `show-systematics`.
- Regression: the legacy configuration reproduces the reference spline file
  graph by graph and the four reference parameter files and binnings.

## 0.2.0 — 2026-09-29

Stage 2: the fake detector.

- Response tables (`Response`, `TableSource`, `configs/response/*.yaml`): true
  event → detection probability and per-column probabilities from constants,
  inline arrays, ROOT histograms or formulas; the legacy HK mixing matrix is
  `Response::fromLegacyFile` / `configs/response/hk_2025.yaml`.
- `Smearing` with fractional, absolute, legacy-binned, stochastic (a/√E ⊕ b),
  histogram and formula resolutions; decay-electron tagging.
- Reconstruction (`reco::`): legacy-exact lepton smearing, CCQE and
  calorimetric energy estimators, `Theta_rec`; `RecoObservables`.
- Analysis samples (`SampleDefinition`, `Binning` with GUNDAM txt read/write,
  `configs/samples/*.yaml`), `split`/`draw` assignment, per-sample hypotheses and cuts.
- `Oscillator` (vendored NuFast, `third_party/NuFast`) and the NuOscillator grid.
- Stages `recon` and `select` (`Experiment::reconstruct/select`, CLI), the
  FakeSelect tree with `ProductionIndex`, `isAppearance`, `SampleProb` and both
  estimators, the rate table (`RateTable`, csv + TH2), `show-response`, `show-samples`.
- Regression: the legacy configuration reproduces the reference `_fakerecon.root`
  bit for bit and the reference `fakeselect_*.root` entry by entry.

## 0.1.1 — 2026-09-29

- Rename the generator-run concept from "MC sample" to **production**:
  `McSample` → `Production`, `SampleLayout` → `Layout`, `SampleContext` →
  `ProductionContext`, `Experiment::addMcSample/mcSamples/...` →
  `addProduction/productions/...`, YAML `mc_samples:` → `productions:`, CLI
  `list-samples` → `list-productions`. "Sample" now means a reconstructed and
  selected analysis sample only.

## 0.1.0 — 2026-09-29

Stages 0 and 1 of the implementation plan.

- Build system with the ROOT dictionary (`libZaku_rdict.pcm` + `libZaku.rootmap`),
  install layout, `ZakuConfig.cmake`, generated `setup.sh`, PyROOT shim
  (`import zaku`), the `zaku` command.
- Core types: `NuSpecies`, `BeamMode`, `Target`, `Detector`, `Units` (GeV
  everywhere; MeV inputs declared and converted), `version()`, `schemaVersion()`.
- `Experiment`: the single user-facing object with its YAML builder,
  `describe`, `status`, `validate`, `potWeight`.
- Flux layer: `FluxSpec`, `FluxHistograms`, `FluxPreparer::makeUniform`,
  `FluxCovarianceInput`.
- Event layer: `FlatEvent`, `TruthVariables`, `FlatTreeCursor`,
  `T2KRWWeightCursor`, `EventJoin`.
- `Manifest`, `Provenance`, `RunSettings`; shipped data and configs; docs on
  units, configuration and the pipeline file contracts.
- Tests: unit, Python, Cling, CLI, and regression against the reference
  production (uniform flux reproduction, flux integrals, the flat/weight join,
  the legacy POT weight).
