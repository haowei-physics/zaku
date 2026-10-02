# Configuration reference

A Zaku experiment is one YAML file, loaded with `zaku::Experiment::fromConfig()`
(C++), `zaku.Experiment.fromConfig()` (Python) or `zaku <file> <stage>` on the
command line. Every key maps one-to-one onto a setter of `zaku::Experiment`, so
anything the file can express can also be built in code, and
`Experiment::toYaml()` writes the file back.

Conventions:

- Paths are resolved relative to the directory of the YAML file unless they are
  absolute. `${VAR}` is replaced by the environment variable `VAR`; `setup.sh`
  sets `ZAKU_DATA` to the directory holding the shipped `data/` and `configs/`.
- Energies are GeV everywhere ([units](units.md)). A histogram file in MeV must
  say `units: MeV` next to the file name.
- Unknown keys are errors, so typos do not pass silently.
- Sections this build does not interpret yet (`generator`, `gundam`) are kept
  verbatim and written back by `toYaml()`.
- A block can be `{ include: other.yaml }` (paths relative to the including
  file) where noted; the included file holds the block's content.

## Keys

| Key | Setter | Meaning |
|---|---|---|
| `name` | `Experiment(name)` / `setName` | Experiment name, used in reports and provenance. |
| `baseline_km` | `setBaseline` | Source-to-detector distance in km. |
| `matter.density_g_cm3`, `matter.electron_fraction` | `setMatter` | Matter density and electron fraction for oscillations (`rhoYe()` = product). Defaults 2.6 and 0.5. |
| `detector` | `setDetector` | See below. |
| `flux_units.pot_norm`, `flux_units.bin_norm_GeV` | `setFluxUnits` | Default normalisation of every flux histogram: neutrinos / cm² / `bin_norm_GeV` GeV / `pot_norm` POT. Defaults 1e21 and 0.05 (T2K release). `bin_norm_GeV: 0` means the contents are per-bin integrals. |
| `beam_modes.<FHC\|RHC>.pot` | `setPOT` | Exposure of that beam mode in protons on target. |
| `beam_modes.<FHC\|RHC>.flux` | `setFlux` | See below. |
| `flux_covariance` | `setFluxCovariance` | See below. |
| `productions` | `addProduction` | List of generator productions, see below. The list order defines `SampleIndex`. |
| `response` | `setResponse`, `setAssignment` | The response tables ([response.md](response.md)), inline, `{include: file}` or `{legacy_file: response.root}`; plus `assignment: split \| draw`. |
| `smearing` | `setSmearing` | Reconstruction resolutions and decay-electron tagging, see below. |
| `energy` | `setEnergyEstimator` | `ccqe` (T2K-like) or `calorimetric` (NOvA-like): which estimator fills `Enu_rec`. |
| `binding_energy` | `setBindingEnergy` | Binding energy of the CCQE formula, GeV (default 0.027). |
| `samples` | `addSample` | The analysis samples, see below; list order defines `EventType`. Inline list or `{include: file}`. |
| `oscillation` | `setOscillation` | `{preset: asimov_A}` and/or `s12sq s13sq s23sq delta_cp dmsq21 dmsq32 ordering: normal\|inverted`. |
| `systematics` | `setSystematics`, `setDialCatalog` | The dial catalogue and the four parameter sets, see below ([systematics.md](systematics.md)). |
| `shifts` | `addShift` | Energy/momentum scale and offset systematics, see below. |
| `run.work_dir` | `setRun` | Directory under which every stage reads and writes. |
| `run.seed` | `setRun` | Master random seed; per-production seeds are derived from it and the production name unless a production sets its own. |
| `run.legacy_rng_order` | `setLegacyRngOrder` | Reproduce the legacy macros' random-number sequence exactly (regression use). |

### `detector`

Either a material or an explicit composition:

```yaml
detector: { material: water, fiducial_mass_kt: 188 }
detector:
  name: my_target
  fiducial_mass_kt: 40
  composition: [ { Z: 18, A: 40, mass_fraction: 1.0 } ]
  n_nucleons: 2.4e34            # optional: assert the nucleon count instead of deriving it
  neut_target: { num_nucleons: 40, num_bound_protons: 18, num_free_protons: 0 }   # optional override
```

Materials: `water` (O-16 16/18, H-1 2/18 by mass), `argon`, `scintillator` (CH),
`iron`. The nucleon count is `mass[g] × N_A` (the composition cancels): 188 kt
of water gives 1.13216e35; the legacy chain asserted 1.132e35, which
`n_nucleons` reproduces.

### `flux`

```yaml
flux:
  file: ${ZAKU_DATA}/data/flux/t2kflux_2020_plus250kA_nominal.root
  t2k_detector: sk              # shorthand for the T2K release names enu_sk_{numu,numub,nue,nueb}
  histograms: { numu: enu_sk_numu, numubar: enu_sk_numub }   # or an explicit species -> name map
  units: GeV                    # GeV (default) or MeV
  pot_norm: 1.0e21              # override flux_units for this file
  bin_norm_GeV: 0.05
```

Species names accept `numu`, `numubar`, `nue`, `nuebar`, `nutau`, `nutaubar`, the
T2K suffixes (`numub`, `nueb`), `anti-numu` style spellings, or PDG codes.

### `flux_covariance`

```yaml
flux_covariance:
  file: ${ZAKU_DATA}/data/flux/t2kflux_2020_covariance.root
  matrix: total_flux_cov        # a TMatrixT<double>, fractional
  far_detector_prefix: sk       # TAxis objects <prefix>_<numode|anumode>_<species>_bins define the blocks
  near_detector_prefix: nd280
  nd_constraint: 0.8            # 0 = unconstrained; s in C - s B^T A^-1 B otherwise
  units: GeV                    # unit of the TAxis edges
```

### `productions`

```yaml
productions:
  - { beam: FHC, flux: numu, interaction: numu, n: 200000 }          # disappearance production
  - { beam: FHC, flux: numu, interaction: nue,  n: 50000, seed: 7 }  # appearance: thrown on the numu flux
```

A production is one generator run; it is not an analysis sample (samples are
the reconstructed and selected categories defined under `samples`).
`interaction` defaults to `flux`. The production name is
`<flux>_x_<interaction>_<N>_<beam>` with `N` as `200K` for multiples of 1000
(`numu_x_nue_50K_FHC`); its files live in `<work_dir>/<flux>/` (see
[pipeline.md](pipeline.md)).

### `smearing`

```yaml
smearing:
  file: "${ZAKU_DATA}/data/hk/smearing.root"     # optional: the legacy h_dp_vs_p / h_dcos_vs_cos histograms
  lepton_momentum: { type: legacy_fractional, sigma: 0.024 }   # fractional | absolute | legacy_fractional | stochastic {a, b} | {hist: file:name} | {formula: ...} | none
  lepton_angle:    { type: gaussian_deg, sigma: 1.8 }          # sigma of cos(theta) = sin(theta) x sigma_theta (binned as the legacy file)
  hadronic_energy: { type: stochastic, a: 0.20, b: 0.15 }      # sigma/E = a/sqrt(E) (+) b
  electron_energy: { type: stochastic, a: 0.08, b: 0.03 }
  muon_energy:     { type: fractional, sigma: 0.035 }
  neutron_visible_fraction: 0.0
  decay_electron_tag: { mu_plus: 0.96, mu_minus: 0.83, pi_plus: 0.45 }
```

See [physics.md](physics.md) for the formulas.

### `samples`

```yaml
samples:
  - name: FHC_1Rmu               # identifier; index in the list = EventType
    title: "FHC 1R#mu"           # display title (ROOT latex)
    beam: FHC                    # FHC | RHC | any
    column: 1Rmu                 # response-table column feeding it (default: name)
    hypothesis: muon             # muon | electron | true: lepton assumed by the energy estimator
    binning: { file: "${ZAKU_DATA}/data/binning/binning_T2K_FHC_MuLike.txt" }   # GUNDAM txt binning
    cut: "Enu_rec < 30"          # optional reco-level cut (TFormula), also written into the GUNDAM config
  - name: FHC_numuCC
    beam: FHC
    column: numuCC
    binning:                     # or an inline Cartesian product (last axis fastest)
      variables: [Enu_rec, Theta_rec]
      edges: { Enu_rec: { min: 0, max: 3, step: 0.05, extra: [3.5, 4, 5, 10, 30] }, Theta_rec: [0, 20, 40, 60, 80, 100, 180] }
```

`configs/samples/t2k_like.yaml` and `configs/samples/nova_like.yaml` are shipped
sample sets to include.

### `systematics`

```yaml
systematics:
  xsec:
    dials: { include: dials/neut_default.yaml }   # the dial catalogue (or an inline list of {name, tweaks, fit, note})
    exclude: [MDLSPiEj, NucleonFSI_PB, SRCFrac_C, RadCorr_rate]   # not fitted, on top of the catalogue's fit: false
    write_unused_dials: true      # spline files carry every dial of the weight file (false: fitted only)
    enabled: true
  flux: { enabled: true }         # far-detector block of flux_covariance; nd_constraint is set there
  detector:
    mode: normalisation           # normalisation (one per sample) | shifts (from the shift parameters)
    sigmas: { FHC_1Rmu: 0.021, FHC_1Re0de: 0.031, FHC_1Re1de: 0.132, RHC_1Rmu: 0.019, RHC_1Re0de: 0.039 }
    names:  { FHC_1Re0de: DET_NORM_FHC_1RE_0DE }   # optional; default DET_NORM_<SAMPLE NAME UPPER-CASED>
    correlation: 0.0              # rho between the two samples of each pair
    pairs: [[FHC_1Rmu, RHC_1Rmu], [FHC_1Re0de, RHC_1Re0de]]
    range: [0.0, 2.0]             # GUNDAM parametersRange
  oscillation:                    # priors are the experiment's oscillation point
    sigmas: { s12sq: 0.013, s13sq: 0.0007, s23sq: 0.021, dmsq21: 1.8e-6, dmsq32: 3.4e-5, delta_cp: 1.17, ordering: 0.1 }
    fixed: [PMNS_SIGN_MASS_SQUARED_32]
  grid: { bins: 2000, min: 0.0, max: 10.0 }   # true-energy grid of the tabulated oscillation dials, GeV
```

Every sample needs an entry in `detector.sigmas` when the mode is
`normalisation` (`validate` checks). The dial catalogue file:

```yaml
generator: neut
dials:
  - { name: MaCCQE, tweaks: [-2, -1, 1, 2], note: "CCQE axial mass" }
  - { name: MDLSPiEj, tweaks: [-3, -2], fit: false, note: "model switch" }
```

### `shifts`

```yaml
shifts:
  - { name: lepton_momentum_scale, target: lepton_momentum, kind: scale,  sigma: 0.024, apply_to: all }
  - { name: em_energy_offset,      target: lepton_energy,   kind: offset, sigma: 0.01,  apply_to: electron }
```

`target`: `lepton_momentum | lepton_angle | lepton_energy | hadronic_energy | reco_energy`;
`kind`: `scale` (sigma as a fraction) or `offset` (GeV, or degrees for the
angle); `apply_to`: `all | muon | electron` (by the sample's hypothesis).

## Example

`configs/hk_t2k_like.yaml` is the HK-on-T2K-beam experiment; `configs/legacy_hk_mock.yaml`
reproduces the legacy production bit for bit; `configs/tiny_test.yaml` is the
same experiment with 2000 events per production.
