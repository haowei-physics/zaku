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
- Sections this build does not interpret yet (`samples`, `response`, `smearing`,
  `energy`, `generator`, `shifts`, `oscillation`, `systematics`, `gundam`) are
  kept verbatim and written back by `toYaml()`.

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

## Example

`configs/hk_t2k_like.yaml` is the HK-on-T2K-beam experiment; `configs/legacy_hk_mock.yaml`
reproduces the legacy production bit for bit; `configs/tiny_test.yaml` is the
same experiment with 2000 events per production.
