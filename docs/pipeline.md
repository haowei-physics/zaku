# Pipeline and file contracts

Zaku's stages read and write ROOT files below `run.work_dir`. Zaku never
executes the external tools: the *prepare* stages write their inputs and a
manifest of command lines, you run them, and the *processing* stages consume
the results.

```
prepare-generation ─┐   (you run)   ┌─ recon ─ splines ─ select ─ covariances ─ gundam-config
prepare-weights    ─┼──► manifests ─┤
prepare-flatten    ─┘                └─ summary
```

## Directory layout

A *production* is one generator run (a beam mode, the flux species thrown, the
species that interacts, an event count); a *sample* is a reconstructed and
selected event category. For the production `numu_x_nue_50K_FHC`:

| File | Written by | Content |
|---|---|---|
| `numu/numu_x_nue_50K_FHC.root` | the generator (you) | NEUT `neuttree` (or `.ghep.root` for GENIE) |
| `numu/numu_x_nue_50K_FHC_flattree.root` | `nuisflat` (you) | NUISANCE GenericVectors flat tree `FlatTree_VARS` + `FlatTree_FLUX`, `FlatTree_EVT` |
| `numu/numu_x_nue_50K_FHC_Weight.flat.root` | T2KReWeight `EventWeightGen` (you) | `T2KRWFlatTree`, `TweakMetaData` |
| `numu/numu_x_nue_50K_FHC_rw_<dial>.root` | GENIE `grwght1p` (you) | one tree per dial |
| `numu/numu_x_nue_50K_FHC_fakerecon.root` | `recon` | flat tree clone + response, smearing and energy branches |
| `numu/numu_x_nue_50K_FHC_spline.root` | `splines` | recon tree + one `TGraph` branch per dial |
| `fakeselect_FHC.root`, `fakeselect_FHC_rates.csv` | `select` | the `FakeSelect` tree GUNDAM reads, and the rate table |
| `parameters/{oscProb,flux,det,xsec}/…` | `covariances` | parameter-set files |
| `manifests/<stage>.txt`, `.json` | prepare stages | command lines for you to run |
| `gen/<production>/` | prepare-generation | per-job working directory (cards, seed files) |
| `flux/uniform_<beam>.root` | prepare-generation | generator-ready uniform-bin flux copies |

`zaku::Layout` is the one place these names are defined.

## The flat tree

The NUISANCE GenericVectors tree is the generator-neutral interchange format.
Zaku reads it with `zaku::FlatTreeCursor` (classic branch addresses, so an
output tree cloned with `cloneTo()` copies every input branch on `Fill()`).
Required branches: `Mode cc PDGnu Enu_true PDGLep ELep CosLep Q2 fScaleFactor
nfsp pdg px py pz E flagCCQE`; everything else is optional. Energies are GeV,
`E[]` includes the rest mass, and the final-state list contains the outgoing
lepton (for NC, the outgoing neutrino). `fScaleFactor` is ⟨σ⟩/N_gen in
cm²/nucleon.

Events with `nfsp == 0` (NEUT single-pion/gamma modes whose decay produced no
particles: 1505 of the 200000 reference events) have no weight row and are
dropped by the recon stage.

`zaku::TruthVariables` derives the quantities response tables and cuts can use:
`PLep_true`, `Ehad_true` (Σ KE of protons, Σ E of charged pions and kaons, Σ E
of π⁰/γ/η, a configurable fraction of neutron KE; the primary lepton and all
neutrinos excluded), `Evis_true`, `W2`, the final-state multiplicities
`nFS_pip nFS_pim nFS_pi0 nFS_p nFS_n nFS_gamma nFS_mup nFS_mum nFS_ep nFS_em
nFS_kp nFS_km nFS_k0 nFS_other nFS_pi`, and the production context `isRHC`,
`flux_pdg`, `isAppearance`.

## The weight file and the join

T2KReWeight's `EventWeightGen` writes `T2KRWFlatTree` with one
`vector<double>` per dial (the weight at each knot), the nominal weight
`cv_weight`, and the event kinematics; `TweakMetaData` holds, per dial, the
knot values (`<dial>`) and the generated value (`<dial>_from`).
`zaku::T2KRWWeightCursor` reads both (using the highest key cycle; the files
hold two).

The weight tree has no event index and skips the empty-final-state events, so
`zaku::EventJoin` matches rows to flat-tree entries with a forward cursor on
(`Mode`, `PDGnu`, neutrino energy to 1e-6 relative). Lepton quantities are not
part of the key: in 13 multi-lepton DIS events of the reference production
T2KReWeight picked a different final-state lepton than NUISANCE. The join is
strict: a weight row that does not match the current event, or a non-empty
event without a row, aborts with the two file names. Sources with an explicit
event number (GENIE) are matched on it directly.

## Normalisation

```
POTWeight = fScaleFactor × Φ/POT × POT × N_nucleons
Φ/POT     = Integral("width") / bin_norm_GeV / pot_norm      (neutrinos / cm² / POT)
```

with `Experiment::potWeight()`. For the reference production
(`fScaleFactor` = 4.4568e-44, Φ/POT = 1.20736e-14 for the FHC ν_μ flux,
6.75e21 POT, 1.132e35 nucleons) this is 0.411157.

## Provenance

Every file Zaku writes carries a `zaku/` directory with `version`, `schema`,
`stage` and `provenance` (`TNamed` objects). Readers refuse files written with a
newer schema than the library.
