# Systematics

Zaku produces the four GUNDAM parameter sets of the legacy chain, the
per-event cross-section response graphs they need, and a summary of what each
parameter does to each sample. Everything is driven by the `systematics:` and
`shifts:` blocks of the experiment ([config_reference.md](config_reference.md))
and written under `run.work_dir` ([pipeline.md](pipeline.md)).

```
recon files + weight files ──► splines ──► spline files ──► select ──► FakeSelect (with graphs)
dial catalogue + weight file knots ─┐                                       │
flux covariance file ───────────────┼─► covariances ──► parameters/…        ├──► summary ──► summary/…
samples + detector sigmas ──────────┤                                       │
oscillation point + sigmas ─────────┘                                       │
```

## The dial catalogue

`systematics.xsec.dials` lists the reweighting dials of the generator
(`configs/dials/neut_default.yaml`): name, tweak points in units of the dial's
one sigma, whether the dial is fitted, a note. It is the single list of
cross-section systematics: the reweight stage writes it into the engine's card,
the splines stage builds one response graph per dial, the covariances stage
turns the fitted dials into the cross-section parameter set. `exclude:` takes
dials out of the fit on top of `fit: false` (the legacy fit excluded
`MDLSPiEj`, `NucleonFSI_PB`, `SRCFrac_C`, `RadCorr_rate`). `zaku <exp> list-dials`
prints it.

## Splines (`Experiment::makeSplines`, `zaku <exp> splines [--production …]`)

For each production the recon file and the weight file are joined event by
event ([pipeline.md](pipeline.md), "The weight file and the join") and the recon
tree is copied with one unsplit `TGraph` branch per dial. GUNDAM reads the
graph as a `TObject` from the tree it takes events from, so the graphs must
live in that tree. The rules are those of the legacy `MakeSplines.C`:

- Knots are the dial's tweak points sorted ascending, with the nominal inserted
  at weight 1 unless it already is a tweak point (MaCCQE: 0.818, 1.014, *1.21*,
  1.406, 1.602). Without it the interpolated response at the prior would not be
  exactly 1.
- Every stored weight exactly 1: the dial does not touch this event; an empty
  graph is written (GUNDAM treats an empty point list as "no dial", at a
  fraction of the file size).
- Every stored weight 0 (whole events where T2KReWeight's `cv_weight` is 0), a
  negative or non-finite weight, or a weight at a tweak equal to the nominal
  that differs from 1 by more than 1e-6: the response is unusable and an empty
  graph is written. A graph flat at 0 must never be written: GUNDAM would read
  it as a genuine rescaling by 0 and silently delete the event.
- A response flat at a value other than 1 is a real constant rescaling and is
  written in full.

`write_unused_dials: true` (default) writes every dial of the weight file, as
the legacy chain did; `false` writes the fitted dials only. The spline file
carries the knots of every dial in `zaku/dials`, so the covariances stage can
run without the weight files. The stage prints, per dial, how many events had
no response and how many were rejected.

## Covariances (`Experiment::makeCovariances`, `zaku <exp> covariances`)

Every parameter set is a `zaku::ParameterSet` and is written as the legacy
macros wrote it, so the GUNDAM configuration of Stage 4 and the legacy one
read the same objects:

| Object | Type | Content |
|---|---|---|
| `<prefix>_param_names` | `TObjArray` of `TObjString` (single key) | parameter names, in matrix order |
| `<prefix>_param_priors` | `TVectorD` | priors |
| `<prefix>_param_cov` | `TMatrixD` | covariance (never `TMatrixDSym`: GUNDAM's reader prefix-matches the class name) |
| `<prefix>_param_lower`, `_upper` | `TVectorD` | bounds, only when a parameter has any (cross section) |

with `<prefix>` = `osc`, `det`, `xsec`, `flux`. Where a binning file goes with
a set, GUNDAM pairs dial *i* with parameter *i* by position, so both come out of
the same loop.

**Oscillation** (`parameters/oscProb/oscCov.root`): the seven PMNS parameters
GUNDAM's tabulated NuOscillator dials take, `PMNS_SIN_SQUARED_12/13/23`,
`PMNS_DELTA_MASS_SQUARED_21/32`, `PMNS_DELTA_CP`, `PMNS_SIGN_MASS_SQUARED_32`,
priors at the experiment's `oscillation` point (sign ±1 by ordering), diagonal
with the widths of `systematics.oscillation.sigmas` (defaults: the legacy table,
PDG-like for the solar and reactor parameters, loose for what the experiment
measures). `fixed:` names parameters the fit must not move (the ordering, by
default). `parameters/oscProb/grid.root` holds the true-energy grid
(`energyBinning`, `systematics.grid`, GeV) the tabulated dials are filled on.

**Detector** (`parameters/det/detCov.root`, `detBinning.txt`), mode
`normalisation`: one normalisation per sample, prior 1, sigma from
`systematics.detector.sigmas` (the HK numbers are arXiv:2505.15019 Tab. 3,
detector row, T2K-systematics column), an optional `correlation` between the
samples of each entry of `pairs` (0.98 is what the FHC/RHC 1Re ratio column of
that table implies), positive-definiteness checked. Parameter names default to
`DET_NORM_<SAMPLE NAME UPPER-CASED>`; `names:` overrides them. The binning has
one point row `EventType` per parameter, in sample order. This scales samples;
it cannot tilt a spectrum. Mode `shifts` builds the detector covariance from the
shift systematics instead (below).

**Cross section** (`parameters/xsec/xsecCov.root`): one parameter per fitted
dial, prior = the value the events were generated at (`<dial>_from` in the weight
file), sigma recovered from the knots: every knot must give the same
`(knot − nominal) / tweak multiple` to 1e-3, otherwise the catalogue and the
weight file disagree and the stage stops. Bounds are the outer knots (the
nominal included), beyond which the spline has no response. Diagonal. The knots
come from the first production whose weight file exists, else from the
`zaku/dials` metadata of a spline file.

**Flux** (`parameters/flux/fluxCov.root`, `fluxCov_ND.root`, `fluxBinning.txt`):
the far-detector block of the released fractional covariance
(`flux_covariance`), priors 1, names `FLUX_<FHC|RHC>_<NUMU|NUMUB|NUE|NUEB>_<bin>`,
one row per (beam mode, species, true-energy bin) in the release's block
order; the binning file has `isRHC` and `flux_pdg` as point variables and
`Enu_true` as a range. Two matrices are written: `C`, the far block on its own,
and `C − s·BᵀA⁻¹B` with `s = flux_covariance.nd_constraint`, `A` the
near-detector block and `B` the cross block: the Gaussian posterior on the far
block after a near-detector measurement of strength `s ∈ [0, 1]`, a convex
combination of the unconstrained matrix and the perfect-near-detector limit,
positive definite for any `s`. At `s = 0.8` the HK normalisations run 2.1–7.4%
against 4.7–15.2% unconstrained. This is not a near-detector fit (no
cross-section degeneracies, no off-axis flux); the difference between fits with
the two files brackets what the flux systematic costs.

## Shift systematics

`shifts:` lists detector systematics that move reconstructed quantities rather
than scaling counts: a `scale` (sigma as a fraction) or an `offset` (sigma in
GeV, or degrees for the angle) of `lepton_momentum`, `lepton_angle`,
`lepton_energy` (the calorimetric lepton term), `hadronic_energy` or
`reco_energy`, applied to `all` samples or to the `muon` / `electron`
hypothesis ones. `ShiftInjector::apply` moves the quantity and re-derives
everything downstream (`Theta_rec`, both CCQE estimators, the calorimetric
sums) so a shifted universe is one consistent reconstruction.

The summary evaluates every shift at ±1σ. With `systematics.detector.mode: shifts`
the detector covariance is built from these universes: for every (sample,
`Enu_rec` bin) the fractional half-difference δ = (N₊ − N₋) / 2N, summed as
Σ δδᵀ over the shifts (fractional, priors 1, a 1e-8 floor on the diagonal so the
matrix inverts), with the binning `EventType Enu_rec Enu_rec` (the sample's
`Enu_rec` edges when it has them, else a default set). This is a linearised,
bin-normalisation approximation of an event migration; GUNDAM's binned graph
dials would be the exact representation and are a later option.

## Summary (`Experiment::summariseSystematics`, `zaku <exp> summary`)

Reads the FakeSelect files (which must carry the response graphs, i.e. `select`
must have run on spline files) and computes, per sample, the nominal spectrum
in the sample's first binning variable (GUNDAM weight
`POTWeight × OscWeight × SelEffWeight × SampleProb`; only events inside the
binning range count, as in the fit, so the nominal sits slightly below the
rate table when some events have `Enu_rec = −1`) and the spectrum at ±1σ of
every parameter:

| Set | ±1σ universe |
|---|---|
| xsec | the event's graph evaluated (linearly) at prior ± σ, clamped to the knot range; empty graph = 1 |
| flux | 1 ± σ on the events of the parameter's (beam mode, species, `Enu_true`) bin |
| det | 1 ± σ on the events of the sample |
| osc | the oscillation weight recomputed with the parameter moved (CC events; `P'/P`) |
| shift | the reconstruction re-derived with the shift applied; events move between bins |

Output: `summary/systematics_summary.csv`
(`set,parameter,sample,nominal,plus,minus,plus_frac,minus_frac`), the same as
`SummaryTable` in memory (`text()`, `top(n)`), and
`summary/systematics_summary.root` with one directory per sample holding
`nominal` and `<parameter>_plus/_minus` `TH1D`s (GeV axes).
`examples/macros/plot_systematics.C` draws them.

## Verification

`regression_syst_reference` (needs the reference production) checks that the
splines stage reproduces every graph of the reference `_spline.root` of
`numu_x_numu_200K_FHC` (198495 entries × 19 dials, point counts and values to
1e-12), that the covariances stage with `configs/legacy_hk_mock.yaml`
reproduces the names, priors, sigmas, bounds and matrices of the four reference
parameter files and the two binning files to 1e-9, and that `select` on the
reference spline files followed by `summary` gives 77 parameters × 5 samples
with finite effects, the detector normalisations moving their own sample by
exactly ±σ.
