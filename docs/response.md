# Response tables

Zaku does not simulate a detector. The detector design enters as a **response**:
a map from true-event variables to the probability that the event is detected
and the probabilities that it lands in each output column. The legacy
`mock_lbl_inputs` mixing matrix (nine interaction categories × five
neutrino-energy bins) is one instance of it; the same machinery accepts any
number of categories, any truth variables as table axes, and constants, inline
arrays, ROOT histograms or formulas as sources.

## Concepts

- **Category** (`RespCat`): an ordered list of named selectors over the truth
  variables; the first one that is true owns the event. Selectors are ROOT
  `TFormula` expressions (`cc && PDGnu == 14 && flagCCQE`).
- **Table**: for one category (or an inline `select`), optionally one beam mode,
  a set of axes (truth variables with bin edges, GeV where kinematic), an
  `efficiency` x and one probability y_c per **column**.
- **Efficiency** x is applied as a weight (`SelEffWeight`), never drawn. The
  legacy HK tables use it as a tuning factor and it exceeds 1 there.
- **Columns** are matched to the experiment's samples by name and beam mode: a
  sample `{name: FHC_1Rmu, beam: FHC, column: 1Rmu}` takes column `1Rmu` in
  neutrino mode. A column with no sample in the current beam mode is rejected
  probability mass (the legacy RHC veto of `1Re1de` falls out of this).
- **Assignment**: `split` writes one FakeSelect entry per column with y_c > 0,
  weighted by y_c (deterministic, no extra MC noise); `draw` makes one random
  draw per event against the running sum of the columns (the legacy behaviour).
- Rows whose probabilities sum above 1 are an error unless
  `normalise_probabilities: true` (the legacy tables were renormalised when built).

## YAML

```yaml
categories:
  - { name: numuCCQE,    select: "cc && PDGnu == 14 && flagCCQE" }
  - { name: numuCCnonQE, select: "cc && PDGnu == 14 && !flagCCQE" }
  - { name: NC,          select: "!cc" }
tables:
  - category: numuCCQE                 # or select: "<formula>"; optional beam: FHC | RHC
    axes: { Enu_true: [0, 0.4, 0.6, 1.2, 3.0, 30] }            # 1-3 truth variables, GeV
    efficiency: [1.08303, 1.19671, 1.12235, 0.775, 0.40954]  # x per bin (last axis fastest)
    probabilities:                                          # y_c per bin, in draw order
      1Rmu:   [0.74241, 0.90496, 0.93853, 0.94951, 0.94999]
      1Re0de: [0.00525, 0.00239, 0.00098, 0.00010, 0.00001]
      1Re1de: [0.00061, 0.00024, 0.00009, 0.00001, 0.0]
  - category: NC
    axes: { Enu_true: [0, 0.4, 0.6, 1.2, 3.0, 30], W: [0, 1.2, 1.6, 30] }   # 2-D example
    efficiency:   { hist: "data/hk/response.root:h_eff_NC", units: GeV }    # TH1/TH2/TH3 over the axes
    probabilities: { 1Rmu: { hist: "data/hk/response.root:h_mix_NC_1Rmu" } }
  - select: "cc && abs(PDGnu) == 12 && Q2 > 0.2"
    efficiency:   "0.9 * (1 - exp(-Enu_true / 0.4))"                        # formula over truth variables
    probabilities: { 1Re0de: "0.85", 1Re1de: 0.05 }
normalise_probabilities: false
```

Sources: a number is a constant; a list is an inline array over the axes; a
string is a formula; `{hist: "file.root:name"}` is a histogram whose axes are
the table axes in order (lookups clamp to the histogram range, like the legacy
`FindBin` clamp). Paths are relative to the YAML file; `${ZAKU_DATA}` works.
A file holding only `legacy_file: response.root` loads the legacy histograms
directly.

In the experiment file: `response: { include: response/hk_2025.yaml, assignment: split }`
or the tables inline. `configs/response/hk_2025.yaml` is the legacy HK response
exported by `zaku::Response::fromLegacyFile(...).toYaml()`;
`configs/response/example_calorimeter.yaml` is an illustrative NOvA-like one.

## Truth variables available to selectors, axes and formulas

Flat-tree scalars: `Mode cc PDGnu Enu_true tgt tgta tgtz PDGLep ELep CosLep Q2 q0 q3 Enu_QE
Q2_QE W bjorken_x bjorken_y fScaleFactor nfsp flagCCINC flagNCINC flagCCQE flagCC0pi
flagCCQELike flagNCEL flagNC0pi flagCCcoh flagNCcoh flagCC1pip flagNC1pip flagCC1pim
flagNC1pim flagCC1pi0 flagNC1pi0` (the Bjorken `x`, `y` are renamed because TFormula
reserves those letters). Derived: `PLep_true Ehad_true Evis_true W2 nFS_pip nFS_pim
nFS_pi0 nFS_p nFS_n nFS_gamma nFS_mup nFS_mum nFS_ep nFS_em nFS_kp nFS_km nFS_k0
nFS_other nFS_pi isRHC flux_pdg isAppearance`.

## C++ and Python

```cpp
zaku::Response r = zaku::Response::fromFile("configs/response/hk_2025.yaml");
zaku::VariableSet vars;               // the variable dictionary
r.prepare(vars);                      // compile selectors, bind tables
vars.loadTruth(event, zaku::TruthVariables::compute(event, ctx));
zaku::ResponseResult res = r.evaluate(vars, zaku::BeamMode::FHC);
// res.category (RespCat), res.efficiency (SelEffWeight), res.probabilities (per column)
```

`Experiment::setResponse(r)` installs it; `zaku <exp.yaml> show-response [--yaml]` prints it.
