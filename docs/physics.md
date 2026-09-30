# Physics conventions and formulas

All energies, momenta and masses are GeV; `Theta_rec` is in degrees
([units](units.md)). Where a formula reproduces the legacy `mock_lbl_inputs`
macros, the floating-point operation order is kept so that the reference
production is reproduced bit for bit.

## Normalisation

```
POTWeight = fScaleFactor × (Φ/POT) × POT × N_nucleons
Φ/POT     = Integral("width") / bin_norm_GeV / pot_norm      [ν / cm² / POT]
```

`fScaleFactor` is NUISANCE's ⟨σ⟩/N_gen in cm² per nucleon; `N_nucleons` is
`mass[g] × N_A` (188 kt of water: 1.13216e35; the legacy chain asserted
1.132e35, reproduced with `detector.n_nucleons`).

## Oscillation weight

NuFast (P. Denton, vendored under `third_party/NuFast`) with the matter effect,
`N_Newton = 0`, the experiment's baseline and ρ·Y_e (HK: 295 km, 2.6 × 0.5).
`Δm²31 = Δm²32 + Δm²21` (negative for inverted ordering). Rules per event:

- CC: `OscWeight = P(flux flavour → interaction flavour)` at the event's true
  energy, with the sign of the flux species selecting antineutrinos.
- NC on a disappearance production: 1 (the NC rate is flavour blind).
- NC on an appearance production: 0 (already counted by the disappearance
  production of the same flux).

Default point: the T2K "Asimov A" set (sin²θ12 = 0.307, sin²θ13 = 0.0218,
sin²θ23 = 0.528, δCP = −1.601, Δm²21 = 7.53e-5 eV², Δm²32 = 2.509e-3 eV², normal ordering).

## Response

`RespCat` is the first matching category, `SelEffWeight` the table's efficiency
x, and `SampleProb[c]` the column probabilities y_c ([response.md](response.md)).
Lookups clamp to the table range.

## Smearing and decay electrons

The outgoing-lepton momentum `p = sqrt(ELep² − m²)` (`ELep²` in float precision,
as the legacy macro; `m` by the true lepton PDG, 0 for the outgoing NC neutrino)
and `cos θ` are smeared with Gaussians of width given by the `smearing` block:

| Model | σ |
|---|---|
| `fractional` | f × value |
| `absolute` | constant |
| `legacy_fractional` | histogram of f × bin centre, 50 bins on 0–5 GeV, clamped (the legacy `h_dp_vs_p`) |
| `gaussian_deg` | histogram of sin θ × σ_θ, 40 bins on −1..1 in cos θ, clamped (the legacy `h_dcos_vs_cos`) |
| `stochastic` | sqrt((a √E)² + (b E)²), i.e. σ/E = a/√E ⊕ b |
| `hist` | any TH1 binned in the true value |
| `formula` | any expression of `value` and the truth variables |

`p` is clamped at 0, `cos θ` to [−1, 1]. Decay electrons: each μ⁺, μ⁻ and π⁺
among the first 100 final-state particles is a parent (`true_ME`) and is tagged
with probability 0.96, 0.83 (capture included) and 0.45 (`recon_ME`).

Random numbers: one `TRandom3` per production ("legacy stream", seed 1234 in
legacy mode, else derived from the run seed and the production name) draws, per
kept event, Gaus(p), Gaus(cos θ) and one Uniform per decay-electron parent, in
that order, exactly as the legacy macro. A second stream draws everything the
legacy chain did not have (hadronic and calorimetric lepton energy smearing),
so the legacy quantities stay reproducible whatever else is enabled.

## Energy estimators

Both are computed for every event and written to the recon tree; the sample's
`hypothesis` (`muon`, `electron`, or `true` = the true lepton flavour) and the
experiment's `energy` key decide what `Enu_rec` is.

**CCQE (T2K-like)**, under lepton mass m_l, struck nucleon n (neutrino mode) or
p (antineutrino mode), binding energy E_b = 0.027 GeV:

```
E_l   = sqrt(p² + m_l²)            m_eff = m_struck − E_b
E_ν   = (m_out² − m_eff² − m_l² + 2 m_eff E_l) / (2 (m_eff − E_l + p cos θ))
Q²    = −m_l² + 2 E_ν (E_l − p cos θ)
```

with `p²` and `p cos θ` in float precision and `Q²` using the float-rounded
`E_ν`, as the legacy macro. No positive solution ⇒ `Enu_rec = Q2_rec = −1`.
`Theta_rec = acos(cos θ)` in float precision, in degrees.

**Calorimetric (NOvA-like)**: `Enu_rec = Elep_rec + Ehad_rec`, where
`Elep_rec = sqrt(p² + m_l²)` optionally smeared (`muon_energy`,
`electron_energy`), 0 for NC events, and `Ehad_rec` is the smeared
(`hadronic_energy`) visible hadronic energy

```
Ehad_true = Σ KE(p) + Σ E(π±, K±) + Σ E(π⁰, γ, η) + f_n Σ KE(n)   (+ KE of secondary μ±, E of e±)
```

with the primary lepton (or outgoing NC neutrino) and all neutrinos excluded
and `f_n = neutron_visible_fraction` (default 0). `EhadFrac_rec = Ehad_rec / Enu_rec`.

## Sample assignment

For each event the recon stage stores the column probabilities y_c of its
table; the select stage maps columns to the samples of the beam mode and
either splits the event (`split`: one entry per column with sample and y_c > 0,
`SampleProb = y_c`, plus a rejected entry carrying the remainder) or draws once
(`draw`: `r = Uniform()`, the first column with r below the running sum wins;
a column without a sample, or a failed sample `cut`, is rejected; `SampleProb = 1`).
The GUNDAM weight is `POTWeight × SelEffWeight × SampleProb × (cc || !isAppearance)`.
Rejected entries are reconstructed under the electron hypothesis (legacy).
The rate table (`fakeselect_<mode>_rates.csv`, `rates` TH2 in the select file)
sums `POTWeight × OscWeight × SelEffWeight × SampleProb` per sample and production:
the expected event counts at the exposure.
