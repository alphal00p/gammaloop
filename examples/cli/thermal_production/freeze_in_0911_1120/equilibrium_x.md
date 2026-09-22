# Equilibrium LU production benchmark

Equilibrium LU residue selection now retains the mixed-sign scattering cuts.
The native forward-X calculation agrees with the complete Bose production
source at all three tested momenta, within 1.4 Monte Carlo standard errors.
The paired zero-temperature equilibrium and vacuum integrals are identical.

The model remains `model.json`: four distinct real scalars with interaction
`lambda X B1 B2 B3`, massless bath particles, and `MX = T = lambda = 1`.
The saved `equilibrium_x.toml` generates the inclusive forward sunset through
`generate xs X > B1 B2 B3 FI^2==2 [{{2}}]`, with X's incoming momentum fixed
at each requested production momentum. It uses native LU integration throughout.

## Results and normalization

The native integral is the spectral rate `Gamma_rho`. The production source is

\[
S(p)=\left.\frac{df_X(p)}{dt}\right|_{f_X=0}
=\frac{\Gamma_\rho(p)}{\exp(E_X/T)-1},\qquad
\frac{d\gamma}{d^3p}=\frac{S(p)}{(2\pi)^3}.
\]

The factor at E_X is the KMS relation for the bath correlator, not an assumed
population of X. GammaLoop already includes the one-particle flux `1/(2 E_X)`.
The table compares one million samples per process and momentum against the
independent Bose phase-space quadrature documented in `inclusive_x.md` and
preserved in `inclusive_comparison.json`.

| p_X/T | Native S(p) | Bose reference, scattering + inverse decay | Pull |
|---:|---:|---:|---:|
| 0 | 0.0042083491 +/- 0.000138 | 0.0042290850517 | -0.15 |
| 1 | 0.0014833704 +/- 2.5e-05 | 0.0015012405272 | -0.72 |
| 3 | 6.7977695e-05 +/- 1.38e-06 | 6.9892337357e-05 | -1.39 |

At rest, the unconverted native spectral rate is
`0.007231129847 +/- 0.000237683356`, versus `0.007266759995`.
Before this change, native LU returned only the decay contribution,
`0.002544820884 +/- 0.000024397616`.
Detailed counters record one NaN evaluation at p=0 and one at p=1, out of
three million process evaluations at each momentum; p=3 has none. These counts
also occurred in the single-root pilot, and their cause has not been isolated.
The console rounds this incidence to 0.00%. The one-million-sample
cap was reached before the card's 0.1% target accuracy; quoted uncertainties
are the actual Monte Carlo errors. Vacuum rates agree with
`lambda^2 MX^2/(512 pi^3 E_X)` within 0.6 standard errors.

[Hall et al., section 6.3](https://arxiv.org/pdf/0911.1120) uses
Maxwell-Boltzmann occupations and counts the three scattering channels. Its
corresponding differential source and integrated reaction density are

\[
S_{2\to2}^{\rm MB}(p)=\frac{3\lambda^2 T^2}{64\pi^3 E_X}e^{-E_X/T},
\qquad
\gamma_{2\to2}^{\rm MB}=\frac{3\lambda^2 m_X T^3}{128\pi^5}K_1(m_X/T)
=4.6098997578035\times10^{-5}.
\]

The native calculation uses Bose occupations and includes inverse three-body
decay as well as scattering. Its comparison therefore uses the full Bose
reference. The derivation and separate Bose/MB channel values are recorded in
`inclusive_x.md` and the comparison JSON files. The paper's MB integrated number
above is an analytic reference, not a native MB run.

## What changed and what was checked

The generated thermal CFF is exactly equal to the previously saved
`inclusive_cff.json`. The first changed boundary is LU residue selection:
the standalone representation now has seven signed cut groups, selecting
orientation IDs `[7, 3, 4, 5, 2, 6, 1]` once each. IDs 3, 5, and 6 supply the
three scattering channels; ID 7 supplies decay. IDs 1, 2, and 4 have no physical
root for this model and incoming timelike X momentum.

The shared surface selector handles E- and H-surface aliases, raised poles,
numerator cancellation, and denominator signs while preserving thermal weights.
Mixed-sign cuts use all positive radial crossings. The shared solver compactifies
`t = u/(1-u)` and uses derivative enclosures to isolate crossings before
bracketed Newton refinement. An unresolved or tangent crossing is reported as
an error. Each root contributes with an absolute phase-space Jacobian; the
canonical signed derivatives remain in the raised-residue calculation.

The two-root regression uses `2 sqrt(t^2+1) - 1.9 t - 1 = 0`, whose two positive
roots are `(3.8 +/- sqrt(9.76))/0.78`. A single Newton solve misses one branch.
Other regressions cover the massless analytic roots, a boosted scattering root,
quad precision, absent roots, raised H-surface residues, and existing vacuum LU
residue factorization. All 24 selected tests passed. Formatting, cargo check,
and clippy passed without warnings. The known unrelated failing thermal tests
were not run or changed.

`equilibrium_events.json` records the inspection at loop momenta `(3,0,0)` and
`(0,4,0)`: one decay event and three scattering events. Absorbed bath legs are
incoming; four-momentum is conserved in all four events, and their weights sum
to the returned integrand. Thermal threshold and IR subtraction are outside
this change; both are unnecessary for this benchmark. Saved states must be
regenerated because the signed cut-group layout uses state manifest version 7.

## Reproduce

From the repository root, build the CLI and run:

```sh
cargo build -p gammaloop-api --bin gammaloop --features ufo_support --profile dev-optim
./gammaloop --read-only-state -s /tmp/gammaloop-freeze-in-equilibrium-state \
  -n -L off \
  examples/cli/thermal_production/freeze_in_0911_1120/equilibrium_x.toml \
  run generate export integrate p1 p3
```

The targeted regression command was:

```sh
cargo nextest run -p gammalooprs -p three-dimensional-reps --lib \
  --cargo-profile dev-optim -P test_gammaloop \
  -E 'test(equilibrium_lu_sunset) or test(signed_lu_sunset) or test(signed_root_isolation) or test(raised_hsurface_residue) or test(lu_residue_selection) or test(indexed_residue_selection) or test(scalar_cut_residue) or test(cut_groups_include_all_physical) or test(safeguarded_newton) or test(exact_lu_cut_matches) or test(exact_raised_lu_cut_matches) or test(exact_cff_lu_residue_factorizes) or test(exact_raised_lu_residue_factorizes) or test(exact_cff_uncancelled_powered_denominator) or test(exact_cff_quartic_lu)'
```

For the retained-event check, use the same card with `run generate inspect`.
The recorded scan ran its three independent momentum blocks concurrently,
with two workers per momentum, separate read-only scratch states, and the
saved seeds. The card uses the interpreted summed function-map evaluator with
UV and threshold subtraction disabled.

Artifacts beside this report:

- `equilibrium_x.toml`: run card; `model.json`: custom model.
- `equilibrium_results.json`: all raw native results and numerical-stability statistics.
- `equilibrium_comparison.json`: converted sources, reference values, pulls,
  selector/event checks, and executable hash.
- `equilibrium_standalone.json`: all selected thermal LU residues.
- `equilibrium_events.json`: retained events at the analytic inspection point.

The earlier `inclusive_*` artifacts are retained as the before-change evidence.
