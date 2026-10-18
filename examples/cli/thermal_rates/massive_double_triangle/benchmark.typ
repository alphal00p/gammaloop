= Massive scalar double triangle

This restores the infrared-finite virtual and real thermal benchmark recorded on
`thermal_rate_investigation` at `a0b2d4191cad17d8dc2f20891bb0f245f0df18b7`.
The model, labelled graph, independent quadrature, and historical reference numbers
are unchanged. Historical validation is retained separately from new checkout results.

The interaction is `-lam/2 * (X*H^2 + L*H^2 + X*L^2)`, with `lam=1`,
`MH=0.6`, `ML=0.25`, and zero chemical potentials. The graph multiplier is `1/2`.
This is one labelled two-loop contribution, not the complete order-`lam^4` correction.
The observable is the spectral relaxation kernel `Gamma_rho=-Im(Pi_R)/E`.
An initially empty probe's production source is `n_B(E)*Gamma_rho`, with another
`(2*pi)^-3` for the rate per `d^3p`. No bath-averaging integral is added.

The independent calculation evaluates #link("https://arxiv.org/pdf/hep-ph/0406109")[Blaizot
and Reinosa, Eq. (58)]. The above-threshold vacuum decomposition also follows
#link("https://arxiv.org/pdf/hep-ph/9501201")[Bauberger and Böhm, Eqs. (4)–(12)].
These are new numerical evaluations of published formulas at the specified masses,
not numbers copied from published tables. Neither reference script reads a native result.

== Scope

At `P=(1,0,0,0)`, the `LL` virtual cut is open; the `HH` and all-outgoing `HHL`
cuts are closed. Crossed three-line cuts give a negative interference contribution.
The virtual triangle is below threshold, and the uncut real propagators stay bounded
away from zero. At `P=(1.6,0,0,0)`, `HH` and `HHL` open. The two virtual principal-value
poles cancel in their inclusive sum; the real cuts remain free of uncut on-shell poles.
Neither case requires a thermal threshold counterterm. All bath particles are massive.

UV and threshold subtraction are disabled in this ultraviolet-finite example.
The separate #link("../massive_self_energy/benchmark.typ")[self-energy insertion]
checks the ordinary local and integrated UV machinery. These examples do not provide
thermal infrared subtraction, screening, or resummation.

== Reproduction

From the repository root, using a freshly built CLI:

```sh
./gammaloop -n examples/cli/thermal_rates/massive_double_triangle/run.toml run generate below_threshold inspect export save
./gammaloop -n -s /tmp/gammaloop-thermal-double-triangle-state run below_threshold integrate -D samples=250000 -D workspace=/tmp/thermal-double-triangle-below
./gammaloop -n -s /tmp/gammaloop-thermal-double-triangle-state run above_threshold inspect pole_probe integrate -D samples=250000 -D workspace=/tmp/thermal-double-triangle-above
./gammaloop -n -s /tmp/gammaloop-thermal-double-triangle-state run lower_temperature integrate -D samples=250000 -D workspace=/tmp/thermal-double-triangle-t025
python3 examples/cli/thermal_rates/massive_double_triangle/reference.py --output /tmp/thermal-double-triangle-reference.json
python3 examples/cli/thermal_rates/massive_double_triangle/above_threshold/reference.py --output /tmp/thermal-double-triangle-above-reference.json
```

The `save` block persists the generated processes and command history. Top-level `-n`
disables the CLI’s additional save on exit; it does not disable an explicit `save` block. Subsequent
commands load that state instead of regenerating its seven slots. Use distinct
state/workspace/output destinations for simultaneous runs. `samples`,
`batch`, `cores`, `workspace`, and `output` are card substitutions. The default cap is
250,000 samples per process, with 25,000 samples per iteration and seed `271828`.
The card's historical 0.5% target is a stopping criterion, not a prerequisite for development.
An above-threshold virtual contribution is a small difference of large terms and can
need more statistics. Separate `virtual_hh` and `virtual_ll` processes are for the
finite-distance pole probes; their unsubtracted singular integrals are not separately integrated.
`pole_probe` uses forced Arb precision through `inspect -f`.

== Independent targets and historical evidence

#table(
  columns: 5,
  table.header([Probe mass], [Temperature], [Virtual], [Real], [Inclusive]),
  [1], [0], [`1.0980937552959695e-4`], [0], [`1.0980937552959695e-4`],
  [1], [0.25], [`1.534382715382693e-4`], [`-6.961419335046441e-6`], [`1.4647685220322285e-4`],
  [1], [0.5], [`3.266901870244345e-4`], [`-2.0924992962569012e-5`], [`3.057651940618655e-4`],
  [1.6], [0], [`4.086103384814681e-5`], [`2.423421102578654e-6`], [`4.328445495072546e-5`],
  [1.6], [0.5], [`1.628911955288129e-6`], [`3.653175781627507e-6`], [`5.282087736915636e-6`],
)

The old branch recorded eight million below-threshold evaluations, all eight comparisons
within 1.45 Monte Carlo standard errors, and ten million above-threshold evaluations,
all five comparisons within 0.53 errors. Both runs recorded zero NaNs and zero unresolved
stability failures. These results establish the previous benchmark, not its reproduction
on the current checkout. #link("historical_validation.json")[Historical comparisons and
point records] retain the source revision explicitly. Current native runs must be recorded
separately with executable/revision identity, sample count, errors, and stability diagnostics.

The historical exact endpoint failures at probe masses `0.25`, `0.5`, `0.95`, `1.2`, and
`1.45`, and a boosted near-tangent ray, remain separate geometry regressions. A finite
inclusive physical rate does not make a multiple radial root a regular local coordinate.
The acceptance cases above deliberately avoid those exact endpoints.

== Current pole-cancellation control

With the current CLI, four approaches to the common virtual pole at relative distances
`+/-1e-4` and `+/-1e-8` reproduce the historical joint virtual result within `1.54e-16`
relative. At the closer points, separate terms of magnitude approximately `506` cancel
to `5.53088e-6`. The shared evaluator performs this cancellation before returning its
result. Adding separately rounded binary64 exports loses digits and is retained only
as a diagnostic. All twelve forced-Arb process evaluations are finite and stable; the
physical events conserve four-momentum. #link("validation.json")[The current record]
retains the executable hash, coordinates, joint results, and the full artifact location.

== Current integrated comparison

The current CLI completes the below-threshold 25,000-point controls and a separate
100,000-point above-threshold run. The latter uses seed `314159` and trains on the full
thermal observable; its three slots share correlated cube coordinates.

#table(
  columns: 4,
  table.header([Kinematics / contribution], [Samples], [Native estimate and error], [Independent target]),
  [Below / inclusive], [25,000], [`3.01437563e-4 +/- 6.50100e-6`], [`3.05765194e-4`],
  [Below / virtual], [25,000], [`3.28608432e-4 +/- 2.97661e-6`], [`3.26690187e-4`],
  [Below / real], [25,000], [`-2.71708687e-5 +/- 5.87302e-6`], [`-2.09249930e-5`],
  [Above / inclusive], [100,000], [`4.77230988e-6 +/- 1.37286e-6`], [`5.28208774e-6`],
  [Above / virtual], [100,000], [`1.07829859e-6 +/- 1.24842e-6`], [`1.62891196e-6`],
  [Above / real], [100,000], [`3.69401129e-6 +/- 4.84555e-7`], [`3.65317578e-6`],
)

The below-threshold comparisons, including vacuum and zero-temperature controls, are
within 1.07 reported errors. The three above-threshold comparisons are within 0.45 errors,
and their means obey inclusive = virtual + real to roundoff. All completed runs have zero
final NaNs and zero unresolved stability failures. The above-threshold inclusive estimate
is positive at 3.48 reported errors, but its Monte Carlo error is still 28.8%; the smaller
virtual mean is not resolved relative to its own error. The direct pole control above
separately establishes local cancellation. No 0.5% precision claim is made.

The earlier above-threshold 25,000-point pilot trained on the vacuum slot and gave a
real contribution 2.23 reported errors from its target. It is retained as a statistically
limited pilot; the new-seed run is not combined with it. Exact commands, training slots,
reduced chi-squared, precision-retry statistics, and result hashes are in
#link("validation.json")[the current validation record].
