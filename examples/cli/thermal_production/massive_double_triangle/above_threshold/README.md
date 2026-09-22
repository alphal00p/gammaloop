# Massive double triangle above threshold

The same scalar model and labelled graph as the parent example are evaluated
with `MX=1.6`, `MH=0.6`, `ML=0.25`, `T=0.5`, unit cubic couplings, zero chemical
potentials, and external momentum `(1.6,0,0,0)`. This opens both the HH threshold
at `1.2` and the HHL threshold at `1.45`; LL was already open. The saved run card
reuses the parent model and graph, and sets `MX=1.6` explicitly.

The inclusive LU sum stays finite and agrees with the independent reference
**without an additional threshold subtraction or any implementation changes**.
Individual HH and LL virtual cuts have principal-value poles. Their common
radial pole cancels in the sum, while all real cuts remain free of uncut
propagator poles because the HHL vertices have stable masses.

## Literature comparison

The reference independently evaluates [Blaizot and Reinosa, Eq. (58)](https://arxiv.org/pdf/hep-ph/0406109)
with explicit principal values for both virtual triangles. The vacuum result
also follows the arbitrary-mass four-cut decomposition of
[Bauberger and Böhm, Eqs. (4)–(12)](https://arxiv.org/pdf/hep-ph/9501201).
These are evaluations of published formulas at our chosen masses, not numbers
from a published table. No GammaLoop expression enters the reference.
[reference_derivation.md](reference_derivation.md) gives the prescriptions,
formulas, and quadrature convergence.

Rates below are the isolated graph's `Gamma_rho = -Im(Pi_R)/MX`, with graph
factor `1/2`. They are not the complete scalar theory's order-lam^4 correction.
Native errors are Monte Carlo errors. Each slot has two million samples.

| Contribution | GammaLoop | Reference | Difference / MC error |
|---|---:|---:|---:|
| vacuum | 4.327448326e-05 ± 1.158e-07 | 4.32844549507e-05 | -0.086 |
| thermal | 5.218447424e-06 ± 2.947e-07 | 5.28208773692e-06 | -0.216 |
| virtual | 1.505229551e-06 ± 2.68e-07 | 1.62891195529e-06 | -0.462 |
| real | 3.713217874e-06 ± 1.139e-07 | 3.65317578163e-06 | +0.527 |
| zero_temperature | 4.327448326e-05 ± 1.158e-07 | 4.32844549507e-05 | -0.086 |

All comparisons are within 0.53 MC standard errors. The reference is stable to
`2e-19` between quadrature orders 128 and 512. Vacuum and zero-temperature
equilibrium agree to rounding; the inclusive result equals virtual plus real
to rounding with correlated sampling. Exact comparisons and production-source
conversions are saved in `comparison.json`; the raw native output is
`native_result.json`. `pilot_result.json` retains the initial 100,000-sample run.

The independent thermal virtual cuts are `+7.98677900462e-5` (LL) and
`-7.82388780909e-5` (HH), leaving `1.62891195529e-6`: approximately 98% of the
positive virtual term cancels. The aggregate real contribution is
`+3.65317578163e-6`. Thus the aggregate real and virtual terms are both positive
at this point; cancellations occur between the virtual cuts and between signed
real channels. This differs from the parent's below-threshold aggregate signs.

The production source for an initially empty X population is
`S = n_B(MX) * Gamma_rho`, with a further `(2*pi)^-3` for the rate per `d^3p`.
The reference source is `2.24458980512e-7`; its native value and error are in
`comparison.json`. No bath-averaging script is used.

## Local threshold cancellation

The pole-probe blocks fix the light loop momentum to
`(sqrt(MX^2/4-ML^2),0,0)` and set the heavy loop momentum to
`(0,(1+delta)*sqrt(MX^2/4-MH^2),0)`. Both two-body radial roots coincide at
`delta=0`. Approaching from both signs of delta gives:

| delta | HH virtual cut, f128 | LL virtual cut, f128 | Joint virtual sum, f128 |
|---:|---:|---:|---:|
| -1e-4 | +0.0506497262849 | -0.0506441957258 | 5.53055914752e-6 |
| +1e-4 | -0.0506206790602 | +0.0506262102674 | 5.53120715667e-6 |
| -1e-8 | +506.352038625 | -506.352033094 | 5.53088359849e-6 |
| +1e-8 | -506.352009684 | +506.352015214 | 5.53088366329e-6 |

The complete approach, from `|delta|=1e-2` through `1e-8`, is saved in
`pole_cancellation.json`. The automatic joint evaluation escalates from double
to quadruple precision at `|delta|=1e-6` and `1e-8`, reproducing the forced-f128
result. Sum the cuts inside the native joint evaluation: subtracting separately
rounded f64 cut values can lose the finite remainder. The exact isolated pole
at delta zero is not evaluated.

The exported group counts partition as `thermal 20 = virtual 6 + real 14`, and
`virtual 6 = HH 3 + LL 3`. At the ordinary inspection point, ten physical events
are retained: one LL, one HH, and four on each HHL cut support. All conserve
four-momentum to `2.62e-13`; inclusive minus virtual minus real is `8.47e-22`.
These checks are saved in `point_check.json`.

There were zero NaNs and zero unresolved stability failures across all ten
million integration evaluations. Quadruple precision was used for 0.10958% of
evaluations; arbitrary precision was never needed. This result demonstrates
inclusive cancellation for this graph. It does not establish general
mixed-sign thermal threshold subtraction, which is not exercised here.

## Reproduce

From the repository root with the existing CLI build:

```sh
./gammaloop --read-only-state -s /tmp/gammaloop-massive-double-triangle-above-state -n -L off \
  examples/cli/thermal_production/massive_double_triangle/above_threshold/run.toml \
  run generate export inspect pole_probe pole_probe_quad integrate

python3 examples/cli/thermal_production/massive_double_triangle/above_threshold/reference.py
python3 examples/cli/thermal_production/massive_double_triangle/above_threshold/reference.py \
  --mass 1 --output /tmp/double_triangle_below_check.json
```

The card uses the existing `SummedFunctionMap` evaluator, four integration
workers, two million samples per slot, seed `271828`, and disabled threshold/UV
subtraction. It writes generated diagnostics and the restarted integration
workspace under `/tmp/gammaloop-massive-double-triangle-above-*`. The isolated
HH/LL processes are used only for the pole probes; their unsubtracted integrals
are not separately Monte Carlo integrated.

The NumPy reference uses analytic angular integrals and an explicit symmetric
principal-value pairing before numerical quadrature. Its `--mass 1` run
reproduces the parent below-threshold result, `0.00030576519406186564`, and its
vacuum triangles are independently checked with Feynman parameters. No core
code or existing tests were changed for this experiment.
