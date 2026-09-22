# Finite rates that the current radial solver cannot evaluate

At source revision `21896a0b99f08c94caeb84290aebb4c92919c6e1`, the same massive
scalar double triangle fails at exact external-energy branch points. This is a
numerical limitation: the independent inclusive rate is finite. All bath masses
remain positive and stable, with `MH=0.6`, `ML=0.25`, `T=0.5`, unit cubic
couplings, zero chemical potentials, and external momentum `(M,0,0,0)`.
The graph factor is `1/2`; this is one labelled graph, not the complete model's
order-lam^4 correction. No implementation or existing tests were changed.

## Exact thresholds and controls

The table uses the parent's normalization `Gamma_rho = -Im(Pi_R)/M`.
For an initially empty X population, multiply by `n_B(M)` for the production
source, and by a further `(2*pi)^-3` for the rate per `d^3p`.

| M | Stationary cut energy | Independent Gamma_rho at T=0.5 | Native equilibrium inspection |
|---:|---|---:|---|
| 0.25 | H - H + L | 6.47981779354e-4 | Error: non-simple radial root at zero |
| 0.50 | L + L | 2.27286537633e-4 | Same error |
| 0.95 | H + H - L | 3.16520790648e-4 | Same error |
| 1.20 | H + H | 4.75083945917e-4 | Same error |
| 1.45 | H + H + L | 5.93406524897e-5 | Same error |

All ten neighboring inspection points `M_threshold +/- 0.001` complete.
The fixed LMB point is `(0.4,0.2,0.1,-0.3,0.5,0.2)`. Inspection values in
`native_results.json` are pointwise integrands, not integrated rates.
Integration also exits with an error at `M=0.95` and `M=1.2`, before completing
its first iteration. The runs do not quietly return a zero rate.

Neighboring Monte Carlo integrations also agree with the independent rates:

| M | Native Gamma_rho | Reference | Difference / MC error |
|---:|---:|---:|---:|
| 0.949 | 3.168765953e-4 +/- 1.814e-6 | 3.16864296031e-4 | +0.007 |
| 0.951 | 3.170659817e-4 +/- 1.471e-6 | 3.16187510417e-4 | +0.597 |
| 1.199 | 4.435575255e-4 +/- 2.150e-6 | 4.45072410144e-4 | -0.705 |
| 1.201 | 4.628332516e-4 +/- 2.055e-6 | 4.63218591297e-4 | -0.187 |

The sample counts are respectively 200000, 200000, 125000 and 150000 (the
last two reach the card's 0.5% stopping target). All 675000 evaluations have
zero NaNs and zero unresolved stability failures. Raw results and exact
threshold integration errors are saved in `integration_results.json`.

`M=0.95` is the most useful thermal-specific reproducer. It is the stationary
point of `X + L -> H + H`, with all three bath momenta zero. It is below the
vacuum HH and HHL thresholds, and the original stable masses are unchanged.

| Control at M=0.95 | Result |
|---|---|
| Full thermal cut sum, ordinary inspection | Non-simple radial root at zero |
| Full thermal cut sum, `inspect -f` | Cannot isolate simple radial roots near zero |
| Vacuum mode | Finite inspection result |
| Zero-temperature equilibrium mode | Non-simple radial root at zero |
| Thermal HH + LL supports only | Finite inspection result |

The zero-temperature equilibrium failure is especially diagnostic: absorption
weights eventually vanish, but the implementation attempts to enumerate their
roots before reaching those weights. The vacuum reference is finite,
`Gamma_vacuum=1.08018533191e-4`. At the ordinary HH threshold, `M=1.2`, vacuum
also fails, while removing only the HH support makes inspection succeed.
At `M=1.45`, retaining only the two virtual supports succeeds. These support
restrictions isolate the failing root families; they are not general recipes
for the inclusive rate, since removing thermal supports also removes genuine
positive-radius contributions.

`inspect -f` selects `Precision::Arb` in the current implementation, despite
an outdated f128 help comment. Both `M=1.2` and `M=1.45` still report the exact
non-simple endpoint root in this mode. At `M=0.95`, arbitrary precision reaches
an unresolved interval at compact radius approximately `1.98e-301` instead;
it still fails to produce a result.

## First failing boundary

Native and reference routes share the model, graph normalization, external
kinematics, equilibrium Bose weights and complete allowed cut supports. Their
integration representations differ:

- Native: graph/CFF generation, thermal LU selection, signed radial root
  enumeration, residues and evaluator, cut sum, Monte Carlo integration.
- Reference: published spectral cut representation, analytic angular and
  principal-value integration, endpoint quadrature, inclusive sum.

Generation succeeds. Standalone exports at both `M=0.95` and `M=1.2` contain
32 orientations, 20 cut groups, 44 parameter slots and no counterterms, as
recorded in `native_results.json`. Changing M changes the radial energy surface; the first
observed failure is root enumeration, before the native evaluation can form
the inclusive sum. It does not implicate the hyperbolic evaluator, UV
subtraction, or the previously demonstrated cancellation between virtual cuts.
For a rest-frame ray, the uncompactified cut equation has the form

```
F(r) = sum_i sigma_i sqrt(m_i^2 + r^2 |q_i|^2) - M,
F(0) = sum_i sigma_i m_i - M,
F'(0) = 0.
```

At the listed masses the surface therefore has a multiple endpoint root.
For all-positive signs this is a normal threshold; mixed signs give the
stationary thermal absorption surfaces. Their nonanalytic external-energy
dependence is compatible with a finite rate. For example, the independent
rate has a cusp at `M=1.2`; it does not diverge there.

The ownership chain is:

1. `crates/gammalooprs/src/integrands/process/cross_section/mod.rs:1577`
   calls `positive_radius_roots` and propagates its error.
2. `crates/gammalooprs/src/cff/esurface.rs:800` enumerates roots on the closed
   compact interval `[0,1]`; the strict positive-radius filter runs only after
   enumeration returns successfully.
3. `crates/gammalooprs/src/utils/newton_solver.rs:723` rejects a root whose
   derivative is zero. The unresolved-interval error is at line 773.
4. `evaluate_from_source` also propagates a root-finding error out of its
   current stability level. An ordinary root error therefore does not trigger
   the next precision level automatically.

At a normal opening threshold the new channel has zero phase space, while
other channels still contribute. For mixed signs, other positive radial roots
may remain. Handling the endpoint must preserve these contributions. An
external-energy branch point is not removed by summing inclusive cuts, and
cancelling internal virtual poles does not supply this missing endpoint rule.

## Literature reference

`reference.py` reuses the independent evaluation of
[Blaizot and Reinosa, Eq. (58)](https://arxiv.org/pdf/hep-ph/0406109)
from `../above_threshold/reference.py`. The values are our numerical evaluations
of the published integral representation at these masses, not published
numbers or GammaLoop-generated expressions. `reference.json` saves every
quadrature order through 1024 at all five thresholds and their neighbors.

The wrapper changes endpoint coordinates and cancellation-sensitive algebra,
allows positive M below the old `M>2L` input guard, and gives the exact HHL
opening zero phase space. In particular, at `M=2H` it factors the vanishing
vacuum-triangle discriminant before quadrature. The unmodified reference's
slow convergence at that endpoint is a separate quadrature limitation.
The wrapper checks every source replacement and stops if the shared reference
has changed; it is an experimental extension of the existing oracle.

For the general origin of finite thermal branch points from stationary
emission/absorption energies, see
[Weldon, section III and Appendix A](https://arxiv.org/pdf/hep-ph/0203057).
His equal-mass examples provide the general mechanism; the unequal-mass
stationary energies above follow directly from our cut equation.

## Reproduce

From the repository root, the saved card defaults to the thermal `M=0.95`
failure:

```sh
./gammaloop --read-only-state -s /tmp/gammaloop-thermal-threshold-state -n -L off \
  examples/cli/thermal_production/massive_double_triangle/failure_modes/run.toml \
  run generate inspect
```

Replace the final line to select a control or another failure:

```sh
run generate inspect -D mass=0.949
run generate inspect -D mass=0.951
run generate inspect -D mode=vacuum
run generate inspect -D mode=zero_temperature_equilibrium
run generate inspect -D precision=-f
run generate inspect -D 'cuts={ H H, L L }'
run generate inspect -D mass=1.2
run generate inspect -D mass=1.45
run generate integrate
```

Use `-D output=/tmp/unique-prefix` to separate outputs. The integration block
uses 1000 samples by default; `-D n_start=25000 -D n_max=200000` selects the
neighboring-rate validation settings. Add `export` between `generate` and the
failing block for a standalone structural export. Generation and runtime
threshold subtraction, as well as UV subtraction, are disabled as in the
validated parent benchmark. All generated workspaces are under `/tmp`.

The reference requires NumPy:

```sh
python3 examples/cli/thermal_production/massive_double_triangle/failure_modes/reference.py
```

A separate [boosted inspection example](boosted_tangency.md) demonstrates a
precision-sensitive tangent radial surface away from these external branch
points. Its saved inputs fail in both thermal and vacuum modes at ordinary
precision and succeed with arbitrary precision. Bounded generic boosted and
unstable-mass pilots found no additional demonstrated failure; their large
Monte Carlo errors do not establish accurate rates.
