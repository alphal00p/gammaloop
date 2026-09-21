# Freeze-in scalar benchmark (arXiv:0911.1120, section 6.3)

The requested forward-X construction and momentum scan are documented in
[inclusive_x.md](inclusive_x.md), using [inclusive_x.toml](inclusive_x.toml).
They need no separate bath averaging driver. The current native cut selection
keeps inverse decay but omits the scattering cuts needed for section 6.3.

The earlier two-incoming-particle controls below are retained as diagnostics.
They integrate a generated cross section in both vacuum and
`thermodynamic_equilibrium` modes, then averages each over Maxwell-Boltzmann
incoming bath momenta. The equilibrium result is a diagnostic: its statistical
weight differs from the freeze-in observable in the paper.

## Result

At `MX = T = lambda = 1`, including all three bath-species channels:

| Quantity | Rate density (mass^4) |
|---|---:|
| Paper, equation (6.24) | 4.6098997578035e-5 |
| GammaLoop vacuum cross section + MB bath average | (4.61249 +/- 0.00597)e-5 |
| GammaLoop equilibrium cross section + same MB bath average | (8.26139 +/- 0.01276)e-5 |
| Independent equilibrium-cut reference | 8.2536274083515e-5 |

The vacuum control differs from the paper by 0.0561%, or 0.43 standard errors.
The native equilibrium diagnostic is 79.21% above the paper and agrees with its
own analytic reference within 0.61 standard errors. This is an occupation-factor
mismatch between observables, not evidence of a failed numerical integration.
The existing equilibrium mode is therefore not validated as a freeze-in
production-rate implementation.

The final run used 64 bath nodes and 100,000 Monte Carlo samples per node
(6.4 million samples per process), two CPU workers, and the interpreted summed
function-map evaluator. No unstable samples were recorded. Quoted uncertainties
are Monte Carlo errors; the much smaller separately checked quadrature errors
are saved in `comparison.json`. A few individual nodes fluctuate by up to
3.28 standard errors; the weighted total agrees with the analytic checks.

`comparison.json` stores the final rates and provenance. `bath_node_results.json`
stores every node's kinematics, integration weight, native results, errors, and
analytic references. `quadrature_convergence.json` stores the independent
order-8 through order-96 convergence check. Earlier pilot runs with overlapping
iteration seeds are excluded from these final results.

## Model and reference

The model contains four distinct real scalars, one interaction
`L_int = lambda X B1 B2 B3`, and no other interactions. `X` has mass `MX`;
`B1`, `B2`, and `B3` are massless. All widths and chemical potentials vanish.
The vertex has magnitude `lambda`, hence the tree-level squared amplitude is
`lambda^2`. The overall sign of the vertex does not affect this benchmark.

We compute `B1 B2 > B3 X` and multiply its rate by three for the other two
incoming species pairs. There is no identical-particle symmetry factor. We set
`MX = T = lambda = 1` and report the rate coefficient; a feeble coupling simply
multiplies the answers by `lambda^2`. Rates have mass dimension four.

The paper uses the dilute Maxwell-Boltzmann approximation, omits final-state
Bose enhancement, and sets the produced X occupation to zero. Its equation
(6.24), with the s integral performed as in (6.25), gives

\[
\gamma_{\rm paper}(T)=\frac{3\lambda^2 m_X T^3}{128\pi^5}K_1(m_X/T).
\]

At the chosen point this is `4.6098997578035e-5`. This is the instantaneous
production rate density, not the cosmologically integrated relic yield or
Omega h^2. Source: [Hall et al., section 6.3](https://arxiv.org/pdf/0911.1120).

## What the earlier two-incoming-particle calculation computes

The ordinary cross section is

\[
\sigma_0(s)=\frac{\lambda^2(s-m_X^2)}{16\pi s^2}.
\]

The equilibrium CFF cut adds the weight
`1 + n_B(E_B3) + n_B(E_X)`. In particular, it treats X as an equilibrium boson.
The result is a spectral cut weight, not a freeze-in gain term with `f_X = 0`.
For zero total spatial momentum and `sqrt(s) = 2`, the final-state energies are
`E_B3 = 0.75`, `E_X = 1.25`; the thermal/vacuum ratio is
`2.2968062528953563`. The native paired integrations reproduce that ratio.

The incoming momenta in GammaLoop are fixed external inputs. The native thermal
mode supplies neither their thermal occupation factors nor their integration.
`bath_average.py` performs that additional integration explicitly in the bath
rest frame. It compares the two native outputs using exactly the same incoming
Maxwell-Boltzmann measure. The resulting thermal number is intentionally a
hybrid diagnostic, not a physical freeze-in prediction.

## Direct fixed-kinematics results

With 100,000 shared samples, the native cross sections at `sqrt(s)=2` are:

| Mode | Native integral (mass^-2) | Analytic value |
|---|---:|---:|
| Vacuum | 0.00372322069 +/- 0.00000788061 | 0.00373019398 |
| Thermal equilibrium | 0.00855151656 +/- 0.00001810024 | 0.00856753285 |

Their sample-by-sample ratio is constant, giving the predicted equilibrium
weight to floating-point precision. The full native result is saved as
`fixed_kinematics_result.json`.

## Bath integration and checks

For massless incoming particles and this scalar contact interaction,

\[
\gamma=\frac{3}{32\pi^4}\int_{m_X^2}^\infty ds\,s
 \int_{\sqrt{s}}^\infty dq^0\,
 \sqrt{(q^0)^2-s}\,e^{-q^0/T}\,\sigma(s,q^0;T).
\]

The cross section depends only on the total incoming four-momentum because the
amplitude is constant. The driver realizes that momentum with collinear,
oppositely directed massless inputs of energies `(q0+Q)/2` and `(q0-Q)/2`,
where `Q = sqrt(q0^2-s)`. It does not assume the pair is at rest in the bath.

With `T = MX = 1`, use `sqrt(s)=1+u` and `q0=1+u+v`. The driver uses ordinary
Gauss-Laguerre quadrature in u and generalized Gauss-Laguerre with alpha=1/2
in v. Every quadrature node performs a genuine GammaLoop Monte Carlo
phase-space integral, for both generated processes. Non-overlapping ranges of native iteration seeds allow
quadrature-summed errors to be combined in quadrature. GammaLoop uses
`seed + iteration`, so the driver spaces node seeds by more than the maximum
iteration count; consecutive node seeds would incorrectly correlate the errors. The two columns share
samples within each node and are correlated; their errors must not be treated
as independent when forming a difference or ratio.

The independent equilibrium check analytically averages each Bose distribution
over the outgoing angle. For energies uniform on `[Emin,Emax]`, the average is
`T*log[(1-exp(-Emax/T))/(1-exp(-Emin/T))]/(Emax-Emin)`.
This is used only as a reference; it never replaces a native integration value.

For the analytic cross sections, increasing the quadrature order from 8 to 96
changes the vacuum rate from `4.6100890268499e-5` to `4.6098997578027e-5`,
and the equilibrium diagnostic from `8.2534828177315e-5` to
`8.2536274083515e-5`. The order-8 errors are respectively 0.00411% and 0.00175%,
well below the Monte Carlo uncertainty. The fully converged equilibrium
reference is about 1.79041 times the paper's value.

No UV or threshold counterterms are generated. Threshold subtraction is also
disabled at runtime. There are no virtual propagators or loop corrections to
the scattering amplitude. The ordinary final-state opening at `s=MX^2` is a
phase-space boundary, not a threshold singularity requiring subtraction.

The bath card fixes `--loop-momentum-bases GL0=3`, selecting the massless B3
edge. At zero loop momentum its cut energy is `sqrt(MX^2+Q^2) < q0`, so the
radial origin is inside the cut surface. The automatically chosen X-edge basis
can instead have `MX+Q > q0`; a pilot run spent its time in precision-escalated
radial root solving and was stopped. This was resolved using the existing basis
option, without changing the solver or disabling stability checks.

## Reproduce

Requires a built `./gammaloop` binary and Python with NumPy. Run from the
repository root:

```sh
./gammaloop -s /tmp/gammaloop-freeze-in-fixed-state -n -L off \
  examples/cli/thermal_production/freeze_in_0911_1120/fixed_kinematics.toml \
  run generate integrate

python3 examples/cli/thermal_production/freeze_in_0911_1120/bath_average.py \
  --order 8 --samples 100000 \
  --work-dir /tmp/gammaloop-freeze-in-bath-n8-independent
```

`fixed_kinematics.toml` is the direct native integration card.
`bath_average.toml` supplies generation and one parameterized bath-node block.
`bath_average_executed.toml` preserves the exact card used for the final run,
including all 64 node commands and their non-overlapping seed ranges.
The driver also writes `executed_run_card.toml`, containing the exact node commands,
and runs it with `run generate bath_nodes`. Logs, checkpoints, individual
`integration_result.json` files, and the invocation live in the chosen /tmp
workspace. `--skip-run` recomputes the summary from that workspace using the
same order and sample limit. The cards use `--restart` for the experiment's
integration workspaces and `-n` to avoid saving session state.

The benchmark used the existing dev-optim binary (GammaLoop 0.1.0), built on
2026-09-19, with checkout `e5f4155eb4ff28199578e0ea25acdf24da8f4aed`.
No Rust implementation or existing tests were changed.
