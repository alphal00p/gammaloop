# Massive scalar self-energy insertion: raised-cut investigation

The original investigation at commit `37f2497da` found that the exact
repeated-line graph evaluated in vacuum and zero-temperature equilibrium, but
finite-temperature generation failed while constructing its raised-cut
evaluator. The original evidence is retained below. A subsequent temporary
workaround is described under [Hyperbolic dual workaround](#hyperbolic-dual-workaround).

The subsequent [UV-renormalized integration](uv_integration.md) agrees with
an independent literature-based spectral reference, with both UV subtractions
enabled and medium vacuum subtraction explicitly false.

## Inputs and scope

The model is the existing [massive scalar model](../massive_double_triangle/model.json)
with cubic XHH and HHL interactions, coupling 1, and zero chemical potentials.
We use `MX=1.6`, `MH=0.6`, `ML=0.25`, `T=0.5`, and external momentum
`(1.6,0,0,0)`. The [graph](self_energy.dot) is an X self-energy bubble with an
HL bubble inserted on one H line. Its propagator structure is

```text
D_H(Q)^2 D_H(P-Q) D_H(K) D_L(Q-K),   D_m(Q) = 1/(Q^2-m^2+i0).
```

Edges 0 and 3 are the serial copies carrying Q. The graph has overall factor
1: this is a labelled graph diagnostic, not the complete order-coupling^4 rate
of the model. The selected cut supports are HH and HHL, including equilibrium
thermal crossings.

The repeated H shell is outside the inserted bubble's spectral support:
`(MH-ML)^2 = 0.1225 < MH^2 = 0.36 < (MH+ML)^2 = 0.7225`.
Thus the bubble has neither a decay cut nor a thermal Landau cut on that shell.
The external energy also avoids the HH and HHL endpoints, 1.2 and 1.45.

The original investigation and workaround results below are **bare pointwise integrands**, at momentum-space input
`[0.4,0.2,0.1,-0.3,0.5,0.2]`. UV and threshold subtraction are disabled to
isolate raised-cut handling. Unlike the double triangle, this graph has a
logarithmically UV-divergent bubble subgraph. An integrated rate requires its
ordinary vacuum mass counterterm. Subtracting the vacuum result from the
thermal result does not by itself remove this subdivergence: the divergent
vacuum subgraph also multiplies the thermal outer loop.

## Original native results

| Mode / cut selection | Generation | Real integrand |
|---|---|---:|
| Vacuum, HH + HHL | succeeds | 1.7841968873755623e-5 |
| Zero-temperature equilibrium, HH + HHL | succeeds | 1.7841968873755620e-5 |
| T=0.5 equilibrium, HH + HHL | panics | unavailable |
| T=0.5 equilibrium, HH only | same panic | unavailable |
| T=0.5 equilibrium, HHL only | succeeds | 1.7586444075576102e-7 |
| Vacuum, edge 3 mass changed to 0.601 | succeeds | 1.7924492278743268e-5 |
| Zero-temperature equilibrium, edge 3 mass 0.601 | succeeds | 1.7924492278743268e-5 |
| T=0.5 equilibrium, edge 3 mass 0.601 | succeeds | 2.1070414566177320e-5 |

These entries use double precision. Exact vacuum and zero-temperature raised
cuts also agree in arbitrary precision, with returned value
`1.7841968873755768e-5` in both modes. In this CLI, `inspect -f` selects `Arb`;
the saved evaluation metadata confirms the precision actually used.

The finite-temperature failure is:

```text
called `Result::unwrap()` on an `Err` value:
"To vectorize the function 'symbolica::coth', the symbol 'symbolica::coth_v'
must be defined that takes the vector index as an additional tag"
Location: crates/gammalooprs/src/integrands/process/evaluators.rs:1447
```

The vacuum export has 14 orientations and two cut groups, with LU residue
orders `[[1,2],[1]]`. Zero-temperature equilibrium has 32 orientations and
12 groups: one `[1,2]` group and eleven `[1]` groups. Both contain the same
order-two derivative shape `[[0],[1]]`. The successful finite-temperature
HHL-only and split-mass controls contain only simple LU cut residues.

## Original first failing boundary

The routes share the graph, masses, external momentum, evaluator settings,
and momentum-space sample. Their stages are:

1. Import the same graph and model.
2. Generate CFF orientations and medium-dependent distribution functions.
   This is the first stage with different inputs to subsequent shared code.
3. Group cut surfaces and construct Laurent coefficients for each residue
   order. Both vacuum and equilibrium routes can request order two.
4. Build a hyperdual evaluator for the raised residue. Finite temperature
   fails here when vectorizing `coth`.
5. Search radial roots, evaluate residues, and integrate. The failing route
   never reaches these runtime stages.

The existing owners are:

- `crates/gammalooprs/src/processes/cross_section.rs`,
  `build_derivative_structure_atom`: shared higher-order Laurent construction.
- `crates/gammalooprs/src/integrands/process/cross_section/mod.rs:772`:
  derive the dual shape from the cut index and build the evaluators.
- `crates/gammalooprs/src/integrands/process/evaluators.rs:1444`:
  `tree.vectorize(&dualizer).unwrap()`, the observed fatal boundary.
- `crates/gammalooprs/src/graph/feynman_graph.rs:268`:
  thermal distribution bodies and their explicit energy derivatives.
- `crates/three-dimensional-reps/src/cff_recursion.rs:277`:
  repeated-chain distribution derivatives.

These paths are relative to the repository root. The existing comment in
`param_builder.rs:2005` already distinguishes native hyperbolic symbolic
differentiation from unsupported evaluator hyperdual vectorization.

At finite temperature the Bose weights contain
`(sign+coth(beta*E/2))/2`; their first energy derivative is
`-beta*csch(beta*E/2)^2/4`. These explicit distribution derivatives are already
registered, and derivative symbols occur in the generated zero-temperature
CFF structure, where their values vanish. They do not replace the additional
radial derivative required by a raised LU cut. That derivative must also pass
through the energy dependence of the thermal functions. Vacuum and
zero-temperature equilibrium at zero chemical potential reduce the weights
to constants, so they avoid this unsupported operation.

The controls exclude an exact external threshold, failure of the radial root
solver, and ordinary evaluation of thermal weights as explanations of this
panic. They identify a concrete evaluator-construction limitation; they do
**not** certify the full finite-temperature raised-cut algebra. No existing
configuration-only bypass was found. Turning compilation off does not skip
this shared evaluator-construction step.

## Split-mass consistency check

Only the mass of edge 3 is varied; the bubble and spectator masses remain fixed.
For squared masses a and b,

```text
D_a(Q) D_b(Q) = [D_a(Q)-D_b(Q)]/(a-b) -> dD_a(Q)/da = D_a(Q)^2.
```

This identity holds before the thermal sums. Consequently the simple-cut
split graph supplies a native consistency check on the coincident-mass limit.
Changing the global MH parameter would also change other propagators and
would not isolate this limit.

[scan.py](scan.py) evaluates the split graph for edge-3 masses `0.6 +/- step`
and averages the two results. All 20 evaluations succeed using arbitrary
precision. The symmetric means are:

| Mass step | Vacuum | T=0.5 equilibrium |
|---:|---:|---:|
| 1e-2 | 1.785206562362233e-5 | 2.0949512144908932e-5 |
| 1e-3 | 1.784206989505487e-5 | 2.0918414625421167e-5 |
| 1e-4 | 1.7841969921013616e-5 | 2.0918103661502923e-5 |
| 1e-5 | 1.7841964244575660e-5 | 2.0918090889135240e-5 |
| 1e-6 | 1.7842309134594377e-5 | 2.0918803472582180e-5 |

At step `1e-4` the vacuum mean agrees with the exact raised result within
`5.9e-8` relative. The thermal sequence supports a finite limit near
`2.09181e-5` at this point. It is not an independent literature rate or a
high-precision determination of the exact thermal raised result.

Convergence deteriorates at the smallest steps despite the Arb rotation
stability check reporting a stable result. This remains a separate numerical
issue to investigate. A candidate input boundary is that literal edge masses
are converted to f64 for root kinematics (`graph/edge.rs:226`), while evaluator
energies retain the symbolic mass atom (`graph/feynman_graph.rs:235`). A small
mismatch can be amplified near coincident cuts, but its responsibility for
this scan has not been established.

## Thermal derivative reference

Raised cuts require derivatives of thermal occupation factors as well as the
vacuum kinematics. A simple independent identity illustrates the extra term:

```text
T sum_n 1/(omega_n^2+E^2)^2
  = (1+2*n_B(E))/(4*E^3) + n_B(E)*(1+n_B(E))/(2*T*E^2).
```

It follows by differentiating the simple bosonic Matsubara sum with respect
to E squared. Merely multiplying the vacuum raised result by a Bose weight
misses the second term. Thermal mass-derivative cutting relations are discussed
in [Blazek and Matak, section 3, Eqs. (17)-(21)](https://link.springer.com/article/10.1140/epjc/s10052-022-10165-8).
The vacuum LU treatment of raised propagators is described in
[Capatti, Hirschi and Ruijl, section 2](https://arxiv.org/pdf/2203.11038).
These support the required derivative structure; their numerical observables
are not this scalar graph. The double-triangle formula used by the previous
benchmark is not a reference formula for this different topology.

## Hyperbolic dual workaround

Until upgrading the pinned Symbolica to version 3, dual evaluators locally
register exponential definitions through the existing `FunctionMap`. For
`q=exp(-2*x)`, these are `coth(x)=(1+q)/(1-q)`,
`csch(x)=2*exp(-x)/(1-q)`, `tanh(x)=(1-q)/(1+q)`, and
`sech(x)=2*exp(-x)/(1+q)`. Symbolica can vectorize these arithmetic and
exponential operations. Function-map resolution also reaches native
hyperbolics nested inside shared thermal function bodies, so the complete LU
residue retains its thermal derivatives.

The override is confined to `GenericEvaluator::new_from_raw_params` when a
dual shape is requested. Ordinary evaluators retain their native hyperbolic
functions. No dependency version, public setting, root solver, or thermal
cutting prescription is changed. The four temporary registrations can be
removed after upgrading Symbolica.

This representation is appropriate for the positive real thermal arguments
in this massive, zero-chemical-potential example. Decaying exponentials avoid
UV overflow on that domain. Large negative real arguments can overflow, and
near-zero Bose arguments lose precision through `1-exp(-2*x)`. It is not a
replacement for robust native hyperbolic vectorization across all regimes.
The added regression checks values and first hyperdual derivatives against
the independent MPC oracle, including a moderate complex argument and positive
UV arguments up to 1000. Existing native-hyperbolic tests are unchanged.

With the workaround, all six native checks generate and evaluate successfully:

| Exact graph / cut selection | Precision | Real integrand |
|---|---|---:|
| T=0.5 equilibrium, inclusive | Double | 2.0918100443146220e-5 |
| T=0.5 equilibrium, inclusive | Arb | 2.0918100443146576e-5 |
| Vacuum, inclusive | Arb | 1.7841968873755768e-5 |
| Zero-temperature equilibrium, inclusive | Arb | 1.7841968873755768e-5 |
| T=0.5 equilibrium, HH | Arb | 2.0742236002390790e-5 |
| T=0.5 equilibrium, HHL | Double | 1.7586444075576102e-7 |

The exact thermal result agrees with the earlier symmetric split-mass mean at
step `1e-4` within `1.54e-7` relative. Double and Arb agree within `1.71e-14`
relative. HH plus HHL closes onto the inclusive result within `3e-20` absolute.
Vacuum and zero-temperature Arb values, and the ordinary thermal HHL value,
are unchanged from the original investigation. These remain pointwise checks;
no UV-renormalized rate is claimed.

[workaround_results.json](workaround_results.json) retains the commands and
full native evaluations; [workaround_comparison.json](workaround_comparison.json)
records the comparison arithmetic. Formatting, cargo check, the CLI build,
and clippy passed. Both the new MPC derivative regression and the existing
parameterless-function hyperdual regression passed:

```sh
cargo nextest run -p gammalooprs --lib --cargo-profile dev-optim \
  -P test_gammaloop \
  -E 'test(=integrands::process::evaluators::hyperbolic_tests::hyperbolic_hyperdual_thermal_workaround) or test(=integrands::process::evaluators::tests::parameterless_function_body_preserves_multiparameter_hyperdual_derivatives)'
```

## Reproduce and retained evidence

Run from the repository root after rebuilding the CLI. This command reproduces
the original generation failure before the workaround and evaluates the exact
thermal raised graph after it:

```sh
./gammaloop --read-only-state -s /tmp/self-energy-thermal-state -n -L off \
  examples/cli/thermal_production/massive_self_energy/run.toml \
  run generate export inspect -D output=/tmp/self-energy-thermal
```

Use distinct state/output paths for each control:

```sh
./gammaloop --read-only-state -s /tmp/self-energy-vacuum-state -n -L off \
  examples/cli/thermal_production/massive_self_energy/run.toml \
  run generate export inspect -D mode=vacuum -D output=/tmp/self-energy-vacuum

./gammaloop --read-only-state -s /tmp/self-energy-zero-state -n -L off \
  examples/cli/thermal_production/massive_self_energy/run.toml \
  run generate export inspect -D mode=zero_temperature_equilibrium \
  -D output=/tmp/self-energy-zero

./gammaloop --read-only-state -s /tmp/self-energy-real-state -n -L off \
  examples/cli/thermal_production/massive_self_energy/run.toml \
  run generate export inspect -D 'cuts=H H L' -D output=/tmp/self-energy-real

./gammaloop --read-only-state -s /tmp/self-energy-split-state -n -L off \
  examples/cli/thermal_production/massive_self_energy/split_mass.toml \
  run generate export inspect -D output=/tmp/self-energy-split

python3 examples/cli/thermal_production/massive_self_energy/scan.py
```

To isolate the HH support that originally failed, use `-D 'cuts=H H'` with the exact card.
To use arbitrary precision there, add `-D precision=-f`. The split card uses
arbitrary precision by default and fixes only edge 3 to mass 0.601. The scan
writes each literal-mass graph under its `--work` directory. It overwrites its
`--output` file, which defaults to the adjacent `mass_limit.json`.

- [initial_results.json](initial_results.json): exact modes, full native point
  results, and failure log.
- [controls.json](controls.json): cut-support and one-sided mass-split controls.
- [precision_controls.json](precision_controls.json): exact vacuum and
  zero-temperature Arb evaluations.
- [structure.json](structure.json): cut-index/dual-shape inventories from
  successful standalone exports; full export paths refer to temporary files.
- [mass_limit.json](mass_limit.json): all symmetric split-mass evaluations,
  precision/stability metadata, and means.

The original investigation changed no Rust code or existing tests and stopped
at the first fatal finite-temperature boundary. The subsequent workaround
changes evaluator construction and adds a separate derivative regression; it
does not change existing test expectations. The later renormalized integral is
reported separately in [uv_integration.md](uv_integration.md).
