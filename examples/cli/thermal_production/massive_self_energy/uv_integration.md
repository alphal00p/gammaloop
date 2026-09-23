# UV-renormalized integration of the massive self-energy insertion

The exact raised-line diagram now integrates successfully with both local and
integrated UV subtraction. The full equilibrium result, including its vacuum
part, agrees with an independent spectral calculation. No further core changes
were needed beyond the previously recorded hyperbolic dual workaround.

## Setup and normalization

The [run card](uv_run.toml) uses the existing [graph](self_energy.dot) and scalar
model: `M=1.6`, `H=0.6`, `L=0.25`, unit cubic coupling, graph factor 1,
`P=(1.6,0,0,0)`, and `T=0.5`. The observable is the graph contribution to
`Gamma_rho=-Im(Pi_R)/M`. It is not the complete order-coupling^4 rate of the model.

The card explicitly sets

```toml
[cli_settings.global.generation.medium]
vacuum_subtraction = false

[cli_settings.global.generation.uv]
subtract_uv = true
generate_integrated = true

[default_runtime_settings.general]
mu_r = 1.6
m_uv = 1.0
renormalization_localization_scale = 1.0
```

The default logarithmic UV prescription, MUV with Vakint MSbar normalization,
combines the local subtraction and finite integrated addback into MSbar
renormalization. Only the inserted HL bubble is superficially UV divergent;
the full graph has degree -2. Its subtraction leaves the renormalized bubble

```text
B_MSbar(s) = -1/(16*pi^2) int_0^1 dx
  log([x*H^2+(1-x)*L^2-x*(1-x)*s-i0]/mu_r^2).
```

The Taylor bubble's integrated finite coefficient is proportional to
`log(mu_r^2/m_uv^2)`, which is nonzero here. The standalone export contains
this scale logarithm and the integrated localization kernel; compact evidence
is retained in [uv_structure.json](uv_structure.json). The physical scale is
mu_r. The auxiliary m_uv and localization scale must cancel after integration.
Threshold subtraction remains disabled. These masses keep the repeated H
shell outside the inserted bubble's spectral support, and M is away from
external threshold endpoints.

## Results

The reference uses the vacuum mass-derivative and dispersion relations in
[Bauberger et al., hep-ph/9409388, Eqs. (3), (13)-(15)](https://arxiv.org/pdf/hep-ph/9409388)
and thermal bubble ingredients in
[Nishikawa, Morimatsu and Hidaka, hep-ph/0302098, section III](https://arxiv.org/pdf/hep-ph/0302098).
The final thermal insertion reduction is derived in
[reference_derivation.md](reference_derivation.md). **These are our independent
numerical evaluations of the formulas, not numbers tabulated in the papers.**
No native results enter either reference script or its normalization.

| Mode / contribution | Native samples | GammaLoop +/- MC error | Independent reference | Difference / MC error |
|---|---:|---:|---:|---:|
| thermal | 2,000,000 | 1.4041890505e-04 +/- 7.137e-07 | 1.40113504288e-04 | +0.43 |
| vacuum | 1,000,000 | 6.9621305384e-05 +/- 4.515e-07 | 6.91115327439e-05 | +1.13 |
| zero_temperature | 300,000 | 6.9459787485e-05 +/- 9.746e-07 | 6.91115327439e-05 | +0.36 |
| thermal_aux | 2,000,000 | 1.4045877446e-04 +/- 7.382e-07 | 1.40113504288e-04 | +0.47 |
| thermal_hh | 2,000,000 | -8.8119106070e-05 +/- 6.769e-07 | -8.82587226166e-05 | +0.21 |
| thermal_hhl | 200,000 | 2.2725037156e-04 +/- 1.099e-06 | 2.28372226905e-04 | -1.02 |

`thermal_aux` changes `m_uv` from 1 to 1.4 and the localization scale from 1
to 0.7, keeping mu_r=1.6, all physical parameters, and medium vacuum subtraction
unchanged. Its central value differs from the baseline by only `3.99e-8`,
much less than either run's MC error. Vacuum and zero-temperature equilibrium
also agree with the same reference.

The negative thermal HH raised-cut term cancels 38.65% of the positive HHL
continuum in the reference. Both pieces separately agree with the independent
calculation. The native inclusive-minus-HH-minus-HHL difference is `1.29e-6`,
comparable to the MC errors of these separately adapted integrations; no
sample-by-sample closure is implied by their integrated estimates.

All six runs, totaling 7.5 million evaluations, have zero NaNs and zero
unresolved stability failures. The baseline thermal run uses double precision
for 99.8235% of samples and quadruple precision for 0.1765%, with no arbitrary
precision escalation. HHL stops at 200,000 samples after reaching the requested
0.5% relative MC error; the other runs reach their specified sample caps.

Gauss orders 64-512 agree on the reference at approximately 14 significant
figures. Complex mass-derivative steps 1e-8 and 1e-20 agree. A separate vacuum
dispersion calculation and an analytic thermal derivative independently check
the normalization and the raised HH contribution.

The production source for an initially empty X mode is
`S=n_B(M)*Gamma_rho`. At T=0.5, the native result is
`(5.96701264e-06 +/- 3.033e-08)`, compared with
`5.954034823894e-06`. The rate per `d^3p` is S divided by
`(2*pi)^3`; both conversions are saved in [uv_comparison.json](uv_comparison.json).

## Reproduce

From the repository root, with the CLI containing the hyperbolic workaround:

```sh
./gammaloop --dev-optim --read-only-state -s /tmp/self-energy-uv-state -n -L off \
  examples/cli/thermal_production/massive_self_energy/uv_run.toml \
  run generate export inspect integrate \
  -D output=/tmp/self-energy-uv -D samples=2000000

python3 examples/cli/thermal_production/massive_self_energy/vacuum_reference.py
python3 examples/cli/thermal_production/massive_self_energy/thermal_reference.py
```

Use distinct state/output prefixes for each control. Vacuum uses
`-D mode=vacuum -D samples=1000000`; zero-temperature equilibrium uses
`-D mode=zero_temperature_equilibrium -D samples=300000`. The auxiliary-scale
control adds `-D m_uv=1.4 -D localization=0.7 -D samples=2000000`.
HH and HHL use `-D 'cuts=H H' -D samples=2000000` and
`-D 'cuts=H H L' -D samples=1000000`, respectively. The exact commands for
these controls are retained in their result JSON. The reference scripts need
NumPy and overwrite their adjacent JSON output by default.

The native integrations use two workers, eager SummedFunctionMap evaluators,
summed orientations, 50,000 samples per iteration, seed 271828, and a target
relative error of 0.005. Workspaces and full standalone exports are under
`/tmp`; the durable result files are:

- [uv_native_result.json](uv_native_result.json): full baseline thermal result.
- [uv_control_results.json](uv_control_results.json): vacuum, zero-temperature,
  and auxiliary-scale controls, including commands and point evaluations.
- [uv_cut_results.json](uv_cut_results.json): HH and HHL integrations.
- [uv_point.json](uv_point.json), [uv_structure.json](uv_structure.json):
  UV-subtracted point and compact export evidence.
- [uv_comparison.json](uv_comparison.json): all comparisons and production factors.
- [vacuum_reference.py](vacuum_reference.py), [vacuum_reference.json](vacuum_reference.json):
  independent vacuum calculation and convergence.
- [thermal_reference.py](thermal_reference.py), [thermal_reference.json](thermal_reference.json):
  independent thermal calculation, channel decomposition, and convergence.

This continuation adds cards, reference scripts, and results. It does not alter
core implementation or existing test expectations. Python formatting/lint,
JSON/TOML parsing, native integrations, and independent reference checks pass.
