= BNL Benchmark Investigation

#quote(block: true)[
#strong[Status:] Archived, non-reproducible performance investigation

This record preserves measurements captured on 2026-05-25 at evidence revision
`4fdbf430b29edd24b9e1292c07ab24a151427dc7`. The BNL card, DOT graph, and ratio script remain,
but the generated `examples/cli/BNL/states/` directory, named workspaces, and result logs do not.
The command blocks named below are also absent from the current card. Treat the values as
historical evidence; reproducing them requires reconstructing the evidence-era setup.
]

This note records the BNL setup, generated states, and
runtime/integration scans used for the 2 -> 5 amplitude graph in
`examples/cli/BNL/BNL.dot`. All heavy generation and integration runs were
executed behind a 30 GiB RSS guard. The final selected setup is CFF,
`SingleParametric`, spherical coordinates with linear mapping, summed graphs,
summed orientations, no LMB multichanneling, one Euler-angle stability rotation,
and Double+Quad stability levels.

== Kinematics and helicities

The final scan uses the same massless 2 -> 5 phase-space point and helicity
assignment for all four masses; only `MT` and `ymt` are changed.

```toml
e_cm = 1000.0
helicities = [-1, 1, -1, 1, 1, -1, 1]
momenta = [
  [500.0, 0.0, 0.0, 500.0],
  [500.0, 0.0, 0.0, -500.0],
  [163.5141597587927, -155.7860143642699, -12.32488779805292, 48.12167194606676],
  [146.5267245560251, -18.94527369806406, 33.44069886802352, -141.3961713514227],
  [162.3143263141399, 4.025183527597962, -151.4903823484999, 58.14122874765845],
  [296.8259181161096, 175.6763342487531, 137.9119270904214, -195.5089553360614],
  [230.8188712549325, -4.970229714017088, -7.53735581189216, 230.6422259937584],
]
```

The kinematic point was selected from the random 2 -> 5 all-massless samples
tested during the scan. The selected helicities flip the two initial leptons and
the outgoing muon pair relative to the original trial point; this avoided a
suppressed below-threshold central value and gave substantially smaller relative
errors.

== Threshold localization knobs

Two independent threshold-localization widths matter:

- Integrated threshold CT: `subtraction.integrated_ct_settings.range.h_function_settings.sigma`.
- Local threshold CT: `subtraction.local_ct_settings.uv_localisation.gaussian_width`.

No new setting was needed for the local width. The local counterterm uses
`evaluate_uv_damper`, with `gaussian_width` as the width and `sliver_width` as
the hard support gate. Event-level inspect checks confirmed that changing
`gaussian_width` changes individual `threshold_counterterm:*` additional
weights while leaving the original contribution unchanged.

The selected production values are:

```toml
subtraction.integrated_ct_settings.range.h_function_settings.sigma = 1.0
subtraction.local_ct_settings.uv_localisation.gaussian_width = 1.0
subtraction.local_ct_settings.uv_localisation.sliver_width = 10.0
```

Important scan observations:

- `sigma = 0.25` was unstable for the selected threshold point.
- `sigma = 4.0` sometimes reduced variance but changed the finite-statistics
  central value strongly in small samples, so it was not used for the final
  four-mass comparison.
- `gaussian_width` in the range `0.5..4.0` was viable; `1.0` was retained as a
  stable default and is the value used in the final integrations.
- `cpe_iterations = 10000` stayed within generation budget but was slower at
  runtime than the optimized `cpe_iterations = 5` state.
- A `summed_function_map = true` state was generated, but the corresponding
  `SummedFunctionMap` runtime path failed with an invalid helicity-type error in
  this branch, so it was not used for the final numbers.

== Generated states

At capture time, the clean final states were kept under `examples/cli/BNL/states/`.
That generated directory is not part of the retained repository evidence.

#table(
  columns: 6,
  table.header([*Key*], [*Threshold CTs*], [*Optimized*], [*Pieces*], [*Guard time*], [*Peak RSS*]),
  [N], [no], [no], [1], [23.325 s], [3.477 GiB],
  [N-opt], [no], [yes], [1], [26.375 s], [4.343 GiB],
  [T], [yes], [no], [4], [28.413 s], [3.715 GiB],
  [T-opt], [yes], [yes], [4], [32.411 s], [3.383 GiB],
)

The evidence-era state names and generation invocations were:

- N: `final_no_threshold_reference`; `run generate_bnl_no_threshold_reference -D mt=510.0 -D ymt=510.0`
- N-opt: `final_no_threshold_reference_optimized`; `run generate_bnl_no_threshold_reference_optimized -D mt=510.0 -D ymt=510.0`
- T: `final_threshold_reference`; `run generate_bnl_threshold_reference -D mt=400.0 -D ymt=400.0`
- T-opt: `final_threshold_reference_optimized`; `run generate_bnl_threshold_reference_optimized -D mt=400.0 -D ymt=400.0`

Recorded generation command:

```bash
./target/release/gammaloop --clean-state --no-save-state \
  -s ./examples/cli/BNL/states/final_threshold_reference_optimized \
  ./examples/cli/BNL/BNL.toml \
  run generate_bnl_threshold_reference_optimized -D mt=400.0 -D ymt=400.0
```

The generation command blocks call `save state -o` explicitly. `--no-save-state`
only prevents an extra save attempt at process exit.

== Inspect benchmark timings

All timings below used `inspect --bench 5 --n_batches 10 --minimal-integrand`.
The capture-time warmup-based auto-sizing undershot the requested wall time on this
point, but the reported uncertainty is still the standard error over 10 batches.

#table(
  columns: 5,
  table.header([*State key*], [*Threshold runtime*], [*Total / sample*], [*Evaluator / sample*], [*Samples*]),
  [N], [off], [`1.9658e-3 +/- 1.48e-5 s`], [`1.8846e-3 +/- 1.47e-5 s`], [743],
  [N-opt], [off], [`4.8728e-4 +/- 1.40e-6 s`], [`4.1889e-4 +/- 1.48e-6 s`], [961],
  [T], [on], [`2.5149e-3 +/- 3.88e-6 s`], [`2.4001e-3 +/- 3.70e-6 s`], [716],
  [T-opt], [on], [`7.2629e-4 +/- 2.00e-6 s`], [`6.2314e-4 +/- 2.01e-6 s`], [862],
)

Recorded benchmark commands:

```bash
./target/release/gammaloop --read-only-state \
  -s ./examples/cli/BNL/states/final_no_threshold_reference_optimized \
  ./examples/cli/BNL/BNL.toml \
  run bench_bnl_reference \
  -D mt=510.0 -D ymt=510.0 \
  -D disable_threshold_subtraction=true \
  -D output_file=final_no_threshold_reference_optimized_bench.json

./target/release/gammaloop --read-only-state \
  -s ./examples/cli/BNL/states/final_threshold_reference_optimized \
  ./examples/cli/BNL/BNL.toml \
  run bench_bnl_reference \
  -D mt=400.0 -D ymt=400.0 \
  -D disable_threshold_subtraction=false \
  -D integrated_h_sigma=1.0 \
  -D local_gaussian_width=1.0 \
  -D output_file=final_threshold_reference_optimized_bench.json
```

== Final four-mass integrations

The final integration scan uses only the optimized setup. R1 and R2 are above
threshold and keep threshold subtraction active. R3 and R4 are below threshold
and disable threshold subtraction.

#table(
  columns: 6,
  table.header([*Label*], [*`MT = ymt`*], [*Threshold*], [*Samples*], [*Time / sample / core*], [*Unstable*]),
  [R1], [400.0], [on], [2,000,000], [1.5759 ms], [0.01855%],
  [R2], [450.0], [on], [2,000,000], [1.4609 ms], [0.01930%],
  [R3], [510.0], [off], [3,000,000], [0.9246 ms], [0.04133%],
  [R4], [550.0], [off], [3,000,000], [0.9088 ms], [0.05173%],
)

#table(
  columns: 3,
  table.header([*Label*], [*Real part*], [*Imaginary part*]),
  [R1], [`+1.076(99)e-16`], [`-1.9963(96)e-15`],
  [R2], [`+3.308(73)e-16`], [`-1.9622(74)e-15`],
  [R3], [`+1.4173(26)e-15`], [`-5.336(27)e-16`],
  [R4], [`+7.683(19)e-16`], [`-2.541(19)e-16`],
)

Recorded integration commands:

```bash
./target/release/gammaloop --read-only-state \
  -s ./examples/cli/BNL/states/final_threshold_reference_optimized \
  ./examples/cli/BNL/BNL.toml \
  run integrate_bnl_reference \
  -D mt=400.0 -D ymt=400.0 \
  -D disable_threshold_subtraction=false \
  -D integrated_h_sigma=1.0 \
  -D local_gaussian_width=1.0 \
  -D workspace=final_bnl_R1_mt400_threshold_2M_5c \
  -D n_start=200000 -D n_max=2000000 -D seed=990401 -D n_cores=5

./target/release/gammaloop --read-only-state \
  -s ./examples/cli/BNL/states/final_no_threshold_reference_optimized \
  ./examples/cli/BNL/BNL.toml \
  run integrate_bnl_reference \
  -D mt=510.0 -D ymt=510.0 \
  -D disable_threshold_subtraction=true \
  -D workspace=final_bnl_R3_mt510_no_threshold_3M_5c \
  -D n_start=300000 -D n_max=3000000 -D seed=990403 -D n_cores=5
```

The other two final commands are the same with `mt=450.0`, `seed=990402`,
`workspace=final_bnl_R2_mt450_threshold_2M_5c`, threshold subtraction enabled;
and `mt=550.0`, `seed=990404`,
`workspace=final_bnl_R4_mt550_no_threshold_3M_5c`, threshold subtraction
disabled.

== Ratios

At capture time, `BNL.py` read the four final `integration_result.json` files by default and
computed the ratios with first-order propagated independent real/imaginary uncertainties. Those
default workspaces are absent from this checkout. The recorded invocation was:

```bash
python3 ./examples/cli/BNL/BNL.py \
  --json-output ./examples/cli/BNL/outputs/BNL_final_ratios.json
```

Final ratios:

#table(
  columns: 3,
  table.header([*ratio*], [*real part*], [*imaginary part*]),
  [R2/R1], [`+0.9889(60)`], [`+0.1124(61)`],
  [R3/R1], [`+0.3046(40)`], [`+0.6936(39)`],
  [R4/R1], [`+0.1476(22)`], [`+0.3769(22)`],
)

== Extended precision and corroboration update

The old final workspaces could not be resumed after the later branch changes:
the integration workspace manifests contain generated-integrand fingerprints
that no longer match the current loaded generated states. The resume guard
correctly rejected this with:

```text
Workspace integrand fingerprints do not match the current generated integrands
for BNL@default.
```

For that reason, the high-statistics extension below was run as a fresh
current-branch workspace, not by bypassing the fingerprint check.

#table(
  columns: 4,
  table.header([*Label*], [*Samples*], [*Time / sample / core*], [*Unstable*]),
  [R1 extended], [10,000,000], [1.5650 ms], [0.01541%],
  [R2 reference], [2,000,000], [1.4609 ms], [0.01930%],
  [R3 reference], [3,000,000], [0.9246 ms], [0.04133%],
  [R4 reference], [3,000,000], [0.9088 ms], [0.05173%],
)

#table(
  columns: 3,
  table.header([*Label*], [*Real part*], [*Imaginary part*]),
  [R1 extended], [`+1.002(48)e-16`], [`-1.9973(46)e-15`],
  [R2 reference], [`+3.308(73)e-16`], [`-1.9622(74)e-15`],
  [R3 reference], [`+1.4173(26)e-15`], [`-5.336(27)e-16`],
  [R4 reference], [`+7.683(19)e-16`], [`-2.541(19)e-16`],
)

The corresponding workspaces were `current_bnl_R1_mt400_threshold_10M_5c`,
`final_bnl_R2_mt450_threshold_2M_5c`, `final_bnl_R3_mt510_no_threshold_3M_5c`, and
`final_bnl_R4_mt550_no_threshold_3M_5c`. None is retained in the repository.

The best independent-error ratios available in this record are:

#table(
  columns: 3,
  table.header([*ratio*], [*real part*], [*imaginary part*]),
  [R2/R1], [`+0.9882(44)`], [`+0.1160(43)`],
  [R3/R1], [`+0.3020(23)`], [`+0.6945(22)`],
  [R4/R1], [`+0.1462(14)`], [`+0.3773(13)`],
)

This improves the previous ratios but does not reach the requested 0.25%
relative uncertainty for every real and imaginary component. The bottleneck is
the small real part of R1. After 10M samples its uncertainty is still
`4.76e-18`; propagating independent errors gives a 3.74% relative uncertainty
on `Im(R2/R1)`. Reaching 0.25% by brute-force independent sampling would require
roughly `(3.74 / 0.25)^2 ~= 224` times more R1 statistics, i.e. about
`2.2e9` R1 samples at the observed variance. At the measured runtime this is an
order-100-hour computation on five cores, so the ratio target is not practical
with the independent-result estimator used in this record alone. A correlated multi-slot
ratio estimator, or a variance-reduced estimator specifically for the small R1
real component, is needed for that precision target.

Corroboration runs were also made with independent seeds and safe runtime
variations. The threshold-subtracted checks keep the UV scale and threshold
localization fixed and only vary the spherical mapping parameter, because
changing `m_uv` or threshold-localization widths changes the low-statistics
central values in this local-subtraction-only setup and is not a clean
corroboration knob.

#table(
  columns: 3,
  table.header([*Label*], [*Variation*], [*Samples*]),
  [R1 map], [`b=0.8`, new seed], [500,000],
  [R2 map], [`b=1.2`, new seed], [500,000],
  [R3 map], [`b=0.8`, new seed], [500,000],
  [R4 MC], [graph/orientation Monte Carlo, new seed], [500,000],
)

#table(
  columns: 3,
  table.header([*Label*], [*Real part*], [*Imaginary part*]),
  [R1 map], [`+5.7(1.9)e-17`], [`-1.998(19)e-15`],
  [R2 map], [`+3.17(16)e-16`], [`-1.941(16)e-15`],
  [R3 map], [`+1.4046(65)e-15`], [`-5.312(68)e-16`],
  [R4 MC], [`+7.78(31)e-16`], [`-2.51(31)e-16`],
)

The missing workspace names were `corrob_map_R1_mt400_threshold_b0p8_500k_5c`,
`corrob_map_R2_mt450_threshold_b1p2_500k_5c`,
`corrob_phys_R3_mt510_no_threshold_b0p8_500k_5c`, and
`corrob_phys_R4_mt550_no_threshold_graph_orientation_monte_carlo_500k_5c`.
