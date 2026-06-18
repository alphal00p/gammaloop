# qq_hhh_2L Integration Debug Report

Date: 2026-06-16

## Summary

I did not proceed to a production 1% integration because inspect calls show that the threshold-subtracted integrand currently produces pervasive NaNs and extreme finite threshold-counterterm weights at ordinary fixed points in the unit hypercube. This is not limited to one endpoint point, not fixed by `SingleParametric`, not fixed by regenerating with `check_esurface_at_generation=false`, and not cured by enabling all LMB channels with inverse-Jacobian multichannel weights.

The strongest isolation test is:

- With threshold subtraction enabled, `random_bulk = [0.37, 0.42, 0.61, 0.28, 0.74, 0.53]`, group 0, LMB channel 0 gives `nan=true`, event weights `NaN`, and threshold CTs up to about `1e277`.
- With the same regenerated diagnostic state, same point/group/LMB channel, and `subtraction.disable_threshold_subtraction=true`, the inspect result is finite and stable:
  - evaluation: `+1.9635608626661010e-5 -2.4892454384933456e-5 i`
  - integrand result before parameterization Jacobian: `+4.9045290575357402e-25 -6.2175697308675424e-25 i`
  - metadata `nan=false`
  - all six event weights finite, and only `original` additional weights are present.

This points to GammaLoop threshold CT generation/evaluation for this exact-xi auxiliary topology as the current blocker, not the imported DOT file itself and not simply the compiled `SummedFunctionMap` backend.

I also tested the tempting generation-time filter `assume_positive_external_energies=true`. It does reduce the generated threshold surface set, but it is not usable for this exact-xi DOT: runtime inspect aborts because a trimmed threshold surface is still required by the configured kinematics. The steering file therefore intentionally keeps `assume_positive_external_energies=false`.

## Card Changes

The steering card now sets `global.generation.threshold_subtraction.check_esurface_at_generation=false` in all relevant generation command blocks and global defaults. It also explicitly keeps `global.generation.threshold_subtraction.assume_positive_external_energies=false`; setting it true trims exact-xi helper surfaces needed by runtime validation. The DOT was not edited.

I also added `inspect_nan_points.py`, a task-local inspect matrix driver that can target alternate state folders and runtime blocks. It now records fatal GammaLoop inspect errors explicitly in `summary.tsv` / `summary.json`.

## Generated States

- `state`: original compiled state, with `SummedFunctionMap` / assembly output.
- `state_no_esurface_check`: attempted full compiled regeneration with `check_esurface_at_generation=false`; interrupted because it remained in heavy evaluator construction for too long.
- `state_no_esurface_check_light`: successful diagnostic regeneration with:
  - `check_esurface_at_generation=false`
  - `compile=false`
  - `summed_function_map=false`
  - `summed=false`
  - `store_atom=true`
- `state_assume_positive_light`: diagnostic regeneration with `assume_positive_external_energies=true`, `check_esurface_at_generation=false`, `compile=false`, and `store_atom=true`. Generation succeeded with fewer threshold CTs, but inspect is invalid for the physical runtime settings because GammaLoop reports that a trimmed threshold counterterm is required.
- `state_all_lmbs_light`: diagnostic regeneration with `override_lmb_heuristics=true`, `check_esurface_at_generation=false`, `assume_positive_external_energies=false`, `compile=false`, and `store_atom=true`. Generation succeeded with peak RAM about `1.05 GiB`, active f64 backend `eager`, `.bin size 225.09 MiB`, `.so size 0 B`, 19 LMB channels in group 0 and 21 in group 1.

The light diagnostic state generated successfully with peak RAM about `1019.80 MiB`, active f64 backend `eager`, `.bin size 460.00 KiB`, `.so size 0 B`.

## Original Integration Attempts

Using the compiled `SummedFunctionMap` state:

- `integrate_smoke` with 128 samples returned `~1e292` central values with infinite errors and `nan_or_unstable=87.5%`.
- `integrate_lmb_summed` returned a first iteration around `1e23` but with `nan_or_unstable=95.90%`.
- `integrate_no_lmb_mc` reached `1e296` scale weights and panicked during grid update.
- `integrate_graph_summed_lmb_summed` returned all zero while reporting `nan_or_unstable=100%`.

Those integrations are not meaningful as central-value estimates because the sampling is dominated by NaN/unstable threshold CT evaluations.

## Inspect Findings

Important logs:

- Compiled-state initial smoke/max-weight inspect logs:
  - `logs/04_inspect_smoke_max_graph1_lmb0.log`
  - `logs/06_inspect_lmb_summed_max_graph1.log`
- Endpoint scan:
  - `logs/nan_endpoint_scan/nan_endpoint_scan_summary.tsv`
- Broader compiled-state matrix:
  - `logs/nan_point_matrix/summary.tsv`
- Corrected-generation light-state matrix:
  - `logs/nan_point_matrix_no_esurface_check_light/summary.tsv`
- `assume_positive_external_energies=true` validation check:
  - `logs/nan_point_matrix_assume_positive_light_check/summary.tsv`
- All-LMB light-state matrices:
  - `logs/nan_point_matrix_all_lmbs_light_g0/summary.tsv`
  - `logs/nan_point_matrix_all_lmbs_light_g1/summary.tsv`

The broad matrix showed NaNs not only near the suspected endpoint but also at bulk points such as `[0.5]*6` and `[0.37, 0.42, 0.61, 0.28, 0.74, 0.53]`.

For the corrected light state, the focused matrix gave:

- `lmb_summed_max`: all inspected group/LMB combinations are `nan=true`.
- `unit_mid`: group 0 all LMB channels are `nan=true`; group 1 has finite but enormous CT-dominated values in some channels.
- `random_bulk`: group 0 all LMB channels are `nan=true`; group 1 channel 0 is finite but unstable/huge, channels 1 and 2 are `nan=true`.

Representative corrected-state bad point:

`logs/nan_point_matrix_no_esurface_check_light/random_bulk_group0_lmb0.log`

Examples from its event rows:

- `threshold_counterterm:14:0 ~ -4.23e276 -8.80e276 i`
- `threshold_counterterm:6:0 ~ -3.66e277 -4.85e277 i`
- `threshold_counterterm:30:0 ~ -4.11e277 -9.97e276 i`
- several CT entries are explicitly `NaN`.

All-LMB diagnostic:

- Group 0, `random_bulk`, LMB channels 0 through 18: every channel has `nan=true`. The largest finite exponents remain in the `1e274` to `1e297` range.
- Group 1, `random_bulk`, LMB channels 0 through 20: some channels avoid literal NaNs, but the finite CT-dominated values are still around `1e292` to `1e294`, and many channels still have `nan=true`.

Positive-external-energy diagnostic:

- `state_assume_positive_light` generated successfully, reducing the active threshold surface list.
- Inspect then fails before numerical evaluation with:
  `Amplitude integrand 'subtracted' was generated with specialized threshold-subtraction assumptions, but the current runtime model parameters require a trimmed threshold counterterm ...`
- This is consistent with the exact-xi helper legs being part of the external signature used by generation filtering.

## What Was Ruled Out

- **DOT mismatch:** the DOT is byte-identical to the pygloop file according to the existing context and prior verification.
- **Compiled backend only:** `SingleParametric` inspect still has NaNs.
- **One endpoint only:** bulk and midpoint samples also produce NaNs.
- **`check_esurface_at_generation=true` only:** regenerating with this setting false still produces NaNs and huge CT weights.
- **Base/IR-subtracted DOT alone:** disabling runtime threshold subtraction at the same bad point gives finite stable values.
- **Positive-external-energy threshold filtering:** `assume_positive_external_energies=true` generates a smaller state but makes runtime inspect invalid because required threshold CTs were trimmed.
- **LMB heuristic coverage:** `override_lmb_heuristics=true` activates all LMB channels, but group 0 remains `nan=true` for every channel at a bulk point.

## Current Blocker

The threshold counterterms generated/evaluated for this exact-xi auxiliary amplitude are pathological. They create NaN event weights and finite values around `1e293-1e297` in fixed inspect calls, before integration has any chance to converge.

Because threshold subtraction is required for the intended integration, the next useful work is to debug or alter the threshold CT construction/evaluation, not to continue VEGAS/Cuhre tuning. Multichanneling or max-weight monitoring can help after the CTs are pointwise finite; at present they cannot cure the pointwise NaNs.

## Suggested Next Debug Steps

1. Add a way to inspect raw threshold CT stacks at the actual sampled momentum point, not only event-level reconstructed additional weights.
2. Use the successful `store_atom=true` light state for standalone stack debugging, but avoid waiting on `rust-script` unless its dependency cache is ready; it did not produce output within about one minute.
3. Compare generated threshold surfaces between `state_no_esurface_check_light`, `state_assume_positive_light`, and `state_all_lmbs_light`; the positive-energy filter removes surfaces that runtime still needs, while the all-LMB state keeps the problem pointwise visible in every group-0 channel.
4. Add or use a generation mode that excludes exact-xi helper external legs from the positive-energy threshold filter while still generating every threshold CT that `get_existing_esurfaces` can require at runtime.
5. Only after fixed-point inspect is finite with threshold subtraction enabled, return to `SummedFunctionMap` generation and integration with 3 cores.

## Follow-up Raw-CT Trace

The follow-up debug copy `../qq_hhh_2L_experiment` contains a focused no-compile `SingleParametric` inspect trace at `x = [0.37, 0.42, 0.61, 0.28, 0.74, 0.53]`, `group = 0`, `lmb = 0`. The detailed write-up is in `../qq_hhh_2L_experiment/qq_hhh_2L.md`, and the main raw logs are:

- `../qq_hhh_2L_experiment/debug_logs/inspect_threshold_paramcheck.stdout`
- `../qq_hhh_2L_experiment/state_threshold_trace_paramcheck/logs/gammalog-threshold-paramcheck.jsonl`

That focused trace narrows the first nonfinite value to the raw threshold CT evaluator output itself: r-star solving, h factors, UV damping, pass-one parameters, original graph evaluation, and multichannel prefactors are finite. Out of 108 threshold CT evaluations, 24 raw pass-one threshold evaluator calls return `NaN + NaN i`, with representative affected esurfaces including `GL05:{10,33,14,29,28}`, `GL05_isr_p1_ct:{29,6,30}`, `GL06:{37,21,34}`, and `GL06_isr_p1_ct/GL06_isr_p2_ct:{35,38}`. This points to an internal undefined operation in generated threshold CT expressions at finite r-star points for the exact-xi auxiliary topology.

The later raw-evaluator tape trace in `../qq_hhh_2L_experiment/qq_hhh_2L.md` identifies the first concrete undefined operation: for `GL05`, esurface `10`, `left_th_1`, orientation `1`, the f64 evaluator inverts
`OSE(9) + OSE(12) + OSE(13) - Q(0,0) - Q(1,0)`. A follow-up xi variation shows that this is a pass-on/tree denominator that becomes degenerate with the active helper threshold surface only because the original helper vector `xi=(1000,0,0,0)` equals `Q(0)+Q(1)` at `sqrt(s)=1000 GeV`. Changing all helper vectors to `xi=(2000,0,0,0)` removes this specific raw-CT NaN in the focused inspect.

The most explicit runtime trace is in
`../qq_hhh_2L_experiment/debug_logs/rstar_coincident_xi1000_ids/inspect_display.stdout`
and
`../qq_hhh_2L_experiment_xi2000/debug_logs/rstar_coincident_xi2000_ids/inspect_display.stdout`.
For xi1000, `GL05` esurfaces `8`, `10`, `16`, and `33` all vanish on each
other's selected r-star. For xi2000, `10` and `33` remain degenerate with each
other, but they no longer coincide with the physical incoming-energy
denominators `8` and `16`. This is why xi2000 removes the observed raw CT NaN
without proving that the threshold evaluator is robust to all coincident-surface
limits.

The xi2000 experiment copy also reran the earlier fixed-point matrix over five
points, both graph groups, and LMB channels `0,1,2`. All 30 inspect calls are
`nan=false` and stable; the summary is in
`../qq_hhh_2L_experiment_xi2000/debug_logs/nan_point_matrix_xi2000/summary.tsv`.

## xi=2000 Integration Probe

I then tried integration with the xi2000 debug state:

```text
../qq_hhh_2L_experiment_xi2000/state_threshold_trace_esurface_xi2000
```

The state is uncompiled `SingleParametric`, with `check_esurface_at_generation=false`,
threshold subtraction enabled, graph Monte Carlo sampling, orientations summed,
and LMB multichanneling enabled.

Two runs were made with 3 cores:

- 128 samples:
  `../qq_hhh_2L_experiment_xi2000/workspaces/integrate_smoke_single_parametric/integration_result.json`
  reports `nan_percentage = 0.0` and `nan_or_unstable_percentage = 0.0`.
- 2048 samples:
  `../qq_hhh_2L_experiment_xi2000/workspaces/integrate_2048_single_parametric/integration_result.json`
  reports `result = -4.268903936148298e-05 - 1.0356238568615117e-03 i`,
  `error = 1.5903700965315063e-03 + 1.793993148690794e-03 i`, and
  `nan_percentage = nan_or_unstable_percentage = 0.1953125`.

So xi2000 is integrable in the sense that the previous catastrophic
`1e292-1e299` threshold-CT weights disappear, but this uncompiled f64 setup is
not production-clean: 4 out of 2048 sampled evaluations were still marked
nonfinite, and the statistical error is nowhere near 1%.

The reported max-weight points from the 2048 run all inspect as finite and
stable. The largest recorded evaluations were only `O(1)`:

```text
re+ graph 1: +1.4708489536031286
re- graph 0: -1.0049571860639608
im+ graph 0: +1.2230554187507177
im- graph 0: -2.1077816907847247
```

A replay of the same seed with 512 samples captured one nonfinite point:

```text
graph group = 1
x = [
  0.04243809892765449,
  0.6501739569562547,
  0.6345802747849795,
  0.9999282161172643,
  0.3850124895740257,
  0.3432911025875355,
]
```

The logs are:

```text
../qq_hhh_2L_experiment_xi2000/state_threshold_trace_esurface_xi2000/logs/gammalog-xi2000-integrate-nonfinite-512.jsonl
../qq_hhh_2L_experiment_xi2000/state_threshold_trace_esurface_xi2000/logs/gammalog-xi2000-bad-sample-group1-f64.jsonl
../qq_hhh_2L_experiment_xi2000/debug_logs/parametric_nonfinite_bad_sample_xi2000/
```

This remaining xi2000 NaN is not the old raw threshold-CT `1/0` problem. At the
bad point, all 588 threshold `amplitude_threshold_pass_one_raw` records have
`raw_nonfinite=false`; the threshold CT sums are finite. The first nonfinite
quantity is the pre-threshold bare/original CFF evaluator, e.g.

```text
GL12 bare_cff = NaN + NaN i
cts = finite, about 1e-100 for the first failing record
diff = NaN + NaN i
```

Enabling the parametric nonfinite dump shows that the f64 original evaluator,
and also the f128 upcast evaluator, return `NaN`, while the arb evaluator
returns finite tiny values, typically around `1e-45` to `1e-60` for the dumped
orientation outputs. There are 7984 dumped nonfinite orientation evaluations
across the graph group, all with `f128_nonfinite_count=1` and
`arb_nonfinite_count=0`.

This points to a second issue in GammaLoop's stability path: nonfinite f64/f128
orientation results are not forcing escalation to arb in the integration path.
The final f64 sample is accepted as a NaN metadata event and zeroed in the
reported integrand result. For xi2000, the next code-side fix to try is to make
the stability check treat any nonfinite rotation/orientation output as unstable
and escalate to the next precision level. That is separate from the original
xi1000 threshold-CT degeneracy.

For the user's specific question about the remaining `10`/`33` coincidence:
the xi2000 traces show `10` and `33` still coincide in helper-helper clusters,
but this is not currently producing the observed `1/0`. The old `1/0` required
the helper surface to also coincide with a physical/pass-on denominator such as
`OSE(9)+OSE(12)+OSE(13)-Q(0,0)-Q(1,0)`. With xi2000 that physical/helper
coincidence is gone in the focused `GL05` trace; the remaining `10`/`33`
coincidence is still conceptually something a robust threshold implementation
should handle as a grouped limit, but empirically it was not the source of the
remaining xi2000 integration NaNs.

## Post-Rescue xi3000 Integration Update

I added a stability-rescue patch in
`crates/gammalooprs/src/integrands/process/mod.rs`: if a nonfinal stability
level produces any nonfinite rotated/orientation result, that level is now
treated as unstable and the evaluation escalates. The previously captured bad
xi2000 integration point now escalates f64 -> f128 -> arb automatically and
returns the same finite value as explicit arb:

```text
-1.9680617192257578e-9 + 1.7812184200189852e-9 i
```

This confirms that the remaining xi2000 NaNs after the threshold-CT degeneracy
fix were a stability-rescue acceptance problem in the base/original CFF
evaluator, not another raw threshold-CT `1/0`.

The xi-choice scan found:

```text
xi=(1000,0,0,0): raw_nonfinite=16, GL05 cluster {8,10,16,33}
xi=(2000,0,0,0): raw_nonfinite=0,  GL05 clusters {5,8,16,31}, {10,33}
xi=(3000,0,0,0): raw_nonfinite=0,  GL05 clusters {5,31}, {8,16}, {10,33}
```

So xi3000 reduces coincidences further than xi2000, but does not completely
lift `10`/`33`. The remaining `10`/`33` pair did not cause a raw CT NaN in the
focused inspect. The attempted spatial vector `xi=(2000,0,0,1000)` was
inconclusive because the debug trace was too slow and only wrote a partial log.

I then tried several 3-core integrations with orientations summed and threshold
subtraction enabled, all using the uncompiled `SingleParametric` xi2000
experiment state. All post-rescue runs had `nan_percentage=0`.

The best short probe was xi3000 with common-radial LMB-summed sampling:

```text
512 samples:
re = -0.002632496494663313 +/- 0.0016044710966468282
im =  0.002562385132968542 +/- 0.0020449130799561816
arb = 0%
```

The largest completed run was:

```text
xi = (3000,0,0,0)
sampling = spherical_common_radial, LMB channels summed, inverse-Jacobian weights
graphs = monte_carlo
orientations = summed
neval = 16384
cores = 3
```

Result:

```text
re =  0.014629629118902268 +/- 0.012149741876449939  (83.05%)
im = -0.025686244113766153 +/- 0.02684835496656811   (104.52%)
nan_percentage = 0
nan_or_unstable_percentage = 0
f64 = 99.9755859375%
f128 = 0.01220703125%
arb = 0.01220703125%
real chi^2 = 8.003461698959141
imag chi^2 = 2.590717782587208
```

This is not percent-level accurate. The max-weight points from this run inspect
as stable in f64 and agree with explicit arb, for example:

```text
graph 1 max re+/im+:
default = 76.163074005184839 + 45.880456999624741 i
arb     = 76.163074004555384 + 45.880456999245574 i

graph 0 max im-:
default = 3.8839380452326684 - 18.265505503870678 i
arb     = 3.8839380452003183 - 18.265505503718334 i
```

Therefore the current blocker to a 1% central value is no longer pointwise
NaNs, and raising precision will not solve the variance. The remaining problem
is real sampling variance from stable near-boundary enhancements, especially
points with `x[3]` close to one. The useful next debug step is to scan around
the 16k max-weight points and build a map/channel targeted to that boundary
structure.

I did not use the existing compiled state under `examples/cli/qq_hhh_2L/state`
for the final result. It contains `SummedFunctionMap` artifacts, but it was
generated with the old `xi=(1000,0,0,0)` defaults and
`check_esurface_at_generation=true`. It should be regenerated from the current
`qq_hhh_2L.toml` before being used for production statistics.

## Fresh Compiled-State Follow-up

I then regenerated a clean compiled `SummedFunctionMap` state from the current
production TOML into:

```text
../qq_hhh_2L_experiment_xi2000/state_compiled_no_esurface_xi2000
```

This used all-LMB generation (`override_lmb_heuristics=true`),
`check_esurface_at_generation=false`, threshold subtraction enabled, and
`summed_function_map=true`. Generation completed with:

```text
1080 generated C++ sources
peak RAM = 11.12 GiB
```

Compiled xi3000 common-radial LMB-summed, 2048 samples:

```text
re =  0.006433623626125279 +/- 0.004639988169733452  (72.12%)
im = -0.004373117159518832 +/- 0.00363038765515719   (83.02%)
nan_percentage = 0
arb_percentage = 0
average_total_time = 0.141591464 s/sample
```

Compiled xi3000 common-radial LMB-channel Monte Carlo is much faster:

```text
4096 samples:
re = -0.004131466377325436 +/- 0.004325668075877935  (104.70%)
im = -0.01073900630204229  +/- 0.010986384085259385  (102.30%)
average_total_time = 0.007632997 s/sample

32768 samples:
re =  0.00047410859333351074 +/- 0.0024297372578875242  (512.49%)
im =  0.014239183449377572   +/- 0.012576294983522572   (88.32%)
nan_percentage = 0
arb_percentage = 0.006103515625
average_total_time = 0.007050956 s/sample
```

The 32768-sample dominant max point is graph 0, LMB channel 14:

```text
x = [0.30784578180460542, 0.90170572955766448,
     0.71240971451878965, 0.75722185209864434,
     0.36452634459630484, 0.265760901346222]
```

Inspect with event/additional-weight output is stable:

```text
default = -1.3627481340754974 + 10.782175613902147 i
arb     = -1.3627481340613690 + 10.782175613790022 i
nan = false
```

The visible event-level decomposition is finite. One representative large term
is `threshold_counterterm:11:0 ~ -1.93e-17 + 2.27e-16 i` multiplied by the full
factor `4.5066e16`, so the max point is a real stable contribution, not an
undetected precision failure.

Conclusion: the cleaned-up threshold setup and stability rescue make the
integrand pointwise usable, and compiled SFM with LMB MC gives practical
throughput. It still does not get close to 1% with the tested maps. The
dominant remaining work is to study/remap stable LMB-channel enhancements,
especially graph 0 / LMB channel 14 around `x[1] ~= 0.90` and `x[3] ~= 0.76`.
