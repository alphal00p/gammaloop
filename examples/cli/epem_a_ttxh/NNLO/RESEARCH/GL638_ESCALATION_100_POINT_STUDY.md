# GL638: 100 actual precision escalations

Study date: 2026-09-15. Reference: `7d580bd01` on `advanced_sampling`.

## Result and production qualification

`min_abs_wgt_for_escalation=1e-12` had **not yet demonstrated a reduction in
production escalations** at the time of this diagnostic. The then-current
implementation required a previous signed integral estimate with at most 1%
relative Monte Carlo error. None of the inspected production checkpoints met
that requirement. The implementation has now been changed to use the
absolute real-component integral estimate and a 10% relative-error gate.

A separate experiment nevertheless provides strong evidence that the floor
would remove most expensive rescues once its reference becomes available:
96 of the first 100 actual escalations have fully weighted f64 real probes
below `1e-12 * 2e-5`. These include all 62 points that reached Arb. Their
largest real-weight correction after rescue is `1.8501e-18`. This is a
conditional replay result, not a measured production suppression rate.

The diagnostic therefore records the old behavior and its counterfactual
counts. The current waiver uses the absolute real-component accumulator,
`integral |Re|`, and accepts it once its relative error is at most 10%.

| Archived optimized-LMB checkpoint | Value |
| --- | ---: |
| Samples | 81,920 |
| Re | -2.26835e-5 +/- 2.51375e-5 |
| Integral of \|Re\| | 2.89294e-4 +/- 2.51173e-5 |
| Relative error of Re | 110.82% |
| Relative error of integral \|Re\| | 8.68% |
| Final f128 fraction | 3.2690% |
| Final Arb fraction | 3.7830% |
| Total escalated fraction | 7.0520% |

Source anchors at the reference commit:

- `crates/gammalooprs/src/integrate/mod.rs:2870`: construction of signed
  `current_integral_estimates`; the active real/imaginary accumulator is chosen
  at line 2876. Absolute accumulators are updated at lines 396–400.
- `crates/gammalooprs/src/integrands/process/mod.rs:1948`: finite reference,
  nonzero central value, at most 10% relative error, and small-weight
  comparison. The norm variant repeats this gate at line 2202.

The archived checkpoint is in
`../workspaces/GL638_optimized_lmbs_workspace_before_quiet_20260915/results/`.
At 17:42 UTC the fresh advanced run had completed iteration 1; the fresh
optimized run had not. Both were live. The archived values above are labelled
explicitly so they are not confused with new completed iterations.
With the revised absolute-real, 10% criterion, the archived optimized-LMB
checkpoint's 8.68% absolute-real error is sufficiently precise to provide a
waiver reference; the fresh optimized run had not yet produced such a
checkpoint.

## Independent experiment

A copy of the optimized-LMB state and a separate integration workspace were
created under `../local_state_GL638_escalation_diagnostic/`. The capture run
used 4,096 points on eight workers, all orientations, six optimized-LMB
channels sampled by Monte Carlo, `sampling_channel_weight=map_density`,
`E_cm=600 GeV`, `mu_r=91.188 GeV`, and `m_uv=50 GeV`. The floor was explicitly
zero. The production states/workspaces were not used for diagnostic evaluation.

A temporary, bounded diagnostic hook recorded the first 100 distinct samples
that actually retried at higher precision. It retained original hypercube
coordinates, discrete choices, outer integrator weight, native probe values,
failure reasons at every attempted level, and final native results. The hook
was removed after the experiment and the production executable restored.
Ordinary loop-momentum debug messages were not counted as escalations.

The Double relative target was `1e-6`, Quad `1e-10`, and Arb `1e-12`.
Stability followed the real integration phase, including its absolute
observable; imaginary values remained available for monitoring.

| Complete 4,096-point diagnostic run | Value |
| --- | ---: |
| Re | -6.18844e-5 +/- 4.49343e-5 |
| Integral of \|Re\| | 2.06430e-4 +/- 4.48288e-5 |
| Im | -6.67016e-6 +/- 4.39329e-5 |
| Integral of \|Im\| | 3.18457e-4 +/- 4.36502e-5 |
| Final Double | 3,796 / 92.6758% |
| Final Quad | 123 / 3.0029% |
| Final Arb | 177 / 4.3213% |
| Unstable or NaN | 0 |

The first 100 captures are a sample conditioned on escalation, not a random
100-point integrand sample or an estimate of the remaining integration error.
All had `integral_estimate=null`, as expected during the first iteration.

## Floor counterfactual

For each captured point, `MAX` is the maximum absolute fully weighted real
probe at Double, considering both signed and absolute-observable checks.
Sampling Jacobians and channel factors are already present in the probes.
The captured outer Havana weight, equal to six for every point in this first
iteration, is applied exactly once.

All 100 first failures were `ErrorThreshold`, for both checks. No captured
failure originated from sampling reconstruction or a threshold-counterterm
error. Final precisions were 38 Quad and 62 Arb. `MAX` ranged from
`6.15e-273` to `7.79e-7`, with median `3.73e-41`.

The following counts are historical counterfactuals from the old signed,
1%-gate implementation. The current implementation instead uses the absolute
real-component estimate with a 10% gate, so these counts should not be read as
the measured waiver rate after the change.

| Relative floor | Waived / 100 with reference 2e-5 | Waived / 100 with reference 2e-4 |
| --- | ---: | ---: |
| 1e-12 | 96 | 97 |
| 1e-10 | 97 | 98 |
| 1e-8 | 98 | 98 |
| 1e-6 | 98 | 99 |
| 1e-4 | 99 | 99 |

At `1e-12` and reference `2e-5`, the four retained points are #19, #49,
#53, and #94. All 62 Arb rescues are waived. The largest discrepancy between
the Double probe mean and the rescued real result among waived points is
`1.8501e-18` (#31). Raising the floor far above `1e-12` has little additional
benefit in this captured sample. No distribution-wide bias bound follows
from these 100 points; tiny inaccurate Double probes could in principle hide
a larger exact result, although that did not occur here.

## Nature of the failures

Eighty-two captures have a Double rotation mismatch of at least 1%; 77 have a
Double-to-rescued real correction larger than the rescued real value. The
relative discrepancies are real. The inefficient aspect is paying to resolve
them when the complete real contribution is negligible for the integral.
Calling them all erroneous stability comparisons would be inaccurate.

Seven points were replayed with cut-level diagnostics. Those traces separate
internal complex-evaluator cancellation from cancellation between real cut
totals. A small final real component by itself does not establish either
mechanism or an exponentially suppressed LU h-function.

| Point | Double MAX, fully weighted Re | Rescued Re | Final precision | Evidence |
| --- | ---: | ---: | --- | --- |
| 1 | 5.7927e-18 | 5.7031e-18 | Quad | Small CT real component and numerical bare real residual; active-cut h is order one. |
| 2 | 2.3322e-25 | 9.0117e-45 | Arb | Dominant cut(2,6,10) bare real residual shrinks with precision while its imaginary part persists; CT real component also needs rescue. |
| 19 | 2.1287e-17 | -2.1286e-17 | Quad | Small real contribution beside Im=-1.0044e-5; all cut h values are about 1.37–1.65. |
| 31 | 2.9361e-18 | -9.2314e-33 | Arb | Spurious real residual beside a finite imaginary component; negligible cancellation between rescued real cut totals. |
| 49 | 5.0531e-11 | -5.0532e-11 | Quad | Genuine cancellation between real cut CTs, condition about 1.23e6; cut gluon13 energy about 0.0179 GeV. |
| 53 | 7.7929e-7 | 7.7929e-7 | Quad | Genuine cancellation between real cut CTs, condition about 1.75e4. |
| 94 | 1.3029e-14 | 1.3028e-14 | Quad | Dominant cut(2,6,10) CT changes by a few parts per million; no large real inter-cut cancellation. |

For #31, the unweighted bare contribution of cut(2,4,12) has real part
`-2.55e-55` at Double and about `-1.23e-340` at Arb, while its imaginary part
remains `2.453e-39`. This is consistent with roundoff leakage into an almost
purely imaginary bare contribution. The rescued real result comes from the
CT. The relevant cut has `h=1.31`, so LU-tail suppression does not explain
that failure. Two other cuts have h around `1e-3290`, but those spectator
tails do not dominate the unstable real component.

The condition quoted for #49 and #53 is the sum of absolute real cut totals
divided by the absolute real sum, evaluated within one precision/rotation.
It is not a condition estimate of every internal numerator operation.
Gluon energies were reconstructed from the appropriate edge ray momentum
and solved cut rescaling, not inferred from a canonical loop norm. In this
graph edges 13 and 14 are gluons; canonical LMB `[3,4,7,10]` consists of top
lines and cannot by itself diagnose gluon softness.

These observations support retaining real-component checks and a small
absolute-weight waiver, rather than replacing the checks by a complex norm
or describing every failure as harmless exponential underflow. In 37 of
the 100 captures the rescued absolute imaginary weight exceeds `1e-12`.
The waiver counts refer specifically to the real-only production objective.

## Replay validation and artifacts

All 100 points were replayed with ordinary `inspect`, and again with
`inspect --use_arb_prec`, each batch loading the isolated state once.
The largest absolute discrepancy from captured native results was
`7.33e-23` for Re and `4.17e-22` for Im. For the 88 normal-range real
references, ordinary replay's maximum relative discrepancy was `1.40e-16`.
Forced-Arb replay's largest relative real discrepancy was `1.69e-10` at
point #72, whose real weight is only about `2.75e-33`.

Twelve real references are subnormal or below Double range. Six round to
zero at the final CLI boundary. Inspect returns a Double value before the
outer weight is reapplied offline, which explains subnormal rounding such
as point #24. Native captured values were retained for floor comparisons.
Forced-Arb replay leaves all 100 displayed imaginary values unchanged; its
largest real change from the normal-stack replay is `6.18e-37`.

Event-rich inspect serialization failed for non-string additional-weight
keys, and for an individual Arb event component below Double range. Scalar
replay therefore explicitly disabled event retention and additional event
weights. This avoids those output-boundary issues without changing the
physical integrand evaluation. The issues were not patched in this study.

Durable artifacts:

- [GL638_ESCALATION_100_POINTS.json](GL638_ESCALATION_100_POINTS.json): exact
  source coordinates/discrete choices, per-level reasons and maxima, initial
  weighted probes, native final values, and replay comparisons.
- [GL638_ESCALATION_SELECTED_CUTS.txt](GL638_ESCALATION_SELECTED_CUTS.txt):
  compact cut-level trace for the seven inspected points. These numbers are
  raw physical contributions, before the complete sampling weight, unlike
  the point table above. Evaluations occur in precision/rotation order.
- [GL638_ESCALATION_SELECTED_CUTS.json](GL638_ESCALATION_SELECTED_CUTS.json):
  native first/accepted cut contributions and reconstructed gluon energies.

The copied integrand binary SHA-256 is
`a33d3408acb7a2c5b553bd9fb9669061169f23272b9745009d02a28f65d92c7c`.
Raw logs, generated replay cards and copied state remain in the ignored local
diagnostic directory. They are not required to read the committed dataset.

To replay any recorded point using the production card and an available
matching state, issue `inspect -p epem_a_tth -i NNLO -x <coordinates> -d
<discrete>` after `run set_model_parameters configure_runtime`,
`set process -p epem_a_tth -i NNLO defaults`, and disabling event retention.
Coordinate and discrete arguments are space-separated. Multiply the returned
x-space weighted result by the recorded outer integrator weight when
comparing with `final_weighted_result`.

## Production log audit

Structured logfile tracing was already disabled. Verbose terminal iteration
reports were about 48–54 KB per iteration, not a current runaway file, but
could accumulate gigabytes at the 10-billion-point cap. Commit `7d580bd01`
removed optional per-bin/max/grid display flags from both production cards.
Per-iteration JSON, including maximum points, remains enabled.

The old workspaces and terminal logs were archived. Resume was rejected by
the saved integrand fingerprint check, so fresh runs were started with the
same seed and quiet cards; this was a restart, not a successful continuation
of the archived samples. Both production runs retain 50 workers and the
configured `1e-12` floor. The isolated diagnostic did not replace their
states, steering or integrator checkpoints.
