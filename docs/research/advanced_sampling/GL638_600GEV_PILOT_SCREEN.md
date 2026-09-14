# GL638 at 600 GeV: first-grid pilot screen

Recorded 14 September 2026. Twelve exploratory runs compare channels, UV mass,
localization and threshold weights. Each uses **4096 draws, seed 60101, one
initial-grid iteration and 20 workers**, with SymJIT O3 compression. These are
GL638 graph contributions to `epem_a_tth@NNLO`, not the complete process cross
section. All 936 orientations are summed, all six Cutkosky cuts and nineteen
threshold variants remain active, with 3D local UV and integrated UV CTs.
`sqrt(s)=600 GeV`, incoming momenta are `(300,0,0,±300) GeV`,
`mu_r=91.188 GeV`; `m_uv=50 GeV` except the two indicated comparisons.

All twelve integrations completed their requested counts with finite output
and no NaNs, but **retain failed/ineligible build11 status** because of the
native instability counts below. The separate [native-precision investigation](GL638_NATIVE_PRECISION_ACCEPTANCE.md#build11-follow-up-suppressed-real-component-instability)
explains the sixteen optimized-baseline flags. This note does not relabel
historical results or claim a passing run after that repair.

## Runtime choices and channel definitions

The real-component stability stack uses `check_on_norm=false`, Double/Quad/Arb
relative targets `1e-6/1e-10/1e-12`, and the rotation `(0.1,0.2,0.3)`.
The Double large-weight escalation threshold is `.9`; it is disabled in the
higher lanes. Exact zeros do not force escalation. Native sampling uses the
validated Fixed256 source for mixed surface catalogues and FixedQuad for the
optimized LMB baseline; evaluators use the configured stability stack.

Sampling uses `sampling_multichanneling=true`, `sampling_channels="monte_carlo"`,
`sampling_channel_weight="map_density"`, spherical coordinates, linear mapping,
`b=.3`, `power=2`; orientation summation is explicit. The grid has 16 bins,
minimum 64 samples per update, continuous/discrete learning rates `.25`, and
`train_on_avg=false`. There is no accuracy-based early stop. Adaptation occurs
only after this pilot's sole iteration, so its measurements do not assess the
trained-grid performance. Common seeds match initial random inputs, but
changing channel catalogues generally changes their physical interpretation.

- **C:** six channels `lu_cut_0` through `lu_cut_5`, targeting respectively
  `(2,6,12,13)`, `(2,6,10)`, `(2,6,7,13,14)`, `(2,4,12)`,
  `(2,4,10,13)`, `(2,4,7,14)`.
- **J:** `cut1_joint_HZ`, Cut 1 followed by joint normals
  `surface(2,4,12)` and `surface(3,10,13)`, parent LMB `[3,6,7,10]`.
- **R:** `cut3_right_5_10`, Cut 3 followed by right threshold
  `surface(5,10)`, parent LMB `[4,5,7,12]`.
- **S:** `soft_6_12`, explicit full-volume `lmb(6,12,13,14)`.
  The optimized baseline resolves six automatic LMB channels plus S; one
  automatic channel uses this same LMB, so the recorded baseline has seven
  entries. The reference `C+J+R+S` has nine.

Cut radial proposals use the existing `lu_h` log-logistic approximation with
broad fraction `.02`. The physical Cutkosky h-function remains
`poly_exponential`, sigma 1, power 3. The default local threshold envelope has
Gaussian width 1, sliver width 10, `smooth_sliver=false`, `dynamic_width=false`;
no dampers are forced to one. Both envelope variations remain symmetric about
the threshold radius. The bias family changes neither this envelope nor the
sampling maps: both WH and WF are transformed by
`T_c(w)=w/(w+c*(1-w))`, using exact paired complements. See the
[weight-bias prescription and controls](GL638_WEIGHT_BIAS_CONTROLS.md).

## Native results

Means, errors and maximum absolute real weights are in pb. Errors are the
native one-iteration Monte Carlo errors, without pooling or postprocessing.

| Setup | Re integral ± error | Integral of absolute Re ± error | Maximum absolute Re weight | Unstable / 4096 |
|---|---:|---:|---:|---:|
| Optimized LMBs + S | -6.38132e-05 ± 4.27e-05 | 0.000204836 ± 4.259e-05 | 0.134112 | 16 |
| C + S | 1.51441e-05 ± 5.135e-05 | 0.000211416 ± 5.125e-05 | 0.104672 | 1 |
| C + J + S | 9.10384e-06 ± 5.703e-05 | 0.000221577 ± 5.693e-05 | 0.119625 | 1 |
| C + J + R + S (reference) | -3.59679e-05 ± 5.29e-05 | 0.000230721 ± 5.278e-05 | 0.134419 | 3 |
| C + J + S; m_uv=91.188 | 8.43764e-06 ± 5.714e-05 | 0.000221093 ± 5.703e-05 | 0.119625 | 1 |
| C + J + R + S; m_uv=91.188 | -3.69751e-05 ± 5.289e-05 | 0.000229955 ± 5.277e-05 | 0.134419 | 3 |
| (C without Cut 1) + J + S | -4.93863e-06 ± 5.418e-05 | 0.000208077 ± 5.408e-05 | 0.104688 | 1 |
| (C without Cut 1) + J + R + S | 9.30025e-05 ± 0.0001103 | 0.000337487 ± 0.0001102 | 0.389058 | 3 |
| Reference; Gaussian width .5 | -3.59637e-05 ± 5.265e-05 | 0.000229427 ± 5.253e-05 | 0.134463 | 3 |
| Reference; Gaussian/sliver .5/.5, smooth | -3.70907e-05 ± 5.255e-05 | 0.000231201 ± 5.243e-05 | 0.134501 | 3 |
| Reference; c=.25 | -0.000106459 ± 0.0001038 | 0.000346945 ± 0.0001037 | 0.355389 | 3 |
| Reference; c=4 | -1.78574e-05 ± 3.727e-05 | 0.000177916 ± 3.716e-05 | 0.079034 | 3 |

The nine-channel reference is retained for its known Cut 3/right coverage,
not because this pilot proves it optimal. Changing localization reduced the
reported errors by only .47–.65%, with slightly larger maxima; there is no
useful evidence to tune it. Changing m_uv from 50 to 91.188 likewise offers no
clear improvement. The replacement catalogue with R encountered a `.389058`
maximum. Its precise replay attributes most of it to Cut 1, with integrated
counterterms A=(7,8) and Z=(3,10,13) contributing about +24.31 and -23.91.
Their projected radial roots differ by .04648 GeV. The sampled point is far
from H/Z and the gluons have finite energies; this is a projected A/Z near
intersection with stable dual cancellation. A surviving asymptotic singularity
has not been established. The reference maximum instead combines sizable
opposite-sign contributions from several cuts in S. All four distinct extrema
agree with forced Arb within the configured real tolerance. The
[maximum attribution](/tmp/gl638-final-physics-screen/maxima/results-diagnostic-o3-build11/MAXIMUM_ATTRIBUTION.md)
retains the exact cut/CT and geometry evidence.

The c=4 pilot reduced both real and absolute-real errors by **29.6%**, and its
observed maximum by **41.2%**, relative to c=1. This warrants independent
validation. The c=.25 pilot approximately doubled the errors; retain its data
but do not prioritize it for the next round. The c=1 and c=.25 maxima occur at
the same complete Sample in S; c=4 moves the maximum to R. Replays must include
CT multipliers at their projected points, not only at the sampled point.

For a fixed subtraction prescription, the absolute observable takes the
absolute value after summing physical graph/cut/CT contributions at the same
point, and before summing contributions from different sampling-channel
points. Explicit channel summation therefore needs one covariance-preserving
update of that full per-draw absolute estimator. It must not use the absolute
value of the signed sum over different points. Across different c or envelope
settings, the signed integral is the invariance target; the absolute integral
may change because the subtraction representation changes. No such change
alone proves improved sampling. Aggregate errors do not provide the covariance
needed for a formal paired-difference test. No final cross-section estimate or
global maximum bound follows from these first-grid pilots.

Timing was recorded but shared host activity prevents ranking settings by
wall time. Reference native sampling/evaluator/remaining-overhead timers are
8.81/56.97/5.27 ms per outer draw; c=4 gives 8.65/55.48/5.21 ms with nearly
identical precision counts. The localization runs overlap diagnostics and
link activity and show much larger times, despite similar rescue counts.
Neither the wall-time residual nor these differences measure cloning cost.

## Next validation stage — prepared, not executed here

Keep c=1/default localization/m_uv=50 for two strict 50-worker throughput
calibrations: MC seed 61401 with two iterations of 16384 draws, and summed
channels seed 61403 with two iterations of 2048 draws. These schedule 32768
and 36864 channel points respectively, so they are throughput measurements,
not an exact equal-work variance comparison. Use their actual timings to
freeze multi-iteration counts before four independent, approximately
15-minute runs on 50 workers:

| Seed | Setup |
|---:|---|
| 61501 | Reference C+J+R+S, c=1, MC channels |
| 61503 | Joint replaces standalone Cut 1, retaining R and S, c=1, MC channels |
| 61507 | Reference catalogue and physics, explicitly summed channels |
| 61511 | Reference C+J+R+S, c=4, MC channels |

Each uses a fresh workspace and at least four completed native iterations.
Reuse c=1 MC throughput to budget c=4 because its maps are identical and its
pilot evaluator cost and rescue mix are nearly identical; the resulting
15-minute duration remains approximate. Inspect each run's exact saved real
maxima, precision rescue, cut/CT attribution and later-iteration behavior.
Native failures remain failures. Only proceed to an hour-scale run if this
evidence makes a real relative error below 10% plausible. If the signed result
is near zero through cancellation, assess absolute-real convergence and quote
the signed absolute uncertainty instead. Do not assume square-root error
scaling when maxima dominate. Reserve seed 60317 for that later run.

## Reproducibility

The following are local execution artifacts, not committed state payloads.
Each summary embeds the authenticated request, effective settings, inventories,
native output and saved-state before/after hashes. All input states remained
unchanged. Full maximum coordinates and precision/timing records are in the
[combined analysis](/tmp/gl638-final-physics-screen/analysis-all12-build11/analysis.json).

| Native summary | SHA256 |
|---|---|
| [Eight channel/UV cases](/tmp/gl638-final-physics-screen/results-prune8-candidates-o3-build11/summary.json) | `9d4b7aad82de1224b08ca01e362423ce24510b810bba64eb13f26c0d2b9e2c64` |
| [Two localization cases](/tmp/gl638-final-physics-screen/compact-runtime-screen-build11/results-localization-o3-build11/summary.json) | `0368c1482d518eb57b2546e169b4b3dcb65d1d272eb8fcda5ba877043377fe6d` |
| [Two metadata-bias cases](/tmp/gl638-final-physics-screen/compact-runtime-screen-build11/results-metadata_bias-o3-build11/summary.json) | `f8a5ec87458ff0f207ef040f64c37eef0ae0dd3190c959f5e530864f68ae1db7` |

The [source build manifest](/tmp/hosted-joint-optimized-build11-manifest.json)
has SHA256 `d5a2701bc0caeddff526272cc37a71a239d5f44d11706b527ba856a4de4108bf`;
the [candidate-client link record](/tmp/gl638-final-physics-screen/final-screen-candidates-build11.build.json)
has SHA256 `dd7db67cb4da41c3463ff40a2cb9844bb6b8895ba6da2cf791e825c5dd7c5d48`.
The client allows continued exploration only after a finite, count-complete
screen with recorded instability; it preserves each failed status.

[Original-state hashes](/tmp/gl638-hosted-joint-gate/provenance/state_sha256_expected.json)
authenticate the 35-file state; [parameterized-state hashes](/tmp/gl638-final-physics-screen/weight_variants/state-sha256.json)
authenticate the 13-file state required for c=.25/4. The
[O3 configuration](/tmp/gl638-symjit-next/configs/symjit-o3-compress-on.toml),
[reference runtime card](/tmp/gl638-final-physics-screen/component_cards/joint_cut3_soft_muv50.toml),
[calibration preflight](/tmp/gl638-final-physics-screen/calibration50-build12/preflight.json)
and [c=4 validation request](/tmp/gl638-final-physics-screen/validation15-build12/request-validation15-retained9-c4-monte_carlo-o3-build12-template.json)
retain the exact choices and hashes. The future requests still require the
actual build12 link records and calibrated counts; they are not run results.

The completed [c=4 maximum replay](/tmp/gl638-final-physics-screen/maxima/results-diagnostic-c4-o3-build11/MAXIMUM_ATTRIBUTION_C4.md)
passes scalar and native-geometry audits; ordinary/Arb real differences are
below 7.2e-14 relative. Its largest weight, +0.079033984 in R, is at exactly
the c=1 positive-extremum point (+0.081350612). Only Cut 1 changes there,
from +0.195405938 to +0.193089310. Offline reweighting of the captured bare
CTs predicts that the former negative maximum shrinks from -0.134419040 to
-0.058801695; this prediction is not a new physical evaluation. The same
algebra reproduces the directly replayed common positive point to 3.7e-300
relative. The actual new negative maximum is -0.062649818 at a different
S point, dominated by stable Cut-1 integrated A/Z cancellation with hard
gluons (171.7 and 157.1 GeV). Thus the maximum changes location through
finite CT coefficients, without evidence of a new asymptotic singularity.
The [capture](/tmp/gl638-final-physics-screen/maxima/results-diagnostic-c4-o3-build11/summary.json)
and [scalar audit](/tmp/gl638-final-physics-screen/maxima/results-diagnostic-c4-o3-build11/scalar-audit.json)
preserve the failed source-pilot status and production ineligibility.
