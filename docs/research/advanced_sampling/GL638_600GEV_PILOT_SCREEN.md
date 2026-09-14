# GL638 at 600 GeV: pilot screen and 50-worker validation

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
explains the sixteen optimized-baseline flags. The table retains those historical
statuses; later strict build12 calibrations and two validations are
reported separately below.

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

## Completed 50-worker calibrations and first validation

All three runs below retain the reference nine-channel catalogue, c=1/default
localization/m_uv=50, all 936 orientations/six cuts/nineteen CT variants and
SymJIT O3 compression. They used the strict build12 client, fresh workspaces,
real-component stability checks and no E_cm-relative near-zero allowance.

| Run / seed | Iterations × outer draws | Scheduled points | Re ± SEM [pb] | Absolute Re ± SEM [pb] | Status |
|---|---:|---:|---:|---:|---|
| MC calibration / 61401 | 2 × 16384 | 32768 | -6.301250e-5 ± 4.752757e-5 | 2.999858e-4 ± 4.749994e-5 | Failed: 15 unstable |
| SUM calibration / 61403 | 2 × 2048 | 36864 | -6.122516e-5 ± 8.577555e-5 | 3.568572e-4 ± 8.565422e-5 | Passed: zero unstable/NaN |
| C SUM validation / 61507 | 4 × 8192 | 294912 | -2.866179e-5 ± 3.029388e-5 | 2.631879e-4 ± 3.031614e-5 | Passed: zero unstable/NaN |

The MC calibration's first eight flagged points have finite real residuals
around 2.75e-308–7.39e-306. Its historical failure remains recorded while an
explicit E_cm-relative component tolerance is implemented and tested. The
other seven flags are not relabeled by that diagnosis. The SUM calibration
scheduled 12.5% more physical points than MC and had larger reported errors;
it was a throughput calibration, not an accepted equal-work ranking against
the failed MC result.

C completed in **843.678 s of native integration** (14.06 minutes), versus the
862.849 s forecast, and 942.954 s total client time. It loaded for 49.097 s,
activated evaluators for 10.867 s and peaked at 156.49 GiB RSS. Inventory and
saved-state hashes were unchanged. Every native estimate below is cumulative;
no postprocessed estimator replaces the production result.

| C completed iterations | Re ± SEM [pb] | Absolute Re ± SEM [pb] | Absolute RSE | Largest absolute signed weight |
|---:|---:|---:|---:|---:|
| 1 | -6.032115e-5 ± 2.553699e-5 | 2.012040e-4 ± 2.555594e-5 | 12.70% | .105039 |
| 2 | -6.119655e-5 ± 2.011012e-5 | 2.322593e-4 ± 2.027240e-5 | 8.73% | .117910 |
| 3 | -2.076620e-5 ± 3.920173e-5 | 2.662285e-4 ± 3.923251e-5 | 14.74% | .871842 |
| 4 | -2.866179e-5 ± 3.029388e-5 | 2.631879e-4 ± 3.031614e-5 | 11.52% | .871842 |

A new positive maximum in iteration three supplies **77.14% of the final
signed raw second moment** and 92.83% of the signed mean's magnitude. At the
same cube the absolute-estimator maximum is .871916595, supplying **76.86% of
its raw second moment** and 10.11% of its mean. Distinct channel points can
cancel in the signed sum but not in that absolute estimator. The negative
signed maximum is -.153427736. The completed [native maximum replay](/tmp/gl638-final-physics-screen/maxima/results-validation15-C-summed-o3-build12/MAXIMUM_ATTRIBUTION.md)
passes ordinary/Arb scalar checks and all eighteen channel geometry contexts.
The positive maximum is almost entirely the named soft-LMB channel evaluated
at hard momenta: Cut-1 gluon energies are 264/369 GeV, and the dominant
integrated CTs are A +.53835, P +.27889 and Z +.05605. A's multiplier is
W_H=.23920. This point is neither a soft/HZ endpoint nor a near-coincident
projected A/Z configuration; the native radial roots remain separated. The
ordinary/Arb real discrepancy is 1.81e-14. The negative maximum is also
dominated by the soft-LMB channel at hard momenta. Thus these are stable
finite CT weights; this replay does not establish a global bound or justify
an hour-scale run. With its maps and outer factor fixed, offline c=4
reweighting predicts +.49877 at the positive point and -.15365 at the
negative point, not new observed maxima.

Successive disjoint blocks have absolute means/errors (in 1e-4 pb)
`2.0120±.2556`, `2.6331±.3147`, `3.3417±1.1049`, `2.5407±.2920`.
Their Quad/Arb counts are `33/7`, `27/9`, `29/1`, `33/4`; all blocks have zero
failures. These diagnostics subtract pooled sums and
`Q=N*mean²+N*(N-1)*SEM²` between checkpoints, with output-f64 rounding;
adapted blocks are not identically distributed replicate integrations.
Native table chi-square per degree of freedom is .29645 for Re and .94396
for absolute Re (the raw accumulated values are 1.18581 and 3.77585). C's eightfold count increase over the SUM calibration happens to reduce
its error by approximately sqrt(8), while its maximum grows by 3.50 times.
The dominant new tail point prevents interpreting that aggregate scaling as
established convergence. The signed mean is only .95 SEM from zero, so its
105.69% relative error is presently uninformative; absolute convergence needs
further independent evidence.

| Run | Native integration wall [s] | Sampling / remaining overhead / evaluators [ms per outer draw] |
|---|---:|---:|
| Failed MC calibration | 379.167 | 8.790 / 6.398 / 71.528 |
| Passed SUM calibration | 328.038 | 88.697 / 44.034 / 287.952 |
| Passed C validation | 843.678 | 77.976 / 39.840 / 282.849 |

For C, dividing by K=9 gives `8.664 / 4.427 / 31.428 ms` per scheduled
channel point. These are worker elapsed times including rescues, not process
CPU time or wall time per sample. Support checks and rotations make scheduled
points different from actual body calls. The overhead row includes integrand
and orchestration; for C, 38.680 ms is within the integrand timer and 1.160 ms
outside it. Sampling is 19.46% of measured outer-evaluation time. The native
integration loop clones workers each iteration; unassigned wall time also
includes persistence and imbalance and is not a measured cloning profile.

The full [C assessment](/tmp/gl638-final-physics-screen/validation15-build12/analysis-C15/C15_ASSESSMENT.md)
and [checkpoint moments](/tmp/gl638-final-physics-screen/validation15-build12/analysis-C15/checkpoint_analysis.json)
retain exact extrema, hashes, disjoint block moments and timing. The
[comparison](/tmp/gl638-final-physics-screen/validation15-build12/analysis-C15/analysis.json)
keeps failed MC and cheap-pilot results ineligible; their means are compatible
with C within the reported errors, without establishing equal-work superiority.

## Completed c=4 SUM validation: lower errors, one unresolved flag

A second fresh-seed SUM run (61611) completed the same four iterations of
8192 outer draws on 50 workers: K=9, 294912 scheduled points, all 936
orientations/six cuts/nineteen variants, 600 GeV, mu_r=91.188 GeV, m_uv=50 GeV,
default localization and O3 compression. It changes only the common threshold
bias to c=4 using the validated parameterized state. Build12 still uses pure
component-relative checks. **It remains failed/ineligible:** one unstable
sample appeared in iteration three; no NaNs occurred. Its cause is not yet
established, and passing saved-extremum replays do not diagnose that sample.

| Completed iterations | Re ± SEM [pb] | Absolute Re ± SEM [pb] | Absolute RSE | Largest absolute signed weight | Cumulative unstable |
|---:|---:|---:|---:|---:|---:|
| 1 | -1.015209e-4 ± 4.838387e-5 | 2.690222e-4 ± 4.851057e-5 | 18.03% | .238863 | 0 |
| 2 | -1.003237e-4 ± 2.983734e-5 | 2.392922e-4 ± 2.990955e-5 | 12.50% | .238863 | 0 |
| 3 | -5.521255e-5 ± 2.890250e-5 | 2.520454e-4 ± 2.892908e-5 | 11.48% | .468984 | 1 |
| 4 | -6.037950e-5 ± 2.253297e-5 | 2.427317e-4 ± 2.255516e-5 | 9.29% | .468984 | 1 |

The final errors are 25.6% smaller than C1, and the largest absolute signed
weight is 46.2% smaller. The signed means differ by .84 combined reported
SEM. This is encouraging provisional evidence, not a production ranking:
the source run fails stability, and the absolute target itself changes with
c. The new positive maximum still supplies **40.34% / 40.12%** of the signed /
absolute raw second moments. Its contribution is 23.70% / 5.90% of the mean
magnitudes. Absolute RSE below 10% alone therefore does not establish tail
convergence. Signed RSE remains 37.32%. The native table chi-square-per-dof is
.95769 Re / .47147 absolute; the undivided values are 3.83076 / 1.88589.

Disjoint absolute blocks (in 1e-4 pb) are `2.6902±.4851`, `2.0956±.3500`,
`2.7755±.6288`, `2.1479±.2465`; corresponding Quad/Arb/unstable counts are
`31/4/0`, `34/6/0`, `33/4/1`, `26/4/0`. The checkpoint-moment interpretation
and covariance rules above apply unchanged.

Both extrema pass ordinary/Arb scalar checks and all eighteen native channel
geometry contexts. The positive **+.468984346** comes primarily from
`lu_cut_1` (+.468294521), with integrated P=(3,12), variant 7, **+1.524224867**,
Z=(3,10,13), variant 6, **-1.025989671**, and A=(7,8), variant 9,
**-.029718315**. P/Z multipliers are one; A's transformed W_H is .037661263.
They share the two-loop subspace [3,7]. Cut-1 gluon energies are 363/402 GeV,
H/Z residuals are 618/557 GeV, and P/Z/A radial roots are 429/400/299 GeV.
Their finite derivatives and separated roots identify a hard CT weight, not
the HZ corner, a soft endpoint or demonstrated projected-root coalescence.

The negative **-.190742759** comes mainly from S, with Cut-1 total -.174284885
and Cut-3 total -.018464318. Its leading CT is Cut-1 A variant 9,
-.146528141 with W_H=.502694969, alongside Z -.038470531 and P +.023202554.
Cut-3 integrated and mixed CTs contribute with both signs. Cut-1 gluons have
energies 70/107 GeV; H/Z residuals are 70.5/-7.05 GeV. This is closer to Z,
but neither a joint HZ nor soft endpoint. Ordinary/Arb relative discrepancies
are 5.38e-15 and 2.96e-15 for the two extrema. Their successful replays do not
change the unresolved run status or demonstrate a global bound.

Native integration took **830.912 s** (13.85 min); total client time was
925.806 s, including 45.942 s load and 10.765 s activation, with 156.61 GiB
peak RSS. Sampling / remaining overhead / evaluator timers were
**83.699 / 42.704 / 270.856 ms per outer draw**, or **9.300 / 4.745 / 30.095 ms**
per K=9 scheduled point. These are worker elapsed timers including rescues,
not CPU times or measured wall-per-point throughput; the modest difference
from C1 is not a controlled timing comparison. Sampling is 21.07% of measured
outer evaluation. Final Double/Quad/Arb counts are 32626/124/18. The input
state remained unchanged.

The [C4 assessment](/tmp/gl638-final-physics-screen/validation15-build12/analysis-C4-15/C4_ASSESSMENT.md)
and [checkpoint moments](/tmp/gl638-final-physics-screen/validation15-build12/analysis-C4-15/checkpoint_analysis.json)
retain exact values, source hashes and the provisional comparison. Continue
with the flagged-sample diagnosis, stable independent runs and planned OSE
versus inverse-density controls before committing to an hour-scale run.

## Remaining independent validation

| Seed | Setup | Status |
|---:|---|---|
| 61501 | Reference C+J+R+S, c=1, MC channels | Prepared; await E_cm-relative policy and build13 validation |
| 61503 | Joint replaces standalone Cut 1, retaining R and S, c=1, MC channels | Prepared; same requirement |
| 61507 | Reference catalogue and physics, explicitly summed channels | Completed above; maximum replay passed |
| 61511 | Reference C+J+R+S, c=4, MC channels | Prepared; same requirement |
| 61611 | Reference C+J+R+S, c=4, explicitly summed channels | Completed above; failed with one unresolved unstable sample, extrema replays passed |

Remaining runs use fresh workspaces, 50 workers and at least four complete
iterations, with approximately 15-minute budgets calibrated from measured
50-worker costs. Proposed GL638 cards declare integrated energy dimension
`-2` and enable dimensionless `ecm_relative_tolerance_for_re=1e-100` only in the
final Arb level; Im and earlier levels remain zero. These cards are prepared,
not results, and require the new owner/tests before execution. The completed
C1 and C4 SUM runs did not use this policy. Native failures remain failures.

Inspect each run's exact real maxima, precision rescue, cut/CT attribution
and later iterations. Only proceed to an hour-scale run when this evidence
makes a real relative error below 10% plausible. If the signed result is near
zero through cancellation, assess absolute-real convergence and quote the
signed absolute uncertainty. Do not assume square-root error scaling when
maxima dominate. Reserve seed 60317 for that later run.

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
retain the exact historical choices and hashes. The completed C request is
authenticated by the actual build12 records. Remaining
[build13 request templates](/tmp/gl638-final-physics-screen/validation15-build13/manifest.json)
require the new build hashes and frozen MC counts; they are not run results.

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

The later completed summaries preserve each native status and authenticated
build12 settings:

| Native summary | SHA256 |
|---|---|
| [Failed MC50 calibration](/tmp/gl638-final-physics-screen/calibration50-build12/results-calibrate-retained9-monte_carlo-o3-build12/summary.json) | `c5c7e2d787bdfae070c467ae49fd3ad22990142b784be2fb96c3e6c4b7b25b37` |
| [Passed SUM50 calibration](/tmp/gl638-final-physics-screen/calibration50-build12/results-calibrate-retained9-summed-o3-build12/summary.json) | `57b15a350f056003bc3be978ddd972d449f07facea83fcccc93455b27631617b` |
| [Passed C validation](/tmp/gl638-final-physics-screen/validation15-build12/results-validation15-retained9-c1-summed-o3-build12/summary.json) | `1a6809041c70b9b4145fb1c03e59a8a0fbf342c82cf03103c5e8d3dfc31b48a6` |

The [build12 source manifest](/tmp/hosted-joint-optimized-build12-manifest.json)
has SHA256 `df8b1f614730003b61628f66ad94c0fd5a66e790cd0388f7c0212b19401277b3`;
the [strict client link record](/tmp/gl638-final-physics-screen/final-screen-build12.build.json)
has SHA256 `7aea7113f6860b005364f9d3cbf346f0dea8f35faf34fb81de99ed8957503cc3`.
The [frozen C request](/tmp/gl638-final-physics-screen/validation15-build12/request-validation15-retained9-c1-summed-o3-build12.json)
has SHA256 `8d066f77ab43615562c65cf5660a7ed18fb5390bbee47f8c50a6c3696ac9dfca`.

The [C extrema capture](/tmp/gl638-final-physics-screen/maxima/results-validation15-C-summed-o3-build12/summary.json)
has SHA256 `5cf6d094821bed8d28f70e1ecb43600ad8d5772f7b47e0e3f0e194cf82e88e04`;
its passing [scalar audit](/tmp/gl638-final-physics-screen/maxima/results-validation15-C-summed-o3-build12/scalar-audit.json)
has SHA256 `8eb997bfb485e9f143e94719bdfb13b023bdd1cd252e06cd92dad3b819e62aa9`;
the [all-channel geometry audit](/tmp/gl638-final-physics-screen/maxima/results-validation15-C-summed-o3-build12/geometry/summary.json)
has SHA256 `3e8fcc6ae0dd05f656ac7b0fef9c9a2456f81a8f8e5dd50b33d5531eb397a677`.

The later [failed C4 SUM validation summary](/tmp/gl638-final-physics-screen/validation15-build12/alternatives-summed/results-validation15-retained9-c4-summed-o3-build12/summary.json)
has SHA256 `504330d1fb823af11f635c0ff176d0ff093ebc8d83f2d3334c3b347e181b6e16`.
Its [extrema capture](/tmp/gl638-final-physics-screen/maxima/results-validation15-C4-summed-o3-build12/summary.json)
has SHA256 `75f2bf441bdf0de0ad339abc75cc7c63fde6e58928fa07bc8bb43d9332be4f59`;
the passing [scalar audit](/tmp/gl638-final-physics-screen/maxima/results-validation15-C4-summed-o3-build12/scalar-audit.json)
has SHA256 `6c779c300044b94c2bbeb881ccffc46b64fbad19776b27d2b50b42636a340284`;
the completed [all-channel geometry audit](/tmp/gl638-final-physics-screen/maxima/results-validation15-C4-summed-o3-build12/geometry/summary.json)
has SHA256 `20373b61e897c645ede32871e2afcdbeb759d8d564dc38a2473c9b99976baae5`.
