# GL638 threshold-subtraction validation

The Higgs-aware WH/WF metadata restores the expected soft powers on the studied
representative rays. Full physical Standard Model generation succeeds for all
936 orientations with direct three-dimensional local UV, orientation
localization, and integrated UV. The prescription still has two distinct
limitations: coincident A/P radial residues need a stable combined evaluation,
and a hard H/Z corner has measured inverse-radius scaling on the studied rays.
This is consistent with a locally integrable first moment but potentially
logarithmic variance under the present smooth proposal.

The executable prescription is in
[GL638.dot](../../../examples/cli/epem_a_ttxh/NNLO/graphs/GL638.dot), with generation
and integration controls in
[the card](../../../examples/cli/epem_a_ttxh/NNLO/epem_a_tth_NNLO_test_GL638.toml).
[The implementation plan](../../../IR_SAFE_THRESHOL_UPDATE.md) records the
grouping invariants and extended research history. Automatic construction of
threshold metadata remains deferred.

[Portable replay commands and inputs](REPLAY.md) cover the soft precision
failure, coincident/nearby A/P roots, and hard/smooth sliver boundary controls.
Regenerate saved integrands after changing DOT metadata or generated exact-pi
coefficients; localization and Euler stability settings apply at runtime to
loaded states.

[The runtime follow-up](runtime-followup.md) checks the remaining mechanisms
with 54 A/P stability comparisons, three 30,000-sample pilots on 20 cores,
and current-binary H/Z and maximum-weight replays. Generic Euler rescues the
tested nearby A/P points; exact coincidence remains undefined. Adding the
already generated LMB `[2,10,13,14]` improves the short pilot. A density aligned
with the threshold normals targets the remaining hard-corner variance problem.

## Prescription

Use Q=1000 GeV, mt=173 GeV, mh=125 GeV. In generation LMB [3,4,7,10],
p=q3=q12, s=q7=q8, t=q10=q5, a=q2=q4−q3. The gluons are q13=p−t and q14=s−p.
Only physical cut 1, (2,6,10), has explicit variants, with parent [3,6,7,10].
Subspace [7] varies s; [3,7] varies p and s at fixed a,t. Other cut sides retain
their maximal one-loop defaults. No additional group IDs are assigned.

At each variant's own cut-preserving star point, define

```
H = eta(2,4,12), P = eta(3,12), Z = eta(3,10,13)
G0 = eta(2,6,12,13), G4 = eta(2,4,10,13)
WH = H² / [H² + (P Z/Q)²]
WF = G0² G4² / [G0² G4² + P² Z²]
```

| Threshold on cut 1 | Shared [7] weight | Native [3,7] weight |
|---|---:|---:|
| A (7,8), U (8,12,14) | 1−WH | WH |
| F (8,10,13,14) | 1−WF | WF |
| V (3,7,14), P (3,12), Z (3,10,13) | absent | 1 |

The exact common zero in WH uses the supplied complementary lazy conditional,
WH=0 and 1−WH=1. This assigns the point; it does not establish continuity.
There is no denominator regulator. The full resolved catalogue has nineteen
variants. Every physical cut first solves its own LU rescaling; shared threshold
groups use the resulting distinct external data. Overlaps, centers and
e-surface multichanneling remain independent for each compatible cycle group.

## Numerical evidence

The historical soft scans used the physical numerator and local UV with
integrated UV disabled before its generation fix. Lambda is a gluon's spatial
momentum radius in GeV. Entries are fitted powers of the momentum-space density;
single-soft and fixed-ratio double-soft measures are lambda² d lambda and
lambda⁵ d lambda. These are representative-ray results, not uniform bounds.

| Treatment | Soft q13 | Soft q14 | Both soft |
|---|---:|---:|---:|
| CTs off, individual cuts | −3 | −1 | −4 |
| CTs off, sum over cuts | −1 | −1 | −2 |
| CTs on, default metadata | −3 | −1 | −3 |
| WH/WF candidate, sum over cuts | approximately −1 | −1 | −2 |

The current matched full-UV comparison uses two freshly generated states with
all 936 orientations: the current WH/WF DOT and a copy with only the
`threshold_counterterms` attribute removed. Topology, generation LMB,
factorized numerator, model parameters, direct-3D local UV with orientation
localization, integrated UV, and orientation catalogues match exactly.
Both states use release SHA256
`57ccd168da44ad2b75b1e1b998b0680a8ccf8fc406ce3deb775356fa3567845a`
with the clean local Symbolica backport
`a277a8dadbaae5a8d6605f2f381991a601e1f345`; this is the rebuilt release with the
integrated-profile, overlap-witness and nonfinite-stability fixes.

There are 24 frozen raw points: q13, q14 and both soft, both ray signs, and
Lambda=1, 0.1, 0.01, 0.001 GeV. Lambda is the original total spatial soft
radius before each cut's own LU rescaling; in the double-soft family each
gluon has radius Lambda/sqrt(2). Each point is evaluated with CTs off,
default CTs and WH/WF CTs, in normal and forced Arb precision. Another 24
forced-Arb CT-off evaluations in the candidate state check the bare integrand.
All 168 outputs are finite and finally Stable. The bare totals and all six
individual cuts agree bit-for-bit between states at every point, isolating the
observed difference to threshold handling.

The table gives Arb powers fitted over the last decade, Lambda=0.01 to
0.001 GeV; intervals span the two signs. They describe finite-range trends,
not exact exponents or global bounds. CT-off Re is an arithmetic residual
and is not assigned a physical power.

| Treatment | Soft configuration | Summed Re power | Summed Im power |
|---|---|---:|---:|
| CTs off | q13 | residual | −0.99679 to −0.99657 |
| Default CTs | q13 | −3.00004 to −2.99994 | −3.00005 to −2.99991 |
| WH/WF CTs | q13 | −0.99842 to −0.98957 | −1.00105 to −0.99872 |
| CTs off | q14 | residual | −1.00002 to −0.99996 |
| Default CTs | q14 | −1.00003 to −0.99986 | −0.98109 to −0.97980 |
| WH/WF CTs | q14 | −1.00003 to −0.99987 | −1.03413 to −0.99519 |
| CTs off | Both | residual | −2.00001 to −1.99836 |
| Default CTs | Both | −3.00001 to −2.99998 | −3.00002 to −2.99997 |
| WH/WF CTs | Both | −1.00147 to −1.00013 | −2.05646 to −2.03057 |

Before summing, the dominant CT-off imaginary cut powers are approximately
−3 for q13 (cuts 0, 1, 3, 4), −2 for q14 (cuts 0, 2, 3, 5), and −5 for
both (all six cuts). Default CTs leave uncompensated q13 Re terms of order
Lambda^-3 in cuts 0, 3 and 4, while cut 1 remains approximately constant.
With WH/WF, cut 1 participates in that cancellation and the sum approaches
Lambda^-1. Including the single-soft measure Lambda² d Lambda, the default
q13 behavior is approximately d Lambda/Lambda, whereas WH/WF gives
Lambda d Lambda. The former indicates logarithmic radial nonintegrability
if its nonzero leading coefficient persists over an angular neighborhood;
the two rays alone do not prove that angular statement. Both q14 treatments
are consistent with Lambda d Lambda. With the double-soft measure
Lambda⁵ d Lambda, default CTs give approximately Lambda² d Lambda;
WH/WF gives Lambda⁴ d Lambda in Re and Lambda³ d Lambda in Im. These
fixed-ratio radial first moments are consistent with integrability, while
angular and hierarchical limits remain unproved.

Wider fit windows retain relevant sign changes. Positive-q14 WH/WF Im has
a tail-three slope of −1.38052, versus −1.03413 over the last decade.
The WH/WF double-soft Im changes sign in wider windows; the last-decade
slopes are −2.03057/−2.05646. The changed individual-cut powers relative to
the historical table cannot be attributed uniquely to integrated UV because
orientation coverage also changed.

Every matched evaluation resets the full runtime profile and uses norm
checking with Euler(0.1,0.2,0.3), the Double/Quad/Arb ladder, and the broad
local envelope (Gaussian width 1, sliver width 10, fixed width, smooth cutoff
disabled). The integrated profile is unchanged. Normal CT-off/default/WH-WF
runs accept Double/Quad/Arb in counts 11/11/2, 12/12/0 and 6/16/2.
Their worst relative complex normal-versus-Arb discrepancies are respectively
7.2531e-7, 7.2355e-8 and 1.7805e-6 for totals, and 7.0054e-11, 4.4447e-8
and 2.2000e-8 for individual cuts. A Stable norm flag is a heuristic rather
than a strict componentwise error bound.

Earlier full-UV scans using the special pi/2 z probe exposed a precision
blind spot: positive double-soft Lambda=0.01 with CTs off gave Double Im
+3.098e-26 versus Arb −4.799e-31, while individual cuts were about 1e-15.
The common rounding error survived that probe. Sixteen controls showed that
the generic Euler probe catches all four selected CT-off/on failures, while
separate Re/Im checks with the same z probe still miss three. The current card
therefore exposes the Euler angles as `stability_alpha`, `stability_beta` and
`stability_gamma`; the matched results above all use that corrected policy.

[The 168 recorded rows](ordinary-soft-matched-values.csv) contain raw totals,
all six cuts, flags, precision histories, measure powers and source hashes.
[The 24 unchanged input points](ordinary-soft-points.json) include the exact
executed model/runtime reset commands and generation/source provenance.
The CSV values are the recorded JSON numbers, including f64 exports of Arb
results; no sampling Jacobian or radial measure has been folded into them.
The local executed protocol is archived under
`twenty-core/default-metadata-soft/{PLAN.json,run.py,execution}` and requires
both generated states. [The portable replay instructions](REPLAY.md) illustrate
the same raw-inspect/reset mechanism, but their script consumes a different
JSON schema and does not directly run this 168-case panel.

Current complete-graph runs include integrated UV, all 936 orientations, all
six optimized LMB channels, and the full factorized numerator. The original
runs and three-profile pilots use the pre-normalization-fix release binary with
SHA256 `ce79cac84b5dca80413441fa07fbbb4ab3d3ec60bdaba9dcb73624d6c48a4c5a`.
Their default integrated profile is mathematically unchanged by the fix;
these records are not post-fix binary replays.
The smooth-cutoff follow-up below identifies its rebuilt binary separately.
Explicit final orientation summation retains direct-3D localization during
generation; it reduces repeated evaluator calls and fits five workers in about
24.4 GB. The earlier selector-sum representation exceeded 30 GB at five workers.

| Proposal | Completed samples | Nominal Re [pb] | Nominal Im [pb] | Absolute relative errors Re / Im | Finite unstable / NaN |
|---|---:|---:|---:|---:|---:|
| Linear, OSE | 50,000 | (−1.826 ± 1.135)e−4 | (0.903 ± 2.274)e−4 | 14.65% / 17.00% | 2 / 0 |
| Linear, inverse Jacobian | 52,500 | (1.379 ± 1.158)e−4 | (−0.401 ± 1.773)e−4 | 13.88% / 13.82% | 0 / 0 |
| Power 2, inverse Jacobian | 60,000 | (−2.828 ± 1.730)e−4 | (−3.019 ± 3.190)e−4 | 18.40% / 19.75% | 19 / 0 |

All used seed 1337, five pinned CPUs, real-component training and at most ten
iterations. Wall times were 693.24 s, 788.21 s and 904.93 s. Inverse Jacobian
stopped during its eighth iteration because a watchdog process query timed out;
power 2 stopped cleanly at its 900 s budget during iteration nine. Only complete
iterations contribute. The ordinary linear run is the best balanced tested
nominal option, not an established optimum. At matched 52,500 samples power 2
had 1.59/2.03 times larger Re/Im errors.

The monitors measure integrals of |Re f| and |Im f| separately. Finiteness of
both would imply complex absolute integrability through
`max(|Re f|, |Im f|) <= |f| <= |Re f| + |Im f|`.
Their estimates
are compatible across proposals and support a finite first moment, but late
large samples still move the errors appreciably. Finite unstable final-precision
values are retained in totals and grid training. Error bars therefore neither
bound their numerical errors nor establish finite sampling variance. Complete
iteration statistics are in [reference](reference-iterations.csv),
[inverse Jacobian](inverse-jacobian-iterations.csv), and
[power 2](power2-iterations.csv). Snapshot paths in those files identify local
research records rather than files distributed with this repository.

Both preferred-run extrema agree between normal and Arb precision within
1.35e−12 per phase. One forced-Arb replay narrowly misses its 1e−12 rotation
tolerance. Cut 1 dominates both. At the imaginary extremum the native two-loop
U CT, with multiplier 0.561917, supplies the largest local term, partly canceled
by Z and V. The historical adaptive PDF is unavailable for those replays.

### Longer runs and maximum-weight diagnostics

A 20-core baseline completed 300,000 samples in 2,715.50 seconds, with four
75,000-sample iterations, linear b=1 and inverse-Jacobian channel weighting.
Peak process-tree memory was 70.08 GB. The 2,700-second budget interrupted an
unfinished fifth iteration, which contributes no samples to the results.
The absolute integrals were (0.83692 ± 0.07330)e−3 pb and
(1.43250 ± 0.15476)e−3 pb for Re and Im. The largest-weight impacts were
6.29% and 5.22%. Seven finite unstable values were retained;
there were no NaNs. The absolute estimates stayed compatible over the four
checkpoints, supporting a finite first moment without certifying a finite
variance. See [iteration records](twenty-core-baseline-iterations.csv).
Here maximum impact means `|w_max|/(N*|estimate|)`, using the corresponding
signed or absolute estimate. It is a descriptive statistic, not the exact
contribution of that point to the estimate combined across adaptive iterations.

The new maxima identify two different regions. The real maximum is dominated
by native two-loop A: its local and integrated terms are +3.40429i and +2.95777
after the replay's parameterization and channel factors. Its own-star
multiplier is WH=0.175659. Actual debug tracing, validated against release Arb,
finds a common center close to p=s=t. At the A star, H*=−0.90172 GeV and
Z*=−12.78656 GeV, whereas sampled H=+13.46493 and Z=+27.06685 GeV.
Sampling only the base H surface therefore need not concentrate on the
counterterm's projected geometry.

The largest imaginary weight is instead single-soft: |q13|=0.776 GeV in cut 1.
Shared-one-loop A contributes +679.781i and native-two-loop A contributes
+75.1248i, canceled strongly by A terms −481.659i, −228.000i and −56.3559i
on cuts 0, 4 and 3. The replayed total is −2.89318−9.82744i. Normal and Arb
agree to 1.8e−10 and 2.3e−11 relative per component. Saved Monte Carlo maximum
magnitudes are 15.78085 and 22.41442; replay values need not equal these because
the historical adaptive PDF is absent. Multipliers of different variants are
evaluated at different stars and need not sum to one at a sampled point.

A scratch normalized proposal mixes 25% H-surface sampling with the six
ordinary channels. A two-sided quadratic map around the analytic H root gives
q_H proportional to |H|^−1/2. On the saved hard H/Z rays, the measured radial
second-moment power changes from approximately −1 to −0.54, approaching the
integrable −1/2 limit. This controls that transverse family; weights still
grow as R^−1/2. In a 256-draw paired full-graph pilot, no global efficiency
gain was established: nominal absolute errors remain 44–86%, and the largest
samples occur away from the targeted surface. All 457 unique evaluations were
finite and stable. Four selected extrema agree with Arb within 1.85e−12 in
complex relative difference. No failures were dropped or values replaced.

The completed controlled follow-up changes only b to 0.3 and uses the same
300,000 samples and approximately 45 minutes on 20 cores:

| Radial scale | Absolute Re integral [pb] | Absolute Im integral [pb] | Largest-weight impact Re / Im | Finite unstable / NaN |
|---|---:|---:|---:|---:|
| b=1 | (0.83692 ± 0.07330)e−3 | (1.43250 ± 0.15476)e−3 | 6.29% / 5.22% | 7 / 0 |
| b=0.3 | (0.74662 ± 0.05147)e−3 | (1.33760 ± 0.10871)e−3 | 3.55% / 3.51% | 12 / 0 |
| b=0.3, Euler probe, before overlap fix | (0.74662 ± 0.05147)e−3 | (1.33760 ± 0.10871)e−3 | 3.55% / 3.51% | 23 / 0 |
| b=0.3, Euler probe, overlap fix | (0.74662 ± 0.05147)e−3 | (1.33760 ± 0.10871)e−3 | 3.55% / 3.51% | 16 / 0 |

Thus the nominal errors fall by approximately 30% in both components, with
nominal estimator variances about 0.493 times the baseline. This is a useful
finite-run improvement; it changes neither the singularity's asymptotic power
nor the outstanding numerical limitations. The absolute estimates remain
compatible. The b=0.3 run took 2,715.15 seconds and peaked at 69.91 GB.
[Its iteration records](twenty-core-b03-iterations.csv) preserve the full
completed prefix. The same seed does not imply identical adaptive grids.

The corrected Euler-probe repeat, before the overlap fix, completed all four
iterations in 2,372.36 s,
with 69.89 GB peak process-tree RSS. Its signed Re/Im estimates are
(-2.07668 ± 5.14918)e−5 pb and (6.13720 ± 10.87401)e−5 pb. Absolute central
values differ from the earlier b=0.3 run by only -3.18e−11 and -4.85e−11 pb;
nominal errors and recorded maximum weights are effectively unchanged.
The stronger probe flags 23 finite unstable samples, versus 12 previously,
and reaches Double/Quad/Arb for 297,973/2,004/23 samples. All are retained.
This supports the finite-run comparison under the improved policy; it does
not certify those unresolved samples or finite variance. The elapsed times
are not matched throughput benchmarks: the earlier runs reached their time
caps while attempting a fifth iteration. See the
[Euler iteration records](twenty-core-b03-euler-iterations.csv).
Both extrema have exactly the same recorded coordinates and discrete channels
as in the earlier b=0.3 run. Fresh Euler-policy replays are finite and Stable
in Double and forced Arb; the relative complex differences are 3.21e−13
and 1.83e−15. Replayed Jacobians agree within 6.44e−16. These checks validate
the saved extrema, not the 23 unidentified unstable integration samples.

After the overlap fix, the matched repeat completed 300,000 samples in
2,320.02 s with 70.08 GB peak RSS. Double/Quad/Arb counts became
297,989/1,995/16; finite unstable evaluations fell from 23 to 16, with no NaNs.
Signed Re/Im estimates are (-2.07668 ± 5.14918)e−5 pb and
(6.13721 ± 10.87401)e−5 pb. Their shifts are -2.99e−11 and +1.30e−10 pb;
absolute shifts are +2.89e−11 and +6.27e−11 pb. Nominal errors are unchanged
within 2.49e−10 relative. All sixteen per-iteration extrema retain identical
coordinates, and the four final signed weights are bit-identical. The fix
reduces flagged evaluations in this run without establishing a further
variance gain or validating the sixteen remaining unstable samples. See the
[post-overlap iteration records](twenty-core-b03-euler-socp-fixed-iterations.csv).
This repeat precedes the separate nonfinite-status correction below.

A further scratch proposal targets H at the native A star, using the ideal
fixed-complement origin or (t,t) center and the full radial-projection Jacobian.
The measured SOCP center differs from (t,t) by about 1e-6 GeV at the saved
maximum: proposal normalization is exact, but runtime-star alignment is
approximate. Symmetric
Gaussian/sliver local-CT profiles were compared on identical points.
These investigations preserve the threshold catalogue, common-center treatment
and the physical integral; changing localization changes its local representation.

A larger 4,096-draw pilot per proposal compares the A-star H map with either
a rational radial law or the normalized integrated-CT radial profile. Each
mixture retains 75% ordinary sampling and includes the full projection
Jacobian. All 6,128 unique evaluations are finite; one rational-proposal point
remains unstable and is retained. Empirical Re/Im variances relative to ordinary
sampling are 0.313/0.534 for the rational law and 0.424/2.129 for the integrated
profile. The latter therefore offers no measured global improvement. These
panels are dominated by very few draws: their largest Im weights supply
70.7% and 83.0% of the absolute estimates, versus 72.7% for ordinary sampling.
The apparent rational-law gain needs a larger independent experiment.

All twenty-six Arb replays are finite. Among the selected Monte Carlo points,
agreement is within 4.58e-12 in relative complex difference; diagnostic rays
reach 4.46e-7. All precision flags remain reported.
Before the overlap fix, the two A-star maxima narrowly missed the requested
Arb rotation tolerance; another low-weight point retained a 9.40% rotation
discrepancy even at Arb. Repeating its Arb value exactly did not resolve that
discrepancy. Matching
the integrated radial profile helps some fixed points and worsens others;
it does not by itself address the remaining local CTs or other projections.

A fresh paired Arb replay of that point and its exact z rotation reproduces
the discrepancy: the total differs by 17.75% in complex relative difference
(9.40% under the two-probe norm metric), chiefly through cut 3. With threshold
CTs disabled, every cut and the total agree exactly at the exported precision.
This localizes the discrepancy to threshold subtraction. A completed debug
trace reproduces the release primary total and all six cuts exactly. It finds
the first divergence in the shared C7 overlap structure: the identity probe
builds seven centers/subsets, while the z probe finds one center containing
all seven instances. The latter has at least 70.9 GeV interior margin and
therefore certifies that the full overlap exists. Both probes have identical
catalogues; fixed complements rotate exactly. The C5 centers agree within
1.54e−12 GeV and the native two-loop origin is identical.

Capturing and replaying the exact SOCP matrices isolates a status-handling bug:
identity terminates with `NumericalError` after 23 iterations, whereas z
terminates with `InsufficientProgress` after 24. Their candidate centers agree
after inverse rotation within 2.62e−10 GeV, and both lie at least 70.9 GeV
inside all seven surfaces. The previous status whitelist accepted only the
second candidate. Three one-ULP changes in fixed coordinates suffice to change
the termination status; the catalogue, matrix ordering and cone ordering are
otherwise identical to the reduced reproducer.

Both overlap implementations now recompute every physical E-surface at the
candidate before interpreting the solver status. A finite, strictly interior
center certifies overlap without requiring an optimal epigraph solution. The
existing margin test, masses and native complements are unchanged. An invalid
candidate yields no overlap only for `PrimalInfeasible`; other unresolved
statuses propagate an error instead of silently removing intersections. The
permanent GL638 regression retains all seven instances and checks one full
overlap under identity, z and Euler rotations. All 21 focused overlap and
localization tests, the workspace check and Clippy pass. The correction is
published as `cb4e82ddf`.

The repaired release (SHA256
`2de024dd8dac62d9688bb67b55e197a0ddccd9b2f381a0ea6c8d4f5989fb5b2e`)
completed 71 full-state replays in 99.78 s (5.35 GB peak
process-tree memory). At the repaired point, explicit z and global Euler Arb
rotations differ from identity by 3.01e−12 and 3.34e−12 in relative complex
value; the largest per-cut discrepancy is 3.62e−12. The old z result is
unchanged, while identity moves by 17.745% to agree with it. Normal evaluation
now returns Double Stable and agrees with Arb within 1.08e−11. Forced Arb
still reports finite instability against the strict 1e−12 target: the internal
Euler norm estimate is 5.20e−11. Explicit Euler inputs have also passed through
f64 coordinate rounding, unlike an internal Arb rotation.

Both integration maxima (normal and Arb), all 48 ordinary-soft Euler controls,
and all 11 finite portable controls are bit-identical to their frozen pre-fix
results. Both known exact A/P invalid cases persist, with no new invalid
cases. Thus the false overlap decision is repaired, while the separate
coincident A/P residue problem remains. The matched 300,000-sample, 20-core
b=0.3 Euler repeat completed with the statistics reported above.

A separate stability fix rejects NaN/Inf probes before either accuracy metric
or the single-probe shortcut, including at final precision. Release SHA256
`57ccd168da44ad2b75b1e1b998b0680a8ccf8fc406ce3deb775356fa3567845a`
completed 19 replays in 93.27 s with the source state unchanged: all 13 portable
controls, the repaired point under normal/Arb Euler checks, and both maxima
under normal/Arb checks. All 17 finite controls retain bit-identical numerical
payloads and precision flags. Both exact A/P cases now finish Arb Unstable with
no accuracy estimate, retaining `is_nan=true`, null cut-1 weights and sanitized
zero totals. Their finite intermediate Quad rotation error, when present,
remains unchanged. This corrects failure reporting; the coincident-root failure
and the known nearby-A/P false acceptance remain. See [the replay notes](REPLAY.md).

The self-contained regression requires no generated SM state:

```bash
cargo test --locked -p gammalooprs --lib subtraction::overlap_subspace::regression_tests::gl638_c7_keeps_a_certified_full_overlap_after_solver_failure -- --exact --nocapture
```

With the structured SOCP log enabled, its identity and z matrices match the
captured live problems exactly and reproduce their respective solver statuses
on the tested platform. The test asserts physical geometry and covariance,
not a platform-dependent termination status.

### Symmetric local-CT localization

Existing runtime controls implement
`D(rho)=exp[-((rho-r_star)/(g Q))²] theta(L Q-|rho-r_star|)`.
The positive and negative radial mirror terms give the local smearing kernel
zero signed-radial PV contribution; this does not set the integrated CT's real
component to zero. The envelope is even about r-star and equals one at the pole.
The integrated
threshold profile is separate and remains unchanged in this comparison.

| Local profile | g / L | Empirical variance / baseline Re / Im | Maximum magnitude / baseline Re / Im |
|---|---:|---:|---:|
| Default | 1 / 10 | 1 / 1 | 1 / 1 |
| Narrow Gaussian | 0.1 / 10 | 1.087 / 3.090 | 1.024 / 2.951 |
| Hard sliver times broad Gaussian | 1 / 0.1 | 0.983 / 4.287 | 0.985 / 3.452 |
| Smooth compact sliver times broad Gaussian | 1 / 0.1 | 1.143 / 6.157 | 1.048 / 4.267 |

These are paired 1,024-draw ordinary-proposal pilots at b=0.3, not long adapted
integrations. Each profile's broader comparison evaluated 1,548 unique points,
all finite and finally stable under the normal precision policy. The sliver
removes the imaginary 1/R growth on the two studied A-star H/Z rays; the real
integrated contribution still grows as 1/R. The Gaussian reduces that imaginary
coefficient to about 0.345 of baseline without changing its power. Both changes
create larger imaginary weights elsewhere in this pilot.

An apparent sliver improvement in the A-star mixture pilot came from replacing
the worst ordinary draw, rather than improving its density. At that same point
the ordinary Im weight is −0.83078, whereas the A-star mixtures would give
approximately −1.107. This is why that small pilot does not establish a global
gain. No sample was removed from an estimate because of its value.

Selected Arb comparisons for the original three profiles preserve an accuracy
caveat: all twelve rays and two
older maxima are finite but fail the stricter 1e−12 rotation tolerance. Ray
normal/Arb values differ by at most 4.46e−7 relatively, with matching 1/R
exponents; the two newer b=0.3 maxima pass Arb. This supports the measured
power without certifying the requested Arb accuracy. Current GL638 cut orders
are one. Differentiating a moving hard boundary for raised cuts requires a
separate investigation; the simple-cut result does not settle that question.

The hard boundary has a separate, demonstrated problem even for simple cuts.
Twelve full-graph points approach a soft q13 configuration while the common
limiting A radial displacement approaches the 100 GeV sliver boundary. Set the
s-ray offset from its limiting A root to `100 GeV + nu*lambda` and let lambda
range from 0.1 to 0.0001 GeV. Individual finite-lambda cut displacements differ
by O(lambda); their CT support has not been traced directly. Individual cuts
0/1/3/4 scale as lambda^-3. The baseline cancels these
terms to approximately lambda^-1 for all three tested values of nu. With the
hard sliver, the Im sum remains lambda^-1 at nu=-1 but scales as lambda^-3 at
nu=0,+1; forced-Arb tail exponents are -2.99997 and -3.00003. The same-point
sliver-minus-baseline difference confirms this is a localization effect.

If the surviving term persists across a transverse strip of width O(lambda),
the measure lambda² d lambda du still gives a finite absolute first moment.
However, even a soft proposal q proportional to lambda^-2 then has a
logarithmic second moment. This is a conditional local argument, not a global
bound. Normal precision also gives misleading Re powers on the straddling
rays: the configured norm-based stability test does not constrain a much
smaller component separately. At lambda=0.0001 GeV, nu=+1, the accepted Double
Re is +5.17448e-27 versus Arb -1.33310e-30, while Im is about 1e-19. Arb
restores the approximately lambda^-1 Re behavior. The existing runtime option
`stability.check_on_norm=false` instead checks the components separately;
the card exposes it as `-D stability_check_on_norm=false`.

The new opt-in runtime field `smooth_sliver=true` multiplies the Gaussian by
`exp[-delta²/(W²-delta²)]` inside `|delta|<W` and returns zero at and outside
the boundary, where `delta=rho-r_star` and `W=L*Q` (or `L*r_star` for dynamic
width). It remains even, equals one at the pole, and is flat to all orders at
the boundary. Scalar and dual paths implement the same profile. Four new tests
cover settings serialization, pole/symmetry, dynamic-width derivatives and
boundary behavior, and the generated helper's signed radial PV with its
negative-radius mirror. Nine focused tests, workspace check and clippy pass.
The option, tests and schemas are published in commit `163687912`.

The rebuilt release has SHA256
`5309304eceafcd6b27fcca2a5087ae983d96c16a08fef4399ada3bb3e6bff84b`.
Fresh baseline Arb values at the twelve boundary points are identical to the
earlier binary. With the smooth cutoff, the straddling families return to
approximately lambda^-1 in Im (tail exponents -0.99977 and -1.00000). Four
additional soft points approach a 50 GeV displacement, where the limiting
envelope is 0.7147422: their total Re/Im tail exponents are -0.99436/-0.98509,
while individual cuts remain approximately lambda^-3. All sixteen values are
finite; all forced-Arb results miss the requested rotation tolerance. The last
interior Im point also drifts by about 7% in lambda*f, so these are qualified
finite-ray fits rather than certified asymptotic bounds.

The smooth profile consequently removes the demonstrated boundary-strip
problem, but it worsens the pilot variance. It reduces the old hard-sliver
maximum's Im weight from -0.83078 to +0.07756 under the same PDF and creates a
larger weight elsewhere. It remains an experimental opt-in control; the
default profile is unchanged.

The new maximum is a hard configuration: q13=59.7–81.6 GeV and
q14=182–250 GeV across cuts. Its ordinary-PDF Im weight changes from -0.09218
to -1.02700, chiefly through cut 1. It agrees with Arb within 5.50e-14 in
relative complex difference, so numerical error does not explain this
deterioration. Of twenty-one selected Arb replays, five are Stable and sixteen
remain finite Unstable; all values are retained. The A-star rays still have
real 1/R growth, while their imaginary part tends to a finite value.

The card exposes `ct_gaussian_width`, `ct_sliver_width`, `ct_smooth_sliver`, `ct_dynamic_width`,
`radial_scale`, `integration_cores`, and `training_phase`. For a fresh matched
300,000-sample run after generation:

```text
gammaloop examples/cli/epem_a_ttxh/NNLO/epem_a_tth_NNLO_test_GL638.toml run integrate -D radial_scale=0.3 -D integration_cores=20 -D n_start=75000 -D n_max=300000 -D workspace=./GL638_b03_300k
```

Use a different workspace when comparing localization settings. The card sets
sample and worker counts; the research runner separately enforced wall time
and process-tree memory limits.
The original two 300,000-sample runs used the earlier pi/2 z probe. The
completed 20-core b=0.3 repeat above uses the corrected Euler probe and rebuilt
release; its precision flags are reported separately.

## Remaining mechanisms

For A/P, exact equal radial roots give `is_nan=true` in both normal and Arb
evaluations of the full sum; they now correctly finish Unstable with no accuracy
estimate. At eta_P=1e−6 GeV and eta_A/eta_P=1.0001, the historical z-probe
Quad evaluation passes its rotation check but its imaginary value differs from
Arb by a factor of 1120. The current generic Euler probe rejects that Quad
value and recovers the direct Arb reference, as documented in
[the runtime follow-up](runtime-followup.md). This is a threshold energy
offset, not a soft radius. Defining the exact coincident limit requires a
common divided-difference/confluent evaluation preserving
Jacobian, damping, multiplier and group factors. More qualifiers do not define
the missing removable limit; no same-amplitude iterated CT is proposed.

For a smooth reduced radial kernel `g(r)/[(r-a)(r-b)]`, the integrated residue
sum is `[g(a)-g(b)]/(a-b)`, tending to `g'(a)` when the roots coalesce. The
separate terms are singular even though their completed sum is regular.
Finite precision leaves errors proportional to machine precision divided by
the root gap; at exact equality, increasing precision cannot define either
individual term. The combined evaluator must retain the full smooth kernel
and actual denominator co-occurrence within an amplitude. A shared SOCP group
alone does not identify such a residue identity. This is an algebraic
conditioning problem; changing sampling alone cannot supply the missing value.

At the hard H/Z corner, eight Arb points with
R²=H²+(P0 Z/Q)², R=0.2 down to 0.0002 GeV, retain a nonzero R f on two rays.
The last-three-point powers are approximately −1 in both components, after the
complete cut/orientation sum; all gluons remain hard. Under local persistence
of this coefficient over an angular/tangential neighborhood, the first moment
is finite, integral R dR/R, while the smooth-proposal second moment contains
integral dR/R. Two rays provide strong numerical evidence, not a proof over
every neighborhood or hierarchy.

A dedicated normal-coordinate density q proportional to 1/R would control
this transverse degree, provided its coordinate Jacobian, inverse branches,
support and global mixture normalization are accounted for. Ordinary physical
edge-radius power maps stay smooth at this hard point. Sharpening WH alone
does not establish a cure. Further sampler experiments are documented in the
plan; A/P evaluation remains a separate issue.

More generally, in this local two-dimensional normal plane, a Cartesian
density `q proportional to R^-a` is normalizable for `a<2`; the singular
contribution to the second moment scales as `integral R^(a-1) dR` and is
finite for `a>0`. Taking `a=1` also makes the leading `f/q` weight bounded,
whereas `0<a<1` leaves unbounded but square-integrable weights. On a disk of
radius Rmax, a uniform angle and `R=Rmax*u^(1/(2-a))` give the normalized
density `(2-a)/(2*pi*Rmax^(2-a))*R^-a`. This local chart calculation assumes
regular tangential coordinates and a bounded angular coefficient; implementing
it for the physical integrand still requires the full chart Jacobian and
global mixture described above.

## Implementation fixes and reproducibility

The integrated-UV fix returns the existing typed zero for an empty selected
reduced-cograph projection, preserving cut-order support. It is also published
on the phase-fix branch as `cbcffb4aed35b19da5b626a0a706564f9ac385df`.
LU guesses skip constant radial energies; every cut's root is checked before
building kinematics. The root solver handles exhausted representable brackets
through its existing residual and local-consistency diagnostics. Twenty-three
focused tests, workspace check and clippy pass. Card settings checks pass.

Commit `11bc3107d` separately fixes the integrated threshold profile's radial
normalization in scalar and dual paths. The helper's measure leaves dr, so its
density must be `h(r/r_star)/r_star`. The old inverted argument gave the inverse
second moment of h: for power zero it was 1/sigma² rather than one. The default
power-zero, sigma-one profile used by all runs above is inversion symmetric and
mathematically unchanged. Two new tests verify normalization across profiles
and scales and analytic first/second derivatives. Both pass, as do workspace
check and clippy. This fix changes no local-CT PV envelope.

The [pure Symbolica reproducer and patch](symbolica-mre/README.md) isolates an
independent external-constant index-remapping error during evaluator merging.
It fails against the old pin and passes with the supplied patch, preserving
builtin pi. No Symjit defect was established.

The equivalent correction is now published upstream as
[Symbolica ba4958c](https://github.com/symbolica-dev/symbolica/commit/ba4958cb43b9370bfe9230b7f541df8ed47b5968).
The minimal backport
[`a277a8d`](https://github.com/alphal00p/symbolica/commit/a277a8dadbaae5a8d6605f2f381991a601e1f345) changes
only eight lines on the existing `4d0a833e` pin and requires no GammaLoop API
changes. A fresh remote fetch on 2026-09-12 confirmed that phase-fix head
`cbcffb4ae` also still pins `4d0a833e`, which lacks this correction.
Updating to upstream `ba4958c` also requires dependency and API migrations (GMP/MPFR
features, Numerica/Graphica 3.0, algebraic-module naming and LicenseManager
location, and the removal of matching's `level_range` API). The numerical runs
retain the tested minimal backport. It is published on the Symbolica fork's
`codex/gl638-evaluator-constant-merge` branch, and GammaLoop pins Symbolica,
Numerica and Graphica to that exact public commit. The dependent threshold
prefactor helper retains builtin exact pi instead of replacing it with an f64
constant. No local dependency override is required. The recorded integrations
used the identical source tree in a local checkout; no upstream API migration
has been performed. The standalone MRE deliberately retains the old pin as its
negative control.

The current card also passed a fresh generation-and-pilot smoke with release
`57ccd168…` and that clean backport. Generation took 304.23 s and verified all
936 distinct orientations, the exact saved DOT metadata, model parameters and
direct-3D plus integrated-UV settings. The separate read-only default pilot
completed three 1,000-sample iterations on five workers in 116.05 s, with
Double/Quad/Arb counts 2,943/57/0 and no NaN or final unstable samples. The
pipeline peaked at 24.14 GB; state content and timestamps remained unchanged
by the pilot, as did the source/card/DOT hashes throughout. This validates the
historical local workflow without establishing convergence or curing A/P.

Publication validation used the public Git pin above and release SHA256
`029d4ecee30c4e038889ef62725ccfcdbf6e88f5d9e7dbfea08b627d238712a3`.
Formatting, unchanged Hakari output, locked workspace/all-target check, twelve
focused tests, Clippy with warnings denied, and the release build pass. The
fetched dependency matches all 159 tested source files byte-for-byte; the lock
changes only the three Git source revisions. Fresh generation verified all 936
orientations and the exact metadata/UV/model settings in 304.86 s. Its read-only
three-by-1,000 pilot took 119.02 s, with no NaN or final unstable samples.
All thirteen fresh-state diagnostic evaluations and nineteen prior-state
controls reproduce the previous numerical payloads and precision flags exactly.
Both panels retain the two known invalid A/P cases; their finite unstable counts
are four and five. The complete smoke took 610.68 s with a 24.44 GB peak, and
source inputs and read-only states remained unchanged. This validates publication
reproducibility, including the documented failures, rather than a full GL638 cure.

The supplied private notes `ttH_defo.pdf` (especially p. 4) and
`thesis_zeno.pdf` (printed pp. 96–97, 111–115 and 135–148) motivate the sector,
factorization and radial-subtraction analysis; they are not redistributed here.
The detailed local research archive is `/common/dev/gl638_integration_validation/`.
