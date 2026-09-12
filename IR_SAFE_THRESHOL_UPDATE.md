# IR-safe threshold update and rebase plan

## GL638 implementation and research status (2026-09-12)

The direct three-dimensional local-UV route is retained, as requested, including
orientation localization before final summation. Full generation now succeeds
with integrated UV for all 936 source orientations. The current threshold
metadata is a tested research candidate, **not a complete cure**: exact A/P
coincident residues still need a stable combined evaluation, and the complete
integrand has a measured hard H/Z `1/R` corner that threatens logarithmic
sampling variance. These are distinct from the fixed UV and LU defects.

The publication checkpoint pins the public Symbolica backport `a277a8d` and
records the investigation in [the GL638 report](docs/research/gl638/README.md).
The public-pin release passes the build/check/Clippy gates and twelve focused
tests. Fresh full generation, a 3,000-sample pilot, and 32 diagnostic evaluations
across two panels pass their reproducibility checks; the panels retain the
known invalid and unstable cases. Follow-on work should address stable
coincident-residue evaluation and H/Z sampling, with the current counterexamples
as controls. This closes the current investigation/publication stage, not the
unresolved physics and numerical questions.

The accepted work is to fix the direct-3D integrated-UV boundary in both
branches, provide the pure Symbolica merge MRE/patch, fix the massive-at-rest LU
`0/0`, validate representative and complete orientations, study maxima and signed
and absolute integrals, compare release sampling on five CPUs within about
15 minutes and ten iterations, and update the GL638 DOT/card with actual
results and limitations. No same-amplitude iterated CT is proposed. Automatic
construction of threshold metadata remains deferred.

Detailed evidence is under `/common/dev/gl638_integration_validation/`:

- `GL638_INTEGRATION_REPORT.md`: current implementation, runs and limitations.
- `GL638_CURRENT_PHYSICS_ASSESSMENT.md`: exact cut/threshold inventory, weights,
  soft powers, measures and links to the supplied PDFs and earlier research.
- `full-explicit/hard-hz-check/RESULTS.md`: current full-sum H/Z scans and CT
  decomposition; `HZ_HARD_CORNER_VARIANCE.md` gives the proposed sampler.
- `full-explicit/ap-check/RESULTS.md` and `AP_COMMON_ROOT_DECISION.md`: current
  full-sum coincident/nearby-root failures and the required confluent operation.
- `symbolica-merge-mre-and-patch.zip`: portable upstream reproducer and patch.

The superseded chronological checkpoints from this section are preserved in
`PLAN_CHECKPOINT_HISTORY.md` in that directory. Earlier prescriptions and raw
soft studies remain in `/common/dev/gl638_threshold_research/`.

### Current metadata and its scope

Only physical cut 1 `(2,6,10)` has explicit variants, all with full
`parent_lmb=[3,6,7,10]`:

| Threshold | Shared `[7]` variant | Native `[3,7]` variant |
|---|---|---|
| A `(7,8)`, U `(8,12,14)` | `1-W_H` | `W_H` |
| F `(8,10,13,14)` | `1-W_F` | `W_F` |
| V `(3,7,14)`, P `(3,12)`, Z `(3,10,13)` | none | 1 |

At each variant's own cut-preserving star point,

```
H = eta(2,4,12), P = eta(3,12), Z = eta(3,10,13), Q = 1000 GeV
G0 = eta(2,6,12,13), G4 = eta(2,4,10,13)
W_H = H^2 / (H^2 + (P*Z/Q)^2)
W_F = G0^2*G4^2 / (G0^2*G4^2 + P^2*Z^2)
```

The exact common zero in the first denominator has the supplied complementary
lazy conditional (`W_H=0`, `1-W_H=1`), not an epsilon regulator or a claim of
continuity. Other cuts retain their maximal native one-loop defaults. Thus
left cuts 0/1/3/4 share cycle C7, right cuts 3/5 share C5, and the native two-loop
variants form their own group. No extra group IDs are assigned. The DOT has
nine explicit variants; the complete resolved catalogue has nineteen.

Every cut's LU root is solved and checked before any cut kinematics are built.
The common SOCP uses each cut's independently prepared external data, including
its own t-star rescaling. Groups independently own overlaps, centers and
e-surface multichanneling; incompatible signed-cycle subspaces are never merged.
LMB importance sampling is a separate operation and does not alter those groups.

Historical representative scans restore generic cut-summed soft powers near
`lambda^-1` for either gluon and `lambda^-2` when both are soft. Their measures
are `lambda^2 d lambda` and `lambda^5 d lambda`. Those scans used the physical
numerator and local UV, but disabled integrated UV before its generator fix;
they are not retroactively certified as complete full-UV tests. The current
full-UV generation, pilots and hard-corner tests are separately recorded.
Shrinking F surfaces must remain present: generation and runtime existence
tolerances are `1e-14`; the earlier `1e-7` setting could falsely improve an IR
power by dropping the constraint.

### Implemented fixes and validation

The absent projected-UV branch returns the existing typed zero while preserving
cut-order keys. It is published on the phase-fix branch as
`cbcffb4aed35b19da5b626a0a706564f9ac385df`; PR #103 remains draft. Its unchanged
physical GL638 regression passes, as does the focused unit test. The same fix
is present in this working tree.

The pure Symbolica MRE fails unpatched and passes with the owning constant-index
merge correction, for factorized and expanded helper expressions. The main
workspace pins the published patched revision; threshold helpers retain
builtin pi at active precision. No Symjit defect was established. The archive
contains no license key, GammaLoop dependency or application initialization.
The published fix is commit `a277a8dadbaae5a8d6605f2f381991a601e1f345`, directly
on the existing `4d0a833e` dependency pin: only the copied destination constant
descriptor is remapped, preserving incoming indices while they are scanned.
It requires no GammaLoop API changes. A fresh remote fetch on 2026-09-12
confirmed that phase-fix head `cbcffb4ae` still pins `4d0a833e`; it does not
yet include the equivalent upstream Symbolica fix `ba4958cb`.

LU radius guesses now skip constant radial energies, handle absent radial
growth cleanly, and use existing safeguarded root diagnostics. Root accuracy is
relative to the dimensionless LU root, without mixing it with an energy scale.
An exhausted representable bracket can pass only with valid finite values,
positive root/derivative, sign/slope probes, and residual bounded by
`epsilon*tolerance*Q + |f'|*bracket_width`. This is a numerical consistency
criterion under continuity, not a proof for arbitrary discontinuous callbacks.
All 23 focused tests, workspace check and clippy pass. Normal/Arb replays of the
previously failing saved maxima now succeed.

The initial validated release uses opt-level 3 and a process-only thin-LTO override,
SHA256 `ce79cac84b5dca80413441fa07fbbb4ab3d3ec60bdaba9dcb73624d6c48a4c5a`.
The complete explicit-sum state generated/exported in 335.41 s, with an 8.32 GB
peak and 2.20 GB saved state. All 936 source keys, 19 active variants and the
component registry match the selector-based reference. At one saved hard
maximum, all six cuts and 50 CT occurrences agree; meaningful component
differences are at most 1.53e-14 relative. Evaluation takes 0.067 s versus 0.810 s
there. This is a point comparison, not an exhaustive symbolic equivalence proof.

The original selector-sum state exceeded 30 GB during five-worker cloning before
completing a calibration iteration. Explicit summation fits: its 25-sample
calibration peaked at 24.11 GB. Each production trial pins the whole process tree
to five logical CPUs and limits auxiliary threads. Earlier debug-pilot times
were not pinned and are not five-CPU benchmarks.

### Complete-graph integrations

All runs retain full UV and physical factorized SM numerators. Training uses
the real component; monitors record both signed components and separately
`integral |Re I|` and `integral |Im I|` in pb. The latter are not the complex
modulus integral. Only completed iterations enter the table; the runner's
shutdown allowance and any operational interruption are recorded separately.

| Sampling | Completed N / iterations | Wall time | Nominal Re [pb] | Nominal Im [pb] | Finite unstable / NaN |
|---|---|---|---|---|---|
| Linear, OSE | 50,000 / 10 | 693.24 s | `(-1.826 +/- 1.135)e-4` | `(0.903 +/- 2.274)e-4` | 2 / 0 |
| Linear, inverse Jacobian | 52,500 / 7 | stopped at 788.21 s | `(1.379 +/- 1.158)e-4` | `(-0.401 +/- 1.773)e-4` | 0 / 0 |
| Power 2, inverse Jacobian | 60,000 / 8 | 904.93 s | `(-2.828 +/- 1.730)e-4` | `(-3.019 +/- 3.190)e-4` | 19 / 0 |

The inverse-Jacobian run stopped because a watchdog process query timed out;
no memory-limit or integrand failure was reported. Its incomplete eighth
iteration is excluded. Absolute relative errors are 14.65%/17.00% for OSE and
13.88%/13.82% for inverse Jacobian. At matched 45k, inverse Jacobian's real
error is 8% larger and imaginary error 18% smaller: no uniform improvement is
established. The power-2 trial stopped cleanly at the 900 s budget, excluding its
unfinished ninth iteration. Its absolute relative errors are 18.40%/19.75%; at
matched 52,500 samples its signed errors are 1.59/2.03 times larger than linear
inverse Jacobian, with fourteen finite unstable samples versus zero.

Linear inverse Jacobian is the best balanced tested nominal setting. Its
absolute integrals are `(8.340 +/- 1.158)e-4` and `(1.282 +/- 0.177)e-3 pb`,
with estimated estimator variances `1.3399e-8` and `3.1397e-8 pb^2`. These
finite runs do not establish an optimum or certified convergence; the hard
H/Z analysis below independently limits what their error estimates mean.

Finite unstable values are included in both the estimator and grid training;
they are not discarded. Their coordinates are not retained in the aggregate
stability records. Thus the table's statistical errors do not bound those
numerical errors. The largest OSE maximum replays stably in Double and Arb
within 1.1e-12 relative, but cannot be identified as either flagged sample.
Both preferred-run absolute extrema were also replayed: normal/Arb phase
differences are at most 1.35e-12. One forced-Arb replay narrowly fails its
rotation tolerance (1.06554e-12 versus 1e-12); that flag is retained. The other
is stable in both modes, with cut 1's native two-loop U CT dominating its
imaginary part. Historical adaptive PDFs are not reconstructed by these checks.

### Remaining failures and supported next steps

1. **A/P evaluation:** exact equal roots at `eta_P=1 GeV` give `is_nan=true`
   even in the complete orientation sum. Returned zeros are sanitization, not
   physical values. At nearby `eta_P=1e-6 GeV`, `eta_A/eta_P=1.0001`, Quad passes
   its rotation check while its imaginary value is 1120 times the Arb result.
   Here eta_P is a threshold offset, not a soft-gluon radius. The remedy is a
   common, conditioning-aware divided-difference/confluent residue evaluation,
   preserving actual denominator co-occurrence, Jacobians, damping, multiplier
   jets and group weights. More threshold qualifiers or changed centers alone
   do not define the missing limit. No same-amplitude iterated CT is needed.
2. **Hard H/Z variance:** eight full-sum Arb points over three decades give
   nonzero `R*Re I` and `R*Im I` on both tested rays, where
   `R^2=H^2+(P0*Z/Q)^2`. All gluons remain hard; normal controls agree. Cuts 1/3
   carry the growth, including native Z and shared A/U CTs and the existing
   opposite-side CT products. Under local persistence over an angular and
   tangential neighborhood, the first moment is finite but the second moment
   is logarithmic with a smooth proposal. This is strong numerical evidence,
   not an all-direction analytic proof. Ordinary edge-radius power tuning
   remains smooth here. A dedicated `q(u,v)~1/R` normal-coordinate channel,
   including its exact Jacobian and inverse branches, is a concrete proposed
   remedy; it is not implemented. Raising W_H powers or inserting an epsilon
   is not an established cure. A more elaborate coherent metadata construction
   has not been ruled out or demonstrated.

The tracked GL638 DOT retains the candidate, and its card exposes shared
generation/runtime blocks and named controls. Full explicit summation defaults
on; representative orientation aliases force it off for diagnostic access.
The card defaults to linear mapping with inverse-Jacobian LMB weights, five
workers, 7,500 samples per iteration and at most ten iterations; its pilot uses
1,000 samples for three iterations. The earlier seven settings checks and the
final edited-card checks pass. Both unresolved limitations are documented in
the card. The bounded research runs additionally enforced five-CPU affinity,
a 900 s budget and a 30 GB tree-memory cap; the card alone sets worker count,
not a wall-time or whole-process CPU-affinity limit.

The UV/LU fixes and agent guidance are committed and pushed on this branch as
`b04fd885e`. The independently verified integrated-threshold profile correction
is also pushed as `11bc3107d`: the normalized radial density is
`h(r/r_star)/r_star`, in both scalar and dual paths. The previous inverse
argument changed the integrated residue under nondefault profile settings.
All completed GL638 runs use the mathematically unchanged inversion-symmetric
default. New normalization and derivative tests, workspace check and clippy
pass. The phase-fix UV correction is also committed and pushed. The
remaining GL638 changes now pin the published Symbolica backport
`a277a8dadbaae5a8d6605f2f381991a601e1f345` on the existing fork's
`codex/gl638-evaluator-constant-merge` branch. All three companion crates use
this public revision; the local override used during the recorded research is
no longer required. The independent MRE, patch and standalone Cargo
setup are committed and pushed as `c86090ac7` under
`docs/research/gl638/symbolica-mre/`; the broader research summary is prepared
alongside it. Its old dependency pin remains an intentional negative control.
No Symjit patch is required by the evidence collected here.

### Follow-on sampling experiments

The active goal now covers completing the GL638 validation and improvement:
publish the tested changes reproducibly, investigate maximum weights and
integrable singularities with sampling and local-CT localization controls,
fix demonstrated implementation defects, and report quantitative results and
any remaining obstacles. New controls, if needed, belong in runtime integrand
TOML. All studies preserve direct-3D local UV, the full orientation catalogue,
common-subspace/common-center dual cancellation, the local envelope's zero
signed-radial principal-value contribution, and the normalized integrated
residue. Automatic metadata construction remains deferred.

The authorized study uses twenty physical CPU cores and longer runs.
The explicit-state calibration completed 1,000 samples without unstable or
nonfinite evaluations in 110.28 s, with a 69.98 GB peak. The host has ample
available memory and no cgroup limit; a scratch copy of the existing watchdog
uses a 100 GB guard for these experiments. The production watchdog is unchanged.

The b=1 baseline completed 300,000 samples in 2,715.50 seconds with 70.08 GB
peak memory. Absolute Re/Im estimates were (0.83692 ± 0.07330)e−3 pb and
(1.43250 ± 0.15476)e−3 pb; maximum impacts were 6.29%/5.22%, with seven finite
unstable evaluations and no NaNs. The matched b=0.3 control also completed
300,000 samples in 2,715.15 seconds: its nominal errors are about 30% smaller
in both phases, maximum impacts fall to 3.55%/3.51%, and twelve finite unstable
evaluations remain, with no NaNs. This improves finite-run efficiency without
changing the hard-corner power. The portable research summary records the
iteration data and maximum replays.

The localization study uses the existing
`subtraction.local_ct_settings.uv_localisation` controls first. Its Gaussian
times hard sliver is even in the signed radial displacement from r-star;
both radial mirror terms are required for the PV-zero identity. Compare a
narrower Gaussian with a narrower hard sliver at fixed Q, retaining unit value
at the pole and the existing integrated-CT normalization. This is separate
from the global LU h function and from the integrated-CT sampling profile.
The current GL638 cuts have simple order one. The derivative implementation's
handling of moving hard boundaries for raised cuts needs a separate audit;
do not infer its validity from these simple-cut experiments.

The paired 1,024-draw localization pilot finds no global improvement from
either g=0.1 or L=0.1: their ordinary-proposal Im variances are 3.09 and 4.29
times the baseline. The narrow hard sliver removes the local Im H/Z cone on
the studied rays but leaves the integrated Re cone. A separate twelve-point
soft/boundary experiment confirms a new boundary strip: baseline cut sums
scale approximately lambda^-1, while a hard sliver cutting between the
nearby cut projections leaves Im lambda^-3 on two straddling families.
Forced Arb confirms those powers. An opt-in `smooth_sliver` runtime setting
now retains evenness, pole value one, both radial mirrors, and flat boundary
derivatives. Nine focused tests, workspace check and clippy pass. Full-graph
controls restore approximately lambda^-1 on the boundary families and retain
cancellation inside the support where damping is nonzero, although the strict
Arb rotation criteria remain unmet. The shared-point pilot finds higher
ordinary-proposal Re/Im variances, 1.143/6.157 times baseline; the new option
therefore fixes the demonstrated boundary issue without establishing a global
sampling improvement. The default remains unchanged. The updated research
summary records the rebuilt binary and the direct baseline comparison.
The setting, tests and schemas are committed and pushed as `163687912`.

A full-UV repeat of the three ordinary soft families also exposed a separate
precision-policy issue: the card's exact pi/2 z rotation can preserve a large
cut-cancellation error and accept a wrong Double sum. This occurs even with
threshold CTs disabled. Generic Euler controls catch all four selected failures,
escalate to Quad and agree with Arb within 4.02e-12; phasewise checking alone
with the same permutation misses three. The card now restores the existing
global Euler(0.1,0.2,0.3) default with named angle controls. Eight edited-card
checks pass. All 48 Euler ordinary-soft replays are finite and Stable, with
worst complex Arb agreement 1.78e-6 and fitted-power agreement within 8.76e-7.
A controlled 20-core b=0.3 integration with the corrected probe completed
300,000 samples in four iterations, 2,372.36 s, with a 69.89 GB peak. Absolute
Re/Im estimates are (0.74662 ± 0.05147)e−3 and (1.33760 ± 0.10871)e−3 pb,
within 4.86e−11 pb of the earlier b=0.3 central values. Nominal errors and
maxima are effectively unchanged, while finite unstable samples rise from
12 to 23, with no NaNs. Double/Quad/Arb counts are 297,973/2,004/23.
The old precision flags alone do not certify accuracy, nor does this repeat
resolve the retained unstable values. All results and precision records are
retained; the tracked research report includes the completed iteration CSV.
Both final maxima also pass fresh Double/Arb replays with complex differences
below 3.21e-13. Thirteen portable raw-point comparisons reproduce the recorded
results, including the unresolved A/P NaNs and finite unstable boundary values.
An additional exact-z rotation comparison isolates the persistent pilot-point
discrepancy to threshold CTs: every CT-off cut agrees exactly, whereas cut 3
dominates the CT-on difference. The completed debug trace matches all release
primary values and exposes a false-negative full-overlap decision: identity C7
has seven centers/subsets, but the rotated probe finds one center for all
seven surfaces with at least 70.9 GeV interior margin. Catalogues and fixed
complements agree under rotation. Exact SOCP capture and standalone replay now
identify the cause: identity returns `NumericalError`, while z returns
`InsufficientProgress`; both candidate centers independently certify all seven
surfaces and agree within 2.62e-10 GeV after inverse rotation. The former status
was incorrectly treated as absence of overlap. Both overlap owners now validate
the physical candidate before interpreting status, and propagate unresolved
solves as errors rather than deleting intersections. The unchanged strict
interior test uses each instance's own kinematics, masses and native complement.
A permanent seven-instance GL638 regression checks identity, z and Euler
rotations. All 21 focused overlap/localization tests, workspace check and Clippy
pass. The fix is committed and pushed as `cb4e82ddf`. The rebuilt release
completed all 71 full-state replays: explicit z/Euler complex discrepancies
fall to 3.01e-12/3.34e-12, with a largest per-cut discrepancy of 3.62e-12.
Normal evaluation agrees with Arb within 1.08e-11; forced Arb still retains
finite instability against its strict 1e-12 target. Both maxima, all 48
ordinary-soft controls and all 11 finite portable controls are bit-identical
to frozen references. The two exact A/P invalid cases remain. The matched
300,000-sample, 20-core b=0.3 Euler repeat completed in 2,320.02 s with 70.08 GB
peak RSS. Finite unstable samples fell from 23 to 16, with zero NaNs and
Double/Quad/Arb counts 297,989/1,995/16. Signed shifts are -2.99e-11/+1.30e-10 pb;
absolute shifts are +2.89e-11/+6.27e-11 pb. Nominal errors are unchanged within
2.49e-10 relative, and all final signed extrema have identical weights and
coordinates. This supports the repair without establishing another variance
improvement or certifying the remaining unstable samples.

The subsequent stability-reporting fix is pushed as `569f03097`. Both existing
stability checks now reject nonfinite probes at every precision before the
single-sample shortcut. The primary value and existing sanitization remain
unchanged. Four focused tests under the curated repository profile, workspace
check and Clippy pass. A broader selection also encountered the unchanged,
already quarantined zero-duration renderer fixture; no expectations changed.
All nineteen full-state replays pass: seventeen finite controls retain identical
numerical payloads and precision flags, while both exact A/P failures now end
Arb Unstable with absent accuracy instead of spurious Stable/zero accuracy.
Their raw invalid weights and sanitized totals persist. This fixes reporting,
without supplying the missing coincident-residue limit or changing the nearby
A/P finite false-acceptance counterexample.

The current GL638 card then passed a fresh end-to-end smoke with the same
release and clean local Symbolica backport. Generation completed in 304.23 s,
with all 936 distinct orientations, exact saved DOT metadata, intended model
parameters, direct-3D local UV and integrated UV. The separate read-only default
pilot completed three 1,000-sample iterations on five workers in 116.05 s:
Double/Quad/Arb counts were 2,943/57/0, with no NaN or final unstable samples.
The pipeline peaked at 24.14 GB; saved state content and timestamps were
unchanged by the pilot. Source, card and DOT hashes remained unchanged. This
certifies the tested local workflow, not convergence or the A/P remedy.

A matched full-UV baseline now closes the earlier orientation/UV comparison
gap. Removing only the DOT's threshold-metadata attribute and using the same
current card/binary generates all 936 identical source signatures in 285.01 s,
with identical model and generation settings. The 168 raw evaluations use the
same 24 frozen momenta and Euler normal/Arb policies. All are finite and Stable;
the largest normal/Arb complex differences are 7.25e-7 (CT-off), 7.24e-8
(default thresholds) and 1.78e-6 (candidate), so the rotation check remains a
heuristic. All 24 forced-Arb bare totals and their six cuts are exactly identical
between states, excluding a change to the underlying integrand.

| Current full-UV treatment | q13 soft sum | q14 soft sum | Both soft sum |
|---|---|---|---|
| CTs off | Im approximately lambda^-1 | Im approximately lambda^-1 | Im approximately lambda^-2 |
| Default threshold metadata | Re/Im approximately lambda^-3 | Re/Im approximately lambda^-1 | Re/Im approximately lambda^-3 |
| WH/WF candidate | Re/Im approximately lambda^-1 | Re/Im approximately lambda^-1 | Re approximately lambda^-1; Im approximately lambda^-2 |

With lambda^2 d lambda, the default q13 behavior gives approximately d lambda
/ lambda on both tested rays, while the candidate gives lambda d lambda.
The default last-decade q13 slopes are within 9e-5 of -3. This is evidence of
the problematic radial power; concluding a divergent full integral additionally
requires its nonzero coefficient to persist over an angular neighborhood.
For both gluons soft, the lambda^5 measure makes the measured default power
radially integrable (lambda^2 d lambda), despite its poorer cancellation.
All-four, tail-three and last-decade fits, individual cuts, sign changes and
raw flags are retained under `twenty-core/default-metadata-soft/execution/`.
CT-off total Re is an arithmetic residual. The q14 and candidate double-soft
crossovers remain explicit; these fits are not general angular or hierarchical
bounds. Replay took 279.63 s, both states and frozen inputs stayed unchanged,
and the complete pipeline peaked at 5.78 GB.

The completed controlled comparison retained linear/inverse-Jacobian sampling
and real-component training, requesting 75,000 samples per iteration, at most
ten iterations, and a 45-minute cap. Its follow-up changed only the radial
scale to b=0.3. Imaginary-component training remains an untested option, not
evidence of an improvement. At the three saved maxima,
b=0.3 increases the unadapted six-channel proposal density by factors about
2.81, 1.39 and 4.15. These are proposal-density comparisons, not reconstructed
historical adaptive weights or measured integration improvements.

A completed scratch pilot additionally tested focusing on the single H=0 surface.
At fixed original a,t, cut 1's t-star is independent of p,s. In the prepared
frame H(p)=E(a+p)+E(p)−E(a+t)−E(t) has its minimum at p=−a/2, and its positive
radial root is analytic. A normalized two-sided power map around that radius
can produce q_H proportional to |H|^(-beta), with 0<beta<1. For a bounded
1/R coefficient in two transverse normals this makes the local second moment
proportional to integral R^(beta−1) dR, which is finite. The degenerate zero-radius
contact is outside this argument. Power-map and finite-difference Jacobian
controls passed, but the 256-draw full-graph pilot established no global gain.
Floating mapping failures must be recorded without silently resampling or clipping.

The focused proposal was mixed with the existing soft-aligned channels,
retaining global support and exact normalization. Raw full-graph momentum
evaluations bypass the old LMB partition; estimates use their value divided by
the full mixture density. Follow-up native-A-star proposals compared rational
and integrated-CT radial laws at 4,096 draws per proposal. Their apparent
gains remain dominated by extreme draws and need independent confirmation;
the integrated radial law worsened Im variance. The tracked research report
records forward/inverse and Jacobian checks, saved maxima, H/Z rays, and Arb
replays. All direct-3D UV and threshold group invariants remain unchanged.

## Goal

The sections below preserve the earlier accepted grouping/rebase plan and its
historical checkpoints. The current GL638 status above supersedes older pending
decisions, generation failures and validation limits where explicitly updated.

Rebase `ir_safe_threshold_subtraction` onto
`origin/codex/raised_energy_cff_wip_optimized_phase_fix`, retain the target
branch's raised-energy CFF, three-dimensional-representation, archive, and CLI
architecture, and correct IR-safe threshold grouping.

The runtime must solve only compatible thresholds together. For every solve
group, maximal overlap discovery, E-surface classification, SOCP construction,
center assignment, radial-root construction, and threshold multichanneling are
performed independently. Groups are then recombined only when forming the final
counterterm value.

## Group identity and directives

For a fully resolved threshold subspace, define its canonical solve-space
signature from the selected defining edges and their normalized oriented
fundamental-cycle coefficients over the full graph:

```
solve_signature = sorted((defining_edge, sorted((edge, signed_coefficient), ...)), ...)
```

Each cycle is normalized along its defining edge. Parent identity, local loop
slots, and parent ordering do not enter this key. Different parents may give
identical keys; the same defining edge can give different cycles. Comparisons
retain every signed cycle separately, not just the union of active edges.

Every explicitly supplied CT variant must include the complete `parent_lmb`.
`subspace` selects active defining edges in that parent. Existing fields suffice;
cycles are derived rather than supplied redundantly. Missing metadata continues
to synthesize maximal cut-side defaults with internally resolved parents.

For grouped amplitudes, all supplied member metadata is interpreted against the master graph
and its topology/LMBs. The user is responsible for aligned edge IDs and meaningful
interpretation. Member residue evaluators and masses remain member-specific;
metadata does not trigger automatic discovery of graph correspondences or of
an IR-safe duplication/weight strategy.
Groups with only implicit defaults retain the existing shared full-space solve,
using each member's native propagator routing and masses. Their graph-local edge
IDs need not match the master's topology.
In mixed groups, implicit full-space instances retain those native equations;
supplied instances retain master-graph equations. Their shared full-space center
uses the exact master generation frame, including loop order. The presence of
an implicit neighbor must not change an explicit instance's interpretation.
Non-master graph exports preserve supplied metadata and omit autogenerated
documents, so reimport does not turn native defaults into master directives.

Cross-section graph groups retain graph-local threshold treatment. Explicit
threshold metadata on a non-master cross-section graph is rejected, as agreed;
sharing threshold groups across distinct cross-section supergraphs is deferred.
This does not restrict sharing across the different cuts of one master graph.

Automatic grouping uses the solve signature. An explicit `group_id` is an
additional subdivision label:

- equal signature and equal group ID: one group;
- equal signature and different group IDs: separate groups;
- different signatures and different group IDs: separate groups;
- different signatures and the same explicit group ID: a clean resolution error,
  never an implicit merge.

The error must be deterministic and exhaustive. It will report the graph, side,
explicit ID, threshold variant names/IDs, physical cut and cut-group IDs,
threshold/cut edge sets, requested subspace and parent LMB, and the conflicting
canonical `(defining edge, signed cycle)` signatures.

Corresponding variants combined into one degenerate higher-power residue must
agree on their group IDs and physical subspaces, as well as their remaining
subtraction directives. Reject inconsistent constituents before generation and
list every constituent's cut, threshold, variant and grouping definition.
Separate duplicated variants may still occupy different groups.

Group IDs are serialized on threshold-counterterm variants and propagated
through parsing, validation, resolution, materialization, runtime metadata,
archives, display, and Python APIs. Cross-section IDs are graph-local;
grouped-amplitude IDs use the master graph namespace. Omitted IDs use automatic grouping and never join
an explicitly tagged group. IDs start at zero and are validated consistently in
each namespace.

When no threshold metadata is present in a loaded DOT graph, preserve the
legacy default: all left thresholds of a cut use the maximal available left
loop subspace and therefore form one left group; all right thresholds use the
maximal available right loop subspace and therefore form one right group.

## Cross-section kinematics and shared groups

Cross-section evaluation must not construct a shared SOCP before the cut
kinematics exist. The graph evaluation will be split into explicit phases:

1. For every active Cutkosky cut group, solve its kinematics, including all
   required `t^\star` rescalings and derivative/residue data, and retain the
   resulting representative sample and external data.
2. Classify the threshold instances using each cut's solved kinematics and
   cut-local reversed-edge rules.
3. Group compatible instances across cut groups and side labels. Each
   group carries the per-instance solved kinematic/external data needed to
   construct its E-surface constraints, while sharing only the compatible
   solve-space coordinates and center problem.
4. Build one group-local maximal-overlap/SOCP structure and one threshold
   multichanneling structure for each group. A group with several overlap
   centers gets its own center channels and PDFs.
5. Retain one complete catalogue per solve group through center, root, and
   multichannel evaluation. Local cut/member IDs dispatch only the owned CTs;
   foreign-cut surfaces remain in every relevant multichannel complement.
6. Transform complete cut samples into a representative parent for each group.
   Copy active center coordinates by defining edge into each native frame and
   preserve its sampled complement. Radial norms use these same independent
   coordinates; unimodular parent changes add no Jacobian. Preserve the existing
   per-variant radial powers and each cut's LU Jacobians.
7. Evaluate every complement at the target active root using its own fixed data
   and requested precision. Keep target-cut LU derivative jets live; another
   cut's already-solved data are constants under the target auxiliary derivative.
   Apply internal multichannel factors before raised-threshold derivative
   extraction, while user multipliers retain their opaque semantics. Use an even
   suppression power sufficient for every raised pole in the group.
8. Merge actual iterated pairs by their commuting physical cycle displacements,
   preserving both native frames. No common parent is required across unrelated
   variants, and generation must never silently replace a selected cycle.
   This linear commutativity is necessary but insufficient for independent
   nonlinear threshold projections. The GL638 validation below exposes a coupled
   pair requiring an additional subtraction design or a clear unsupported-case
   diagnostic; the choice is pending with the user.

Left and right labels retain residue ownership. Different solve signatures
remain independent; a shared explicit group ID on incompatible signatures errors.

## Rebase procedure

- Fetch and verify the remote target and common ancestor.
- Preserve a recovery ref for the pre-rebase branch.
- Replay all existing branch commits with `git rebase --onto` from the common
  ancestor onto the target tip.
- Resolve shared files by preserving target CFF/3D-representation/archive and
  command-template changes, then restoring threshold directive, metadata,
  multiplier, and runtime behavior around them.
- Keep target-only files and dependencies intact.
- Run formatting and compilation after the core threshold replay and after the
  complete commit series.

## Implementation areas

- Add `group_id` to the graph-level threshold variant schema and all resolved
  threshold draft/metadata structures.
- Add a canonical solve-signature formatter/comparator on `SubspaceData` using
  selected defining edges and signed cycles, with diagnostics suitable for
  incompatible explicit IDs.
- Validate explicit group-ID compatibility after all parent-LMB and subspace
  resolution/rebasing is complete, before threshold generation or runtime setup.
- Partition amplitude variants before existence classification and overlap
  construction; give every group independent overlap, SOCP, centers, radial
  roots, multiplier workspace, and multichanneling.
- Enumerate complete maximal-overlap candidates. Coverage of every pair by
  previously found larger overlaps does not justify dropping a candidate;
  discard it only if one already-found overlap contains the whole candidate.
- Add cross-section group metadata/cache construction after all cut kinematics
  are solved, with per-instance external data and local-ID remapping.
- Keep evaluator/archive indexing and legacy no-directive paths stable except
  where version markers must change for new serialized metadata.

## Tests and physics validation

Add focused tests for schema round trips, automatic signatures, same-cycle /
different-parent invariance and same-edge / different-cycle rejection, explicit-ID incompatibility diagnostics, omitted IDs,
master-graph namespaces, and independent multichanneling. Test two maximal sets
`{A,B}` and `{A,C}` with foreign-cut complements: for eta_B=1 and eta_C=3,
weights are 9/10 and 1/10. Independent groups or `group_id` refinements must not
change these weights. Include nonzero centers, reordered parents, full precision,
raised/mixed derivatives, and iterated products.

Validate GL297 first, then GL638. Extend the existing GL297 and GL638 tth@NNLO fixtures and command cards. For
both graphs, verify a soft-limit cancellation in which threshold counterterms
from different cuts share one compatible solve signature and are solved in one
common group, while thresholds with different signatures remain separate.
Verify the cancellation after the per-cut `t^\star` kinematics and external
data have been prepared. Retain the existing forced-cut, full-cut,
serialization, no-directive-equivalence, raised-energy CFF, archive, and overlap
regressions.

Run `cargo fmt`, targeted `cargo check`, focused `cargo nextest` tests, broader
workspace checks, and clippy.

## Validation record, 2026-09-11

Fresh GL297 generation uses the process specification and selected orientation
from the actual NNLO CLI card, with full UV treatment and no `force_cuts`.
All four explicit overrides retain their intended native cycles. Runtime tracing
certifies two groups, each with six instances and one center at the inspected
point:

- Defining edge 5, cycle `+5,+6,+7,-13`: cuts `(2,4,9)`, `(2,11,14)`,
  `(2,4,11,12)`, using parents `(4,5,9,11)`, `(3,5,11,14)`, `(4,5,11,12)`.
- Defining edge 3, cycle `+3,+4,+9,+12`: cuts `(2,6,7)`, `(2,11,14)`,
  `(2,6,11,13)`, using parents `(3,6,7,11)`, `(3,5,11,14)`, `(3,6,11,13)`.

The 162-point f128 scan uses nine magnitudes from 0.001 to 10 GeV on both signs
of each of three rays. Here lambda is the unrescaled total soft spatial radius;
the double-soft ray keeps a fixed ratio between the two gluon momenta. Each
entry below is the fitted integrand power / power including the radial measure.
Single-soft measure is `lambda^2 d lambda`; double-soft measure is
`lambda^5 d lambda`. Powers greater than -1 are radially integrable.

| GL297 treatment | q12 soft | q13 soft | Both soft |
| --- | --- | --- | --- |
| No threshold CTs | -1 / +1 | -1 / +1 | -2 / +3 |
| Default thresholds, no metadata | -3 / -1 | -3 / -1 | -3 / +2 |
| Supplied metadata, common cycle groups | -1 / +1 | -1 / +1 | -2 / +3 |
| Leading individual cut, no threshold CTs | -3 / -1 | -3 / -1 | -6 / -1 |

All points are finite. Metadata-run exponents differ from the displayed integers
by less than 0.0006 on these rays; the no-CT results exactly reproduce the prior
control JSON. Default threshold treatment therefore leaves a logarithmic
single-soft divergence here, while common-group metadata restores the tested
radial cancellations. This is a check of these rays, not a proof over all
angular or simultaneous threshold configurations. Full-precision summed values
are fitted directly, rather than reconstructing cancellations from f64 events.

Reproduction cards, exact commands, raw JSON, fits, binary hashes and group
traces are retained in `/tmp/ir-cycle-validation/GL297/`.

Fresh GL638 generation with the source card's full UV settings fails before
state save with `direct local-3D residue branches cannot be empty`. The same
failure occurs with threshold metadata removed, isolating it from these
directives. To test threshold handling, scratch cards disable integrated UV
generation while retaining local UV and local/integrated threshold CTs. The
source card and GL638 numerical metadata are unchanged; the DOT comment now
records the ordinary-threshold failure below. All GL638 conclusions below carry
this UV limitation.

The GL638 generic scan comprises 162 finite f128 points with the same
dimensionful radial convention and range as GL297. Entries again show fitted
integrand power / power including the radial measure, rounded to integers.

| GL638 treatment | q13 soft | q14 soft | Both soft |
| --- | --- | --- | --- |
| No threshold CTs | -1 / +1 | -1 / +1 | -2 / +3 |
| Default thresholds, no metadata | -3 / -1 | -1 / +1 | -3 / +2 |
| Supplied metadata, common cycle groups | -1 / +1 | -2 / 0 | -4 / +1 |
| Leading individual cut, no threshold CTs | -3 / -1 | -1 / +1 | -4 / +1 |

The supplied metadata restores integrability on these rays but loses one power
in q14 and two powers in the double-soft sum relative to the no-CT baseline.
An additional 54-point correlated q13/two-threshold scan is finite and gives
approximately -1, -3, -1 for no CT, default CT and supplied metadata respectively.
That historical ray uses a dimensionless parameter whose physical soft radius
is `326.2192514 * abs(lambda)` GeV; it is kept separate from the generic scan.

However, the present GL638 prescription fails an ordinary-threshold test away
from soft pinches. On cut `(2,4,12)`, left variant `embedded_2l`, selected edges
`[3,7]` in parent `[3,4,7,10]`, moves energies 12 and 13 in right threshold
`(5,12,13)`, whose native selected edge is 5 in parent `[4,5,7,12]`. Its linear
defining-coordinate projections commute, but the left projection changes the
right threshold equation. The current independently constructed iterated
helper therefore uses inconsistent right-root data.

An 18-point physical-cut scan holds eta_L at 20, -20 or 1 GeV and approaches
eta_R from 1 to 0.00001 GeV. Both gluon spatial norms exceed 339 GeV. Original
plus all right CTs remains finite, while the full all-cut result scales as
`eta_R^-1` (fitted powers -0.999658, -1.000295, -0.993456 respectively).
Postprocessing isolates the entire residual pole to the weighted
`embedded_2l` / right-threshold iterated pair. Its user multiplier is exactly
one, and recording on/off agrees; grouping and detailed output do not explain
the failure. Removing this pair is only a diagnostic, not a valid cure.

Supporting this prescription requires ordered subtraction with the second
root and inverse-map/derivative helper recomputed after the first projection,
including the integrated terms. Merely finding a point on both surfaces is
insufficient. Alternatively, generation can reject pairs whose threshold
energy equations depend on the opposite selected cycles. This restriction
would concern the iterated helper, independently of valid SOCP group identity.
The user has been asked which scope to pursue; neither change is implemented
pending that decision. GL638 is therefore not certified as a successful cure.

The subsequent read-only ordered-subtraction audit narrows that decision. For
the triangular local model `eta_L=u`, `eta_R=w-c*u`, the exact right-first
operator `(1-T_L)(1-T_R)` has a local smoothness certificate with complete
raised derivative jets. Reverse order leaves an image pole: for
`F=1/(u*(w-c*u))`, the remainder is `-1/(u*w)`. The correct double term uses the
projected right equation `w`, including its damping and weight derivatives;
changing only the second root is insufficient. The companion exact-arithmetic
script passes 2,163 rational/Taylor checks. These establish the local algebra,
not the GL638 global subtraction or its integrated addbacks.

There is also a concrete GL638 root-domain obstruction. Keeping the prior
physical cut `(2,4,12)` and setting the pre-LU defining momentum q7 to
`(50,0,0)` GeV preserves that cut and its t-star. The admissible zero-center
left `[3,7]` projection sends `abs(q3)` from 197.1641 to 885.1287 GeV.
For `eta_R(y;p)=E_t(y)+E_t(p)+abs(p-y)-1000`, with `p=q3`, the exact global
minimum is `2*E_t(p)-1000`; it changes from -475.3947 to +803.7536 GeV.
Thus no right root exists anywhere after that fixed left projection, regardless
of the right center. This is a geometric calculation, not a new CLI evaluation,
and it does not rule out choosing another left projection or a broader ordered
construction. The disappearance passes through a soft degenerate configuration,
so an ordered implementation needs a domain/boundary prescription as well as
nested helpers.

The current right-then-left-then-LU CFF residue channels appear reusable for a
certified triangular construction. Helper extraction must retain outer-radial
and LU derivatives through the entire inner helper and all four local/integrated
pieces. Complete-group MC, frozen versus varying centers, absent nested surfaces,
and stage-specific opaque multiplier contexts must be specified consistently.
In particular, literal differentiation of an inner weight is not automatically
compatible with the accepted opaque-multiplier contract. The supplied PDF
motivates weighted splitting and cross-cut shared subspaces but does not define
these iterated rules. No production patch has been made for this unresolved
physics extension.

Proof and code audit: `coupled-threshold-study/ordered_projection_proof.md`,
`coupled-threshold-study/ordered_projection_checks.py`, and
`ordered-subtraction-review/CODE_FEASIBILITY.md`, under the GL638 scratch
directory below. The geometric inputs and evaluated bound are recorded in
`ordered-subtraction-review/nested-root-nonexistence.json`.

GL638 cards, raw JSON, fits and counterexample are in
`/tmp/ir-cycle-validation/GL638/`, including `soft-fits.json` and
`coupled-threshold-study/VERDICT.md`.

Code gates completed so far: workspace all-targets `cargo check`, targeted
all-targets clippy with warnings denied, CLI build, and all 17 initial focused
new or migrated regression tests. The broader suite passed 27 tests directly.
Five fixtures that create symbolic settings before initialization passed
unchanged when run serially in one process after an existing initializer test;
this includes both GL297 fixtures and the GL638 runtime roundtrip. Export/reimport
also passes unchanged with a longer runner limit (118 seconds versus the default
60). The remaining empty-catalogue forced-center failure was fixed in production
by preserving the existing forced-center transformation, without changing its
assertions. All 29 tests in the final focused rerun pass, including that fix, mixed amplitude geometry,
native default export roundtrips, and actual raised iterated evaluation with two
centers and recorded/unrecorded component parity.
The final CLI rebuild also passes. Saved-state manifest version 7 records the
changed threshold payload layout; older saved states must be regenerated.

A direct triangle/box CLI control reproduces the native-default edge-mapping
error before the amplitude correction and evaluates after it. This fixture has
negative incoming energies and therefore requires
`generation.threshold_subtraction.assume_positive_external_energies=false`;
otherwise the existing generation specialization correctly rejects the runtime
request for a trimmed surface. The successful control has two overlap centers.
Cards and before/after logs are retained under
`/tmp/ir-cycle-validation/grouped-amplitude-native-defaults/` and
`/tmp/ir-native-default-*.log`.
