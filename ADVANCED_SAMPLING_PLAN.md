# Advanced sampling channels for GammaLoop

Status: implementation in progress, 2026-09-13. This document records the accepted
design and current implementation status. It is the implementation authority
and applies to arbitrary loop order and topology, to amplitudes and
cross sections, and to both ordinary and threshold-adapted sampling.

Current implementation status: there is one `SamplingChannelCatalogue` and one
`SamplingChannelId` domain. The Symbolica selection parser, explicit parent-LMB
validation, active `subspace_lmb` metadata, exact affine LMB routing, graph/default
selection and canonical map-density bridge are implemented. Generated LMB basis
numbers remain graph-routing metadata; they never define another channel axis.

Summed and Monte-Carlo sampling now use the same canonical bridge for amplitudes
and cross sections. The `LegacyLmb` request and its LMB-only partition/prefactor
implementation have been removed. Summed sampling applies `J_c w_c` once per
channel to the graph result and event weight, with outer Jacobian one. Monte
Carlo also combines its selected map/partition factor with the physical value
at native precision, keeping only grid probability separate. Mapping is
performed before applying each stability rotation. Direct momentum input retains
its supplied raw point. The remaining default single-basis reinterpretation and
obsolete weight-setting aliases still require cleanup; they are not an
alternative multichannel enumeration.

The map kernels include exact eager/dual Symbolica Jacobians, affine LMB maps,
regular and implicit radial maps, bounded products and ordered `then` maps.
Products support disjoint, independently evaluable blocks; ordered maps derive
later centers/root contexts from earlier blocks and their inverse. The X2 source
generalizes the shared traversal to multiple explicit surface blocks; nested
mixed compositions remain rejected until their dependency semantics are covered.
Root residuals must be certified; narrow brackets alone do not suffice.
Inverse maps no longer clip recovered coordinates, and native finiteness checks
retain arbitrary precision. The implicit kernel has a normalized fallback for
a clearly positive origin residual under its increasing-root contract; that
alone does not establish global absence. Numerically uncertain origin signs
raise a precision error, while explicitly classified pinched/absent surfaces
retain their fallback. The amplitude host classifies existence separately. Proper fibers prepare an
interior center from the sampled complement; full-space maps freeze that
preparation at warmup. Two-energy rank-one fibers use their analytic convex
minimum, while general fibers can reuse the existing SOCP center solver. An
SOCP result is only a candidate: its actual native energy residual must be
strictly negative. Failed solves or ambiguous minimum signs produce typed
numerical errors rather than certifying absence. Per-draw errors enter precision
rescue. The X2 source also retains typed frozen-binding failures in the existing
native caches. Warmup requires one configured precision usable for every graph;
structural errors or failure at every configured precision invalidate the epoch.
The focused and broader X2 core/API numerical gates pass.
Soft/collinear primitives and general joint normal/star charts remain unfinished.

A standalone `phase_space(cut(...))` channel now uses the actual graph energy
equation, warmup masses, fixed external momenta and full parent-frame radial
chart. It retains the auxiliary LU scale instead of reducing the integration
dimension. Thus it solves its own directional cut radius before evaluating the
physical integrand. Cut IDs are validated metadata, not channel IDs; identical
energy-edge sets with incompatible energy shifts are rejected. X2 binds
cut/left/right compositions through the same qualified block plan and native
implicit kernels. Exact signed routing must certify that a side displacement
preserves the host cut and that each target depends only on available inputs.
The generated-graph numerical gates pass, including both sides and boosted frames.

The conditional boundary consumes actual preceding raw coordinates, solves the
host's native `t*`, and prepares physical complement data before the side map.
Each foreign inverse performs this preparation at its own supplied point.
The active physical-to-raw affine map retains both `t*` and the external-momentum
shift of the native parent. One registry holds geometry qualified by target,
host/side, parent, active edges and ordered prerequisites. The earlier unused
complete prepared-sample records have been replaced by this existing embedding
boundary; unsampled active coordinates are never passed as a complete LU sample.

The partition supports exact map densities and explicit Symbolica-compiled
`singularity_proxy` scores in complete raw coordinates. Every selected channel
must provide a valid positive score in proxy mode. No proxy is silently inferred
from a map; the selected map always supplies its actual determinant. Runtime
contexts are supplied independently to every denominator score at the same raw
point. Proxy asymptotic certification and automatic score construction remain
part of the unfinished generic machinery.

Production standalone-integrand owners have been removed; a test-only probe
remains for generic UI/integration tests. The bridge acceptance harness tests
normalization, moments, partition sums, inverses and Jacobians, including mixed
named/LMB catalogues. Physical and Gaussian-reference targets now use the same
graph traversal and LMB routing boundary. Summed channels retain their separate
mapped points and weighted moments until per-draw aggregation; statistics average
once over requested draws. Reference mass accounts for graph and visible
orientation selection, and shifted references use the original raw frame under
stability rotations. Invalid contributions fail acceptance. A single explicit
discrete selection is a partition contribution; unit normalization requires the
complete summed or correctly probability-weighted Monte Carlo estimator.

The unified-bridge milestone passed 140 isolated sampling tests and the API
physical summed/MC regression with rotations. This includes actual generated
massive/massless cut charts, independent shell/determinant checks, four-loop
Cartesian inverses and canonical Gaussian acceptance. These gates exposed and
verified corrections to the radial outer-branch `(1-split)` derivative, affine
inverse coordinates and the massless-origin radial right derivative. Formatting,
core/API library checks, core test checking and clippy passed. The following
amplitude/reference/hardening slice has passed all 146 library sampling tests,
including real kite geometry and unequal-group/orientation weighting. Both API
regressions pass: saved/reloaded shifted bubble and six-dimensional kite
normalization/moments, fully subtracted kite summed/MC agreement, and an
independent default-LMB raw-frame check. Formatting, core/API test checking and
clippy pass; warnings introduced in this slice were addressed.

The surface profile now concentrates near the threshold from both sides: its
inner branch powers the distance to the threshold, using a stable native
power-complement evaluation and the same inverse owner as the implicit map.
Nine focused gates pass, including signed-normal scaling, native tiny/subnormal
coordinates, six-dimensional determinants, Gaussian moments and physical kite/
cut charts. This verifies the map law; it does not yet establish physical
variance improvement or remove the precision boundary below.

The component, affine/composition, context, bridge, eager/dual and partition
owners now preserve native Double/Quad/Arb values. Their 122-test numerical gate
covers actual six-dimensional E-surface geometry, native cut data, composed
Jacobians and all foreign inverse densities. Tiny positive scores cannot become
zero support: unrepresentable values raise typed numerical errors. Original
invalid cube endpoints remain structural errors, whereas endpoints produced by
rounding an inverse are retryable. Quad adds mantissa precision with Double's
exponent range; extreme range tests correctly require Arb recovery.
The broader 162-test sampling suite and both API regressions also pass, including
saved-state reference acceptance and physical summed/explicit-channel equality.

The production sampler now binds native geometry and recomputes maps and foreign
densities from the original binary64 draw through the existing stability stack.
Eight focused tests pass, including real kite and cut Double-to-Quad rescue,
Arb external improvement and cache invalidation, a positive sampling-factor
underflow guard, and precise API output outside f64 range. Ordinary reporting
rejects unrepresentable contributions/events; it does not silently zero them.
The broader 169-test core gate and both saved-state/summed API regressions pass;
six integration API tests also pass for precise output, grouped events and
X-space/momentum-space event and histogram weights. Formatting, core/API checks,
the integration API check and clippy pass, with no warnings on changed lines.
A radius rounded onto the threshold
raises a typed precision error, never discards a finite band or returns zero.
For `R=3`, `beta=2`, an independent half-ULP estimate
of the combined radial-coordinate band is about `9.38e-9` at power 2 and
`4.44e-6` at power 3; arithmetic cancellation can enlarge it. Native component
tests alone do not establish production map rescue. The new physical tests
establish native reconstruction; they do not yet certify root uncertainty
against the requested signed distance and density accuracy. The new bounded
gate evaluates inverse density at the supplied raw point and checks selected
forward/inverse density consistency against the requested stability accuracy.
A deterministic approximate directional root can define a valid proposal even
when its focus is displaced from the physical surface. After these consistency
gates and independent physical kite/cut localization checks, proceed with the
bounded power-2 amplitude measurements and report their focusing accuracy.
Do not delay those numerical measurements for a formal global enclosure proof,
or describe them as such a proof. Rigorous true-surface asymptotic claims still
need a justified physical-root enclosure. Native retry now also covers the reference harness. Gaussian-body underflow
and final adaptive-grid range handling remain separate acceptance requirements.
Structural map errors still fail immediately rather than retrying as numerical
instability.
The owner-by-owner migration, root-localization versus density-accuracy criteria,
and original-source rescue gates are specified in
[SAMPLING_PRECISION_RESCUE.md](docs/research/advanced_sampling/SAMPLING_PRECISION_RESCUE.md).
The native host audit also found that the existing `EdgeMass::Evaluator` owner
computes derived mass expressions in f64, even through its generic accessor.
Direct model masses are the original input boundary; native derived-expression
masses require a later extension of that shared physical/sampling owner. Do not
describe promoted derived masses as fully native geometry.

The driver now retains the bridge in the existing `RuntimeCache`/warmup owner
and borrows its catalogue on the hot path. Settings mutation invalidates it;
warmup publishes the complete set transactionally after improved external data
are ready. Worker clones own their mutable eager score buffers, and partition
evaluation borrows those buffers without cloning them per draw. The milestone
passes all 154 library sampling tests, both saved-state/summed API regressions,
core/API test checking and clippy; no warnings remain on changed lines. The bounded frozen-X1 amplitude and GL638 pilots below supply initial physical
variance and evaluation-cost measurements; equal-time convergence studies remain
open.
Runtime cut roots and conditional complements remain
evaluator inputs, not reasons to rebuild a Symbolica program for each point.
The audited ownership and invalidation design is recorded in
[WARMUP_CACHE_DESIGN.md](docs/research/advanced_sampling/WARMUP_CACHE_DESIGN.md).

The amplitude candidate study generated a massive two-loop kite in under a
second, checked its four thresholds and same-orientation intersections, and
identified a massive planar double box as a second inexpensive topology. Explicit
full-rank amplitude surfaces now use actual catalogue equations, masses and
externals in the shared implicit kernel, with rank, frame and ambiguity checks.
The generated kite regression verifies its genuine six-dimensional C/D surfaces;
physical propagator sets remain distinct from the common output-coordinate
frame. Explicit proper subspace/complement registration now works in the master parent
LMB, with active cube axes in canonical parent order. A surface shorthand orders
the complement first; explicit dependent products are rejected in favor of
`then`. X2 resolves alternative complete parents and their requested order by
reordering a clone of the existing generated basis, retaining external routing.
Automatic amplitude discovery remains open.

Named standalone cut channels now accept `radial_profile="lu_h"` and its expert
table. They inherit actual runtime h settings, compile the fitted log-logistic
mixture CDF and dual derivatives once, and invert that full normalized mixture.
The broad component preserves the ordinary raw radial law independently of the
cut root. Each name retains its own proposal; equivalent physical cut groups
contribute their largest raised order before geometry deduplication. This
initial fit does not yet tune a derivative envelope. Physical h, CT weights and
residue differentiation remain separate from sampling. The new independent
normalization tests found and corrected existing `poly_left_right_exponential`
table entries for powers 15 and 16 in both h evaluators; other entries are intact.

The analytic and implicit radial inverses now share a density derivative
evaluated at the supplied radius. The bridge checks selected forward/inverse
density consistency with a budget derived from configured stability precision;
proxy mode needs only the selected inverse. Numerical inconsistency triggers
native rescue. This practical check does not enclose the true physical root or
certify arbitrary sub-epsilon floating-point errors.

The X1 validation milestone passes the 177 selected core tests across broad and
focused runs, both physical/saved-state API regressions, core/API and Python-
feature test checking, formatting and clippy with no warnings on changed lines.
The saved-state reference test includes two LU-h profiles on the same cut and
8192 points per fixture. The detailed evidence and reproduction commands are in
[the X1 ledger](docs/research/advanced_sampling/LU_H_MATCHED_SAMPLING.md#x1-validation-milestone).

The conditional-fiber milestone extends the existing implicit map and native
stability driver. Complement preparation happens once per forward or foreign
inverse, outside radial root iteration. Analytic minima account for fixed-energy
terms, boosts and unequal or zero masses; conservative native sign checks include
large cancelling routing inputs. General centers use the existing SOCP owner and
native recertification. Only certified absence receives a full-support fallback.
Known pinched status is supported by the kernel; ambiguous host geometry still
requires rescue. No separate map engine or channel enumeration is introduced.

Validation for this milestone passes the three focused conditional/native tests,
the broader 139-test core gate, both saved-state/physical API regressions and all
11 differential tests. Core/API and Python-feature test checking pass. The final
lint cleanup replaces two clones of Copy values with copies and uses an iterator
for the finite-difference test matrix. Formatting and clippy pass with no warnings
on changed lines.

Gaussian reference values and their raw-frame moments now share the physical
original-draw precision retry loop. Both observables must pass stability checks;
reporting rejects nonfinite or unrepresentable native results. This does not yet
solve underflow inside the Gaussian body: log-domain weighting and heavier-tailed
normalized probes remain acceptance work.

The frozen X1 all-orientation kite pilot compares 8 LMB channels, one genuine
six-dimensional surface channel and their union, after Gaussian acceptance.
The union reduced the pooled reported real-part variance by 27% and largest
real weight by 48% in three 10,000-draw runs, at about 3% more evaluation time.
The surface alone worsened both variances; this is bounded pilot evidence, not a
bounded-weight or asymptotic-convergence claim. Inputs, settings, maxima and
limitations are retained in
[the kite matrix](docs/research/advanced_sampling/AMPLITUDE_X1_KITE_MATRIX.md). The corresponding
[all-98-orientation double-box matrix](docs/research/advanced_sampling/AMPLITUDE_X1_DOUBLE_BOX_MATRIX.md)
also passes its reference gates and all twelve physical runs. There, the single
surface reduces reported real variance to 0.346 times the LMB value, while the
mixture reduces imaginary variance to 0.793 times baseline. Neither channel list
wins on both components; the reported maxima and cost limits remain part of the
evidence rather than an automatic-selection prescription. The same frozen binary
also completes the boosted-kite matrix with all 18 orientations: the mixture's
reported variance ratios are 0.736 (real) and 0.864 (imaginary), with maxima down
52% and 25%, at 3.1% more evaluation time. The three prepared graph/frame cases
now have twelve profile-specific reference checks and 36 physical pilot runs.
The follow-up uses these existing cases for proper-subspace and equal-time
comparisons; it does not add more amplitude topologies before cross-section work.

The frozen X1 GL638 pilot retained all 936 orientations and full local/integrated
UV. Four proposals each ran three paired 2,048-draw seeds on 20 cores with no
final invalid samples. LU-h reduced one large ordinary-cut excursion, but its
errors worsened in the other two seeds. Cut-only absolute estimates remain
strongly undersampled; no reliable efficiency or convergence gain is established.
[The GL638 pilot](docs/research/advanced_sampling/GL638_X1_PILOT.md) retains exact
settings, maxima, state hashes and the limitations of this comparison.

These source changes and bounded pilots do not establish a GL638 sampling improvement. Conditional
side validation, physical pinched-root classification, automatic channel construction,
the complete saved-state acceptance harness and all-orientation GL638 improvement
remain open. See the independently reviewed
[LU profile study](docs/research/advanced_sampling/LU_H_MATCHED_SAMPLING.md).

The portable research bundle is in
[`docs/research/advanced_sampling/README.md`](docs/research/advanced_sampling/README.md).
Its latest detailed companion is
[`SOPER_AND_SAMPLING_API.md`](docs/research/advanced_sampling/SOPER_AND_SAMPLING_API.md);
the bundle also contains the geometry, API, literature and GL638 notes. The
precedence order is this plan, then the latest Soper/API document, then the
genericity addendum, then the earlier report. Earlier illustrative names such
as `channel_mode` and `channel_partition` are superseded by the names below.
The old research material remains useful evidence, not an alternate contract.

## 1. Goal and execution

### Goal

Add graph-aware sampling channels that can target regular energy surfaces,
intersections, projected threshold-star images, soft/collinear loci, and
Cutkosky-cut/left/right compositions. They must compose with ordinary LMB
channels through the existing multichanneling machinery and use one unbiased
raw-coordinate estimator. The implementation must work for amplitudes as well
as cross sections, with all cut preparation and physical cut sums preserved.

Include generic radial profiles matched to the existing Cutkosky/LU localization
function. A cut-aligned channel should sample its dimensionless LU scale from a
normalized density approximating `lu_h_function`, with its actual Jacobian and
inverse density retained exactly. This extends the implementation goal to
reducing artificial LU radial weight variation as well as the GL638 threshold
tails. The proposal must account for raised-cut derivatives and reference-test
tail coverage; it must never modify the physical localization function.

Validate the generic implementation on multiloop amplitudes as well as cross
sections. Include inexpensive UV-finite scalar two-loop amplitude benchmarks
with existing and intersecting thresholds, testing both map correctness and
variance improvement without any Cutkosky rescaling. Two-loop six-photon graphs
are an optional later stress test, not a prerequisite if their generation cost
would impede the core validation. Choose the scalar topologies and kinematics
from the dedicated amplitude-candidate study and preserve their UV-finiteness
including every loop subgraph.

Keep a future light-cone/collinear-segment chart possible for IR-pinched surfaces.
The dedicated pinched-surface study defines the necessary geometry and interface
guardrails. This is architectural preparation, not a requirement to implement
new pinched sampling in the current milestone.

The intended GL638 outcome is a reproducible improvement of the remaining
regular H/Z and projected-star large-weight regions while retaining ordinary
soft coverage, local 3D UV orientation localization, integrated UV and all 936
orientations. A local bounded-weight argument is not a global guarantee: exact
coincident A/P residues, pinched/tangent strata, UV tails and numerical
conditioning must remain separately diagnosed.

Threshold subtraction metadata and its automatic construction are outside this
sampling project. Sampling channels never change threshold groups, SOCP common
centers, maximal overlaps, multipliers, or dual cancellation. A channel only
partitions and samples the complete already-subtracted physical integrand.

### Non-negotiable execution rules

* Read `CONTRIBUTING.md` before every implementation phase. Preserve comments,
  use the existing owners, and avoid compatibility shims; internal Rust/Python
  and saved-state compatibility is not a goal.
* The first tracked implementation edit is this plan copied to the requested
  implementation branch location when that phase starts. Register the explicit
  implementation goal at that point; do not register it during this planning
  phase.
* The root agent owns interfaces, invariants, integration, commits and pushes.
  Agents work in bounded slices, do not commit or push, and independently
  review one another's code, tests and numerical claims. Consult the portable
  research documents whenever a geometry or estimator choice is uncertain.
  This is a deliberately hard, cross-cutting task: the root agent should act
  primarily as an orchestrator, delegate implementation, audit, testing and
  research slices to multiple agents, and repeatedly re-read the existing
  research (especially `SOPER_AND_SAMPLING_API.md`, `REPORT.md`,
  `API_AND_GENERICITY_ADDENDUM.md`, `SOPER_ANGULAR_REVIEW.md`,
  `LU_H_MATCHED_SAMPLING.md`, `PINCHED_COLLINEAR_SAMPLING.md` and
  `AMPLITUDE_BENCHMARK_CANDIDATES.md`) before
  settling interfaces or numerical claims. Extend that research when the
  implementation exposes a case not covered by those documents.
* Preserve the untracked user directory `TMP_TO_IGNORE/`.
* Commits and push authentication must use the verified repository-local
  identity `ValentinHirschi <valentin.hirschi@gmail.com>`. The verified remote
  is `git@github.com:alphal00p/gammaloop.git`, authenticated as
  `ValentinHirschi`; use SSH `git push`, not an alternate `gh` identity.
* Make coherent milestone commits on `advanced_sampling` and push only after
  the root agent has reviewed each milestone. Do not force-push.
* Do not introduce a parallel LMB/channel enumeration while migrating. Every
  grid, evaluator, diagnostic and API consumer must resolve through the same
  `SamplingChannelCatalogue` and `SamplingChannelId`; legacy LMB fields may be
  consumed only as catalogue inputs until deleted.

The graph-evaluation boundary has one `Option<SamplingChannelId>` field.
The former `SamplingChannelEvaluation` wrapper and its `LegacyLmb` variant are
gone. No separate advanced/legacy channel fields may be reintroduced; mapping
and partitioning belong to the canonical driver outside physical evaluation.

Event metadata and observable axes now use `sampling_channel_id` and
`sampling_channel_edge_ids` (with corresponding `SamplingChannel*` enum
variants).  These values describe the canonical catalogue position and
selected map edges, including surface channels; the remaining LMB-specific
names are restricted to temporary basis/weight internals slated for the same
retirement sequence.

The discrete sample representation is likewise named `SamplingChannel`; the
old `Advanced` variant is gone, so callers cannot mistake mapped channel
samples for a second enumeration.

The discrete settings variant formerly called `DiscreteMultiChanneling` is
now `SamplingMultiChanneling`. It selects one canonical channel with the grid,
whereas the summed mode evaluates all entries of that same catalogue. This
choice of estimator must never introduce a separate channel index type.

The remaining retirement work is explicit: route default single-basis sampling
through the same compiled maps, remove its reinterpretation helper, remove
obsolete OSE/alpha settings and aliases, and replace `LmbMultiChannelingSetup`
where it still owns generic catalogue behavior. The old LMB partition helpers
are already removed. The ordinal `SamplingChannelId` remains only as the
catalogue's stable position, never an independently generated LMB index.

### Milestones

0. **Plan and research snapshot.** Land this plan and the portable research
   bundle only. No GammaLoop behavior changes.
1. **Canonical model.** Add the resolved sampling-channel data model, Symbolica
   parser/pattern validation, renamed runtime settings, graph/default selection,
   inspection and transactional errors. Ordinary sampling remains unchanged.
2. **Map kernel.** Add Symbolica eager forward/inverse residual evaluators and
   dual-number Jacobian evaluators. Implement ordinary LMB and one regular
   full-support radial surface map, with exact raw-coordinate density and
   analytic/nonconstant normalization tests.
3. **Composition.** Add validated product and ordered conditional maps, explicit
   support/fallback branches, soft/collinear primitives, and fixed partition
   estimators. Validate the inexpensive UV-finite two-loop amplitude suite,
   including genuine graph-derived surfaces and variance comparisons; add
   cross-section cut-context tests separately.
4. **Intersections and stars.** Add rank-checked joint normal charts, controlled
   branches, affine CT-star pullbacks where certified, then generic implicit
   star charts. No unsupported map is silently accepted.
5. **Automatic catalogue.** Enumerate candidates, rank/deduplicate them, run
   optional finite pilot tuning, export the frozen catalogue, and validate all
   orientation/cut contexts. Automatic mode does not generate IR metadata.
6. **GL638 validation.** Run the complete all-orientation/full-UV matrix,
   replay old extrema and H/Z/A/U-star controls, compare equal wall-time and
   independent seeds, and record limitations before changing defaults.

Each milestone must pass focused tests, `cargo fmt`, `cargo check` and the
repository's applicable lint/test gates before it is merged into the branch.

### Coordinated amplitude and cross-section execution

The milestones above describe capabilities, not a requirement to finish every
amplitude feature before working on cross sections. Following the 13 September
plan review, use the execution order below. Native precision, conditional
subspaces, interior centers, partitions and the acceptance driver are shared
work. Their gates must include real amplitude **and** real cut geometry.
Amplitude benchmarks provide inexpensive independent checks; the completion
target remains a measured GL638 improvement with the full physical calculation.

For rapid development, generate a single orientation when it materially reduces
fixture cost. Record the selected orientation explicitly and use that state for
focused map, rescue and diagnostic tests. This is a development shortcut only:
repeat end-to-end acceptance and physical comparisons on unrestricted states,
including all amplitude orientations and all 936 GL638 orientations with the
complete cut/counterterm sums. A selected-orientation result never satisfies an
all-orientation gate. Prefer the existing compatible full state when reusing it
is faster than generating a reduced one.

| Work | Concrete deliverable and acceptance | Parallel work / stopping point |
| --- | --- | --- |
| Shared numerical foundation | Finish the warmed catalogue/two-sided-profile milestone, then native component, eager, inverse-density and original-source rescue. A selected point that fails in Double must recover in Quad/Arb on both a physical kite surface and a physical cut, without changing IDs or dropping contributions. | Three bounded source slices: maps/bridge; eager/partition; physical bindings/stability. Agree types before editing. This is the immediate shared dependency, not an amplitude-only project. |
| First amplitude measurements | Saved-state kite at rest/boost and massive double box at rest: physical surfaces, independent determinants, inverses, Gaussian moments, summed/MC equality, then a bounded variance matrix. | Run alongside LU h-profile implementation once native mapping passes. Stop adding amplitude topologies after these two; six photons and extra boost families remain optional. |
| Cut radial profile | Implement the generic `lu_h` profile on the already connected standalone cut chart, with normalized broad coverage, actual forward/inverse density, simple/raised-residue checks and saved-state acceptance. | Does not wait for all side maps, joint normals, stars or automatic discovery. Run an early GL638 radial pilot; label any gain as radial optimization, not a cure of the H/Z shape singularity. |
| Shared conditional fibers | Resolve active signed cycles, sampled complements and an interior center; support deterministic Existing/Absent/Pinched branches and all foreign inverse contexts. Generalize beyond one surface child. | Validate the same owner with kite/box proper subspaces and with a cheap cross section. This unlocks amplitude products and cut-preserving side maps together. |
| Physical cut and sides | Connect actual production prepared-sampling handoff; prepare each host's native `t*`, cut momenta and external data before its dependent blocks. Deliver cut plus one side, then both sides, with dimensional/rank checks and cut preservation. | The cross-section owner works while the amplitude owner runs its bounded measurements. Shared maps and context contracts are audited across the two tracks. No complete amplitude feature catalogue is a prerequisite. |
| GL638 direct and projected targets | Resolve direct H on cut 1 even though H is not a generated cut-1 threshold CT. Add the certified A-star affine pullback using the actual native variant/overlap center. Then add a generic rank-two joint chart for H/Z in one p block and its certified star pullbacks. | Replay hard rays and old extrema as each target becomes available. One-normal channels can establish local finite variance but do not imply bounded leading maxima. U-star and rank-changing/tangent cases retain separate diagnostics. |
| Automatic selection and final comparison | Share discovery, rank/support checks, canonicalization, soft coverage and optional frozen pilot tuning between amplitudes and cross sections. Finish user-facing acceptance, runtime settings/CLI/Python migration and legacy retirement. | Automatic mode follows demonstrated manual channels; manual GL638 pilots begin earlier. Final full-UV/all-orientation experiments compare explicit and automatic selections against current optimized-LMB baselines. |

The root maintains an evidence ledger in this plan and the research notes:
implemented, tested, measured and remaining are separate statuses. Every source
slice names one owner and one independent auditor; agents do not simultaneously
edit shared interfaces. After the common native interface lands, allocate one
agent to cross-section host/profile work, one to the bounded amplitude tests and
variance runs, and one to shared conditional/joint-map implementation and audit.
Reassign the amplitude slot to GL638 replay, integration and performance analysis
once its agreed matrix is recorded. Do not keep expanding amplitude examples
while cross-section production connections remain missing.
The detailed cross-section owners, host/block syntax and delivery gates are in
[CROSS_SECTION_DELIVERY_PATH.md](docs/research/advanced_sampling/CROSS_SECTION_DELIVERY_PATH.md).
Direct-H pilots begin after one host-bound conditional block works; they do not
wait for both-side composition or the joint/star infrastructure.

Following the X1 commit `3f9b1a284`, freeze that binary for the first matched
LU-h/ordinary-cut/LMB comparison and the bounded power-2 amplitude matrix. Build
it in an isolated checkout while developing the shared conditional-fiber owner
and native reference retry in the main worktree. Keep the three responsibilities
separate: cross-section build/pilot, amplitude/reference acceptance, and shared
fiber geometry. Root integrates and audits their interfaces before the next
source milestone. Experiments record the frozen source and binary hashes;
subsequent source edits must not silently change a running comparison.

The conditional-amplitude/reference milestone is committed and pushed as
`f89addb17`. X2 now has three source owners: shared AST/block traversal, one
qualified geometry registry and affine embedding; production cut/side binding
with the serial-bubble and boosted-parent oracles; and numeric warmup-failure
caching through the existing precision retry owner, plus amplitude/interface
migration. Root runs the combined gates after these interfaces are coherent.
The bounded frozen-X1 amplitude inventory is complete; do not expand it while
conditional cross-section production and the GL638 direct-H check are pending.

The X2 implementation now passes all 141 selected core/API tests, including
saved-state reference acceptance, summed/explicit physical agreement, the
conditional nine-dimensional raised-cut fixture and frozen-geometry precision
rescue. A real LU-h profile also passes the conditional-channel inverse and
full-Jacobian checks. The shifted bridge reference moment was corrected and
tested; prior frozen-pilot evidence uses an unaffected reference owner. Core/API
and Python-feature checking pass, as do all eleven differential regression tests.
Formatting and clippy pass; new test warnings are resolved, and the existing
1624-byte map-enum size warning is unchanged from the previous milestone.
The detailed checks, tolerances and limitations
are in [the X2 ledger](docs/research/advanced_sampling/CROSS_SECTION_DELIVERY_PATH.md#x2-generated-graph-validation).
The actual GL638 direct-H map now passes reduced-state reference checks at 8192
and 32768 points, native inverse checks on eight saved rays, and an independent
full twelve-dimensional determinant check. The normalization estimates are
0.998664 and 1.010322, with broad deterministic-draw dispersions; these are
mapping checks, not a precision or physical-variance claim. The density's
measured local exponent agrees with the intended `|H|^-1/2` behavior. Inputs,
provenance and limits are recorded in
[GL638_X2_DIRECT_H.md](docs/research/advanced_sampling/GL638_X2_DIRECT_H.md).
The full-state diagnostic first stopped before reference evaluation because it
confused one exposed summed selector with the production orientation count. The
corrected optimized run verifies all 936 production keys, 32768 finite reference
draws, the twelve-dimensional determinant and native ray inverses; its moments
match the reduced-state run at the same cubes. Reference substitution bypasses
physical orientation evaluation. The subsequent
[physical pilot](docs/research/advanced_sampling/GL638_X2_PHYSICAL_PILOT.md)
evaluates all 936 orientations with full UV. All 24,576 draws are valid, and all
48 saved extrema replay exactly with their historical nested sample weights.
Native physical H/Z rays and soft controls pass the precision and unchanged-physics
gates. The power-2 direct-H proposal shows a smaller pooled imaginary error,
dominated by one power-1 excursion, with more Quad rescue and essentially
unchanged real error; this does not establish a reliable global efficiency gain.
The actual fixed-ray weights grow approximately as `R^-0.498`, reaching complex
magnitudes 230 and 265 pb at `R=0.0002 GeV`. Joint normal maps and star images
remain necessary for the intended corner treatment.

The next projected-target boundary is recorded in the
[affine-star implementation audit](docs/research/advanced_sampling/AFFINE_STAR_IMPLEMENTATION_AUDIT.md).
Certify complement-only dependence for the entire physical solve group before
sharing its actual overlap centers with sampling. For a projecting threshold
that has null directions in the sampling block, sampling and subtraction must
also share the existing root owner's complement-only projection scale; a
separate approximate root does not establish alignment with the counterterm.
Represent stable overlap-membership branches in the same canonical catalogue,
with a finite expansion cap and normalized fallback for absent branches. This
is an unimplemented X4 requirement; direct-H validation does not cover it.

The physical alpha foundation now passes its focused gates. The existing root
owner solves on the unnormalized displacement, and both raised geometry paths
derive physical radial quantities from the same alpha packet. An independent
pre-change regression exposed higher-order IFT factorial conversions; these
are corrected and covered in Double, Quad and Arb, including mixed derivatives
and a shifted native parent. Generated raised-component and conditional-cut
regressions pass. The complete sampling/physical cut, center and alpha handoff
and stable star branches remain unimplemented. The six simple GL638 cuts do
not dispatch the higher-order IFT. The subsequent
[all-936-orientation native replay](docs/research/advanced_sampling/GL638_X4_ALPHA_REPLAY.md)
passes at the saved hard H/Z point (ordinary stack and forced Arb) and
soft control (Arb), retaining all six event identities and precision choices.
The hard-point Double total changes by `6.71e-10` relative; native Arb agreement
is below `1e-288` for both points. This supports unchanged tested physics, not
an accuracy or sampling-efficiency gain. All eleven differential regressions
also pass after this projection change.

The bounded direct-H allocation is complete using frozen optimized `f2f64fb17`:
reference gates, physical H/Z and soft replays, then optimized LMB, direct-H
powers 1 and 2, and power 2 mixed with optimized LMB. Each uses three seeds and
2048 draws on 20 cores. The artifact preserves each maximum's full stored
integration `Sample`, including its historical grid weight; printed coordinates
alone are insufficient. Do not extend these runs before the missing geometry
is ready. The next sequence is shared physical CT-star preparation and
complement-only projection scale, followed by the generic rank-two H/Z chart
and its certified star pullbacks. The joint chart must certify the feasible
normal domain, rank, inverse branches and normalized absent/degenerate fallback
before physical comparisons. Its generic kernel/domain work can proceed
alongside CT-star preparation: a direct joint chart needs shared host tau but
does not depend on CT centers or alpha. Only its star pullback adds those
dependencies. Automatic channel discovery and long adaptive runs are not
prerequisites for this explicit-chart implementation.

The [joint-chart implementation slice](docs/research/advanced_sampling/TWO_NORMAL_PROPOSAL.md#generic-shared-energy-component-and-remaining-graph-binding)
now identifies its independent prerequisites: outside-support inverse results
flowing into the existing zero-score partition, eager differentiation of active
coordinates with prepared parameters held fixed, a routed shared-energy pair
matcher and a rigorously certified complement-dependent normal disk. The first
exact primitive is graph-agnostic but accepts a specified energy-equation class;
other pairs remain explicit binding errors until the generic implicit rank-two
component is implemented under the same AST and catalogue. A compact chart
requires full-support ordinary coverage and deterministic normalized fallback.

The active-Jacobian prerequisite now passes 71 native/eager/map/partition
checks and the saved-state acceptance fixture, including loaded amplitude and
conditional-cut maps, ordinary integration and exact workspace resume.
It uses the existing compiled dual program and determinant owner with
ordered active columns. Static identity-zero seed information prevents unused
prepared-parameter singular derivatives from poisoning those columns; the
`u+sqrt(m)` test at `m=0` gives active derivative one in Double, Quad and Arb.
Requested singular derivatives and structurally active singular intermediates
remain typed failures. This is not a general removable-singularity solver.

The represented-geometry joint component and compact-support migration now
pass all 161 selected core checks. One full-circle map handles the supported
shared-energy pair, with eager active Jacobians, density evaluated at the
supplied inverse point, and a directed enclosure certificate for its normal
disk. The existing partition receives certified outside-support as zero, while
requiring an explicitly selected full-support sibling. Tests cover unequal and
massless energies, native Quad/Arb inputs, both circle branches, independent
Cartesian determinants, normalized Gaussian/raw moments, and foreign support
through conditional affine/permuted maps. A separate `F=1/R` oracle converges
to the same finite analytic weight from either selected channel.
This is a component milestone: graph matching and run-card binding are still
missing. The production host must also enforce prerequisite-only compact/fallback
and radius decisions across native retries; later active-point-dependent rescue
cannot silently select another proposal. These dependencies are recorded in
the existing joint and shared-preparation audits. No GL638 bounded-weight or
variance-improvement claim follows from the fixed-context component tests.
The final nine joint/enclosure tests and two saved-state/canonical-channel API
regressions also pass; core/API checking and changed-line clippy gates are clear.

In parallel, the subtraction owner has separated representative overlap
kinematics from the raised derivative packets, using the existing sample and
group types. The foreign-cut/radial-derivative regression and generated raised
component/roundtrip fixture pass. The acceptance owner exposes the existing
reference target through the normal loaded-state integration workflow, retaining
its grids, workspaces, statistics and resume checks. Root audits these boundaries and the eventual
star selector before assigning further map changes. Neither this overlap
extraction nor a reference-only full-state run establishes star alignment or
all-orientation physical improvement.

The ordinary reference-workflow gate passes all 43 focused core/API checks,
including the reloaded bubble, two-loop kite and conditional-cut fixture,
two-worker Monte Carlo integration, exact completed-iteration resume and an
actual trained-grid sample with unequal probabilities. All 53 normal core
integration-owner checks also pass. An exploratory run additionally includes
the pre-existing `failing::target_accuracy_status_treats_zero_relative_reference_as_inactive`
ETA test and fails its unchanged expectation; the repository's curated profiles
already exclude it. Python-feature checks and all-target clippy pass, with
pre-existing warnings recorded separately. All eleven differential regressions
and six physical multi-integrand/model/batching/workspace-resume checks pass.
Building that test target also exposed two stale event channel field accesses;
they now use the existing canonical `sampling_channel_id`.


The shared fiber slice reuses `Esurface`, `SubspaceData`, the current implicit
radial kernel and the existing overlap center solver. The solver may propose a
center in binary64, but native evaluation must certify it strictly interior.
Solver failure is not an absence certificate. The current invariant-margin
classification is not sufficient for every constrained routing; certify its
applicability or use a justified lower bound. A two-varying-energy fiber provides
a generic exact-minimum case with complete existing/absent coverage and an
independent kite/direct-H oracle. More general unresolved cases remain explicit
until certified. Prepare the center and classification once per complement,
before directional root iteration. A prepared host's side is optional, and an
existing surface's radius belongs to the sampled direction, not to a universal
host record. X2 connects these owners to conditional cuts; its combined
generated-graph and subsequent GL638 gates must establish the claimed support.

Production and diagnostic map preparation must use the same dependency graph.
The former complete prepared-cut records had only test callers and could not
represent a block whose active coordinates were not yet sampled. X2 uses the
existing embedding's native preparation callback instead. Its cut data exist independently of an
optional Left/Right designation; a target can be another energy surface evaluated
on that host. Resolve geometry from graph energy equations, not exclusively from
the threshold-CT registry. The shared Symbolica model uses `at_cut` for the
explicit host and `block` for per-child active coordinates.
Current `on_cut` numbers are CutIds, not edge sets; never silently reinterpret
them. A `cut(...)` expression identifies physical edges and is resolved against
the loaded graph. No sampling selector restricts the physical cut sum.

Use at least three independent seeds for the first bounded amplitude matrix:
`auto:optimized_lmb`, `auto:lmb`, a nonduplicated full-rank surface selection,
and its mixture with optimized LMBs. Compare equal evaluations first, then equal
wall time for the baseline and the most promising valid candidate. Include
surface-only channels only when their declared full support is sufficient.
At rest the kite and box C/D pairs coincide numerically, so use one representative
for performance; retain duplicates as a separate partition correctness check.
This is an explicit benchmark selection, not automatic identity deduplication
based on accidental numerical equality. Distinct canonical targets remain.
Hold all 18 kite / 98 box orientations and the fully subtracted integrand fixed.
Report neutral or worse variance honestly and move on; amplitude success means
correct generic sampling plus measured impact, not a required gain for every
already smooth subtraction.
The independent box routing and actual same-orientation intersection witness
are recorded in
[DOUBLE_BOX_GEOMETRY_CHECK.md](docs/research/advanced_sampling/DOUBLE_BOX_GEOMETRY_CHECK.md).

Keep GL638 measurements on the critical path rather than postponing all of them
until automatic channel construction is finished:

1. Refresh the frozen baseline provenance and current-binary soft, H/Z, nearby
   A/P and maximum-weight replays while implementation proceeds. Reuse the
   compatible full state read-only; regenerate only if a real format/generated
   data change requires it. Use
   [the replay inputs](docs/research/gl638/REPLAY.md),
   [the runtime follow-up](docs/research/gl638/runtime-followup.md),
   [the direct/star geometry](docs/research/advanced_sampling/gl638.md) and
   [the two-normal derivation](docs/research/advanced_sampling/TWO_NORMAL_PROPOSAL.md).
2. At each usable channel milestone run a small paired pilot of the **complete**
   calculation: baseline versus cut-h, direct H, A-star H, joint H/Z, and their
   supported unions. First verify reference acceptance and frozen-point physics;
   then compare estimates, absolute moments, second moments, maxima and runtime.
   Representative rays diagnose local powers, not global integration improvement.
   A statistically paired difference requires common cube draws and frozen grids;
   independently adapting runs are matched-budget comparisons, not paired draws.
3. Run longer comparisons on 20 cores only after those gates, with fixed physical
   settings, separate warmup/training/production budgets, independent seeds and
   both equal-evaluation and equal-wall-time results. Record map/root/inverse and
   physical-evaluation costs, rescue/invalid counts and complete extrema inputs.
   Replay and decompose maxima to choose the next useful channel instead of
   adding channels without evidence. Wait economically while jobs run.
4. Final evidence retains six cuts, all 936 orientations, direct 3D local UV,
   integrated UV, WH/WF metadata and ordinary soft/UV coverage. Require compatible
   physical estimates and a reproducible improvement in useful error per wall
   time and/or the previously diagnosed large-weight tail; state exactly which
   improves. The separate exact A/P confluent-evaluation issue is not cured by
   sampling. An early radial/local gain is not final completion, and failed or
   unchanged comparisons remain part of the report.

The shared acceptance driver is a deliverable for both tracks: expose the
existing process reference overlay through the ordinary loaded-state integration
workflow so grids, probabilities, clones, resume and statistics are exercised.
Do not build separate amplitude and cross-section benchmark integrators.
Anisotropic/cross-loop and heavy-tail reference probes extend this same owner.
Generic joint charts and automatic discovery also require witnesses from both
tracks; graph-specific formulas may be recognized geometric specializations,
never dispatches on a name such as GL638.

## 2. Public API, defaults and migration

### Runtime names and selection

Use one channel axis for LMB and graph-aware maps:

```toml
[sampling]
sampling_multichanneling = true
sampling_channels = "monte_carlo"       # or "summed"
sampling_channel_weight = "map_density" # or explicit "singularity_proxy"
default_channel_selection = ["auto:optimized_lmb"]
channel_selection = { GL638 = ["auto:surfaces", "hard_intersection"] }
```

`lmb_multichanneling` becomes `sampling_multichanneling`, `lmb_channels`
becomes `sampling_channels`, and `lmb_channel_weight` becomes
`sampling_channel_weight`. `sampling_channels` controls summed versus Monte
Carlo execution; `sampling_channel_weight` selects the partition density.
Channel identity is selected by the separate lists.

With multichanneling disabled, ordinary sampling follows the existing route.
When it is enabled without a selection, `auto:optimized_lmb` is used. An
explicit graph entry replaces (rather than extends) the default list. Entries
within one list form an additive union. A named definition does not activate
it until selected. `auto:lmb`, `auto:optimized_lmb`, `auto:surfaces` and named
channels may therefore be combined, including redundant combinations.

`auto:lmb` means every admissible generated LMB in the graph catalogue. For an
amplitude this is the full graph spanning-forest catalogue. For a cross section
it is the full initial-state-adapted catalogue. `auto:optimized_lmb` applies
the current massless-combination/cut-side heuristic and then the elementary
soft coverage audit; it is a useful heuristic, not a global optimum.
`auto:surfaces` adds resolved surface/intersection/star/soft/angular channels
and retains enough optimized ordinary coverage for domain and soft factors.

Exact named selections mean exactly the selected channels: no ordinary LMB is
silently appended. A selected set must nevertheless have full declared domain
coverage; an uncovered domain or missing elementary soft coverage is a clear
error. A surface channel's own normalized tail and absent-fiber fallback count
as part of that channel, not hidden LMB additions.

TOML inline tables use `=` and bare keys are strings, so
`{ GL638 = ["auto:surfaces"] }` equals `{ "GL638" = ["auto:surfaces"] }`.
Python accepts the corresponding nested dictionaries. The same resolved model
is used by CLI, TOML, Python and saved settings.

### Definitions and Symbolica syntax

Definitions are graph-scoped and named by the user. `GL638` in an example is an
existing loaded graph name; `HZ` or `hard_intersection` is an arbitrary channel
name. H and Z are not reserved symbols: they are research labels for the
resolved GL638 surfaces `eta(eset(2,4,12))` and `eta(eset(3,10,13))` in the
specified prepared-cut context.

The concise form is parsed by Symbolica and structurally matched, never by a
regular-expression mini-parser:

```toml
[sampling.channel_definitions.GL638.HZ]
around = "intersect(surface(2,4,12), surface(3,10,13))"
subspace_lmb = [3]
parent_lmb = [3, 6, 7, 10]
on_cut = [1] # Current CutId filter; physical host edges are cut(2,6,10).
```

This is the proposed joint geometry in one three-dimensional p block; H and Z
are two scalar normals there. The current `on_cut` filter alone does not prepare
that physical host. Its explicit host/frame binding and the joint-map compiler
remain implementation work; the runtime must continue to reject this example
until those contracts are satisfied.

The Symbolica vocabulary includes `lmb`, `surface`, `cut`, `soft`, `collinear`,
`complement`, `product`, `intersect`, `then`, `phase_space`, `left`, `right`,
`at_cut` and `block`. Parsing a constructor does not imply that every geometry
it can describe has a compiled map; unsupported primitives remain explicit
errors until their implementation and acceptance gates are complete.
`cut` denotes a physical Cutkosky cut; `surface` denotes an energy equation.
The shorthand `cut(...) x cut(...) x complement(...)` is display sugar for a
validated `product`, and never silently means an intersection or conditional
composition. Long structured fields compile to the same canonical expression.

Every explicit channel resolves to a complete parent LMB, signed active cycles,
complement, routing/frame, energy shift, dependencies, support, profile,
branches and target context. `parent_lmb` is mandatory in canonical exported
definitions. Inferred complements are checked against it; ambiguous shifts,
orientations, host cuts, sides, duplicate constraints and invalid ranks produce
exhaustive transactional diagnostics. Exact CT-star selectors name host cut,
threshold variant, native subspace and common-center branch.

### LU localization profiles

Implement the proposal in
[`LU_H_MATCHED_SAMPLING.md`](docs/research/advanced_sampling/LU_H_MATCHED_SAMPLING.md)
as part of the generic cut-channel machinery. For the current multiplicative
LU flow and zero-centered full-parent chart, let `R_c(n)` be the physical cut
radius and sample `t = R_c(n)/r`. A simple cut, including its residue and raw
radial measure, leaves `h(t) dt` times a shape-dependent coefficient. Therefore
target a normalized `q_t(t)` close to the existing LU `h(t)`, with
`dt/du = 1/q_t(t)` and full raw determinant
`J_Omega R_c(n)^D / (t^(D+1) q_t(t))`. The full determinant is not just `1/h`.

Expose a concise named-channel option `radial_profile = "lu_h"` and an expert
inline-table form, for example
`radial_profile = { kind = "lu_h", approximation = "log_logistic", broad_fraction = 0.02 }`.
Both forms resolve to one profile owner. Inherit the actual runtime
`[h_function]` settings, including sigma/power, and invalidate the warmup
proposal when they change. Expert parameters steer only the proposal. Reject
this cut-specific profile on an amplitude-only surface or a chart without a
certified LU scale. Amplitude sampling retains its generic surface profiles.

Start with the algebraically invertible log-logistic approximation, fitted to
the mode and curvature of the log-t density `t h(t)`. Use a family-appropriate
proposal for `exponential`, whose density is nonzero at t=0. Provide normalized
positive piecewise-exponential surrogates or controlled numerical CDF inversion
behind the same profile interface when useful. Approximation quality affects
variance only: always differentiate/invert the actual proposal, never replace
its determinant by the desired `1/h`. Compile the map and density expressions
together using the existing eager/dual Symbolica evaluators.

For a cut of order m, cover the radial packets `t^j h^(j)(t)`, `0 <= j < m`.
Inspect the generated maximum order automatically and support a positive
derivative-envelope approximation with user-steerable coefficients. Proposal
Jacobians and partitions remain outside all residue derivatives. The physics
`h`, local threshold localization symmetry, and PV cancellation stay unchanged.

Require adequate broad coverage, not just positivity. Pure exact-h proposals
can give the Gaussian reference infinite variance at the raw origin; some
high-power log-logistic fits do too. Include a normalized broad component such
as `a/(a+t)^2` with `a=R_c/beta` and `beta=e_cm*b`. This induces the existing
ordinary raw-radius density `beta/(beta+r)^2`, independently of the cut root.
Alternatively, certify equivalent coverage from a selected ordinary channel
with a positive probability. A mixture must use its actual full density and
CDF/inverse (or all inverse branches), rather than an artificial denominator
floor. Users may disable the component only with explicit diagnostics about
the resulting reference-moment assumptions. Use the same deterministic ordinary
fallback with power one in forward and inverse when the cut has no regular
radial root. The first mixture implementation inverts its monotone analytic CDF
with the existing scalar solver, bracketed by its two component quantiles.
Splitting the cube interval between two full-support maps would introduce two
preimages and is not compatible with returning only one branch's determinant.

Each selected host contributes a canonical channel, while every channel still
evaluates the full physical cut/CT sum. Evaluate all foreign densities at the
same raw momentum point with their own `t_c`; never substitute the selected
host's t into another cut's density. Validate radial normalization, tails,
independent Cartesian determinants, simple-residue radial flattening, and the
raised-residue oracle `integral t^j h^(j)(t) dt = (-1)^j j!` where endpoint
terms vanish. For existing tabulated h normalizers, check the derivative-packet
identity relative to the same measured zeroth moment, and test normalization
separately at the table's accuracy. Exercise the loaded-state reference harness with broad coverage.
Then compare all-orientation GL638 absolute/signed moments, maximum weights,
root ratios, precision-rescue rate and time at equal work budgets. This radial
optimization complements the H/Z and CT-star shape channels; its one-dimensional
model improvements do not establish a GL638 improvement.

### Profiles, coverage and channel identity

Regular one-surface maps use a normalized full-support profile. A signed power
map of order two in the existing compactified radial coordinate gives the
integrable `|delta|^-1/2` density. Finite windows, energy/radial distance,
widths, tails and branch probabilities are explicit typed profile settings;
width units are never inferred. Empty or pinched fibers choose a deterministic
normalized fallback from complement data before active coordinates are consumed.
Rejecting samples or silently dropping a failed root is forbidden.

`map_density` is the default exact estimator and requires an available,
correct push-forward density. `singularity_proxy` is an explicit alternative:
nonnegative scores must be evaluated at the same raw point, with support,
branch multiplicity and asymptotic comparability certified. A proxy may replace
foreign Jacobian evaluation in the partition, but never the selected map's
true forward Jacobian or its sampling PDF. There is no silent score fallback.

Deduplicate by resolved map identity: parent routing, signed cycles, target and
star context, composition order, profiles, support and branches. Names and
common edge lists alone are insufficient. Preserve aliases/provenance and use
deterministic order and IDs in inspection.

### Migration

Migrate Rust settings, Python bindings, CLI completion, documentation, examples,
schemas and persisted state together. Remove old selector fields and competing
paths rather than retaining aliases. Migrate `lmb_basis_ids` into named
`lmb(...)` channel definitions or generated catalogue selections; numeric IDs
remain inspection/replay data only. Update all UnitVolume test callers and
settings arms while preserving their validation purpose. No saved-state
compatibility shim is required.

The former `ChannelIndex` and the current LMB-only discrete enumeration are
migration scaffolding, not a second production channel model. There must never
be two channel enumeration techniques in the final implementation. New
catalogue, map-density, grid and evaluator work must use the resolved
`SamplingChannelId` catalogue. Once that driver covers the existing sampling
modes and the summed/Monte-Carlo equivalence gate passes, remove the remaining
LMB-only index/reinterpretation path, its prefactor helpers and its settings
labels rather than expanding both systems. The canonical catalogue is the
only long-term source of channel identity.

## 3. Engine contracts

### Geometry and composition

For fixed complement `y`, a regular positive-energy surface is a convex bounded
fiber in its active routing space. Remove spectator/null directions, choose an
interior center and solve a positive radial root. The radial determinant is
`r_star^d rho^(d-1)` times the angular/complement factors. Center dependence on
earlier complement coordinates is triangular; active dependence must be included
in the full derivative.

`product` requires a direct-sum block and independent fixed-complement data.
`then` is an explicitly ordered acyclic conditional map with a triangular
determinant. Its block inverse uses the same earlier physical output blocks as
the forward condition, so it can traverse declaration order. Ordinary function
composition instead reverses its component inverses. Both forms must follow
their actual dependency graph; master-frame embeddings propagate the outer
context before appending earlier child outputs.
`intersect` uses scalar normals plus tangent coordinates, requires full rank and
controls all inverse branches and chart patches. A
successful root solve does not establish global injectivity; unsupported or
rank-changing maps fail with a capability diagnostic.

Cut/left/right channels first prepare the Cutkosky kinematics, including each
cut's own solved LU root and `t_star` rescaling. A host cut equation after its
own preparation is not a third varying threshold normal; it organizes phase
space and the auxiliary scale. Left and right threshold maps then act in their
preserving subspaces. The complete physical cut sum remains in every event.
Amplitudes use the same primitives with fixed external data and no artificial
physical-cut context.

### Future pinched-surface maps

Preserve the extension described in
[`PINCHED_COLLINEAR_SAMPLING.md`](docs/research/advanced_sampling/PINCHED_COLLINEAR_SAMPLING.md).
A massless two-line pinch is a collinear segment, with one longitudinal and
two transverse coordinates. In a spatial cylindrical chart its measure is
`d^3k = (P/2) dz d(k_perp^2) dphi`; a factor `k_perp^2` appears after using a
logarithmic transverse coordinate. Uniform sampling in that logarithm down to
zero is not normalizable. Use a normalized fractional-power or suitably
logarithmically softened profile, with explicit full-space longitudinal tails
and joint soft-endpoint coverage.

The alternative exact prolate chart has normal coordinate
`delta=(|k-A|+|B-k|)/P-1`, so the energy deficit is exactly `P delta`.
It has an explicit inverse and determinant and can enhance both the segment
interior and endpoints. It is a candidate for a later focused implementation,
not a mandatory new production map now. Neither construction replaces a
general multiloop pinch classifier: simplex interiors, massive point pinches,
soft faces and overlapping normal spaces require their own rank/support data.

Keep `Existing`, `Pinched` and `Absent` distinct. A pinched chart needs its
geometry, normal rank and frame rather than a fictitious zero-radius regular
shell. Reuse signed-cycle subspace metadata and certify that a conditional
active block leaves the host cut momenta and segment total momentum fixed.
Otherwise reconstruct the changed cut and its full pullback determinant.
Do not overcount shell shape coordinates when adding the auxiliary LU scale.

Endpoint and axis vectors must be typed frame data, or be rebuilt in the
target frame and precision. They must not be placed in the current opaque
runtime vector and silently copied as rotational scalars. Stored f64 prepared
cut records are not a high-precision geometry cache. Forward and all foreign
inverse/proxy evaluations must reconstruct matching native contexts at the
same raw point. These are requirements before enabling conditional physical
maps, even if the pinched primitive itself remains deferred.

### Symbolica evaluators and Jacobians

Regular surface profiles must enhance both signs of the intended normal
distance. For example, with `s=R/(R+beta)` and `p>1`, use
`r=R[1-(1-u/s)^p]` for `u<s`, and retain the exterior distance-power branch.
The inner derivative is `(R p/s)(1-u/s)^(p-1)` and its inverse is
`u=s[1-(1-r/R)^(1/p)]`. Both sides then have radial density proportional to
`|r-R|^(1/p-1)`. Implement stable small-coordinate evaluation, independently
check determinants and inverses, test both signed asymptotic approaches, and
repeat Gaussian/moment acceptance. Retain separate ordinary soft coverage where
needed. This concerns proposal densities and never changes the physical CT
localization function or its PV symmetry.

During warmup, build specialized eager evaluators for forward maps, inverse
residuals, density factors and support predicates. Seed independent map inputs
with first-order dual numbers and extract the output derivative matrix. Use a
stable determinant/factorization and exact block products where certified.
Implicit roots are differentiated from their defining equations (`G_z dz =
-G_u`), not by differentiating an arbitrary solver iteration. Absolute-value,
piecewise and branch boundaries receive explicit one-sided/stratum handling.

The evaluator is specialized to routing, topology, profile and expression shape;
runtime kinematics, roots, centers and widths remain inputs. Dual evaluation
does not prove normalization, support, global invertibility or optimizer-active
set stability, so those are separate contracts.

Retain compiled catalogues/programs in the existing runtime cache owner, with
explicit invalidation on channel definitions, routing, masses and relevant
kinematics/profile changes. Do not add a second cached channel enumeration.
Test that repeated samples reuse the programs and that changed settings rebuild
the relevant geometry without stale roots or masses. Record compilation/warmup
cost separately from per-point root, inverse-density and physical evaluation
cost in the amplitude and GL638 variance benchmarks.

### Estimator and adaptive grids

Every selected chart maps unit-cube coordinates into the same raw/master graph
frame. For exact map densities use the selected forward Jacobian and all branch
probabilities. For proxy partitions use scores `rho_a(K)` at the common raw point
and `w_a=rho_a/sum(rho)`. Apply each outer adaptive/discrete sampling weight
exactly once. Do not reuse the old per-cut LU inverse-Jacobian partition as a
new raw proposal density without proving the coordinate relation.

Fixed channel geometry and support are frozen before production. Existing
discrete/continuous grid adaptation may tune probabilities and bins, subject to
positive floors. Automatic pilot tuning reports candidate, seed, evaluations,
cost, errors and extrema, but never claims global optimality or bounded weights.

Automatic discovery enumerates physical surfaces, feasible intersections,
projected targets, soft/collinear loci and cut/side combinations across all
requested orientations. It removes impossible mass-gap/rank cases, canonicalizes
duplicates, caps combinatorial expansion and exports every decision. It never
creates or edits threshold-subtraction metadata.

## 4. Acceptance harness

The old standalone `UnitVolume`, `UnitSurface` and `HFunctionTest` machinery is
deprecated completely. Migrate its useful volume, angular/coarea, profile,
monitoring, status, resume and multislot validation purposes onto the new
process-level sampling acceptance path and focused map/profile tests; do not
retain a second standalone integrand owner. The new path loads a real saved
state in an isolated workspace and substitutes a known raw-frame reference
function at the physical-integrand boundary while retaining the actual map,
Jacobian, channel partition, branch decisions and outer weights.

The first ordinary-workflow interface is:

```text
integrate -p MyProcess -i MyIntegrand --reference-gaussian '{"width":300}' -w reference_check --n-cores 20 --show-phase both --restart
```

An optional `center` array supplies a shift in the generation-parent raw
spatial coordinates; omission means zero at the resolved loop dimension.
Python's `integrate(..., reference_gaussian=(300.0, []))` selects the same
centered descriptor; a nonempty second tuple item supplies the shifted center.
Each selected process slot currently requires one common nonzero loop count
across its loaded graphs, even if a graph filter selects a smaller subset.
The integration accumulators report normalization in Re and the raw
momentum-squared moment divided by its known expectation in Im, both with
target one. These are explicitly labeled acceptance observables. This lets
ordinary adaptation, Monte Carlo channels, worker clones and completed-iteration
statistics exercise the real maps. Resume repeats the same descriptor; omission
requests physics and cannot silently resume an acceptance workspace. The
workspace manifest is versioned for these new semantics; prior integration
workspaces need `--restart`, while generated states keep their binary layout.
This initial interface does not yet supply anisotropic or heavy-tail probes,
and the documented Gaussian-body underflow limit remains.

Required probes include normalized Gaussian, shifted anisotropic Gaussian with
cross-loop covariance, several widths, known first/second moments and a
full-support heavy-tail reference. Test summed and Monte Carlo channels,
unequal probabilities, trained/frozen grids, exact/proxy partitions, duplicate
channels, absent/pinched fibers, all inverse branches and nonconstant moments.

Structural failures (support, branch, rank, inverse residual, partition sum,
wrong target frame or Jacobian) fail regardless of a central estimate. Statistical
agreement reports estimate, standard error, pull, seed and evaluation budget.
Unit-volume normalization alone is insufficient: a matched wrong Jacobian can
pass a round trip and a finite run can miss a branch.

Tests cover at least one amplitude and one cross section, including a host cut
with left/right maps and correctly prepared `t_star`. They verify unchanged
threshold groups/centers, complete cut sums, event weights, summed/MC
equivalence and raw-frame invariance. All Symbolica evaluators have direct
finite-difference determinant comparisons at regular points.

## 5. Multiloop amplitude validation

Use the reproducible candidates and kinematics in
[`AMPLITUDE_BENCHMARK_CANDIDATES.md`](docs/research/advanced_sampling/AMPLITUDE_BENCHMARK_CANDIDATES.md).
The initial fixture is the five-propagator massive scalar two-loop kite, an
off-shell two-point amplitude obtained
from `tests/resources/graphs/double_triangle.dot` with every internal mass
positive. Its overall UV degree is -2 and every loop subgraph is UV finite.
At masses one and rest-frame energy five, it has both two-particle and
three-particle thresholds. A boosted configuration supplies a regular
intersection of the genuinely coupled three-particle surfaces. The massive
planar double box supplies a secondary four-point scattering topology;
two-loop six-photon production is an
optional later stress test.

Resolve the actual generated thresholds, orientations and requested parent LMB;
do not substitute a catalogue basis number for the parent frame. Test existing,
absent and intersecting surfaces, spectator complements, coupled six-dimensional
maps and runtime changes of external kinematics. No Cutkosky root or rescaling
is introduced in this amplitude path. The acceptance overlay must use the actual
loaded amplitude's compiled maps and cover normalization, nonconstant moments,
independent Jacobians, inverse recovery and summed/MC equivalence.

For physical variance comparisons hold the fully threshold-subtracted amplitude,
threshold groups/centers and all orientations fixed. Compare ordinary LMB,
optimized LMB, explicit surface-only, mixed and automatic channels using several
independent seeds, first equal evaluations and then equal wall time. Record real
and imaginary integrals, estimated variances/second moments, maximum weights,
stability/rescue counts and map costs. Establish statistical agreement before
claiming a gain. Report neutral or worse variance too: a smoothly subtracted
amplitude need not benefit from concentrating on a particular threshold.

## 6. GL638 validation and completion criteria

The first manual channels are the direct H surface in the prepared cut-1 p
block, the certified affine A-star H pullback with its native `[3,7]` context,
and the rank-two H/Z joint chart in the shared `[3]` p block. The H/Z chart uses
two scalar normals plus one tangent; it is not a product of two three-vector
subspaces. Native U-star remains a separately declared coupled chart.

Use the recorded H/Z geometry and exact-star inputs only after resolving them
against the current graph. Test both normal-power and polar profiles, angular
and tangential coordinates, support boundaries, shrinking fibers, center
changes, hierarchical approaches, soft13/soft14/double-soft limits and UV tails.
Re-evaluate the actual A and U overlap centers; never infer a star target from
an edge list alone.

For each candidate compare baseline, optimized LMB, surface union, joint H/Z,
star-target and automatic selections at equal evaluations and wall time, across
independent seeds. Record real/imaginary estimates, absolute-integral behavior,
second moments, maximum weights, unstable/NaN counts, branch visits and cost.
Replay previous extrema at normal and Arb precision. Exact A/P coincidence is
not cured by importance sampling; report it as requiring a confluent divided-
difference evaluator. Nearby distinct roots must retain the existing stability
rescue behavior.

Final GL638 acceptance requires the complete six-cut sum, all 936 orientations,
direct 3D local UV with orientation localization, integrated UV and current
threshold metadata. A local finite-variance or bounded-leading-weight result
for one regular H/Z patch is reported as such, not promoted to a global theorem.
Defaults change only after the full matrix and independent-pilot evidence show
an improvement without loss of ordinary soft/UV coverage or IR cancellation.
