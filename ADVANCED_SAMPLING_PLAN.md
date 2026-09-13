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
channel to the graph result and event weight, with outer Jacobian one; Monte
Carlo keeps its selected-map factor and separate grid probability. Mapping is
performed before applying each stability rotation. Direct momentum input retains
its supplied raw point. The remaining default single-basis reinterpretation and
obsolete weight-setting aliases still require cleanup; they are not an
alternative multichannel enumeration.

The map kernels include exact eager/dual Symbolica Jacobians, affine LMB maps,
regular and implicit radial maps, bounded products and ordered `then` maps.
A product currently supports disjoint blocks with at most one explicit surface;
ordered maps derive later centers/root contexts from earlier blocks and their
inverse. Unresolved numerical roots fail explicitly. The current implicit map
uses a normalized ordinary fallback when its supplied center is outside or on
the surface; this does not prove that the surface is absent. Interior-center
selection and true absent/pinched classification remain required for generic
surface targeting, especially with boosted external kinematics.
Soft/collinear primitives and general joint normal/star charts remain unfinished.

A standalone `phase_space(cut(...))` channel now uses the actual graph energy
equation, warmup masses, fixed external momenta and full parent-frame radial
chart. It retains the auxiliary LU scale instead of reducing the integration
dimension. Thus it solves its own directional cut radius before evaluating the
physical integrand. Cut IDs are validated metadata, not channel IDs; identical
energy-edge sets with incompatible energy shifts are rejected. Cut/left/right
compositions remain guarded until their conditional maps preserve the cut
kinematics and provide consistent forward and foreign inverse densities.

The deferred boundary retains canonical coordinates and a typed prepared cut
context. `PreparedCutSamplingContext::from_lu_sample` uses that sample's positive
`t*` and complete parent frame. `DeferredCrossSectionSamplingState` and
`SamplingChannelRuntimeContexts` carry distinct per-channel cut and evaluator
contexts. `PreparedCrossSectionMapEvaluation<T>` records a mapped sample,
Jacobians and prepared context together. These records are interfaces for the
remaining conditional maps; they do not themselves implement those maps or
replace the high-precision LU solver.

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
named/LMB catalogues. Loaded-process continuous and explicit canonical discrete
coordinate entry points exist, with saved/reloaded scalar-bubble coverage.
The summed reference overlay currently errors explicitly because its one-point
moment record cannot represent several mapped channel points. Moving reference
substitution and moment accumulation through the shared physical evaluation
boundary remains a required harness milestone. Actual physical summed-versus-MC
regressions are separate from this reference-overlay limitation.

The unified-bridge milestone passed 140 isolated sampling tests and the API
physical summed/MC regression with rotations. This includes actual generated
massive/massless cut charts, independent shell/determinant checks, four-loop
Cartesian inverses and canonical Gaussian acceptance. These gates exposed and
verified corrections to the radial outer-branch `(1-split)` derivative, affine
inverse coordinates and the massless-origin radial right derivative. Formatting,
core/API library checks, core test checking and clippy passed; clippy retains
17 pre-existing style warnings. Extreme-tail validation and final numerical root
residual certification remain a separate kernel-hardening slice.

The amplitude candidate study generated a massive two-loop kite in under a
second, checked its four thresholds and same-orientation intersections, and
identified a massive planar double box as a second inexpensive topology. Actual
amplitude surface registration is still missing from the production compile
context; supplied-geometry map tests and LMB fallback are not evidence that
graph-derived amplitude surface maps already work.

LU-h-matched profiles are now part of the goal, with a dedicated independently
reviewed research note and portable one-dimensional checks. They are proposed
profiles, not implemented settings or evidence of a GL638 improvement. General
conditional sides, automatic channel construction, the complete saved-state
acceptance harness and all-orientation GL638 improvement remain open.

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
`sampling_channel_weight`. The latter controls summed versus Monte Carlo
execution only; channel identity is selected by the separate lists.

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
subspace_lmb = [3, 10]
parent_lmb = [3, 6, 7, 10]
on_cut = [2,6,10]
```

Supported function heads are `lmb`, `surface`, `cut`, `soft`, `collinear`,
`complement`, `product`, `intersect`, `then`, `phase_space`, `left` and `right`.
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
as `s/(s+t)^2`, or certify equivalent coverage from a selected ordinary channel
with a positive probability. A mixture must use its actual full density and
CDF/inverse (or all inverse branches), rather than an artificial denominator
floor. Users may disable the component only with explicit diagnostics about
the resulting reference-moment assumptions. Use the same deterministic ordinary
fallback in forward and inverse when the cut has no regular radial root.

Each selected host contributes a canonical channel, while every channel still
evaluates the full physical cut/CT sum. Evaluate all foreign densities at the
same raw momentum point with their own `t_c`; never substitute the selected
host's t into another cut's density. Validate radial normalization, tails,
independent Cartesian determinants, simple-residue radial flattening, and the
raised-residue oracle `integral t^j h^(j)(t) dt = (-1)^j j!` where endpoint
terms vanish. Exercise the loaded-state reference harness with broad coverage.
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
The initial fixture is the five-propagator massive scalar two-loop kite, obtained
from `tests/resources/graphs/double_triangle.dot` with every internal mass
positive. Its overall UV degree is -2 and every loop subgraph is UV finite.
At masses one and rest-frame energy five, it has both two-particle and
three-particle thresholds. A boosted configuration supplies a regular
intersection of the genuinely coupled three-particle surfaces. The massive
planar double box is a secondary topology; two-loop six-photon production is an
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
