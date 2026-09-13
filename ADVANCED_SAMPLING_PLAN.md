# Advanced sampling channels for GammaLoop

Status: implementation in progress, 2026-09-13. This document records the accepted
design before source changes. It is the authority for the implementation that
follows. It applies to arbitrary loop order and topology, to amplitudes and
cross sections, and to both ordinary and threshold-adapted sampling.

Current implementation status: the canonical `SamplingChannelId` catalogue,
Symbolica selection parser, explicit parent-LMB validation, explicit active
`subspace_lmb` metadata for named surface channels, prepared cut/side guards,
exact affine LMB routing, and the discrete map-density bridge are implemented
and pushed on `advanced_sampling`. There is one channel catalogue and one ID
domain: the remaining LMB reinterpretation branch is an evaluation detail of
that catalogue during migration, never a second enumeration or index space.
It is scheduled for removal once all momentum-space consumers use compiled
parent-frame maps. Physical directional E-surface maps, prepared cross-section
`t*` context construction, and final old standalone-integrand removal remain
open milestones.

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
  `API_AND_GENERICITY_ADDENDUM.md` and `SOPER_ANGULAR_REVIEW.md`) before
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

The retirement sequence is explicit: first route direct momentum evaluation
through compiled parent-frame maps for every catalogue entry; then remove the
`DiscreteGraphSample::DiscreteMultiChanneling` compatibility route and the
LMB-only prefactor/reinterpretation helpers; finally rename or remove remaining
LMB-specific labels and API quantities where they describe a generic sampling
channel. The ordinal `SamplingChannelId` remains only as the catalogue's stable
position, never as an independently generated LMB index.

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
   estimators. Add amplitude and cross-section cut-context tests.
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

`ChannelIndex` and the current LMB-only discrete enumeration are migration
scaffolding, not a second production channel model. New catalogue, map-density,
grid and evaluator work must use the resolved `SamplingChannelId` catalogue;
once that driver covers the existing sampling modes, remove the legacy index
and its parallel enumeration path rather than expanding both systems.

## 3. Engine contracts

### Geometry and composition

For fixed complement `y`, a regular positive-energy surface is a convex bounded
fiber in its active routing space. Remove spectator/null directions, choose an
interior center and solve a positive radial root. The radial determinant is
`r_star^d rho^(d-1)` times the angular/complement factors. Center dependence on
earlier complement coordinates is triangular; active dependence must be included
in the full derivative.

`product` requires a direct-sum block and independent fixed-complement data.
`then` is an explicitly ordered acyclic conditional map with reverse inverse and
triangular determinant. `intersect` uses scalar normals plus tangent coordinates,
requires full rank and controls all inverse branches and chart patches. A
successful root solve does not establish global injectivity; unsupported or
rank-changing maps fail with a capability diagnostic.

Cut/left/right channels first prepare the Cutkosky kinematics, including each
cut's own solved LU root and `t_star` rescaling. A host cut equation after its
own preparation is not a third varying threshold normal; it organizes phase
space and the auxiliary scale. Left and right threshold maps then act in their
preserving subspaces. The complete physical cut sum remains in every event.
Amplitudes use the same primitives with fixed external data and no artificial
physical-cut context.

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

## 5. GL638 validation and completion criteria

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
