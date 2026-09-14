# Certified affine CT-star pullbacks: implementation audit

**X4 map proposal; the star chart is unimplemented.** Initially audited after X2
commit `f2f64fb17`, without new numerical runs or literature searches. The first
source extraction separates representative overlap kinematics from raised
derivative packets. The physical projection owner now solves an unnormalized
dimensionless alpha and derives physical radial quantities from it, including
raised and mixed derivatives. Its focused regressions pass in Double, Quad and
Arb. Sharing those actual data with sampling and constructing a star map remain
unimplemented.
This proposes a generic certified affine class, not a GL638 equation engine or
a claim that projected singularities are already covered by direct-H sampling.

## Inputs and existing owners

Research inputs are [the X4 delivery boundary](CROSS_SECTION_DELIVERY_PATH.md),
[GL638 geometry](gl638.md#exact-a-star-targeting-with-another-one-block-ellipsoid),
[the earlier center audit](TWO_NORMAL_PROPOSAL.md#targeting-the-actual-a-star-image),
[the actual prescription](../gl638/README.md#prescription), and
[recorded affine arithmetic](gl638_geometry_checks.json). The earlier arithmetic
supports the geometry; it does not validate the proposed production boundary.

| Existing owner | Relevant responsibility |
| --- | --- |
| [`LUCounterTerm::prepare_shared_overlaps`](../../../crates/gammalooprs/src/subtraction/lu_counterterm.rs), line 1664 | Classify active variants, partition by group ID and complete solve signature, prepare maximal overlaps, and retain local-to-global center identities. |
| Same file, `ThresholdSolveGroup`, `LUSharedOverlaps`, `LUCTKinematicPoint`, lines 1403–1471 | Complete group membership, actual center, distinct per-cut kinematics, and raised LU data. |
| Same file, `new_overlap_builder`, `solve_rstar`, `base_rstar_loop_momenta`, lines 3276–3614 | Copy active center coordinates into each native parent; solve the physical projection and construct the actual star point. |
| [`overlap_subspace`](../../../crates/gammalooprs/src/subtraction/overlap_subspace.rs), lines 158, 527, 574 | Build constraints from fixed complements, validate centers, and construct maximal-overlap membership. |
| [`CrossSectionGraphTerm::bind_sampling_bridge`](../../../crates/gammalooprs/src/integrands/process/cross_section/mod.rs), line 1641 | Resolve physical targets, certify host-cut invariance and bind native conditional geometry. |
| [`SamplingMapContextTransform`](../../../crates/gammalooprs/src/integrands/process/sampling_maps.rs), line 833 | Existing conditional affine frame, inverse and determinant composition. |

Line numbers identify the audited source and may drift; links name the owners.

## Exact dependency certificate

Let W be the sampler's active displacement space and V the selected CT variant's
solve space, with W contained in V. Certify the following from exact signed
routing, the actual complete group and its runtime settings:

1. Every participating cut's momenta are invariant under W. Its LU root and cut
   acceptance decisions therefore depend only on preceding coordinates.
2. After each native-frame transformation, every fixed complement entering
   existence classification and overlap construction is invariant under W.
   This includes foreign-cut members. Group membership, maximal overlaps and
   center selection must consequently be functions of preceding data alone.
3. The projecting threshold has zero momentum coefficients along W. All its
   other required coordinates have already been sampled.
4. Its projection scale is positive, finite and simple; the target fiber meets
   the existing rank, interior-center and inverse-density requirements.

The certificate covers forced-center routing, origin shortcuts and the actual
SOCP input, not merely the projecting threshold's edge list. For current GL638,
native A=(7,8) has V=[3,7], W=[3], parent=[3,6,7,10]. Its six-surface native group
belongs to cut (2,6,10); A depends on s while p is cut-preserving. This is the
documented eligible structure. U=(8,12,14) depends on p and fails condition 3;
its coupled star chart remains outside this bounded implementation.

## Share the actual physical geometry and scale

Physical evaluation currently prepares overlaps after collecting accepted cuts.
The actual centers are stored in f64 `OverlapStructure`, then copied into native
arithmetic. X4 must reuse those exact physical center values, membership and
probe-frame interpretation. A separately solved, nearby SOCP center defines a
normalizable proposal but does not align its density with the actual CT image.

Expose complement-dependent preparation through the existing LU/group owner,
and consume the same prepared result in sampling and subtraction. Keep the
complete per-cut derivative ownership in `LUCTKinematicPoint`; do not fabricate
a complete sample with guessed active coordinates. Group-scoped preparation
should require only the participating cuts whose dependencies were certified.

The physical projection previously normalized the full V-direction, solved r*,
and used r*/r. That ratio is mathematically independent of p for A, but the
normalized numerical solve still saw p. The existing projection/root owner now
solves alpha directly on the unnormalized displacement. Physical reconstruction
uses `c+alpha*Delta`; sampling must consume this same definition when the full
dependency certificate permits it. No second approximate projection solver is
introduced.

The existing `Esurface::compute_self_and_r_derivative_subspace` and radial-guess
algebra already accept an unnormalized displacement, despite their unit-ray
argument names. The existing `RadialRootDiagnostics` budget bounds the energy
residual, so it can solve the dimensionless scale with the same tolerance and
precision diagnostics. Convert its result using `r*=alpha*r` and
`d eta/dr=(d eta/d alpha)/r`, retaining the residual and iteration count. The
physical base star point and both raised geometry paths now use the same
alpha reconstruction. Their common owner derives physical radius derivatives
from the alpha packet; the existing IFT evaluator differentiates the
unnormalized equation at the actual native representative. Mixed threshold
variations remain physical radial variations, inserted as `delta_r/r(t)` in
alpha. Three focused regressions pass, including the independent closed-form
radius coefficients, native null-direction/affine-frame and mixed-derivative
checks, generated raised-component roundtrips and conditional cut sampling.
Quad and Arb derivative coefficients are checked at `1e-27`, with native inputs
perturbed below binary64 resolution. This validates the physical projection
foundation, not the still-missing sampling handoff.

### Existing higher-order derivative regression

Before applying the alpha refactor, an independent analytic case reproduces a
defect in the current `RstarTDependenceEvaluator`. With selected energy `|q0|`,
threshold energy 4, raw normal component `3t`, and a null-coordinate center of
2, the physical radial root is

```
r*(t) = 4 sqrt(9t^2 + 4) / (3t).
```

At t=1 its value and first Taylor coefficient pass, but the existing evaluator
returns second coefficient `2.4463700961728048` instead of
`1.9912314736290266`. This was executed against the unchanged projection owner,
before the proposed normalization correction; the regression remains in the
existing shared-group fixture. Source inspection identified two conversions
now corrected: HyperDual supplies Taylor coefficients while the generated IFT
program expects raw partial derivatives, and its order-n output needs division
by n!, rather than (n-1)!. The corrected implementation retains this independent
physical-radius oracle and passes higher and mixed-order checks.
The frozen full-orientation GL638 state has `raising_power=1` for all six cuts,
as recorded in [the X1 pilot artifact](GL638_X1_PILOT.json). Its generator requests
zero LU derivative orders and stores no IFT evaluator, so those pilots do not
dispatch the defective higher-order calculation. The
[three-point full-orientation replay](GL638_X4_ALPHA_REPLAY.md) passes for the
alpha parameterization with unchanged event identities and precision choices.
It validates the tested zeroth-order physics, without an accuracy or variance
improvement claim.

## Per-draw handoff and stability rotations

The current routes first diverge before subtraction. In
[`evaluate_stability_level_precise` / `evaluate_all_rotations`](../../../crates/gammalooprs/src/integrands/process/mod.rs),
selected-channel mapping and its foreign inverses run in the identity frame,
then the physical sample is rotated. Summed channels replay the identity-frame
map inside each probe before rotating its output. Each physical cross-section
probe subsequently solves its own cuts and prepares its own overlap centers.

The common host-root entry was consolidated as `Esurface::solve_lu_cut`,
preserving the seed policy and existing settings: inside zero, tolerance
factor 1, 2,000 iterations, 64 bracket expansions and native epsilon times Ecm.
The removed sampling-only wrapper's massless-origin right derivative applies
only when radius, actual mass and all three routed momentum components are
exactly zero; a computed energy of zero alone may underflow. The earlier native
endpoint, shifted massive, underflow and unchanged-policy checks passed.

The new geometry source passes focused validation and returns `(EsurfaceRay<T>,
NewtonIterationResult<T>)` from this entry. Its ordered edge occurrences retain
native velocities, offsets, masses and temporal shift; physical raised eta jets
use that same ray. Only explicit LU changes from scale-then-route to
route-then-scale. The unchanged scalar owner is an independent comparison:
binary64 `k=[1e16,-1e16+2]`, `t=0.1` gives routed momentum 0.125 before this
reassociation and 0.2 after it. This is a deliberate represented-equation
difference, not solver equivalence. Amplitude/static CT and fiber scalar
evaluation stay unchanged. One common seed reduction retains their external-only
offset convention; the subspace seed remains untouched. New gates exercise
this discrepancy, repeated energies, eta derivatives and endpoint underflow.
The [bounded GL638 replay](GL638_LU_RAY_REPLAY.md) passes three native calls at
the stored hard `hz_03` and double-soft `soft22` points, retaining all 936
orientations, six event identities/weights and all 35 saved-state hashes.
Intervening changes prevent isolated attribution to LU-ray reassociation;
no precision or variance gain is claimed. Its simple cuts do not exercise
higher raised jets. The complete kite rerun passes in 1550.337 s with the fixed
8192-draw 6% normalization/8% moment bounds and physical/retry/timing checks,
completing 175 selected core tests and three API gates. This unoptimized fixture
duration is not a production sampling-budget measurement.

This consolidation still does not supply identical numerical host data. The
conditional `SamplingMapContextTransform` uses a certified null-direction
representative, discards its returned ray and uses fresh diagnostics; physical
LU uses the complete sample
and the evaluation's persistent diagnostics. The safeguarded solver ignores
the derivative at its inside endpoint. Sharing its entry therefore removes
duplicate routing and policy setup without reconciling those different input
points or diagnostic histories.

Common equation and solve entry points therefore guarantee the same result
only for the same native inputs and diagnostics. `RadialRootDiagnostics` also
indexes repeated calls by occurrence within each precision. Giving map forward,
foreign inverse and physical solves the same text identity would mispair those
occurrences across precision lanes. Actual reuse must prepare each physical
root once for the relevant per-draw dependency and reuse its result; preserve
the existing lower-precision observation owner rather than creating another
retry history.
Independent host solving does not itself bias a normalized proposal with its
own correct inverse density. It limits exact asymptotic-alignment claims about
the physical singularity; current direct-H acceptance is not invalidated by it.

The minimum X4 extension is therefore to the existing cross-section cut
preparation and `LUCounterTerm::prepare_shared_overlaps` owners together:

- Retain a graph-scoped, typed host preparation beside the existing physical
  cut owner: actual cut/group identity, immutable ordered routing/dependency
  plan, declared prior coordinates, routed energy velocities and constant
  offsets, native masses/temporal shift, and the existing
  `NewtonIterationResult<T>`. This is partial physical data, not a guessed full
  sample. The common `Esurface` owner supplies the same routed equation and its
  derivatives. Feed the retained zeroth-order data into physical cut construction;
  keep acceptance, raised packets and residue derivatives on their current
  owners. Prove acceptance is complement-only or reject the chart.
- Extend the existing `LUSharedOverlaps` / solve-group result to retain its
  canonical-frame geometry and complement-only alpha. Pass this transient
  preparation through the existing map/evaluation context for one original
  draw and native precision. Every foreign inverse requests the same physical
  group from its own supplied point's certified dependencies; it cannot reuse
  the selected channel's geometry merely because their names or cut IDs match.
  Reuse requires the same canonical numerical routing plan and per-draw source
  lineage, or a certified conversion to them. Substituting a selected root into
  a foreign map with another routing law can change that map's inverse density.
  Do not key geometry by stringified floats, current center ordinal or counters.
- Retain that preparation for primary physical evaluation and all probes. Warm
  bridge caches contain immutable bindings, not this per-point geometry. The
  current typed map context carries numerical prerequisites and discrete policy
  access, but no native host handoff: this requires
  extending `SamplingChannelRuntimeContexts<T>` with typed physical data and
  transporting it through the map/bridge result, selected
  `DiscreteGraphSample::SamplingChannel` or direct summed branch, and a native
  borrowed `GraphTermEvaluationContext`. The existing physical `lu_solutions`
  loop adopts matching preparations before events and raised packets. Local
  cut IDs from different grouped graphs must remain distinct. This uses the
  existing evaluation flow, without an opaque global callback cache or a
  fabricated complete `MomentumSample` with unknown coordinates.
  The actual primary physical owner must adopt and validate these data against
  the completed mapped point, not merely memoize speculative map preparation.
  In particular, retained native priors do not repair information lost when
  the completed master point undergoes affine cancellation. Check that it
  represents those priors within the supported geometry budget; otherwise
  request native original-draw rescue or report failure.
  Exact group invariance excludes every active or future dependency; construct
  a full physical sample only after its actual coordinates are available.

The agreed first handoff class remains **unimplemented**. It requires the same
ordered resolved parent and canonical host-required prior subset for each
sharing consumer, with every omitted energy coefficient exactly zero. This
restricts certified sampling pullbacks, not user threshold metadata: physical
groups, solve signatures and centers remain unchanged. Ordered permutations
or a general quotient reconciliation require explicit existing-frame proof;
they cannot be assumed by a cache lookup. Native ray/root payloads live only
in the current runtime row and lane. Canonical preparation discards them and
each retry reconstructs them from the original source; metadata retains only
discrete policies and diagnostic history across lanes.

The initial selected forward, or initial selected inverse for direct momentum
input, freezes its authoritative record prefix **before** partition evaluation.
Only this prefix reaches the physical body; repeated selected consumers with
inconsistent priors error rather than create competing authority. Reuse requires
the same plan and exactly equal canonical native source values. If the selected
channel's own inverse or a foreign inverse reconstructs different rounded priors,
it prepares the actual-point ray/density independently and cannot overwrite or
borrow the selected authority through a tolerance. Keep the existing `J*q`
consistency check and the separate completed-point physical-adoption check.
Partition-only work has no authority and each summed generating row is separate.
This is bookkeeping in the existing bridge/context, not another registry,
catalogue, helper engine or retry owner.

Prepare automatic and forced centers once in the identity frame. Convert the
actual stored f64 center to native arithmetic, then rotate its active vectors
once together with the group's fixed data and correctly routed externals.
Preserve the complete membership identities; reuse the scalars tau and alpha,
and certify the rotated cut/projection residuals and every participating
surface's strict interior margin in native arithmetic. A failed certificate
requests original-draw rescue/error, not a different SOCP center for that probe.
Forced centers retain their existing identity-master-to-native routing rule.
This changes the physical owner convention: the current `OverlapGroup.center`
comment explicitly describes a center already solved in the current probe, so
blindly rotating today's result would rotate it twice. No such change is made
by this audit. Higher-precision replay must rebuild native dependencies and
certificates, rather than promote a lower-precision prepared result.

Independent rotated SOCP solves can choose different approximate centers or
near-boundary memberships. For a corner target the induced shift of H(A-star)
can exceed the intended normal distance and appear as failed rotation stability;
f64 center differences are not repaired merely by evaluating residues in Quad
or Arb. This is a source-level risk, not a measured failure here. Sharing the
actual geometry removes that avoidable divergence, but does not prove bitwise
rotation covariance or a rigorous true-root enclosure. Require native residual,
signed-distance/localization and selected/foreign density-budget checks. The
decisive regression compares host tau, complete group membership, rotated actual
center and physical A-star point before comparing final values, including a
foreign-cut raised case, nonidentity-only probes, summed/MC, cache disabled and
native original-draw rescue. It must also count preparations to exclude one
host-root or SOCP solve per probe or foreign inverse of the same certified
dependency data. The common-equation prerequisite additionally needs exact
massless and shifted massive endpoint checks, and a nonzero tiny mass whose
square underflows in Double/Quad but remains representable in Arb.

### Proposal policy across native retries

Implemented and validated in the current supported amplitude class. The complete
generated kite gate passes alongside 170 focused core tests and two API regressions.

For each target block b, define its proposal decision as `P_b(y_b)`, where y_b
contains only that block's ordered raw prerequisites. The decision may select
compact focusing or normalized ordinary fallback, plus a radius recipe and
dyadic exponent. It must exclude the block's active coordinates, the generating
channel and the physical integrand's stability outcome. Retain the discrete
decision through native retries; rebuild coefficients, roots and the radius
recipe from the original source, without promoting derived lower-precision data.

Every foreign inverse needs its own `P_b(y_b)`: reconstruct that target's
prerequisites from the common raw point through its immutable routing and
preceding blocks. Do not reuse the selected map's prerequisites. Graph,
generating-channel/source lineage and resolved target-block identity distinguish
stored records when summed channels map the same cube to different points.
These validate storage, not mathematical policy inputs: the same target and
prerequisites must choose the same law regardless of the generating channel.
Use canonical graph/channel identities and compiled-tree block paths, not
floating-coordinate hashes or mutable occurrence counters.

The implementation uses one fixed preparation phase in
[`evaluate_from_source_precise`](../../../crates/gammalooprs/src/integrands/process/mod.rs),
before norm-based lane selection. It reconstructs participating generator rows
at the existing [`ArbPrec = VarFloat<1000>`](../../../crates/gammalooprs/src/utils/mod.rs)
budget through the same bridge forward and common partition. It reuses
`evaluate_single`/`evaluate_graph_group` traversal, skipping final targets,
physics, rotations, events and the physical shared-overlap pipeline. Existing
map-required roots/centers still run: an implicit surface may call
[`find_center`](../../../crates/gammalooprs/src/cff/esurface.rs) for an approximate
SOCP seed, then validate it natively. A preparation failure aborts explicitly
at this budget; it does not select another proposal or introduce a second retry
loop. Then seal the decisions and enter the existing stability loop.

Only nonnumeric decisions are stored beside
[`EvaluationMetaData`](../../../crates/gammalooprs/src/integrands/evaluation.rs)
root diagnostics. The existing
[`SamplingChannelRuntimeContexts`](../../../crates/gammalooprs/src/integrands/process/sampling_selection.rs)
owner now supplies borrowed policy access, and
[`SamplingMapComponent`](../../../crates/gammalooprs/src/integrands/process/sampling_maps.rs)
forward/inverse methods receive its typed borrowed component view.
Composition and embedding preserve lineage while supplying their existing
ordered numerical prerequisites. Empty-prior amplitude blocks bind once at
warmup through static `Affine`; only nonempty priors use ordered embedding.
There is no mutable warm-bridge cache or
separate policy evaluator: the existing
[`SharedEnergyJointMap::prepare`](../../../crates/gammalooprs/src/integrands/process/sampling_joint.rs)
algorithm resolves the native decision. It compares its **complete choice** with
the sealed canonical entry, including ordinary fallback. Merely certifying that
a retained disk is valid is insufficient: two valid radii can have different
densities. Missing entries, mismatches or unresolved certificates request the
existing native retry/error; they never insert records or change the law.
Diagnostics identify graph, generator, target block, native precision and both
choices.

The source threads this access through `EvaluationSource::build_gamma_sample`,
[`parameterize`](../../../crates/gammalooprs/src/integrands/process/gammaloop_sample.rs),
direct-channel inverses and summed mapping. `loop_norm_sum` and `debug_sample`
only read sealed decisions. A previously unseen foreign policy cannot inherit
a lane chosen by another active map or physics failure; missing-record replay
errors, including a later block newly reached near a support boundary.

For N channels, preparation adds one forward plus N inverse scores in MC, or
N forwards plus N² inverse scores for summed sampling, once per original draw.
It visits participating rows, not every MC generator or every orientation.
Ordinary catalogues need no extra phase; Arb bindings reuse warmed programs.
Measure this overhead. Lazy preparation still needs canonical original-source
reconstruction, not incidental native values.
The 1000-bit source and 2048-bit predicates certify represented inputs, not
exact physical host roots. Original-source/native representability and actual
host authority remain separate gates. An interval containing zero is unresolved;
more precision cannot remove finite-box dependency overestimation. Persistent
disagreement near a true policy boundary must fail explicitly.

A decisive regression uses a compact law uniform in `R in (0,rho)` and the
normalized target `f=q_rho` on that patch, zero outside. Force a physics retry
only when the original radial cube coordinate is below 1/2. If rescue changes
rho to 2rho for those draws, their weight is 2 while the unretried half retains
weight 1: the expectation is 3/2 instead of 1. Freezing the policy gives 1.
That is the compact-generator oracle. With an equally selected full-support
sibling and constant positive support-gated proxies, changing only that
generator gives the complete estimator mean 5/4 instead of 1. Keep this
partition-aware process check distinct from the component oracle.
The added generated-kite source regression uses a fixed native reference density
and controlled Double body failure, then compares original-source Quad replay,
sealed debug access, norm escalation and ordinary-only phase bypass. The bridge
regression exercises separate generating rows and missing foreign records.
The generated triangle separately exercises empty-prior binding with nonzero
offset, compact roundtrip and one production reference evaluation.
These source checks and the all-18-orientation physical consistency and
8192-draw Gaussian/moment gates pass. Nonuniform probabilities, summed acceptance
and fallback retain their existing coverage. Current joint source fixtures
clear rotation probes. These gates are separate from
the kernel's fixed-context normalization and `F=1/R` bounded-weight tests, and
do not complete cross-section host or CT-center/alpha authority. The complete
kite test takes 1537.369 s in the unoptimized test profile; the optimized GL638
10% sampling-cost budget remains unmeasured.

The accompanying timing source passes core/API checks and the full kite rerun. It charges
canonical preparation once and accumulates native source/map/partition,
failed-attempt and norm/debug remapping costs in the existing metadata owner.
Actual target bodies have disjoint inclusive timing across lanes and rotations;
their evaluator/event counters remain subsets. Final metadata refresh includes
post-snapshot debug work. Terminal error-only returns still expose no counters.
This enables measurement; it supplies neither shared-host attribution nor a
passing 10% budget before native transport/adoption is implemented and timed.

## Resolve a star without requiring user threshold directives

A proposed structural descriptor, inside the existing host and active block, is
`at_star(block(lmb(V...),left(surface(A...))),target)`. Here the inner block is
the projecting CT's solve-space descriptor, not another sampled coordinate
block. The surrounding channel supplies the complete parent and the sampler's
active W. Resolve the descriptor against existing runtime thresholds, native
subspaces and exact solve signatures. The left/right qualifier identifies CT
ownership; it does not restrict the physical cut sum. Amplitudes use the same
structural selection without a physical cut host.

The compact alternative `at_star(variant(id),target)` may disambiguate an
explicit metadata variant. It cannot be the only selector: the current
no-directive LU path need not allocate a metadata registry. If the structural
selection matches several user-split groups, report all candidates and require
disambiguation. Neither spelling nor its parser is implemented yet.

Internally, stable overlap branches use sorted complete occurrence identities:
cut group, side, local threshold and variant when present. A variant ID alone
can have several cut associations. Resolve these through the existing owner
and expand the existing canonical sampling catalogue; do not add a second
threshold or channel enumeration. Runtime existence IDs and center-vector
positions are not persistent branch identities.

## Affine composition and center branches

Let o be the native-frame point corresponding to zero master momentum, tau the
host LU root, and c the actual CT center. On the active block,

```text
p_star = c + alpha * [o + tau*(p_raw-o) - c]
p_raw  = o + (c-o)/tau + (p_star-c)/(alpha*tau)
J_active = (alpha*tau)^(-d)
```

Use the existing target fiber in star coordinates and compose these native
affine transformations. The exact certificate makes complement derivatives
off-diagonal in the full triangular map; boosted external shifts remain in o.
For the earlier notation alpha=1/ell, this reproduces `(ell/tau)^d`.

Every containing overlap center requires coverage. Prefer stable membership-set
selectors expanded into the **existing canonical channel catalogue**, with an
explicit finite expansion cap. When a selected membership is not currently a
maximal overlap containing A, use a normalized ordinary fallback determined by
the preceding data. Do not enumerate only centers observed at one sample.
Six native surfaces allow at most 32 A-containing subsets before pruning; this
is a finite candidate bound, not 32 simultaneous physical centers or evidence
that every subset is feasible. Larger groups require an explicit cap/diagnostic.

A dynamic mixture hidden inside one channel is not a free alternative: the
current bridge uses one inverse and checks selected J*q approximately equal to
one. Multiple preimages would require an explicit summed inverse-density
contract. Canonical branch channels avoid adding a second enumeration engine.

## Minimum gates before an efficiency claim

- Resolve the actual GL638 A variant and every containing overlap; compare the
  sampler's H(A-star) with physical star reconstruction, not with base H.
- Vary p at fixed complements: group membership, actual center and alpha must
  stay invariant. Compare the resulting projection with the existing physical
  output across all relevant orientations and raised orders.
- Check full Cartesian determinants independently, including tau, alpha and
  complement derivatives; test a boosted external and alternate ordered parent.
- Test foreign inverses, selected reciprocity, shifted Gaussian normalization
  and raw moments, and physical summed-versus-sampled channel/event equality.
- Exercise missing/changing overlaps, absent or pinched projecting surfaces,
  center/cone boundaries, small/large alpha and tau, and native original-draw
  rescue. Numerical uncertainty is an error/retry, never silent zero or clipping.
- Only then compare full-UV, all-936-orientation GL638 maxima and variance with
  direct H and ordinary controls. No claim here covers U-star or joint H/Z charts.
