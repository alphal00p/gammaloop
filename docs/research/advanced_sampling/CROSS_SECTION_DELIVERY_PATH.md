# Cross-section delivery path alongside amplitude sampling

Delivery sequence and evidence ledger, updated 2026-09-13. This note distinguishes
implemented boundaries from proposed slices under
[the accepted plan](../../../ADVANCED_SAMPLING_PLAN.md). Consult
[GL638 geometry](gl638.md), [LU h matching](LU_H_MATCHED_SAMPLING.md),
[joint normals](TWO_NORMAL_PROPOSAL.md), and
[recorded physical replays](../gl638/runtime-followup.md) at each relevant slice.

## What is implemented, and what is still a contract

The canonical bridge, cached warmup, exact raw-frame partition, real amplitude
surfaces and standalone full-parent `phase_space(cut(...))` are implemented.
The common native-precision migration landed in `882f80eb0`: physical amplitude
and cut tests rescue maps and foreign densities from original inputs through the
existing stability stack, including improved native external data. Root/density
accuracy certification and outer-grid range handling remain separate limits;
native reference retry is now implemented, with Gaussian-body underflow still
requiring separate treatment; the migration alone does not authorize all strongly
concentrated long-run claims. See
[the precision audit](SAMPLING_PRECISION_RESCUE.md).

The X2 source now binds conditional targets through the existing map embedding
and one shared qualified-block plan. It replaces the earlier test-only complete
prepared-sample handoff, which could not faithfully represent unsampled active
coordinates. The combined numerical gates pass. The physical evaluator
still solves every active LU cut before preparing shared threshold overlaps;
its ordering and complete cut sum are unchanged.

The first useful GL638 H target is **another cut equation evaluated on cut 1**,
not a generated threshold CT of cut 1. H must therefore come from the full
physical energy-equation catalogue, independently of the threshold CT registry.
Only an explicit CT-star target requires a CT variant and its overlap centers.

The production threshold owner already copies its left/right equations from the
global CFF surface cache and evaluates them on the solved, rescaled host sample.
Its subspace evaluator splits the same equation into active energies and a shift
containing every fixed/complement energy. For target H and host C, at the host
root this gives `eta_H = sum(H minus C) E - sum(C minus H) E + chi_H - chi_C`,
with multiplicities retained. Cross-section ST candidates built from the same
initial-state cut have identical external-shift lists, so the last difference
vanishes structurally, including in boosted alternate frames. Keep the complete
global equation in the generic map; using only the target-side boundary energies
would omit the negative host-cut energy shift. X2 tests must compare the global
equation, existing subspace evaluator and independent boundary-energy difference
on the actual host root before checking the map determinant.

## One explicit host/frame and per-child block representation

Extend the existing Symbolica AST and resolver with two transparent constructs:

```text
at_cut(cut(edge_ids...), map)
block(lmb(active_edge_ids...), map)
```

`at_cut` selects the physical cut frame of the enclosed target; it does not
select a contribution to the integrand. `block` binds that child's active
coordinates in the complete channel `parent_lmb`; its `lmb(...)` argument is
an ordered coordinate descriptor, not another sampled map. Resolve signed
cycles, complements and affine routing with the existing owners. The existing
top-level `subspace_lmb` remains shorthand for a single unambiguous active block.
Reject conflicting outer and child specifications. No second target resolver
or channel enumeration is introduced.

For a single active block, retain a shorter user-facing form:

```toml
[sampling.channel_definitions.GL638.direct_H]
parent_lmb = [3,6,7,10]
subspace_lmb = [3]
around = "at_cut(cut(2,6,10), surface(2,4,12))"
```

When routing certifies the dependency, infer the complement and generate it
before the active block. Export the explicit ordered form below so the inferred
layout is inspectable. Multiple or ambiguous blocks use explicit scopes.

A concrete GL638 direct-H channel to validate after the cheap generated gates is:

```toml
[sampling.channel_definitions.GL638.direct_H]
parent_lmb = [3,6,7,10]
around = "then(complement(6,7,10), block(lmb(3), at_cut(cut(2,6,10), surface(2,4,12))))"
```

X2 implements this syntax. Raw complements are generated first; exact routing
certifies that the remaining p displacement cannot change any host-cut momentum.
The host can then solve its complete cut kinematics without guessed active data.
The H map operates in prepared p and returns raw p, with the factor `t*^-3`.
The shared parent resolver now accepts the complete native alternative parent
and its requested order through existing generated routing, without a GL exception.

Keep the affine external shift when changing the LU parent. If the native
coordinates are `L = A K + B Q`, while physical LU rescales the generation-parent
coordinates `K` at fixed external data Q, then the prepared native coordinates
are `A(t* K) + B Q = t* L + (1-t*) B Q`. They are not generally `t* L`.
Construct this frame from the actual prepared LU sample through the existing
affine LMB owner. The inverse active volume factor remains `t*^-3` per loop,
but its center includes the shift. A rest-frame GL638 test can hide the omitted
term; add a boosted-external, alternate-parent regression before claiming that
the conditional host works in arbitrary external kinematics.

For a general cut plus both sides, the corresponding structure is
`then(block(lmb(cut_coordinates...), phase_space(cut(...))),
block(lmb(left_coordinates...), left(surface(...))),
block(lmb(right_coordinates...), right(surface(...))))`.
The unique preceding `phase_space` supplies the inherited host. An explicit
`at_cut` resolves ambiguity if more than one host is present. Left/right are
optional ownership qualifiers on children; the host itself has no side.
Represent the prepared host independently of optional side metadata, rather
than assigning H an artificial left/right label.

The cut block retains one auxiliary LU scale and its cut-shape coordinates;
it does not consume the entire parent dimension and then append side variables.
All leaf blocks together consume exactly `3L` raw degrees of freedom. Both
sides must preserve the host cut; mutual dependencies need a joint chart or a
clean diagnostic, not an unsupported triangular determinant.

Keep `on_cut`'s current meaning as a list of numeric `CutId` restrictions.
The plan example `on_cut=[2,6,10]` must not be interpreted as the edges of a
cut. Recorded GL638 cut 1 has edges `(2,6,10)`; resolve and validate that identity
against the loaded graph before exporting an ID. A restriction must agree with
the explicit/inherited host; it never silently changes the target frame.

### Shared implementation boundary for X2

Extend the existing Symbolica AST and compile traversal once. Each block exposes
its qualified target, complete native parent, active edges, ordered preceding
edges and remaining edges to both graph binding and compilation. The geometry
registry key includes that qualified target (host and optional side), the routed
native parent, canonical active edges and ordered preceding edges. The preceding
layout matters because a compiled callback consumes those actual coordinates;
equal targets can appear as the first or second side of different channel
compositions. This prevents collisions between different cut hosts, cycles or
callback input layouts. It migrates the existing
registry; it is not a second catalogue. Named channels remain distinct proposals.

The existing map embedding prepares a conditional affine transformation from
actual preceding raw coordinates. It returns the physical preceding coordinates
and the active physical-to-raw affine map. Forward and every foreign inverse
prepare independently at their supplied point. Multiply forward determinants by
the affine determinant and divide inverse densities by it. Dependence on preceding
coordinates occupies the off-diagonal block of the ordered Jacobian; later
unsampled dependencies must be excluded by exact signed routing, including every
fixed-energy term. Do not assume the right side was sampled before the left.

Wrap the complete compiled native map in the existing native-to-master affine
owner. Generalize that owner to wrap any compiled map rather than retaining an
LMB-only numerical path. One unique cut block/profile per channel is the first
scope cap; its LU-h profile uses the actual host h and raised order. A direct
`at_cut` target can instead infer an ordinary complement block when its routing
certificate permits it. The LU-h profile retains the generation-origin fixed
point (native BQ) and currently requires it to be certified inside the cut. A
nonzero geometric SOCP center cannot replace that point while retaining the
physical interpretation of R/r as the LU scale.

The earlier test-only complete prepared-sample handoff cannot describe unsampled
active coordinates. Migrate its validation responsibilities to the immutable
qualified block and the existing context-affine preparation boundary. Retain or
reshape the prepared-cut record only for actual fixed/cut data needed there;
never fill a supposedly complete physical sample with guessed zeros or keep two
permanent preparation engines. Preserve its native-precision, invalid-t*, host
identity and foreign-context checks in production tests. Proven spectator values
may be represented algebraically without claiming they were sampled.

## Prioritized owner and test matrix

| Slice | Existing owners and bounded deliverable | Decisive tests and scope cap |
| --- | --- | --- |
| X0: reproducible baseline | Existing GL638 card, `ALL_ORIENTATION_VALIDATION.md`, replay JSON, CLI/API reference owners. Freeze current source/state/card/graph hashes and resolve cut/threshold/channel identities. | All 936 orientations, six cuts, 19 variants, 3D local UV and integrated UV. Replay old maxima, H/Z rays and soft controls. This can proceed while native/amplitude code is developed; no new density required. |
| X1: LU-h profile | `settings/runtime.rs`, `sampling_maps.rs`, `sampling_evaluator.rs`, cross-section bridge binding; reuse `lu_h_function`, `h_dual` and existing implicit cut radius. | Normalized proposal, exact forward/inverse/J, simple-cut radial flattening, raised packets, broad-tail reference acceptance. Start with the fitted log-logistic law plus a normalized broad component; defer exact CDFs/optimized derivative envelopes until measured need. |
| X2: one conditional target | `sampling_context.rs`, `sampling_selection.rs`, cross-section `compile_sampling_bridge` and LU-preparation boundary; reuse `SubspaceData`/LMB transforms and `Esurface` classification/root owner. Bind one host and one active block using the AST above. | Actual generated cut graph with a side threshold and changing complement; all-cut physical summed/MC equality. Also GL638 direct H, which is not in the host CT registry. No general joint chart or star solver in this slice. |
| X3: both sides/composition | Same conditional/context owner and existing `SamplingMapComposition::then`/embedding. Resolve independent child frames, inherit host and preserve cut momenta. | A cheap generated graph with loops on both sides: root consistency, independent full Cartesian J, conditional inverses and loaded-state normalization. Reject active blocks that change cut momenta or form unresolved dependency cycles. |
| X4: exact affine star pullback | Existing `LUCounterTerm::prepare_shared_overlaps`, `LUCTKinematicPoint`, threshold variant metadata and map embedding. Expose certified complement-only center data to sampling and physical evaluation. | Actual GL638 A variant `(7,8)`, native `[3,7]`, every containing overlap center; compare runtime H(A-star), affine inverse and independent J. Use general certified affine pullbacks; keep coupled U-star explicitly unsupported pending its own chart. |
| X5: regular joint normals | Shared `Esurface` callbacks, eager/dual normal Jacobians, existing map support/partition owners. One generic rank-two local chart with certified domain and complete inverse branches. | Both coupled amplitude surfaces and GL638 H/Z use the same engine. Analytic GL638 equations are an independent oracle, not the implementation engine. Compare scalar-normal powers and polar normal laws on certified patches; defer global intersection atlases. |
| X6: automatic catalogue | Existing `SamplingChannelCatalogue`/selection owner; pass fully resolved physical targets from amplitude and cross-section hosts. | Initial discovery covers supported single surfaces/cuts with deterministic ordinary soft coverage, then adds already-supported side/joint/star candidates under caps. Stable IDs/aliases, current-context fallback, explicit omissions and no-surface graphs. Do not claim complete auto support while conditional candidates are unimplemented. |
| X7: demonstrated GL638 gain | Existing integration engine, runtime cards and replay/acceptance reports. Compare frozen catalogues at equal evaluations and equal wall time on 20 cores. | Independent seeds, signed/absolute moments, second moments, maxima, rescue/NaN counts, map/physics cost and soft cancellation. Long runs start only after native rescue and map/reference gates. Report neutral or worse outcomes honestly. |

### X2 generated-graph validation

The first complete focused run passes twelve tests in 66.919 seconds. It includes
the actual serial raised-bubble graph, the two-loop amplitude kite, shared parser
and conditional-map checks, frozen-binding rescue, and shifted reference moments.
The final broader core/API run passes all 141 selected tests in 99.301 seconds,
including saved-state acceptance and physical channel sums. Core/API and Python-
feature test checking pass. All eleven differential regressions, formatting and
clippy also pass. New test warnings are resolved; the existing 1624-byte map-enum
size warning is unchanged from the previous milestone.

The serial-bubble fixture has two nine-dimensional proposals: cut, left threshold
and ordinary right coordinates; and cut, left threshold and right threshold.
It exercises parent `[1,4,7]` at rest and native parent `[3,6,9]` with external
momentum `(5,1,0,0)`. The latter requires the nonzero affine external shift.

| Check | Independent comparison and acceptance |
| --- | --- |
| Physical equations | Global/master, native-parent, existing subspace and boundary-energy equations agree within `1e-9`; the simple bubble supplies an analytic oracle. |
| Host preservation | Changing side radii leaves the raw host block unchanged and its LU root unchanged within `1e-10`. |
| Full Jacobian | Nine-dimensional Cartesian finite differences agree within `3e-5` relative, including derivatives of both side maps with respect to cut coordinates. |
| Foreign inverses | Both channels invert the same supplied raw point and recover it within `1e-7`. |
| LU-h composition | An actual runtime LU-h profile on the both-sides channel inherits raised order two, changes the proposal, and passes selected/foreign inverse and full nine-dimensional finite-difference checks in both frames. |
| Normalized reference | 8192 draws per channel and frame, 32768 total; normalization within 6% and raw second moment within 8%. |
| Native range | Tiny raw cut momenta give `t*` of order `1e110`: finite Double affine entries have an unrepresentable determinant, which raises a typed error; Arb returns finite positive reciprocal Jacobians. |
| Physical estimator | Direct momenta with `J*w`, summed channels and explicit canonical channels agree within `1e-8` relative, including retained per-cut events. Controls require valid metadata and a nonzero physical sum. |
| Reference estimator | Summed and explicit-channel reference values/moments agree within `1e-11`/`1e-10`. |
| Diagnostics | A direct hosted target survives removal of its host-side CT association; a side-qualified target still requires that association. Missing dependencies and conflicting blocks fail explicitly. |

The shifted-reference audit corrected the low-level bridge report to accumulate
raw `|K|²`, while its Gaussian density still uses `|K-center|²`. A displaced
three-dimensional regression has expected moment `4.8125`; the old expression
would approach `3`. Existing committed bridge callers had zero centers. The
frozen X1 amplitude reports used the separate, correct process-level reference
owner, and the GL638 pilot used physical integration; their evidence is unaffected.

These generated tests establish conditional-map and estimator behavior. They do
not establish a GL638 efficiency gain, coverage of projected star images, or a
generic solution for dependent/nested joint charts. The serial fixture disables
UV subtraction; the final GL638 comparison still requires all orientations and
the full local 3D and integrated UV calculation.

X1 now passes the resolved runtime LU h through the existing compile owner and
process warmup, together with the largest
`cut_group.related_esurface_group.max_occurence` among equivalent cut groups.
It does not invent defaults or look up settings per draw. Ordinary cache
invalidation covers changes to those runtime settings. The physical and
saved-state gates are recorded in [the X1 evidence ledger](LU_H_MATCHED_SAMPLING.md#x1-validation-milestone).

The [completed frozen-X1 pilot](GL638_X1_PILOT.md) compares four proposals on
all 936 orientations, with three paired seeds and 20 cores. It records no final
invalid samples, but no consistent LU-h improvement across seeds. The next
measurement therefore targets conditional H geometry after X2 correctness gates;
longer runs of the radial profile alone are not the next dependency.

### Reuse the existing serial-bubble fixture for conditional sides

The existing `tests/resources/graphs/ir_safe_thresholds/triple_dotted_bubble.dot`
and `processes/raised_cross_section_tests.rs` provide a cheap X2/X3 host fixture.
Its complete parent is `[1,4,7]`: three serial one-loop bubbles with a raised
line in each. The existing preprocessing assertion identifies a middle cut
group with threshold CTs on both sides and checks that every discovered cut
remains grouped. Resolve the middle cut and both thresholds from those existing
records; do not assume that one of the repeated-line cut IDs is representative.

The left `[1]` and right `[7]` cycles do not change the middle bubble's cut
momenta. Thus a three-dimensional middle-cut block can prepare its LU scale,
followed by three-dimensional left and right maps, for exactly nine variables.
For equal masses and a rest-frame external energy Q, each simple geometric
threshold has prepared radius `sqrt(Q^2/4-m^2)`. Its raw active radius is divided
by the middle cut's `t*`, giving `t*^-3` for each side's volume conversion.
This is an independent oracle for the routing and scale factors, not the
production surface implementation. Vary the middle-cut coordinates when taking
the full nine-dimensional finite-difference determinant, so the off-diagonal
derivatives of both conditional radii are exercised too.

Start with one side, then both; retain the complete physical cut sum and raised
derivatives in the summed/Monte-Carlo comparisons. The existing fixture disables
UV subtraction, which is adequate for map normalization and fixed-point routing
checks, not a substitute for the final full-UV GL638 comparison. It also has
simple always-existing side thresholds at the chosen kinematics; changing-
existence coverage comes from the shared kite fiber and GL638 target gates.

### Baseline comparability

Recheck the three state hashes in
[`runtime-followup-pilot-provenance.json`](../gl638/runtime-followup-pilot-provenance.json)
before loading the existing full state. The 2026-09-13 inventory found them
unchanged; that check alone does not prove that the current binary can load
the state or that its generated payload still has all 936 orientations.
Confirm those properties from the loaded graph, including its six cuts and
19 threshold variants. Keep raw-point replay separate from proposal testing:
the former has no channel or sampling Jacobian and tests the physical sum.

The historical six-channel pilot is **not** the numerical baseline for the new
sampler. Its cut-dependent inverse-Jacobian partition has been removed; the
canonical bridge uses exact densities in the common raw frame. Reusing the
old option spelling does not restore that distribution. Start a fresh control
with `sampling_channel_weight="map_density"` and six named `lmb(...)` channels
whose edge lists are `[6,12,13,14]`, `[4,12,13,14]`, `[6,7,13,14]`,
`[4,7,13,14]`, `[6,10,13,14]`, and `[4,10,13,14]`. Their explicit output
`parent_lmb` is the generation frame `[3,4,7,10]`. Resolve these lists against
the loaded graph instead of assuming that historical numeric basis IDs remain
stable. Compare this control, `auto:optimized_lmb`, and candidate surface
channels with the same physical settings and fresh workspaces. Historical
errors and maxima remain context, not matched-control measurements.

Preserve the original f64 coordinate tokens when replaying recorded maxima and
H/Z/soft points. Inspect finiteness, every precision attempt, each cut, and both
complex components. Exact A/P coincidence is an expected unresolved diagnostic,
not a successful zero or a regression in proposal normalization. Record current
binary/source hashes and all runtime overrides with the new results; a prepared
card or an archived replay cannot count as a current-binary validation.

The first refresh attempt with commit `06d470409` failed while loading the
recorded state, before display or any of the 33 queued numerical cases. The
diagnostic CLI's SHA256 was
`8fa32d1dac7bea40f623b2f6bf95eb295e92498f978c168105ca39b6e2e1415a`.
Cross-section decoding reported an empty string where
`bitvec::order::Lsb0` was expected (`cross_section/mod.rs:384`). The unoptimized
run took 298.2 seconds and peaked at 4,296,572 KiB; these are loading diagnostics,
not integration timings. All three recorded state metadata hashes remained
unchanged. The current saved-state API regressions pass, but this older payload
is not a usable current baseline. Regenerate from the committed card/metadata,
first one orientation for a cheap setup check, then the unrestricted graph.
Do not patch archived bytes or add a second state decoder for this experiment.

Fresh generation at the same `06d470409` checkpoint succeeded with the optimized
CLI SHA256 `583cfd7134b1cf3307b2b22e01dfbfe47bcf4b95fd07cbf4116326fa5234ecf5`.
The single-orientation setup generated in 15.37 seconds and reproduced its
smoke evaluation exactly after a read-only reload. It contains six cut groups
and 19 threshold variants. Unrestricted generation plus live inspection took
304.64 seconds; reload took 45.90 seconds. Its 936 distinct production
orientation signatures match the prior independently exported complete set,
with contiguous IDs 0 through 935. Reload changed no state file hashes.
Both configurations retain 3D local UV, integrated UV and the committed graph
metadata. Native registry inspection confirms all 19 variants are active in the
full calculation, with 54 components and four metadata evaluators. The existing
`soft22_euler` raw control returns `(1.4532236519417696e-48,
-4.798908444994201e-31)` without a sampling Jacobian: Double fails its stability
check, then Quad passes with reported relative accuracy `5.088e-12` and six
accepted cut events. All read-only checks preserve the state hashes.
The [portable baseline inventory](GL638_CHECKPOINT06D_BASELINE.json) records
commands, hashes, exact registry identities and reporting conventions; full
local evidence is under `/common/dev/gl638_checkpoint_06d470409`. These timings
and smoke evaluations establish usable fresh fixtures, not sampling or IR-safety
improvement.

For comparisons across the native-host milestone, compare the complete
estimator or direct `-m` raw-point values. Checkpoint `06d470409` reports the
X-space map factor separately, while the new native host combines it before
reporting and returns a unit unapplied Jacobian. Comparing the unadjusted
`integrand_result` fields across those conventions would be incorrect.

## The conditional correctness boundary

At warmup, certify that every active displacement has zero signed coefficient
in every host-cut edge momentum. Distinct defining edges are not enough.
The required cut/complement variables must be generated before solving `t*`.
Reuse the production LU root owner and preserve its raised-cut derivative data;
do not build a competing approximate cut solver. For fixed cut-preserving
complements, mapping a d-dimensional prepared block back to raw coordinates
adds `t*^-d`; earlier-complement dependence is triangular. Otherwise a full
pullback is required and the bounded first implementation must reject it.

Register the complete equation, external shift and native routing; reject an
unknown selector separately from a known surface that is currently absent.
Classify `Existing`, `Pinched`, and `Absent` from the current complement and
native precision. An existing body needs a certified interior center, not the
assumption that zero is inside. Select a normalized ordinary fallback from
earlier complement data for absent/pinched fibers. Keep the same canonical
channel and reproduce that choice in its inverse, including near transitions.

At the final raw point, each foreign channel reconstructs its own complement,
host LU solution, target frame, branch and density. Reusing the selected host's
`t*` or context for every denominator term is incorrect. The physical host
evaluates the complete cut/CT sum; proposal Jacobians and channel partitions
remain outside all raised-residue derivatives and threshold overlap weights.

For A-star, exact alignment requires the actual common center, not a nearby
SOCP recomputation. Certify which coordinates its complete threshold group and
all participating cuts depend on. Share that prepared geometry with physical
evaluation or account for the full center dependence in the map Jacobian.
Do not change CT metadata, grouping, dual cancellation or localization symmetry.

## Why this sequence reaches GL638 without waiting for everything

X1 can improve auxiliary radial variation but cannot change physical H/Z shape:
`t* K=R(n)n` is invariant under raw global scaling. Its one-dimensional h/q
improvement is not a GL638 result. X2 is the first density with a demonstrated
local finite-variance argument at the regular H/Z corner. Its single-normal
`|H|^-1/2` law leaves leading weights growing as `R^-1/2`. X4 covers measured
integrated-A extrema that direct H misses. X5 targets leading bounded weights
on a regular patch using `q(H,Z)~1/R`; this is not a global boundedness proof.

The amplitude active-subspace/interior-center work is shared with X2, while
its joint-normal engine is shared with X5. Give one agent an explicit X1/X2
cross-section deliverable at each milestone instead of exhausting the complete
amplitude suite first. X0 replay/harness work, X1 profile derivations/tests,
and X4 dependency audits can proceed independently of amplitude source edits.
Cap the second amplitude topology to validation of a new shared capability;
do not require every amplitude benchmark before a first GL638 pilot.

Reuse `ProcessIntegrand::evaluate_reference_*` and its saved-state tests now.
A reusable CLI/state runner and trained-grid/heavy-tail reference cases remain
deliverables of the accepted harness; avoid a new standalone integrand owner.
New maxima artifacts must preserve cube coordinates, canonical definition and
frame, Jacobian, partition, adaptive/discrete probability and raw point.
Historical raw replays lack the original adaptive PDF and cannot reconstruct
their old Monte Carlo weights. Exact A/P coincidence remains a separate
confluent-evaluation failure; importance sampling must not hide or count it as
a successful zero. Ordinary soft/UV coverage and all-orientation gates remain
mandatory even when a hard-corner pilot improves.
