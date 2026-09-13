# Cross-section delivery path alongside amplitude sampling

Read-only source audit, 2026-09-13. This note proposes bounded implementation
slices under [the accepted plan](../../../ADVANCED_SAMPLING_PLAN.md); it adds
no maps or new physical integration evidence. Consult
[GL638 geometry](gl638.md), [LU h matching](LU_H_MATCHED_SAMPLING.md),
[joint normals](TWO_NORMAL_PROPOSAL.md), and
[recorded physical replays](../gl638/runtime-followup.md) at each relevant slice.

## What is implemented, and what is still a contract

The canonical bridge, cached warmup, exact raw-frame partition, real amplitude
surfaces and standalone full-parent `phase_space(cut(...))` are implemented.
The common native-precision migration remains a prerequisite for long runs
with strongly concentrated densities. It must rescue maps and foreign densities
from original inputs, with consistent improved external data in physical and
sampling evaluation, including Arb precision.

`PreparedCutSamplingContext::from_lu_sample`,
`PreparedCrossSectionMapEvaluation::new`, and
`DeferredCrossSectionSamplingState::new` currently have only test callers.
The production cross-section host does not yet bind these sampling contracts.
Its physical evaluator already solves every active LU cut before preparing
shared threshold overlaps; preserve that ordering and complete cut sum.

The first useful GL638 H target is **another cut equation evaluated on cut 1**,
not a generated threshold CT of cut 1. H must therefore come from the full
physical energy-equation catalogue, independently of the threshold CT registry.
Only an explicit CT-star target requires a CT variant and its overlap centers.

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

A concrete future GL638 direct-H channel is:

```toml
[sampling.channel_definitions.GL638.direct_H]
parent_lmb = [3,6,7,10]
around = "then(complement(6,7,10), block(lmb(3), at_cut(cut(2,6,10), surface(2,4,12))))"
```

This is proposed syntax. Raw complements are generated first; exact routing
certifies that the remaining p displacement cannot change any host-cut momentum.
The host can then solve its complete cut kinematics without guessed active data.
The H map operates in prepared p and returns raw p, with the factor `t*^-3`.
Current production still requires the graph's own full parent; supporting this
native alternative parent is a shared frame-binding extension, not a GL exception.

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

X1 needs one additional binding at the existing compile owner: the current
`GraphTerm::compile_sampling_bridge` arguments do not include the runtime LU h.
Pass its resolved settings from process warmup, together with the existing
`cut_group.related_esurface_group.max_occurence`, rather than inventing defaults
or looking up settings per draw. Ordinary cache invalidation already covers
changes to those runtime settings.

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
