# Canonical physical overlap centers

The GL638 Euler instability was caused by independently choosing SOCP centers
for rotated physical probes. Native root solves accurately evaluated different
rays; increasing their precision could not repair that mismatch. The existing
cross-section overlap owner now chooses one center prescription from each
canonical physical source and retains it across rotations and precision rescue.
The [GL638 report](GL638_HOSTED_JOINT_GATE.md#canonical-physical-centers-actual-gl638-failure-resolved)
records the actual failing-point replay, unchanged source factors and limits.

## Source and physical ownership

`EvaluationSource::prepare_draw` supplies an Arb1000 source for selected,
ordinary/default, summed and raw-unselected inputs. Raw-unselected input has no
fictitious channel and needs no sampling map or partition. Before physical lane
selection, the existing Gamma row traversal visits its participating graphs and
prepares cross-section overlaps. Each graph's existing physical LU equation
owner solves the complete accepted-cut representatives, independent of the
generating chart's sampling-host routing. No alternative root algorithm is
introduced. The selector-only event path determines the accepted cut set without
retaining physical events, observables or event counts.

Canonical selector events carry the real optional source sampling-channel ID,
as native events do. Roots, rays and geometry do not depend on that label;
explicit `SamplingChannelId` selectors may intentionally change the accepted
cut set. Channel independence therefore means the same canonical physical
point and accepted cuts, not replacement of explicit user predicates by `None`.

The existing `prepare_shared_overlaps` grouping performs its identity-frame
SOCP/forced-center choice once, preserving complete cross-cut membership and
local-to-global center ordering. Each completed row retains that immutable
authority; clone, materialization and rotation retain it unrotated. A reference
target or disabled threshold subtraction skips this physical preparation.
Amplitude preparation remains on its existing path.

## Native consumption and failure

Physical attempts retain their native fixed complements, per-cut kinematics,
LU/alpha solves and raised derivative packets. They reconstruct the existing
group metadata and validate accepted cuts, complete solve signatures, ordered
parents, membership and prefactor powers against the canonical authority.
Missing or changed membership gives a typed numerical failure before CT
evaluation; it cannot trigger another SOCP choice. An accepted cut with no
existing thresholds is distinct from a rejected cut. The center-choosing branch
requires identity rotation, preventing an accidental rotated prescription from
being rotated a second time.

Both native center consumers promote the stored binary64 components exactly,
then rotate in T and copy active coordinates by their defining edge. They do not
rotate in binary64 and cast afterward, or reconstruct a center by inverse
rotation. The existing forced-center routing remains an explicit binary64
boundary at the canonical identity preparation. Every consuming center is
checked against the actual native group members before its alpha solve.

Canonical physical preparation isolates radial diagnostic history and temporary
event runtime from the ordinary probes, restoring them on error. Its elapsed
work is included once in `integrand_evaluation_time` and exposed separately as
`canonical_physical_preparation_time`; it is physical preparation, not map work.
Native primal solves still occur. Any later performance comparison must show
this subset and avoid crediting repeated preparation as useful physical work;
future sampling/CT-star sharing would need explicit sampling attribution.

## Evidence and remaining scope

Fourteen focused and 27 broader core tests pass, including the complete hosted
8192-draw Gaussian/moment fixture. Tests exercise selected, summed and bare
source authority, no repeated preparation after failure, missing-cut rejection,
real channel-ID selectors, reference/CT-off skips, and center consumption in
Double/Quad/Arb. Nonzero-center native covariance uses a supplied valid overlap;
actual nonzero SOCP behavior is additionally tested by the saved GL638 replay.
The saved failure now passes Euler in Quad and Arb, with matching Pi2Z identity
totals and six cut weights, while canonical point/J/partition remain exact.

These gates do not constitute separate numerical validation of nonmaster graph
or forced-center source cases, and they do not repair amplitude center behavior.
They also do not provide an A-star proposal, coupled U-star chart, confluent A/P
residue evaluation, global bounded weights or improved Monte Carlo variance.
The complete existing group and center prescription is retained for future
star work; no additional sampler, cache, registry or subtraction body is added.
