# Certified affine CT-star pullbacks: implementation audit

**Read-only proposal; X4 is unimplemented.** Audited after X2 commit
`f2f64fb17`, without new numerical runs, source changes or literature searches.
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

There is also a numerical projection boundary: the physical code normalizes
the full V-direction, solves r*, and uses r*/r. That ratio is mathematically
independent of p for A, but the normalized numerical solve still sees p. Extend
the existing projection/root owner to expose a complement-only scale alpha for
certified null directions; physical reconstruction and sampling must use the
same definition. Retain r*=alpha*r and the existing raised-residue derivative
calculation. Do not introduce another approximate projection solver.

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
