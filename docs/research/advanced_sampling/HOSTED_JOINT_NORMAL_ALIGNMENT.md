# Direct hosted joint sampling and original-normal accuracy

Implementation audit, 2026-09-14. All selected generic gates pass; GL638 gates remain pending.

## Implemented behavior

Explicit `intersect(surface(A),surface(B))` now binds as one active 3D block
under an existing physical LU host. The existing routed shared-energy matcher,
compiled Joint program, conditional Embedded transform and whole-parent Affine
map remain the only sampling implementation. The supported class retains the
existing ordered-parent and exact omitted-dependency requirements. Matching
uses signed routing and certified fixed external spatial relations, rather than
edge identity. All masses, repeated energies and original temporal shifts stay
in their declared equations.

Both binding and physical adoption use the same extracted cross-section target
resolver. Unqualified original graph targets need no active threshold-CT
association; `left`/`right` targets still require the corresponding physical
association. The resolver extraction preserves the previous candidate and
ambiguity rules. The existing negative-side regression needed only full error
chain inspection after adding channel context.

The fixed-Arb source already owns each complete raw point, J*w factor and
selected host ray/root. Physical retries materialize that same source. The new
checks neither replay maps nor replace a radius, partition, fallback or root.
Amplitude joint targets now use the same original-normal accuracy owner, with
an exactly zero host defect. Target identity there is selected using original
canonical external data before rotating either evaluated point.

## Mathematical invariants

For parent coordinates L with LU fixed point o=BQ, physical coordinates are
`L_phys = o + tau*(L_raw-o)`. The joint component returns
`x = L_phys_active + c0`, where c0 is its common-energy routing offset at the
same physical complement. The existing affine composition returns
`L_raw_active = x/tau + (1-1/tau)*o_active - c0/tau`, with active determinant
`tau^-3`. Trial radius is `e_cm*b`; the normal scale is isotropic one. The full
ordinary sibling and the joint component's normalized ordinary fallback remain
unchanged.

Original H/Z equations are evaluated at the actual scalar completed points:
rotate the retained canonical point once, then rescale canonical and native
loops with their respective materialized tau through the existing LU owner.
Directed evaluation at radius one includes actual scalar rescale rounding,
original signed routing, fixed energies and temporal shifts. It does not
replace original C1/C2 by identities requiring an exact host root.

The existing EsurfaceRay directed arithmetic/fold is shared by original routed
enclosures and retained-root verification; no second solver or interval
framework was added. At fixed 2048-bit directed precision, let R_lower bound
`sqrt(H_anchor^2+Z_anchor^2)` and epsilon be the existing row sampling budget.
Acceptance requires finite ordered bounds, R_lower>0, and both:

- upper norm(native normals - anchor normals) <= epsilon*R_lower/2;
- upper abs(original host residual at completed native point) <= epsilon*R_lower/2.

Finite unresolved/failed inequalities return typed uncertainty; nonfinite range
failures return typed unrepresentability. The existing physical precision stack handles both without changing the source.
This is a pointwise original-normal/on-cut condition, not a uniform relative
whole-disk bound. All adoption/alignment work, including errors and rotations,
is charged to existing sampling time; physical timing excludes that nested cost.

## Validation and current status

| Gate | Evidence/status |
| --- | --- |
| Shared routed enclosure, native rescale and directed half-budgets; retained LU certificate regression | Three focused tests passed before this audit. |
| Generated triangle joint binding without prerequisites | Passed, including successful materialization followed by a typed rotated Double normal failure, same-anchor Quad success, actual retry and forced-Quad equality. |
| Generated hosted heavy kite, bounded checks | Passed in 4.076 s . |
| Full hosted Gaussian/raw-moment acceptance | Passed in 616.215 s: 8192 draws with unchanged 6% normalization and 8% raw-moment limits. |
| Existing all-18-orientation amplitude acceptance | Passed in 954.918 s, retaining all 18 orientations and the full 8192-draw reference, physical and rescue assertions. |
| Formatting, core/API checks and all-target clippy | Pass; clippy took 63 s, with 50 existing warnings and none on changed lines. |
| Broad core regression | All 214 pass: 213 in the initial run; the expected side-error diagnostic assertion passes separately in 57.603 s after inspecting the full error chain. |

The hosted witness keeps the generated physical cuts, groups, production
orientations and counterterms. Stable edge names resolve indices after external
sewing; no topological surfaces or metadata are injected. It checks original
C1/C2 at a finite off-root tau, surviving original CFF spectator rows,
f64/Quad/Arb component geometry, nonzero common offset, actual full 6D and
active 3D finite-difference determinants, both circle branches and foreign
inverses. Existing source transport checks include immutable host payloads and
physical rotation. Actual selected-X and selected-raw physical totals/events
match the bare same-Double-point control times J*w and w respectively. Expected
raw w comes from the canonical inverse at the exact supplied binary64 point.
A sharp source causes actual Double-to-Quad physical retry and agrees exactly
with forced Quad.

The full acceptance runs last, after those bounded checks, restoring the saved
summed-channel runtime. Every six-dimensional Halton draw sums joint and
ordinary contributions before report moments; it is not a selected-channel
partition contribution. The existing 8-draw finite/cost probe remains. Halton
dispersion is diagnostic, not a randomized confidence interval. A targeted
15-minute timeout preserves the required sample count and tolerances.

Together these cover 216 unique core tests. The single diagnostic-only assertion
correction does not change the physical implementation or numerical criteria.

## Scope and concision audit

The substantial new generated fixture shares one graph and the existing
runtime/reference/source helpers. One local generic body checks all three
native precisions; FD comparisons test distinct full-frame and active-frame
Jacobians. No extra generator, sampler, solver, cache or test framework is
needed. Existing comments were retained, and new comments state original-equation,
actual-rescale and proposal-law boundaries. No correctness blocker was found in
this read-only review. Repeated canonical point/mass preparation is currently
inside sampling timing; defer any optimization until measured GL638 cost.

This enables explicit direct original hosted joint targets. It does not enable
CT-star pullbacks, automatic catalogue discovery, arbitrary WH accuracy bounds,
raised derivative-jet accuracy, global bounded weights or global variance
claims. Surviving physical CFF spectators keep their original equations; there
is no universal cut-energy elimination. No GL638 physical, convergence or gain
claim follows from these generic gates. Actual all-936 GL638 verification and
warm optimized sampling-cost measurement remain later requirements. The user
budget is total sampling cost <=10% of the matched physical integrand cost;
optimize only if needed and stop below that budget.
