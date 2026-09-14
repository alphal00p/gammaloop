# Exact two-normal H/Z proposal and its A-star pullback

Analytic audit begun 12 September 2026; source status updated 14 September.
The bounded generic component and amplitude binding are implemented and validated
for the supported class. Discrete proposal transport passed its original tests,
but a subsequent audit found that continuous map changes during physical rescue
also require a [fixed canonical draw](CANONICAL_DRAW_SOURCE.md); that correction
passes 214 core tests and three API gates. The hosted GL638 channel and A-star
pullback remain proposals. The construction strengthens the one-normal map:
on a regular compact patch, a density proportional to `1/R` bounds a leading
`1/R` weight, whereas `|H|^(-beta)`, `0<beta<1`, gives finite local variance but
leaves the weight unbounded on generic rays.

## Exact conditional momentum chart

Work first in prepared physical cut1 coordinates, with fixed a,t, top mass m,
`Et=E(t)`, `B=E(a+t)+Et`, and

```
h = H = E(a+p)+E(p)-B,
z = Z = E(p)+Et+|p-t|-Q,
u = E(p), S=B+h, D=Q+z-Et.
```

The original unsquared energy signs require
`u>=m`, `S-u>=m`, and `D-u>=0`. Squaring with those conditions gives

```
2 a.p = S^2 - 2 S u - |a|^2,
2 t.p = 2 D u - D^2 - m^2 + |t|^2.
```

For `a x t != 0`, solve these two linear equations in their plane:
`p_parallel=d+e*u`. An orthonormal frame built from a and the part of t
perpendicular to a avoids an unnecessary ill-conditioned Gram-matrix inversion.
Let n0 be the unit normal to that plane. Then

```
p = d + e*u + w*n0,
w^2 = u^2-m^2-|d+e*u|^2.
```

The inverse Jacobian follows directly from the triple product of the gradients:

```
|d^3p / (dh dz du)| = u(S-u)(D-u) / (|a x t| |w|).
```

Both signs of w must be included. A simpler complete parameterization removes
the apparent endpoint singularity. For a physical nonempty H ellipsoid,
`S>sqrt(|a|^2+4m^2)>|a|`; hence `kappa=|e|^2-1>0`. Set

```
uc = -(d.e)/kappa,
delta_u^2 = uc^2 - (m^2+|d|^2)/kappa,
u = uc + delta_u*cos(phi),
w = sqrt(kappa)*delta_u*sin(phi),  0<=phi<2*pi.
```

A full circle includes both inverse branches, including their joining points.
The remaining Jacobian is smooth at `w=0`:

```
Jp = |d^3p / (dh dz dphi)|
   = u(S-u)(D-u) / (|a x t| sqrt(kappa)).
```

Require `delta_u^2>0` and the physical energy signs throughout the sampled
circle. If the interval `[uc-delta_u,uc+delta_u]` must be clipped to
`[m,min(S-m,D)]`, its endpoints determine the allowed phi arcs explicitly via
arccos. Uniform sampling on those arcs uses their exact total length in its
density; it does not reject points or retain only one square-root branch.
The initial compact-patch design should use a complete physical circle with
strict margins, so this additional bookkeeping is unnecessary.

## Normal-plane sampling and normalization

Choose a fixed positive scale `alpha=|P0|/Q` within a patch, with P0 separated
from zero, and use

```
h=R*cos(theta), z=R*sin(theta)/alpha,
|d(h,z)/d(R,theta)|=R/alpha.
```

Uniform R in `(0,Rmax)`, theta in `[0,2*pi)`, and phi in `[0,2*pi)` induce

```
q_p(p | a,t) = alpha*|a x t|*sqrt(kappa)
              / [4*pi^2*R*Rmax*u(S-u)(D-u)].
```

All factors, including Rmax and any restricted angular-domain lengths, belong
in the inverse density. This proposal has `q_p~1/R` on a regular patch. If the
complete integrand is `C(theta,phi,y)/R+O(1)` with bounded coefficient and
smooth nonzero remaining densities/Jacobians, its leading importance weight
is bounded there. This does not prove bounded weights or finite variance in
other regions or near tangencies.

Not every complement a,t admits H=Z=0. A normalized implementation can avoid
both unknown rejection probabilities and missing support in either of two ways:

1. Choose a normalized compact complement patch around a known regular
   intersection. Certify a disk `|h|<=Rmax`, `|z|<=Rmax/alpha` throughout that
   complement support, with strict margins on `|a x t|`, kappa, delta_u squared,
   and the original energy signs. Mix this local proposal with a globally
   supported ordinary proposal. This only addresses the certified patch.
2. Draw complements y from any normalized g(y). Decide from y alone whether
   an admissible disk can be certified. If so, use its explicitly calculated
   Rmax(y); otherwise use a normalized ordinary conditional momentum law.
   The conditional distribution integrates to one in either case. Recompute
   the same deterministic decision when evaluating its inverse density, and
   include the full outer mixture density in every weight.

The explicit polynomial/rational coefficients permit interval bounds for such
a disk. Adaptive reduction of Rmax is acceptable only when the final domain
has a certificate; testing a finite grid does not establish domain coverage.
Computing the entire global feasible H/Z domain is unnecessary. An outer
mixture supplies positive density outside a compact focused chart.

## Targeting the actual A-star image

The same conditional p chart can be used for p-star rather than sampled p.
For the A=(7,8) counterterm, let `RA=sqrt((Q/2)^2-m^2)`, sample
`s_star=RA*n` uniformly on its sphere, and use the same center as its native
two-loop subtraction group:

```
(p_cut,s_cut) = c + ell*((p_star,s_star)-c),  ell>0.
```

For fixed complements and `|cs|<RA`, this covers every `s_cut!=cs` exactly once.
Writing `v=s_cut-cs` and `D_A=RA^2-|cs|^2`, its inverse is

```
ell = [cs.v + sqrt((cs.v)^2+D_A*|v|^2)]/D_A,
(p_star,s_star) = c + ((p_cut,s_cut)-c)/ell.
```

For `cs.v<0`, rationalize the ell expression as
`|v|^2/[sqrt((cs.v)^2+D_A*|v|^2)-cs.v]`. Its six-dimensional Jacobian is

```
J_A = ell^5 * RA^2 * (RA-cs.n) * d^3p_star d(ell) dOmega.
```

Substituting the two-normal p-star chart gives the positive Jacobian factor

```
ell^5 * RA^2 * (RA-cs.n) * (R/alpha) * Jp.
```

With original, unrescaled a,t drawn from a normalized six-dimensional density,
cut1's tau depends only on those complements. Returning the two active
momenta to original coordinates adds **tau^(-6)** to this Jacobian. The density
is the product of the normalized complement, ell, sphere, R, theta and phi laws,
divided by the complete Jacobian. Dependence of tau and c on the complements
is block triangular and adds no other determinant. Drawing prepared a,t freely
instead would require a separate shell/radial parameterization and its Jacobian.

The fixed-center premise is supported for this full GL638 native `[3,7]` group:
its six surfaces are from cut1, and the SOCP's inactive energies/spatial shifts
depend only on the cut-preserving complements. Its objective is unweighted
minimax. If the origin shortcut does not apply, Z and F fix the ideal unique
center `(p,s)=(t,t)`; the actual recorded numerical solver center differs by
about 1e-6 GeV at the examined maximum. Subsets, other metadata or additional
cuts in a group can change this conclusion.

Any chosen fixed approximate center still defines a normalized sampling chart
when its own inverse and Jacobian are used. It does **not** place the singular
density on the exact runtime CT image, and therefore cannot establish its
asymptotic variance regularization. Use the actual common center for that claim,
cover every relevant A-containing overlap center with a normalized mixture,
and maintain positive conditioning margins. The existing full graph and its
threshold-group operations must remain intact.

## Independent checks performed

No GammaLoop execution was needed for this audit. On the existing prepared
hard-max complements, h=0.1 GeV and z=0.2 GeV give
`u in [263.0012,414.4574] GeV`, below the physical upper bound 555.0715 GeV,
with kappa=12.9427. Four circle angles, including the joining point phi=pi,
reproduce the unsquared H/Z equations within 9e-13 GeV. A direct central
finite-difference determinant agrees with Jp within 3.8e-9 relative.

Separately, the six-dimensional A-star cone map with nonzero cp and cs gives
the claimed determinant within 1.5e-9 relative; its inverse recovers ell within
2.3e-16 and p-star within 2.9e-14 GeV. These checks substantiate the formulas
at regular points, not a global domain certificate or an implementation claim.

## Generic shared-energy component and remaining graph binding

The represented-geometry component, support/partition migration, common
`intersect` block and routed-energy matcher have passed their earlier focused
checks. The current source adds amplitude command-card binding and native
proposal-policy transport; their combined numerical gates now pass.
The component is not yet a usable hosted GL638 channel through the run-card API.

The first exact primitive accepts two prepared energy equations in one
three-dimensional active block, rather than graph names or H/Z formulas:

```
f1(x;y) = E0(x) + E1(x+a(y)) - C1(y),
f2(x;y) = E0(x) + E2(x+b(y)) - C2(y),
Ei(v) = sqrt(|v|^2 + mi^2),  mi >= 0.
```

The existing Esurface/routing owner must establish that each target has exactly
these two varying energies; all other energies and external shifts enter the
complement-only constants. Shared energies are identified by equal mass and
exact signed routing, including external shifts, not by a common edge ID.
A common invertible affine change to x is retained by the existing frame map.
The bounded first matcher accepts active coefficients equal up to sign; an
unsupported routing class is a structural diagnostic, not an ordinary fallback.

GL638 supplies `a=q2`, `b=-q10`, masses `(mt,mt,0)`, and the original global
equation constants `C1=Q-Eh(a)` and `C2=Q-Et(t)` in the prepared host frame.
The analytic expression `C1=Et(a+t)+Et(t)` agrees only on the exact host cut;
substituting it at finite LU-root accuracy moves the requested target by the
host residual. Retain the original global equation during binding.

Its shared energy is carried by distinct edges 12 and 3. Graph conservation
gives `q12=q3+Q_in`: their spatial momenta coincide in the specified
centre-of-mass frame, but their full formal external signatures need not match.
The initial matcher required identical complete signed routing and a common mass
expression. The canonical-draw correction adds an exact directed check of this
fixed-external spatial identity, retaining the original rows; a tolerance
comparison or an assumption about boosted kinematics is insufficient. This belongs
to the sampling matcher, not to the structural raised-edge equivalence rule.
The kernel receives none of these graph labels. Amplitudes supply the same
prepared equation class without a cut or projection owner.

The matcher also corrects an existing whole-sign comparison: independently
normalizing internal and external signs falsely identifies `L+Q` with `-L+Q`.
One common sign must relate both parts. The existing raised-edge callers retain
this structural invariant, with dedicated regression coverage. The factory
preserves fixed-energy multiplicities and rejects missing prerequisites,
ambiguous shared occurrences and unsupported energy classes explicitly.

This is a graph-agnostic **exact supported class**, not support for every pair
of routed E-surfaces. The generic machinery is the existing block, context,
support, branch, inverse-density and native Jacobian contract. A later implicit
rank-two component can solve arbitrary routed residuals under that same
contract, after proving local rank, branch uniqueness and its feasible domain.
It must use the same `intersect` AST and canonical registry; no parallel parser
or second catalogue is implied. Until then other equation classes fail clearly
at binding, even when a generic numerical root finder could find one solution.

For requested residuals h,z put S=C1+h, T=C2+z and u=E0(x). The general plane
constraints, retaining the unsquared energy signs, are

```
2 a.x = S^2 - 2 S u + m0^2 - m1^2 - |a|^2,
2 b.x = T^2 - 2 T u + m0^2 - m2^2 - |b|^2.
```

For `|a x b|>0`, solve in an orthonormal plane frame to obtain `x=d+e*u+w*n`.
Define `k=e.e-1`, `B=d.e`, `C=m0^2+d.d`, `D=B^2-k*C`. A regular full circle is

```
uc=-B/k, delta=sqrt(D)/k,
u=uc+delta*cos(phi), w=sqrt(k)*delta*sin(phi), 0<=phi<2*pi,
J = |dx/(dh dz dphi)| = u*(S-u)*(T-u)/(|a x b|*sqrt(k)).
```

The inverse evaluates the original two residuals and u at the supplied point,
then uses `atan2(w/(sqrt(k)*delta), (u-uc)/delta)`. A half-open full circle
includes both signs of w once; its joining points phi=0,pi have finite J.
Compute inverse density at that supplied point, not at a forward reconstruction.
Clipped arcs, tangencies, vanishing energies and additional algebraic branches
are outside this first primitive. The origin of normal polar coordinates is
an explicit singular boundary, not a silently finite density cap.

### A certified disk and a normalized fallback

Use alpha=1 initially, or any fixed positive normal-scale parameter; neither
requires the GL638 physical coefficient P0. For each prepared complement y,
choose a dyadic candidate radius rho and certify the enclosing rectangle
`|h|<=rho`, `|z|<=rho/alpha`. The following sufficient inequalities enforce a
complete physical circle everywhere in that rectangle:

```
|a x b|^2 > 0, k > 0, D > 0,
L0=-B-k*m0 > 0,        L0^2 > D,
L1=k*(S-m1)+B > 0,     L1^2 > D,
L2=k*(T-m2)+B > 0,     L2^2 > D.
```

The last six conditions prove `uc-delta>m0` and
`uc+delta<min(S-m1,T-m2)` without squaring away a sign. These are rational
expressions in the normal coordinates and prepared data; after establishing
positive Gram and k denominators, no transcendental domain proof is needed.
Halving rho is allowed only until these inequalities are certified on the
whole rectangle. A finite point grid or a native epsilon margin is not a proof.

There is currently no rigorous interval owner in the sampling code. The
smallest arithmetic addition is a private bounded enclosure calculation for
these predicates using the existing `rug`/MPFR dependency, with directed
rounding at **every** add, multiply, divide and square. Import the actual native
inputs through their existing numeric wrappers: binary64 exactly, Quad as
its high and low limbs enclosed under directed addition, and Arb from its
stored MPFR value. Do not first narrow to binary64; the existing Quad-to-106-bit
MPFR conversion alone is not a proof of exactness for arbitrarily separated
limbs. Any small exposure of those representations belongs in `utils::FloatLike`,
not a second numeric hierarchy. Keep the enclosure code private to this
certificate until another implemented primitive needs it.

The arithmetic certificate applies to the mathematical energy equations with
those prepared parameters. To certify the exact physical equations rather
than the represented prepared chart, the host must also supply enclosures of
its derived constants/coordinates. In particular it does not certify the
accuracy of the existing LU root or a future CT projection. Native point/J
representability and selected J times inverse-q checks remain separate gates.
Unresolved interval signs must not be called proven absence. Use typed numerical
retry when precision prevents the specified branch decision; a deterministic
complement-only decision to decline the compact chart is a normalized fallback,
not a claim that the physical intersection does not exist.

For a certified disk sample R uniformly on (0,rho), theta and phi uniformly
on their full circles, with `h=R*cos(theta)`, `z=R*sin(theta)/alpha`. The cube
Jacobian is `4*pi^2*rho*R*J/alpha`. If the complement-only certification policy
declines the chart, the **same named channel** uses the existing normalized
ordinary three-dimensional conditional map. Forward and every foreign inverse
must reproduce that policy and rho from the same ordered complement data.

For complements using the compact chart, its inverse returns outside-support
for points outside the disk/image. Keep an ordinary full-support channel in
the same canonical catalogue and use the existing exact-density partition;
this supplies coverage outside the patch. Do not select overlapping compact
and broad laws inside one cube branch and report only the branch Jacobian.
The current fixed-partition estimator is not changed to a mixture-PDF estimator.

The first implementation uses a bounded 2048-bit MPFR enclosure calculation,
with separate enclosures of the mathematical box endpoints. Certified failure
of a conservative lower bound permits dyadic reduction; arithmetic uncertainty
about that endpoint's sign remains a typed error. It tries
`max_radius / 2^j`, `0 <= j < 32`, then uses the ordinary conditional map if
the deliberate domain policy declines every candidate. This arithmetic budget
is not a claim to resolve every native input. In addition to the singular
normal origin and disk boundary, the ordinary fallback retains its existing
open-cube angular endpoints; an exact angular zero there is diagnosed without
clipping or remapping.

Precision replay adds a production gate beyond within-lane forward/inverse
agreement. An active-point-dependent physical retry must not silently change
the compact/fallback choice or dyadic radius and thus select another proposal
law for the same draw. Fixed MPFR certificate precision alone cannot prevent
this when the prepared host geometry changes between native lanes. Distinguish
deliberate box overestimation from arithmetic uncertainty in its bound, and
retain/check proposal decisions through the existing per-draw context owner
for selected and foreign maps. Rebuild native geometric coefficients and
certificates from original inputs; never promote earlier prepared numbers.
The current source implements this transport in the original-source owner:
prepare participating rows at existing 1000-bit Arb precision before norm-based
lane selection, seal discrete decisions, and compare complete native choices
through the same typed borrowed context on every replay. The same map traversal
prepares any required roots/centers; it does not run final bodies, physical
shared-overlap solving, rotations or events. Ordinary catalogues bypass the
phase. Its source regressions and production amplitude acceptance now pass
independently of the earlier fixed-context component tests.

The detailed prerequisite-only policy and its exact bias counterexample are in
[the shared preparation audit](AFFINE_STAR_IMPLEMENTATION_AUDIT.md#proposal-policy-across-native-retries).

### Existing owners and dependency order

1. `sampling_joint.rs` implements the native component under the existing
   `sampling_maps.rs` forward/inverse contract. Outside support is represented
   by `Result<Option<SamplingMapEvaluation<T>>>` and propagated through affine,
   embedding and ordered composition. Positive supported
   densities still undergo existing validation; `None` never represents
   underflow, a failed solve, or a selected forward point's missing inverse.
2. `sampling_selection.rs` passes outside-support to the **existing**
   `SamplingPartition::from_log_scores` None convention: zero channel weight;
   all-None already errors. The migrated map contract separates Full/Restricted
   image coverage from complement dependence. It replaces the previous
   Full/Conditional/Branched enum without an independent domain registry.
   An explicitly selected full-support sibling covers a compact chart.
3. `sampling_evaluator.rs` retains the common eager/dual program. This prerequisite
   is now implemented: select the three cube derivative columns while holding
   prepared parameters fixed, using the existing determinant owner. Static
   identity-zero seed metadata is supplied to Symbolica during vectorization,
   so a prepared-only singular derivative cannot contaminate an active column.
   The `u+sqrt(m)` boundary at `m=0` passes in Double, Quad and Arb; an actually
   requested singular derivative or a singular intermediate depending on active
   inputs still gives a typed numerical error. Expressions remain compiled once
   in the current program cache, with worker-local buffers and native parameters.
4. `cff/esurface.rs` owns the routed-energy matcher and prepared constants;
   `sampling_selection.rs` admits exactly `intersect(surface(...),surface(...))`
   as one existing block/registry key. Resolve both leaves under one host/frame,
   consume three coordinates total, and preserve ordered prior dependencies.
   The amplitude binder now returns this same component with its native affine
   translation: empty-prior blocks bind once at warmup through static `Affine`,
   while nonempty priors use ordered embedding. Cross-section binding awaits the
   physical host gates. No second
   channel enumeration, automatic discovery, arbitrary implicit intersection
   solver or new projected-target resolver belongs in this slice.

A later certified star pullback consumes this identical three-dimensional
component through the existing affine/context embedding when the projection
is affine in its active block. A larger cone chart can compose it as a child.
Actual center/alpha ownership and non-affine sensitivities remain requirements
of that outer map, not assumptions hidden in the joint kernel.

### Bounded validation before physical pilots

- Unequal masses, one massless energy, nonorthogonal shifts and translated or
  reordered parents; both w signs and phi joining points. Check the original
  unsquared equations, direct inverse support and an independent Cartesian
  finite-difference determinant, then native Quad/Arb references.
- Certify complete boxes, deliberately infeasible boxes, rank-zero and nearly
  tangent complements. Include a case where finite probing misses a violating
  interior point, and verify that refinement never labels an uncertain sign
  as absence. Inspect the returned margins/rho, not just sampled points.
- Shifted Gaussian normalization and raw moments with focused plus ordinary
  channels, nonuniform discrete/continuous weights, foreign points outside the
  focused image, and complements switching to normalized ordinary fallback.
  Include full composed determinants with complement-dependent affine frames.
- Reuse actual routed GL638 H/Z records and a routing-qualified amplitude pair
  as graph-backed witnesses; verify their class before inventing a fixture.
  One orientation is suitable during development. Final production acceptance
  must retain every selected amplitude orientation and all 936 GL638 keys,
  six cuts and the complete subtraction/UV configuration.
- On a regular compact patch, test the expected 1/R density and bounded weight
  for a controlled leading C/R integrand. This proves neither global finite
  variance nor star alignment. Benchmark certificate/map costs before a long
  GL638 run; no performance gain follows from the formulas alone. The complete
  sampling pipeline must meet the warmed 10% physical-runtime budget below.

### Implemented component validation

The initial combined run passes 161 core tests, including the existing map,
partition, native solver and raised-cut regressions. After routine explicit
test initialization, all nine joint/enclosure tests pass again. Core and
Python-feature API checks pass; clippy reports no warnings on changed lines.
Two API gates also pass: saved-state reference acceptance and summed versus
selected canonical-channel evaluation. These use the existing loaded-state
and integration owners; no joint graph binding is implied by those regressions.

| Component gate | Evidence |
| --- | --- |
| Circle and native arithmetic | Unequal masses `(2,3,0)`, nonorthogonal shifts `(3,0,0)` and `(1,4,0)`, both circle signs and joining points; independent Cartesian determinant, Quad/Arb perturbations below binary64 resolution and independent worker buffers. |
| Rigorous support | Directed native enclosures, including separated Quad limbs; certified disk, declined chart, explicit singular/boundary errors and certified outside-support. |
| Complete estimator | Compact plus ordinary channels integrate a shifted normalized Gaussian and raw second moment within the fixed 2% and 4% tolerances; compact-only admission is rejected. |
| Foreign channels | Different complement-dependent energy data and affine frames retain their own inverse support through an ordered, permuted six-dimensional composition, in exact and proxy modes. Constant proxies test support only. |
| Intended corner power | In Quad and Arb, `J*w/R` from either selected channel approaches the same independent finite analytic limit over four radial decades; final relative error is below `1e-3`. |

The focused component command is:

```sh
cargo nextest run -p gammalooprs --lib --profile test_gammaloop --retries 0 \
  -E 'test(sampling_joint::) | test(joint_bridge) | test(joint_density_partition) | test(native_mpfr_enclosures)'
```

The final nine-test run takes 103.401 seconds in the test profile; its
Gaussian/partition test accounts for 102.144 seconds. This includes 8192
reference draws and repeated support/certificate evaluations, not a production
map benchmark. Prepared-data reuse and production map cost need measurement
when the graph binder and per-draw context are connected.

The accepted runtime requirement is a conservative
`T_sampling / T_physical <= 0.10` on matched warmed, optimized GL638 samples,
using all 936 orientations, six cuts and the actual UV/threshold/stability
settings on 20 cores. Include the fixed-Arb proposal-policy pass, native maps,
Jacobians, foreign inverses, proxies/partition, certificates and all repeated
sampling work during rescue in `T_sampling`. Count shared preparation once;
conservatively attribute it to sampling and remove it from the physical
denominator, documenting the allocation. The denominator must be the actual
complete physical evaluation, with no redundant preparation added for timing.

Measure representative production draws and stored hard/max-weight samples,
reporting each set's aggregate ratio, per-draw distributions, high quantiles
and extrema with repeated-batch timing uncertainty. Report first compilation,
warmup and training separately. Unoptimized amplitude or component test timings
do not establish this gate. Optimize only enough to meet the 10% bound, then
stop performance work; retain correctness and normalization requirements.
This budget is pending measurement, not a new passing or improvement claim.
The historical optimized
[X2 GL638 pilot](GL638_X2_PHYSICAL_PILOT.md#fixed-budget-pilot) recorded
72.04–79.47 ms total per all-orientation draw, with 0.151–0.170 ms map time for
direct-H selections and 0.075 ms for six optimized LMBs. It predates the joint
chart and fixed-Arb policy phase. These values provide context only; the current
unoptimized amplitude test rate likewise establishes neither success nor
failure against the new matched optimized GL638 budget.

### Graph-binding groundwork validation

The subsequent combined run passes 168 core tests in 141.997 seconds. The
resolver treats the qualified pair as one ordered 3D block, its registry key
retains both equations, and the catalogue shares one compiled joint program.
The generated kite checks the matcher and component in Double, Quad and Arb,
including a parent with nonzero affine translation. Independent evaluations of
the original E-surfaces agree at mapped and unrelated Cartesian points.
A serial-edge graph checks distinct globally reversed shared routes, unequal
partner masses, fixed-energy multiplicities and inconsistent mass diagnostics.
Missing prerequisites/external ports and underflowed fixed energies are rejected.
Existing raised-cut, raised-selection and shared-threshold regressions pass.
These were factory/component checks while the factory was staged under
`cfg(test)`. They precede the current production amplitude binder and policy
transport. The full GL638 external-frame specialization and hosted joint
channel remain unimplemented.

Core checking, Python-enabled API checking and all-target core clippy pass,
with no changed-line warnings. Both saved-state reference-acceptance and
summed/selected canonical-channel API regressions pass in 125.211 seconds.
They validate the existing loaded-state workflow; they do not enable or certify
production joint channels. Formatting and diff checks pass.

### Current source milestone — validated amplitude binding and transport

The amplitude binder accepts the supported shared-energy pair through the
existing `intersect` block and complete parent frame. For the boosted kite the
definition is
`then(complement(4),block(lmb(6),intersect(surface(2,4,6),surface(3,5,6))))`,
with parent `[4,6]` and an explicitly selected ordinary sibling. The first
runtime prescription uses trial radius `e_cm * sampling.b` and normal scale
one. `sampling.power` does not change the joint uniform-R law or introduce
anisotropy. Both target equations and their native translation remain intact.
Empty-prior blocks bind once at warmup through static `Affine`; only nonempty
priors use ordered embedding.

The added tests pass native map comparisons, all-18-orientation Gaussian/value
and moment acceptance, retained choices under selective Double-to-Quad body
retry, norm/debug prepasses, ordinary-phase bypass and generating-row separation.
Actual physical checks establish `X = raw * J*w` and selected-momentum
`raw * w` at the regular point.
`generated_triangle_joint_binds_without_prerequisites` passes the actual
empty-prior binder with nonzero offset, compact roundtrip and one production
reference evaluation.
These source fixtures clear rotation probes; they do not establish behavior
under nonidentity probes.

The combined run passes 170 focused core tests in 110.864 s. The complete kite
test passes in 1537.369 s, including 8192 Gaussian draws satisfying the fixed
6% normalization and 8% moment bounds. Both saved-state/canonical-channel API
regressions pass in 104.496 s; core checking takes 17.15 s and API checking
15.51 s. Clippy passes in 73 s with 52 existing warnings, none on changed lines;
formatting and diff checks pass. The full kite timing comes from the unoptimized
test profile and measures the complete fixture, not production map cost. The
matched optimized GL638 10% runtime gate remains pending.
In the controlled law-switch example, 1 versus 3/2
is the compact-generator conditional oracle; with the equally selected ordinary
sibling and constant support-gated proxies, the complete mean is 1 versus 5/4.

These changes do not supply GL638's exact fixed-external routing specialization
or, at that milestone, authoritative Cutkosky `t*` handoff. CT-star pullbacks
separately need the actual common overlap center and native projection alpha.
The subsequent native handoff implementation below keeps routed equations and
roots in the lane; metadata retains only cross-lane decisions/history. No all-orientation
GL638 bounded-weight or variance improvement follows from amplitude acceptance.

### Geometry and timing milestone — selected gates pass

Explicit LU preparation now retains an ordered native `EsurfaceRay` beside its
root; raised eta jets consume those coefficients. This deliberately routes
before scaling in the LU path. Other scalar CT/fiber evaluators and their
external-only seed convention remain unchanged. The old and prepared routes
can differ under finite-precision cancellation; the new tests preserve that
counterexample rather than asserting general equivalence. This milestone
predates the native handoff below. Hosted joint binding, exact GL638 external-frame
specialization and CT center/alpha handoff stay open.

The existing timing counters now accumulate source/map/partition work and
failed native attempts, norm/debug remaps and one inclusive canonical policy
pass. Actual target bodies are timed separately across all lanes and rotations;
evaluator/event fields are subsets, not extra physical costs. Final snapshots
retain all accumulated time. Focused timing tests pass; exhausted error
returns still emit no timing metadata. The optimized matched GL638 10% budget
remains unmeasured, with shared preparation attribution still to be established
after host adoption. Unclassified work and first native binding costs must be
reported separately from the warmed comparison.

The combined shorter run passes 174 core tests in 124.853 s and three API tests
in 103.571 s. Combined checking passes in 48.73 s; clippy passes in 58.10 s with
68 existing warnings, none on changed lines. The complete kite rerun passes in
1550.337 s, completing 175 selected core tests. Its 8192 draws satisfy the fixed
6% normalization and 8% moment bounds, with repeated physical, body-retry and
timing checks passing. This is an unoptimized complete-fixture duration, not
production map cost; the separate optimized 10% runtime gate remains unmeasured.

The [GL638 LU-ray physical replay](GL638_LU_RAY_REPLAY.md) additionally passes
three native calls: hard `hz_03` in the ordinary stack and forced Arb, and
double-soft `soft22` in forced Arb. All 936 orientations, six event identities
and weights, and the original subtraction/UV settings are retained; all 35
saved-state hashes remain unchanged. Relative complex total changes against
the frozen earlier results are `2.1690e-16`, `1.3995e-294` and `1.3671e-281`,
respectively. The baseline includes intervening commits, so this is a bounded
physical regression rather than isolated attribution or an accuracy gain.
These six simple cuts have no higher raised eta derivatives. Their coverage
comes from generated raised-cut fixtures; no sampling-budget claim follows
from this unoptimized raw-momentum replay.
The linked artifact independently verifies the post-build private tuple type
alias as a semantics-neutral source change; it introduces no new numerical
behavior beyond the tested snapshot.

### Native host handoff — combined gates pass

The existing runtime row retains native ray/root records and freezes the initial
selected forward/direct-inverse prefix before partition work. Exact same-plan
priors may reuse preparation; rounded own/foreign inverse priors prepare their
own density without overwriting authority. Selected/summed transport is native,
rotations retain the original-frame payload, and retries reconstruct the source.
Canonical/norm/debug roots have isolated histories. Physical adoption authenticates
graph/group/parent/prior identity against the selected immutable channel, rotates
the retained ray once, and independently routes the completed sample. The fixed directed 2048-bit certificate bounds only
represented scalar roots, first derivatives and explicit inverse-volume factors.
It does not bound coefficient-vector discrepancies, higher raised jets or the
full conditional/joint density; selected `J*q` reciprocity adds numerical chart
consistency, not such a proof. Hosted H/Z remains guarded until its host sensitivity
and exact external routing are handled. Host-only adoption work is included in
sampling time and subtracted from the enclosing physical-body timer on success
and error.

The final source passes 184 unique core tests: the 176-test focused
gate (118.591 s), seven Gamma sample tests and the complete kite (1504.601 s).
The conditional fixture was rerun after its final test-only update (56.996 s),
with strict 1e-8 nonidentity physical-total agreement. Existing event ownership
buffers/counts only identity probes; surrounding selected/direct event comparisons
remain intact. The complete kite again passes its 8192-draw 6%/8% reference bounds;
its unoptimized duration does not measure production map overhead. Three API
gates pass in 106.024 s, combined checking in 12.56 s and clippy in 57.76 s
(68 existing warnings, none on changed lines). Formatting/diff checks pass.
The optimized GL638 10% budget remains unmeasured; hosted joint binding and
full host-to-H/Z sensitivity propagation remain open.
