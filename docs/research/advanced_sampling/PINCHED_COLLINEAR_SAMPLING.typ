= Sampling pinched collinear configurations
<sampling-pinched-collinear-configurations>
Research note, 13 September 2026. This is a proposal for a later extension of advanced sampling, not a claim that the parsed `collinear(...)` primitive is already executable. No production changes accompany this note.

== Recommendation
<recommendation>
The proposed transverse/light-cone construction is feasible. For a massless one-loop collinear segment, it gives an exactly invertible, normalized map of the full three-dimensional integration space. A useful alternative is an exact prolate chart: its normal coordinate is the energy deficit itself, and its other coordinate resolves the two soft endpoints. Both can fit the same canonical sampling-channel catalogue and the existing conditional-map design.

Three qualifications matter:

+ The factor `k_perp^2` appears in the measure when using a #strong[logarithmic] transverse radius. Sampling uniformly in that logarithm down to zero is not normalizable. A normalized fractional power, or an appropriate logarithmically softened boundary profile, is necessary.
+ A segment is one kind of pinch. Multiple independent collinear loop momenta produce simplices with more transverse directions; massive and overlapping pinches require other local geometries. A generic architecture must represent their rank and geometry without treating every pinch as a regular threshold shell of zero radius.
+ For a cut-dependent segment, its endpoints and transverse frame must come from consistent solved cut kinematics. Changing momenta and retaining an obsolete `t*`, or treating frame vectors as invariant scalar context, invalidates the construction.

For an absolutely integrable residual with a known transverse power, these maps can give locally bounded weights away from soft endpoints and other intersections. This is not a proof of bounded weights for an entire graph. The complete physical cut sum and threshold counterterms remain unchanged. Sampling cannot repair a failure of their local cancellation.

This note complements #link("SOPER_ANGULAR_REVIEW.typ");, #link("TWO_NORMAL_PROPOSAL.typ");, and #link("LU_H_MATCHED_SAMPLING.typ");. The latter optimizes the auxiliary LU rescaling; the present proposal optimizes physical shape variables. They solve different parts of the sampling problem.

== 1. The one-loop segment and the actual integration measure
<1-the-one-loop-segment-and-the-actual-integration-measure>
Fix two spatial points `A` and `B=A+P_vec`, with `P=|P_vec|>0`, and let `n=P_vec/P`. Choose an orthonormal transverse pair `e1,e2`. Write

```
k = A + z P_vec + q_perp,
q_perp = r (cos(phi) e1 + sin(phi) e2),
tau = r^2,
z in R, tau >= 0, phi in [0,2 pi).
```

The momenta of the two massless lines are `k-A` and `B-k`. For an on-shell massless parent, their positive-energy deficit is

```
eta(k) = |k-A| + |B-k| - P
       = sqrt(P^2 z^2+tau) + sqrt(P^2(1-z)^2+tau) - P.
```

The triangle inequality shows that `eta=0` precisely on `tau=0, 0<=z<=1`. There is no interior region with `eta<0` around this pinched surface. This is the massless segment described in the thesis\'s treatment of pinched E-surfaces; it is different from an existing non-pinched ellipsoid.#footnote[Zeno Capatti, #emph[Singularities of Feynman diagrams and their local cancellation in collider cross-sections];, doctoral thesis (2023), #link("https://doi.org/10.3929/ethz-b-000647466")[ETH record];. The user supplied `TMP_TO_IGNORE/thesis_zeno.pdf`. Relevant locations are §2.3 on cycle bases, §5.1.1 on the one-loop segment, §5.1.2.3 on collinear/soft/pointlike pinched surfaces, and §5.2.2 on multiloop collinear and soft limits.]

At fixed `0<z<1`, expanding in the transverse displacement gives

```
eta = tau/[2 P z(1-z)] + O(tau^2/P^3),
```

where the coefficient hidden in the remainder is not uniform at the two endpoints. The expansion requires `r << P min(z,1-z)`. A useful exact form that avoids subtracting nearly equal hard energies is

```
eta = P (|z|+|1-z|-1)
    + tau/[sqrt(P^2 z^2+tau)+P |z|]
    + tau/[sqrt(P^2(1-z)^2+tau)+P |1-z|].
```

At a zero numerator and denominator, use the continuous value, not `0/0`. This formula is a numerical diagnostic and a possible expression for the map\'s geometry; it does not alter the physical evaluator.

The spatial measure follows directly from the cylindrical determinant:

```
d^3 k = P dz r dr dphi = (P/2) dz d tau dphi.
```

Consequently the user\'s proposed factor appears as

```
y = log(r/mu):       d^3 k = P r^2 dz dy dphi,
y_tau=log(tau/mu^2): d^3 k = (P/2) tau dz dy_tau dphi.
```

It is absent when integrating directly in `tau`. This distinction is essential: a uniform probability density over `y in (-infinity,infinity)` does not exist. A Cartesian proposal behaving as `1/tau=1/r^2` all the way to the axis has a logarithmically divergent normalization in the two normal directions.

=== What a literal light-cone variable changes
<what-a-literal-light-cone-variable-changes>
For one positive-energy on-shell momentum, relative to the same fixed axis,

```
k_plus = E+k_parallel,
k_minus = (tau+m^2)/k_plus,
E = (k_plus+k_minus)/2,
k_parallel = (k_plus-k_minus)/2.
```

Direct differentiation at fixed `tau` gives

```
d^3 k       = (E/k_plus) d k_plus d^2 q_perp,
d^3 k/(2 E) = d k_plus d^2 q_perp/(2 k_plus).
```

With `zeta=k_plus/(2P)` for a massless parent,

```
d^3 k       = E/(2 zeta) d zeta d tau dphi,
d^3 k/(2 E) = 1/(4 zeta) d zeta d tau dphi.
```

At the forward collinear locus `zeta=z`; away from it they differ. If the second line has fixed spatial momentum `P_vec-k_vec`, its positive on-shell energy does not generally give plus fraction `1-zeta`. Two non-collinear positive-energy massless daughters cannot sum to a fixed lightlike four-vector. Imposing that extra condition would require a recoil or a different parent kinematics.

Thus the simple spatial chart is sufficient for the current loop integral. A Lorentz-invariant splitting/phase-space chart is also possible, but its recoil map and measure must be included explicitly. The exact phase-space factorization in Catani--Seymour is a useful comparison, not permission to substitute an invariant measure for GammaLoop\'s raw spatial measure.#footnote[S. Catani and M. H. Seymour, #emph[A General Algorithm for Calculating Jet Cross Sections in NLO QCD];, #link("https://arxiv.org/pdf/hep-ph/9605323")[hep-ph/9605323];, §5, particularly the phase-space factorization around equations (5.17)--(5.21).] Existing LTD/on-shell energy factors must not be counted twice.

== 2. An exact normalized cylinder covering the full space
<2-an-exact-normalized-cylinder-covering-the-full-space>
Let `p_z(z)>0` be normalized on the whole real line and let `p_tau(tau|z)>0` be normalized on `tau>0` for every `z`. Sample

```
z   = F_z^-1(u_z),
tau = F_tau^-1(u_tau | z),
phi = 2 pi u_phi.
```

The conditional density in physical Cartesian coordinates and determinant are exactly

```
q_k(k) = p_z(z) p_tau(tau|z)/(pi P),
J_k    = pi P/[p_z(z) p_tau(tau|z)].
```

The inverse obtains `z=(k-A).n/P`, subtracts the longitudinal vector, computes `tau=|q_perp|^2` and the azimuth with `atan2`, then applies the CDFs. The axis and the azimuthal seam are measure-zero coordinate degeneracies, handled with deterministic conventions. There is no missing finite-volume region: allowing `z` outside `[0,1]` is necessary for full support.

For example, a simple longitudinal distribution has positive masses `w_minus+w_middle+w_plus=1` and dimensionless tail scales `ell_minus,ell_plus`:

```
u_z < w_minus:
    z = ell_minus log(u_z/w_minus),
w_minus <= u_z < w_minus+w_middle:
    z = (u_z-w_minus)/w_middle,
otherwise:
    z = 1-ell_plus log((1-u_z)/w_plus).
```

This is normalized, invertible, and allocates explicit probability to the segment. An endpoint-enhanced middle branch can instead use `z=sin^2(pi v/2)`, where `v=(u_z-w_minus)/w_middle`. Its conditional middle density is `1/[pi sqrt(z(1-z))]`. The two exterior branches still matter.

A convenient transverse family is the log-logistic CDF with a momentum scale `b>0` and exponent `a>0`:

```
tau = b^2 [u_tau/(1-u_tau)]^(1/a),
F_tau(tau) = (tau/b^2)^a/[1+(tau/b^2)^a],
p_tau(tau) = (a/b^2) (tau/b^2)^(a-1)/[1+(tau/b^2)^a]^2,
d tau/d u_tau = tau/[a u_tau(1-u_tau)].
```

For `0<a<=1`, this gives `q_k ~ r^(-p)`, where `p=2-2a` ranges from zero to just below two. It therefore exposes an adjustable integrable enhancement of the correct codimension, with a simple exact inverse. Log-domain evaluation avoids overflows in the tails.

The scale may depend on `z` and on previously sampled complement data. The determinant remains triangular: derivatives of `b(z)` enter off-diagonal blocks. The same is true if `A`, `P_vec`, and the transverse frame depend only on earlier complement variables. If they depend on the current active variables, this product determinant is no longer justified. Such dependencies require a different full map and its actual determinant.

These piecewise CDF branches are implementation details of one channel, not a second channel enumeration. If an alternative construction uses overlapping stochastic branches, its density must sum the probabilities of all preimages.

== 3. An exact prolate alternative, including endpoint information
<3-an-exact-prolate-alternative-including-endpoint-information>
Soper\'s elliptic coordinates use the distances to the two endpoints. His sampling organization also keeps the proposal\'s cut label separate from the complete sum of physical cuts. These are directly useful ideas here; his particular endpoint densities and residual powers need not apply to GammaLoop\'s subtraction representation.#footnote[Davison E. Soper, #emph[Choosing integration points for QCD calculations by numerical integration];, #link("https://arxiv.org/pdf/hep-ph/0103262")[hep-ph/0103262];, especially §IV on sampling versus physical cut sums, §VI equations (25)--(35) on elliptical coordinates, and §IX on the collinear region and its endpoints. His NLO examples do not establish the residual powers of GammaLoop\'s higher-loop threshold-subtracted integrands.]

The following formulas rederive that geometry in the notation above. Define

```
xi  = (|k-A|+|B-k|)/P,        xi>=1,
v   = (|k-A|-|B-k|)/P,        -1<=v<=1,
delta = xi-1.
```

Its explicit forward map is

```
z   = (1+xi v)/2,
r   = (P/2) sqrt((xi^2-1)(1-v^2)),
k   = A + z P_vec + r(cos(phi)e1+sin(phi)e2),
eta = P delta,
d^3 k = (P^3/8) (xi^2-v^2) d delta dv dphi.
```

The energy deficit is linear in the normal coordinate #strong[exactly];, not just near the segment. The chart covers the full spatial space; `delta=0` is the segment, while `v=+/-1` also describes its exterior collinear axes. The inverse uses the two distances and the transverse azimuth. The original distance formulas can be replaced by stable rationalized equivalents near the segment without changing the map.

Choose any normalized positive `p_delta` on `(0,infinity)` and `p_v` on `(-1,1)`, with uniform azimuth. Then

```
q_k = p_delta(delta) p_v(v)/[2 pi (P^3/8)(xi^2-v^2)],
J_k = 2 pi (P^3/8)(xi^2-v^2)/[p_delta(delta) p_v(v)].
```

For instance, the same log-logistic quantile for `delta`, with a dimensionless scale and `p_delta ~ delta^(-beta)`, realizes `0<=beta<1`. At a fixed interior position this means `q_k ~ r^(-2 beta)`. A normalized endpoint profile `p_v ~ (1-v^2)^(-gamma)`, `0<=gamma<1`, adds a separately controlled soft enhancement. `gamma=1/2` has the elementary arcsine CDF; general powers can use a normalized beta CDF with a certified inverse.

The joint behavior is informative. Near `A`, write `k-A=rho(cos(theta)n+sin(theta)e_perp)`. At a fixed non-collinear angle,

```
delta      ~ rho(1-cos(theta))/P,
1+v        ~ rho(1+cos(theta))/P,
xi^2-v^2  ~ 4 rho/P,
q_k        ~ rho^(-1-beta-gamma)
             (1-cos(theta))^(-beta) (1+cos(theta))^(-gamma).
```

This provides a joint soft/collinear profile in one exact one-loop chart. Its soft radial power is below three, as required for normalization. The angular limits of this expression must still be studied as coupled limits, not replaced by the fixed-angle power. It is a promising first graph-backed pinch map because it avoids multiplying unrelated endpoint and axis approximations.

The cylinder remains useful when one wants direct control over the physical transverse momentum or when only a fixed local normal plane is available. The prolate option is restricted to this two-focus massless geometry. Neither should be advertised as a universal coordinate system for every multiloop pinch.

== 4. What is sufficient for bounded weights?
<4-what-is-sufficient-for-bounded-weights>
These are local statements about the #strong[complete post-cancellation integrand] in raw coordinates, including its existing physical factors. Suppose `|F| ~ r^(-alpha)` in the two transverse directions, with the longitudinal position in a compact interior part of the segment. Let the actual proposal behave as `q ~ r^(-p)` there.

#figure(
  align(center)[#table(
    columns: 2,
    align: (auto,auto,),
    table.header([Property], [Pure-power criterion],),
    table.hline(),
    [Absolute integrability of F], [`alpha<2`],
    [Normalizability of q near the axis], [`p<2`],
    [Locally bounded F/q], [`p>=alpha`],
    [Finite second moment `integral F^2/q`], [`p>2 alpha-2`],
  )]
  , kind: table
  )

For the cylindrical family, `p=2-2a`. Thus bounded weights require `a<=1-alpha/2`, and finite variance requires `a<2-alpha`, in addition to `a>0`. A residual `1/r` is locally bounded with `a=1/2`; a smooth Cartesian density has logarithmically divergent variance for that same residual.

For normal codimension `c`, replace the number two by `c`: the corresponding conditions are `alpha<c`, `p<c`, `p>=alpha`, and `p>2alpha-c`.

There is no normalized importance density that makes a non-absolutely integrable pure-power singularity absolutely integrable. Borderline forms such as powers times logarithms require their own analysis: for example a normalizable density can approach `r^-c` with a sufficiently strong logarithmic denominator. A universal fixed fractional exponent is not guaranteed to dominate every such integrable boundary case.

Nor does a strongly singular #strong[proxy score] by itself improve a weak sampling map. The true proposal must put sufficient probability near the singularity. With channel probability `a_i`, exact density `q_i`, and a partition `w_i`, the sampled estimator is `F w_i/(a_i q_i)`. One valid balance partition is `w_i=a_i q_i/sum_j(a_j q_j)`, giving `F/sum_j(a_j q_j)`. The current fixed map-density partition instead uses `w_i=q_i/sum_j q_j`, with the outer `1/a_i` accounted for separately; this is also unbiased. A proxy partition must be checked channel by channel; the same cancellation of an arbitrarily chosen score with `q_i` cannot be assumed. The density power criteria above concern the corresponding actual estimator, not a claim that only one of these partitions is valid.

== 5. Soft endpoints and more than one collinear momentum
<5-soft-endpoints-and-more-than-one-collinear-momentum>
The interior quadratic formula ceases to be uniform when `z P` is comparable to `r`. Near `A`, its leading structure is instead

```
eta = sqrt((z P)^2+r^2)-z P + O(r^2/P).
```

On-shell factors may now be soft as well. The prolate endpoint profile above is one way to handle this joint behavior for a pair of massless lines. Additional soft spherical channels centered on `A` or `B` are another useful ingredient. A sum of isolated soft and isolated collinear proposals is not automatically sufficient for every product singularity or every strongly ordered approach. Test paths such as `z~r/P`, `z~(r/P)^2`, and their opposite hierarchy explicitly.

Letting a conditional width `b(z)` collapse at an endpoint can concentrate points in a narrow sector, but its limiting support and inverse must remain defined. Adding a width floor changes its asymptotic enhancement; it is not a proof of endpoint coverage. When `P=0`, the segment collapses to a point and its axis is undefined. Switch to a normalized soft-point chart or an explicit broad fallback rather than dividing by a vanishing `P`.

=== A multiloop massless simplex
<a-multiloop-massless-simplex>
The thesis explains that pinched collinear E-surfaces with several independent loop momenta occupy a simplex, with soft configurations on its faces.#footnote[Zeno Capatti, #emph[Singularities of Feynman diagrams and their local cancellation in collider cross-sections];, doctoral thesis (2023), #link("https://doi.org/10.3929/ethz-b-000647466")[ETH record];. The user supplied `TMP_TO_IGNORE/thesis_zeno.pdf`. Relevant locations are §2.3 on cycle bases, §5.1.1 on the one-loop segment, §5.1.2.3 on collinear/soft/pointlike pinched surfaces, and §5.2.2 on multiloop collinear and soft limits.] The following local metric is derived here to identify a useful map extension.

Take `n` independent momenta and a final momentum fixed by their sum:

```
k_i = z_i P_vec + q_i,                    i=1,...,n,
k_(n+1) = z_(n+1) P_vec - sum_i q_i,
z_(n+1)=1-sum_i z_i,
q_i.P_vec=0.
```

In the simplex interior, every `z_i>0`, including `z_(n+1)`. The pinch has `n` longitudinal directions and `2n` normal directions. Expanding the sum of positive massless energies gives

```
eta = 1/(2P) [sum_i |q_i|^2/z_i + |sum_i q_i|^2/z_(n+1)] + O(q^4),
M = diag(1/z_i) + 1 1^T/z_(n+1),
det M = 1/[product_(i=1)^(n+1) z_i].
```

Here the remainder, with the appropriate momentum dimensions, is again not uniform near soft faces. The matrix is positive definite in the interior. Whiten the two transverse component vectors using `w=M^(1/2) q`. The exact Jacobian of this linear change at fixed longitudinal coordinates is

```
d^(3n) k = P^n [product_(i=1)^(n+1) z_i] d^n z d^(2n) w,
d^(2n) w = R^(2n-1) dR dOmega_(2n-1).
```

The determinant is exact even though the expression for the physical energy deficit is only a local quadratic approximation. Thus a normal-space radial proposal can use the actual codimension `2n`, or independent normal profiles can address hierarchical limits. A simple `1/eta` factor at a two-loop simplex interior behaves as `R^-2` in four normal directions: it is integrable but has borderline infinite variance under a smooth normal density. Actual graph numerators and products determine the relevant power.

This is not yet a full-support global chart: outside the longitudinal simplex this matrix need not be positive definite, and it degenerates at soft faces. A graph-generic implementation would need explicit sectors or another globally defined continuation, with their probabilities and inverse branches accounted for. Direct products of several collinear maps are valid only after verifying independence of their active blocks. Overlapping normal directions require a joint rank analysis.

Pinches with masses can be points rather than segments; mixed soft configurations, degenerate thresholds, and overlapping Landau loci need other geometry. The automatic classification of all of those cases is a separate task from providing an explicit, certified segment map.

== 6. Compatibility with a Cutkosky host and LU h matching
<6-compatibility-with-a-cutkosky-host-and-lu-h-matching>
The segment data for a side of a cross-section cut are conditional physical data. Their construction must follow the host cut\'s kinematics, including `t*`, and retain its orientation, parent frame, masses and external shifts. Two useful cases should be distinguished.

#strong[A block that preserves the host cut momenta.] Let `B_active` embed a subspace variation into full loop coordinates, and let `C_cut` route those coordinates to all cut-edge spatial momenta. The certificate

```
C_cut B_active = 0
```

means that this variation leaves every host cut momentum fixed. For an amplitude-side loop, the pair total momentum can then be fixed by those cut momenta and the earlier complement. Verify that it is fixed under the same variation; do not infer this from two edge labels. In this case one may construct a conditional segment chart without resolving the cut after each local transverse change. Any conversion between physical coordinates and raw LU coordinates still contributes its actual scaling determinant.

#strong[A block that changes the physical cut.] It is not legitimate to change transverse raw momenta while freezing a `t*` previously computed from another point. Use a full chart of physical cut shapes and then the auxiliary flow coordinate, or solve the newly modified cut and account for its full pullback determinant.

More explicitly, let `D=3L` and let `Q(s)` be a genuine parametrization of a regular host shell with `D-1` independent shape coordinates. Under the current multiplicative flow, set

```
K = Q(s)/t.
```

The raw-coordinate Jacobian is

```
|det dK/d(s,t)| = t^(-D-1) |det[partial_s Q, Q]|.
```

The determinant is the `D` by `D` determinant with `D-1` tangent columns and one `Q` column. If `Q=R(n)n`, this reduces to the angular measure times `R(n)^D/t^(D+1)`. Derivatives of `R(n)` cancel against the radial column; derivatives of the actual shape chart do not disappear. An unrestricted three-dimensional segment map cannot simply be appended to `D-1` already independent shell coordinates: the degree-of-freedom count must be correct.

The derivation in #link("LU_H_MATCHED_SAMPLING.typ") then allows the auxiliary density `p_t=h(t)` for a simple residue, with a broad normalized component where necessary for reference-integrand variance. Higher-power cuts require coverage of the `t^j h^(j)(t)` terms, not just the undifferentiated profile. These residue structures follow the normalized LU flow construction.#footnote[Zeno Capatti, Valentin Hirschi and Ben Ruijl, #emph[Local Unitarity: cutting raised propagators and localising renormalisation];, #link("https://arxiv.org/pdf/2203.11038")[arXiv:2203.11038];, §2.1.1 and equation (2.13). The explicit multiplicative-flow derivation and runtime h-family audit are in #link("LU_H_MATCHED_SAMPLING.typ");.] They do not change the collinear chart\'s own Jacobian or cure a remaining physical shape singularity.

For inverse-density evaluation, each channel reconstructs its own host kinematics at the same full raw point. A host cut chosen for sampling never replaces the other physical cuts, their LU solutions, or their contributions. The complete cut and counterterm sum is evaluated at that common point. Threshold grouping, common centers, symmetric CT localization and dual cancellation remain properties of the integrand, independent of proposal selection.

For amplitudes, the same segment geometry can be built from fixed external data and routed loop/complement momenta without a Cutkosky host. An amplitude with an unsubtracted physical IR divergence does not acquire a finite integral from this coordinate change; regularization or the relevant local IR subtraction remains necessary.

== 7. Audited implementation boundaries to preserve
<7-audited-implementation-boundaries-to-preserve>
The following audit describes the working tree on 13 September 2026. Symbols identify the owner because line numbers may move during the current implementation milestone.

#figure(
  align(center)[#table(
    columns: 2,
    align: (auto,auto,),
    table.header([Existing owner], [Contract needed for a future pinch map],),
    table.hline(),
    [#link("../../../crates/gammalooprs/src/integrands/process/sampling_maps.rs")[sampling\_maps.rs];, `SamplingMapDefinition::Collinear` and Symbolica parsing], [Reuse the existing expression language and canonical `SamplingChannelId`. Resolve topology and geometry before enabling the primitive; do not create another channel list.],
    [#link("../../../crates/gammalooprs/src/integrands/process/sampling_selection.rs")[sampling\_selection.rs];, `unsupported_map_error`], [The existing capability error correctly requires routed vectors, a relative angular frame, support/branch rules and a normalized profile. Keep it until a real map supplies those ingredients. A proxy alone does not qualify.],
    [#link("../../../crates/gammalooprs/src/integrands/process/sampling_maps.rs")[sampling\_maps.rs];, `SamplingMapComponent` and `SamplingMapComposition`], [`then` already carries earlier outputs into later contexts and inverts in declaration order. Preserve the triangular dependency proof. Its present callback context and evaluation record are `f64`; generic high-precision runtime maps need a deliberate precision boundary.],
    [#link("../../../crates/gammalooprs/src/integrands/process/sampling_context.rs")[sampling\_context.rs];, `PreparedSurfaceStatus`], [`Existing`, `Pinched` and `Absent` are distinct. A pinch should eventually carry geometry kind, normal rank and a certified frame, not a fictitious positive regular-shell radius. `threshold_radius()` returning `None` is not sufficient to distinguish pinch from absence.],
    [#link("../../../crates/gammalooprs/src/integrands/process/sampling_maps.rs")[sampling\_maps.rs];, `ImplicitSurfaceRadialMap::root_for_direction`], [This kernel assumes a negative origin value before finding a regular enclosing shell. Its no-root return is not a general pinch classifier. A surface with no interior and a zero normal derivative needs its own chart.],
    [#link("../../../crates/gammalooprs/src/momentum/sample.rs")[momentum/sample.rs];, `SubspaceData::solve_signature` and `is_mergable_with`], [Reuse selected defining edges paired with complete signed cycles. Different native parents may represent the same active cycles. Resolve affine shifts and complements explicitly, and certify rank and fixed pair/cut momenta. Bare edge IDs do not determine the subspace.],
    [#link("../../../crates/gammalooprs/src/integrands/process/sampling_context.rs")[sampling\_context.rs];, `PreparedCutSamplingContext`], [Keep graph/cut/orientation/side/full-parent identity together with solved kinematics. Its stored momenta and `t*` are `f64` metadata; they are not a high-precision map cache.],
    [#link("../../../crates/gammalooprs/src/integrands/process/sampling_context.rs")[sampling\_context.rs];, `PreparedCrossSectionMapEvaluation::rotate` and `cast_sample`], [`rotate` rotates the mapped sample and prepared momenta but copies the opaque `runtime_context: Vec<F<T>>` unchanged. Endpoint and axis components cannot be stored there as if they were rotational scalars. Introduce typed vector/frame data or rebuild it in the target frame and precision before conditional maps consume it. Merely casting the opaque values does not recompute geometry accurately.],
    [#link("../../../crates/gammalooprs/src/integrands/process/sampling_selection.rs")[sampling\_selection.rs];, `DeferredCrossSectionSamplingState`], [Forward and all inverse/proxy evaluations must reconstruct their respective contexts at the same common raw point, with matching graph/cut/frame/precision. Never reuse a channel\'s prepared data solely because two channels have the same edge labels.],
  )]
  , kind: table
  )

The rank issue also affects map composition: one segment has one longitudinal and two transverse coordinates, and a shell has one fewer physical degree of freedom before adding the LU flow coordinate. The abstract component dimension methods permit this distinction. Graph routing/coverage checks must not force every intermediate shape block to correspond to three times an independently selected LMB-edge count. Full input and output dimensions must nevertheless remain equal for the complete invertible raw map.

In particular, distinguish #strong[conditional block inversion] from inversion of ordinary function composition. The current `then` map consumes separate coordinate blocks; later blocks depend on earlier output blocks. Its inverse can and does traverse declaration order, since those earlier output blocks are already available in the full point. An ordinary composition `f_m(...f_2(f_1(x)))` instead applies inverse functions in reverse order. The unified requirement is to invert each map in the order dictated by its dependency graph, with exactly the same conditioning values as its forward evaluation; a blanket instruction to reverse all `then` children is wrong.

A deterministic transverse axis can be built by choosing the Cartesian reference direction least aligned with `n`, projecting it orthogonally, and normalizing. Ties define chart seams. For a uniform azimuth the density is independent of that choice; any anisotropic azimuthal profile must define its patch convention and inverse consistently. Under stability rotations, rotating or reconstructing this basis is essential. A covariant physical reference vector can be preferable when one exists and is non-degenerate.

Maps, normalized densities and proxy scores should be built together at warmup using the existing eager Symbolica evaluator direction. Dual evaluators can differentiate the actual algebraic map, including frame and scale dependencies. For a root-defined map, differentiate the defining equation where its implicit derivative is regular, not the history of solver iterations. Pinched zero derivatives are precisely where this regular-root formula must not be used.

Proxy scores must reproduce all intended singular enhancements with at least the required powers and keep their denominator positive. Use a strictly positive, normalized broad component and stable logarithmic sums, with analytic limiting conventions on singular sets. An epsilon added only to a score does not create support in the actual proposal. Geometry changing from existing to pinched or absent must select a deterministic normalized branch; it must not silently remove a channel or truncate a finite region.

== 8. Verification and acceptance tests
<8-verification-and-acceptance-tests>
Before connecting a graph-backed pinch primitive, require:

+ #strong[Map identities.] Forward/inverse and actual determinant checks for interior points, exterior longitudinal tails, small/large transverse radii, translated/rotated frames, conditional widths, and every CDF or axis-patch branch. Compare eager dual derivatives with finite differences away from coordinate seams. Exercise native high precision separately.

+ #strong[Geometry.] A real routed two-massless-line example must reproduce the analytic segment and exact prolate gap. Include a massive point pinch, absent surface, and zero-length segment to check classification/fallback. A solver-origin massless cusp and a pinched energy surface are different cases; both require well-defined diagnostics.

+ #strong[Normalization through production.] Use the saved-state acceptance harness, with a normalized full-space Gaussian and moment checks, through the actual map, inverse density and canonical multichannel bridge. Exercise the pinch channel alone and in mixed catalogues, including conditionally existing/absent surfaces. Preserve support in all branches.

+ #strong[A reference that tests the intended singularity.] Gaussian normalization alone does not establish variance improvement. In `c` normal dimensions, the normalized function

  ```
  f(q) = Gamma(c/2)
         /[pi^(c/2) b^(c-alpha) Gamma((c-alpha)/2)]
         * |q|^(-alpha) exp(-|q|^2/b^2),   0<=alpha<c,
  ```

  multiplied by normalized longitudinal/complement factors, tests known powers. Verify the bounded-weight and second-moment criteria above rather than inferring them from an apparently stable finite Monte Carlo run.

+ #strong[Endpoints and overlaps.] Check single-soft endpoints, strongly ordered soft/collinear paths, two-loop simplex interiors and faces, and overlapping normal directions. Use determinant/rank certificates and normalized references for each supported atlas sector.

+ #strong[Physical invariance.] Compare the complete cut-plus-CT result at common physical/raw points before and after the map. Keep the individual-cut scaling diagnostics separate from the post-cancellation residual. Test host-preserving variations and changing-host variations, along with rotation and precision rescue. Do not change threshold groups to simplify the test.

+ #strong[LU composition.] Combine a shape chart with the h-matched auxiliary map, check both radial endpoints and raised-residue profile coverage, and retain a suitable broad component for finite Gaussian-reference variance.

As a small independent algebra check for this note, a standard-library Python calculation used `P_vec=(2,3,6)`, `A=(0.8,-1.2,0.3)`, a rotated orthonormal transverse frame, exponent `a=0.6`, and symmetric coordinate differences of `10^-6`. Four points exercised longitudinal exterior and interior branches and both moderate and large transverse radii. The maximum relative determinant discrepancies were `1.99e-10` for the cylindrical map and `1.65e-10` for the prolate map. The latter\'s exact energy-gap identity agreed to `2.23e-15`. The cylinder used conditional `b^2=4 sqrt(1+z^2)`; the prolate delta scale was `0.7`. The two-loop metric at fractions `(0.2,0.3,0.5)` gave `det M=33.333333333333336`, matching the formula above. These are finite-dimensional formula checks, not GammaLoop integration tests or evidence of improvement for GL638.

The following standard-library Python snippet reproduces these checks:

```python
import math

P, A, axis, exponent = 7., [.8, -1.2, .3], [2/7, 3/7, 6/7], .6
def norm(v):
    return math.sqrt(sum(x*x for x in v))
def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0]]
e1 = cross(axis, [1., 0., 0.])
e1 = [x/norm(e1) for x in e1]
e2 = cross(axis, e1)

def radial(u, scale):
    value = scale*(u/(1-u))**(1/exponent)
    return value, value/(exponent*u*(1-u))

def chart(u, prolate):
    if prolate:
        delta, derivative = radial(u[0], .7)
        xi, v = 1+delta, 2*u[1]-1
        z = (1+xi*v)/2
        r = P/2*math.sqrt((xi*xi-1)*(1-v*v))
        jacobian = (P/2)**3*(xi*xi-v*v)*derivative*4*math.pi
    else:
        if u[0] < .2:
            z, dz = math.log(u[0]/.2), 1/u[0]
        elif u[0] < .8:
            z, dz = (u[0]-.2)/.6, 1/.6
        else:
            z, dz = 1-math.log((1-u[0])/.2), 1/(1-u[0])
        tau, derivative = radial(u[1], 4*math.sqrt(1+z*z))
        r, jacobian = math.sqrt(tau), math.pi*P*dz*derivative
    phi = 2*math.pi*u[2]
    point = [A[j]+P*z*axis[j]
             + r*(math.cos(phi)*e1[j]+math.sin(phi)*e2[j])
             for j in range(3)]
    return point, jacobian

points = [[.1, .2, .17], [.4, .65, .33], [.9, .45, .79], [.55, .9, .42]]
for prolate in [False, True]:
    errors, gap_errors = [], []
    for u in points:
        point, jacobian = chart(u, prolate)
        columns, h = [], 1.e-6
        for j in range(3):
            plus, minus = u.copy(), u.copy()
            plus[j] += h
            minus[j] -= h
            columns.append([(a-b)/(2*h) for a, b in
                            zip(chart(plus, prolate)[0], chart(minus, prolate)[0])])
        determinant = sum(a*b for a, b in zip(columns[0], cross(columns[1], columns[2])))
        errors.append(abs(abs(determinant)/jacobian-1))
        if prolate:
            gap = norm([point[j]-A[j] for j in range(3)])
            gap += norm([A[j]+P*axis[j]-point[j] for j in range(3)])-P
            gap_errors.append(abs(gap/(P*radial(u[0], .7)[0])-1))
    print("prolate" if prolate else "cylinder", max(errors), gap_errors)
z1, z2, z3 = .2, .3, .5
print("metric determinant", (1/z1+1/z3)*(1/z2+1/z3)-(1/z3)**2,
      "expected", 1/(z1*z2*z3))
```

== 9. Scope and order of implementation
<9-scope-and-order-of-implementation>
No immediate production change is needed solely to keep this option open. The existing canonical catalogue, explicit unsupported-primitive error, conditional composition and distinct pinched status are appropriate starting points. Before enabling conditional physical maps, the frame/precision context contract above must be made real; storing vector components in the current opaque invariant context would be an architectural mistake.

A focused first implementation would be a fixed-segment exact kernel, likely the prolate chart plus elementary normalized profiles, with the cylinder as an alternative when its independent transverse control is useful. This is a moderate map-and-test slice, plausibly several hundred lines using the existing evaluator and map owners. That estimate excludes automatic graph resolution and is not a commitment to a line count.

The next slice should certify graph routing and support a cut-preserving amplitude-side block, including a corresponding amplitude with fixed external data. A fully coupled physical-cut/left/right chart needs explicit shell degrees of freedom and the runtime kinematic/frame contracts; it is more than attaching a profile to the existing string parser. Automatic classification and coverage of arbitrary overlapping multiloop pinches is a subsequent research and implementation milestone.

For GL638, these channels may improve any surviving collinear/soft shape dependence, but they do not by themselves resolve the distinct intersection of threshold and CT-star features studied in #link("TWO_NORMAL_PROPOSAL.typ");. First measure which local normal coordinates and powers dominate its complete residual. Use the corresponding generic map, retain the full cut/CT cancellation, and assess maximum-weight tails and variance alongside the integral. Neither the auxiliary h optimization nor a collinear chart justifies promising globally bounded GL638 weights before those measurements.

== Sources and attribution
<sources-and-attribution>
The coordinate determinants, normalized CDF profiles, power criteria, multiloop normal metric, and proposed software contracts above are derived or audited in this note. The primary sources provide the physical geometry, existing sampling precedent and LU/phase-space context. The local `TMP_TO_IGNORE/ttH_defo.pdf` was also reread for the threshold strategy: its one-/two-loop CT variants and IR/dual-cancellation constraints are retained, not replaced by a sampling prescription.
