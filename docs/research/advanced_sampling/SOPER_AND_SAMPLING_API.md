# Generic angular sampling and a unified channel API

The sampling design should include correlated angular and normal densities,
automatic coverage of massless propagator factors, and conditional maps for
surfaces that appear or disappear with the external data. Both amplitudes and
cross sections use the same primitives. Symbolica expressions provide the
user-facing composition language.

This proposal supersedes details in the earlier API addendum: `auto:surfaces`
now denotes a complete surface-aware strategy, including soft coverage where
needed, and ordered composition uses Symbolica function constructors rather
than the illustrative `x` and `->` syntax. `auto:lmb` selects the unpruned
catalogue, while `auto:optimized_lmb` invokes the existing selection heuristics
and the additional coverage audit. No implementation is implied by the examples
below.

## Angular sampling and the relationship to Soper

Soper includes angular endpoint enhancement and variable shell widths in
Eqs. (39)–(42), correlated soft/threshold sampling in Eqs. (49)–(54), and
angular/radial coupling near secondary surfaces in Eqs. (62)–(65). Sections V
and IX also address soft/collinear configurations beyond regular shell maps.
The quadratic widths in those constructions reflect contour-deformation
behavior.[^soper]

A general surface-normal map with uniform angles does not provide these
capabilities. A generic extension can use local coordinates (y,z,delta), where
y denotes previously generated external/complement data, z denotes coordinates
tangent to a target surface, and delta denotes its normal residual. Its density
has the normalized conditional form

\[
q(y,z,\delta)=q_y(y)q_z(z\mid y)q_\delta(\delta\mid z,y).
\]

Widths in the last factor may depend on geometric distances measured in z:
distance to a soft momentum, to a collinear locus, to another surface, or to a
near-tangent contact. The map remains triangular in this order. Its normalized
conditional PDFs, geometric Jacobian and inverse density are all part of the
resolved channel definition.

For a concrete family, take dimensionless normal distance delta and positive
dimensionless widths A(z), B(z). On -D<delta<D define

\[
K(\delta,z)=\frac{1}{(|\delta|+A)(|\delta|+B)},
\qquad q_\delta=K/I(z),
\]

\[
I(z)=\frac{2}{B-A}
\log\frac{B(D+A)}{A(D+B)}.
\]

At A=B, the continuous limit is
I=2(1/A-1/(D+A)). One may choose A=a omega(z)^p and B=b omega(z)^q for positive
a,b,p,q. To obtain a joint density proportional to a desired envelope
w(z)K(delta,z), the tangential marginal must be proportional to w(z)I(z),
provided its integral is finite. Normalizing the last conditional factor while
ignoring this change to the desired tangential marginal would produce a
different joint density. For two tangential dimensions and widths of orders
omega^2 and omega, I behaves as log(1/omega)/omega; the tangential radial measure
omega d omega makes this local enhancement integrable.

This is an example of a generic density family, not a universal prediction of
the powers in GammaLoop. The relevant widths must be determined from the actual
threshold-subtracted integrand and the geometry of its exact projections.
Adaptive angular grids are useful once the coordinates expose the structure;
they do not by themselves guarantee discovery of arbitrarily narrow correlated
regions.

The geometric inputs are topology-independent: routed momenta, energy
residuals, derivatives, soft loci q_e=0, and physically relevant collinear loci.
Local charts and multiple patches allow the same algorithms at arbitrary loop
order. Automatic enumeration and optimization of all overlapping soft,
collinear, threshold and projected configurations is a separate, difficult
classification problem. A generic library of maps is feasible without claiming
that a finite heuristic discovers the optimum for every graph.

## Three user-facing concepts

1. **Selection:** which named channels or automatic strategies are used for each
   existing graph.
2. **Definitions:** optional user-named channels, described by symbolic geometry
   expressions and ordinary runtime profile settings.
3. **Inspection:** display the resolved targets, coordinates, profiles, coverage,
   branch rules and reasons that an automatic channel was included.

The same settings structure should be used by TOML and the Python API. Graph
names are scoped to the selected integrand. Definitions do not activate channels
until selected. Explicit graph selections replace the default selection;
multiple entries in one selection list are combined and deduplicated.

```toml
[sampling]
sampling_multichanneling = true
sampling_channels = "monte_carlo"

default_channel_selection = ["auto:surfaces"]
channel_selection = { GL638 = ["auto:surfaces", "auto:optimized_lmb", "hard_intersection"], GL297 = ["auto:optimized_lmb"] }

[sampling.channel_definitions.GL638.hard_intersection]
around = "intersect(surface(2,4,12), surface(3,10,13))"
subspace_lmb = [3]
parent_lmb = [3, 4, 7, 10]
on_cut = [2,6,10]
```

`sampling_channels` selects summed versus Monte Carlo evaluation of the channel
axis; it does not select channel identities. The generalized
`sampling_channel_weight` owns the density-partition choice. These replace the
old `lmb_multichanneling`, `lmb_channels` and `lmb_channel_weight` names.

### Names and TOML syntax

`GL638` is the existing loaded graph name, not a newly invented alias in this
table. `hard_intersection` is an arbitrary channel name scoped to that graph.
The previous `HZ` name had the same status, and `G` was only a placeholder for
an actual graph name.

H and Z are research labels for two resolved energy surfaces in the prepared
cut (2,6,10): H=eta(eset(2,4,12)) and Z=eta(eset(3,10,13)). They are not reserved
API tokens, and Z does not signify a Z-boson propagator. The channel name may be
changed without changing those targets. In the refined vocabulary, `surface`
selects an energy surface; `cut` selects a physical Cutkosky cut. Edge lists
must be resolved with the correct energy shift, orientation and context;
ambiguities require explicit qualification and a diagnostic listing candidates.

TOML uses equals signs inside inline tables. Bare `GL638` and quoted `"GL638"`
are equivalent string keys. Colons belong to the Python/JSON dictionary form,
not the TOML form. Graph names containing punctuation outside the bare-key
alphabet must be quoted.[^toml]

```toml
channel_selection = { GL638 = ["auto:surfaces"] }
# Equivalent spelling:
# channel_selection = { "GL638" = ["auto:surfaces"] }
```

The Python dictionary corresponding to that setting is
`{"channel_selection": {"GL638": ["auto:surfaces"]}}`.

## Symbolica expression language

Use Symbolica parsing followed by structural pattern matching and typed argument
validation. Existing threshold-multiplier parsing and orientation matching
already demonstrate the relevant APIs.[^parser]

The proposed basic constructors are:

| Constructor | Meaning |
|---|---|
| `surface(edges...)` | A resolved energy-surface target |
| `cut(edges...)` | A physical Cutkosky-cut selector |
| `soft(edge)` | Focus around a routed massless momentum vanishing |
| `collinear(edge_a,edge_b)` | A physically qualified collinear target |
| `complement(edges...)` | Explicit remaining sampling coordinates |
| `product(maps...)` | Independent maps, validated in the current context |
| `intersect(targets...)` | Joint coordinates for compatible target constraints |
| `then(first,second,...)` | Ordered conditional composition |
| `phase_space(cut(...))` | Construct a physical cut's kinematic context |
| `left(target)`, `right(target)` | Resolve a target in one side of that context |

Each constructor has a typed role. A vector-valued soft constraint cannot
silently be treated as a single scalar normal; ranks and available dimensions
must be checked. Profiles for correlated tangential/normal sampling accompany
these geometric targets. Short forms receive deterministic defaults; the full
resolved specification remains exportable and editable.

Use ordinary, non-symmetric function heads in a dedicated namespace. Their
arguments retain ordering, whereas algebraic products permit reordering and
coalescing repeated factors.[^symbolica] This makes `product(...)` and
`then(...)` preferable to overloading arithmetic `*` or adding a separate
parser for `x` and `->`.

The compiler owns the patterns. It validates the complete expression, then
recursively extracts its ingredients; finding a valid nested subtree is
insufficient. User expressions remain data rather than becoming wildcard
patterns. Allowed heads, arities and integer arguments are checked first;
edge existence, duplicate constraints, unknown references, dimensions and ranks
are checked against the resolved graph and routing. The existing edge-set
diagnostics should be reused.

For a schematic cut C and catalogue threshold labels T_L and T_R, the earlier
cut/left/right example becomes

```text
then(phase_space(C), product(left(T_L), right(T_R)))
```

Here the labels are schematic resolved catalogue selectors, not arbitrary
undeclared symbols; inline edge-based selectors may replace them. The physical
cut data are prepared first. Internal left and right
cycles then preserve those cut momenta. The raw LU auxiliary scale and its
Jacobian complete the map. Selecting a host-cut channel never removes any term
from the physical cut sum. Amplitudes use the same surface/soft/collinear and
composition primitives with their fixed external data, without requiring a
physical-cut context.

## Automatic strategies and absent surfaces

| Selection | Proposed contract |
|---|---|
| `auto:lmb` | All admissible LMB channels in the generated catalogue, without heuristic pruning |
| `auto:optimized_lmb` | Current LMB-selection heuristics, supplemented by the agreed elementary-soft coverage audit |
| `auto:surfaces` | Optimized soft/domain coverage, augmented by useful surface, joint, projected, angular and cut/side channels |
| Explicit channel names | Exactly the selected channels, with domain coverage checked and soft coverage reported |

An ordinary LMB channel can be omitted from `auto:surfaces` when another
selected channel demonstrably supplies its coverage. Automatic mode should
report such decisions. Its ordinary-channel baseline uses the optimized
selection, not the unpruned catalogue. All presets and names combine by union:
adding `auto:optimized_lmb` can restore ordinary channels omitted from a surface
strategy, while adding `auto:lmb` requests the entire admissible catalogue.
Deduplication compares the resolved map, context, density profile, support and
branch rules. Different widths or different projections remain different
channels even when their LMB is identical. Users can choose only new channels through explicit
named selections; no ordinary channel may be silently appended to that exact
selection. Missing domain coverage must be diagnosed. Soft efficiency and
coverage of every compound singularity are stronger properties than basic
domain coverage and must be reported separately.

With no E-surfaces anywhere in the graph, `auto:surfaces` resolves to the
soft/LMB coverage. With no relevant surfaces and no loop-dependent massless
factors, it uses an ordinary full-support map. This behavior applies to
amplitudes as well as cross sections.

### Surfaces that exist only for some complement data

Keep a fixed channel catalogue. In each channel, generate external/cut and
complement variables y before the active variables x. Choose cut-preserving
cycles so the later left/right maps cannot change the prepared host-cut data.
Those data, including the necessary t-star rescaling, determine the energy
equation E_y(x). For an ordinary positive-energy convex residual with a nonzero
coercive active block after separating spectator directions, a strictly
negative minimum supplies a regular fiber with an interior center; a positive
minimum means no surface; a zero minimum identifies a pinched or degenerate
fiber. Arbitrary nonlinear star pullbacks require their separate chart checks
and cannot inherit this classification merely from the original surface.
Reuse the existing existence classifier where applicable, with an interior
certificate for the selected active space.

The initial implementation can use a deterministic conditional map:

\[
T(u;y)=\begin{cases}
T_{\rm shell}(u;y),&\text{regular fiber},\\
T_{\rm base}(u;y),&\text{absent fiber}.
\end{cases}
\]

Both maps must cover the entire active momentum space, not just a finite shell
window. For y=T_y(v), the combined map (v,u) to (y,T(u;y)) is block triangular:

\[
J=|\det D_vT_y|\,|\det D_uT(u;y)|.
\]

Dependence of the center, threshold radius and branch choice on y does not add
a missing diagonal factor. On the regular regions those derivatives occupy
the off-diagonal block; the existence boundary is handled piecewise, without
a distributional delta-function Jacobian. The cut/auxiliary-scale and native-
to-raw coordinate reconstruction factors remain part of their actual map
blocks and are included exactly once.

An explicit full-support shell map illustrates the construction. Compactify
r>=0 with g(r)=r/(s0+r), s0>0. For direction n set a=g(r_star(n,y)), with
0<a<1. Given uniform u in (0,1), choose p>1 and set

\[
s=\begin{cases}
a-a[(a-u)/a]^p,&u<a,\\
a+(1-a)[(u-a)/(1-a)]^p,&u\ge a,
\end{cases}
\qquad r=g^{-1}(s),\quad x=c(y)+rn.
\]

This maps onto the entire radial half-line and focuses the threshold at u=a;
p=2 gives the usual square-root singular density near the shell. The active
volume Jacobian is r^(d-1) (dr/du) J_Omega. Angular derivatives of r_star add
radial components which cancel from this determinant. The absent-fiber branch
uses r=g^{-1}(u), with its ordinary radial Jacobian and full support. Other
normalized profiles can be substituted without changing this architecture.

For a cut followed by left/right threshold maps, four combinations occur:
both regular, only left regular, only right regular, and neither regular.
Each side uses its applicable full-support map; all four combinations remain
valid contributions. Arbitrary projected targets involving other cut frames
may violate the assumed dependency order and must pass the generic conditional
or coupled-chart checks instead.

An optional smooth conditional mixture is

\[
q(y,x)=q_y(y)\left[a(y)q_{\rm surface}(x\mid y)
+(1-a(y))q_{\rm base}(x\mid y)\right].
\]

Both conditional PDFs are normalized. If the fiber is empty, set a(y)=0. Near
a shrinking fiber, a smooth transition can reduce the focused contribution.
Pinched fibers require their appropriate soft/collinear or degenerate chart;
they are not regular shells with arbitrarily small radii. An ordinary fallback
there preserves the reference integral but does not certify efficient sampling
of the physical singularity. Mixture-density weighting must use the actual
mixture. Alternatively an explicit partition over branch estimators may be
used with its corresponding probabilities; the two formulas must not be mixed.

Do not sample a completed point, delete the channels that fail there and
renormalize the remaining catalogue. Physical absence and root-solver failure
are different outcomes; a failed numerical solve cannot silently acquire the
meaning of an empty fiber. A declared deterministic fallback after a geometry
failure can preserve normalization if both generation and density evaluation
reproduce it, but must be reported as degraded sampling. The initial acceptance
implementation should fail such unexpected cases. Conservative activation
margins for sampling must not redefine physical threshold-CT existence.

## Elementary massless-factor coverage

For a loop-dependent massless edge with r=|q_e|, let p_r be the radial marginal
and q_cart the density with respect to d^3 q_e. With uniform solid angle,

\[
q_{\rm cart}=\frac{p_r(r)}{4\pi r^2}.
\]

Matching an isolated factor 1/E_e=1/r therefore requires p_r proportional to r,
or r proportional to sqrt(u) near zero. Ordinary spherical maps with a finite,
nonzero dr/du at zero instead yield q_cart proportional to r^(-2), which makes
the isolated weighted contribution vanish proportionally to r. GammaLoop's
current linear/log spherical maps have this latter property when the actual
edge momentum is the centered sampling coordinate.[^radial]

Automatic construction should enumerate every loop-dependent massless factor
and verify that at least one selected channel exposes its soft locus with a
suitable actual density and nonzero probability. It may add an LMB containing
that edge and complete the basis with massive edges, or add an equivalent
explicit soft map. For projected terms, the audit must use the actual prepared
or star momentum; raw edge zeros are not universally sufficient.

There is a concrete gap to address in the present heuristic: its amplitude
coverage universe is built from massless-edge combinations of length equal to
the entire loop count. If no fully massless basis exists, it falls back to one
basis. This does not establish individual massless-edge coverage when a basis
also needs massive edges.[^lmb]

Coverage of simultaneous soft factors must be assessed jointly.
Separate channels for each edge do not generally bound their product.
Dependent constraints and collinear hierarchies require rank-aware geometry
and correlated profiles. Guaranteeing coverage of individual 1/E factors is
therefore a practical first contract; bounded weights for every compound
singularity need further analysis. A nonintegrable uncancelled singularity
cannot be cured by a normalized proposal density.

The current OSE channel partition is not itself a normalized sampling PDF.
Soft-coverage claims must be checked against the generated map and its
pushforward density, rather than inferred from an OSE weighting exponent.

## Specialized eager evaluators and Jacobians

The resolved channel should construct reusable Symbolica eager evaluators during
warmup for its explicit algebraic maps and implicit residual equations. Graph
routing and expression structure specialize the programs; current kinematics,
widths and solved roots remain runtime inputs. Existing GenericEvaluator and
its dual-evaluator construction provide the appropriate infrastructure.[^dual]
The first-order dual shape has D+1 coefficients for D differentiated inputs;
the Jacobian does not require all mixed higher derivatives.

Seed the independent map coordinates with first-order dual components, evaluate
the output vector, and extract its derivative matrix. The volume Jacobian is
the absolute determinant, computed numerically with a stable matrix
factorization. Block-triangular maps permit a product of smaller determinants.
For y=T_y(v), x=T_x(u;y),

\[
|\det DT|=|\det D_vT_y|\,|\det D_uT_x|.
\]

Consequently, dependence of the active map's center on previously generated y
does not require center derivatives for this determinant. Dependence on the
active coordinates is different and cannot be discarded.

A numerical root is not automatically differentiated merely because its value
is passed to an evaluator. For G(z,u;y)=0, solve for z at the sample point,
evaluate G_z and G_u with dual numbers, and solve

\[
G_z\frac{\partial z}{\partial u}=-G_u.
\]

These derivatives then propagate into later blocks. The scalar case supplies
the derivative of a threshold radius or LU root. This differentiates the
defining equation rather than an arbitrarily truncated solver iteration.

Runtime SOCP centers require the same discipline. Use complement-dependent
centers and triangular determinants where valid. A center that genuinely
depends on differentiated active coordinates requires the appropriate
sensitivities, including branch or active-set handling. An eager evaluator is
not a derivative implementation for a black-box optimizer.

Piecewise profiles and branches should have explicit boundaries. The generated
function vocabulary must be checked against dual-evaluator support; for example,
an abs derivative is not defined at zero. Automatic differentiation supplies
local derivatives. It does not establish density normalization, support,
global injectivity, or discrete branch probabilities.

## Invertibility and geometric-score partitions

Reliable inverses are realistic for the main map families. LMB changes are
linear; radial profiles are monotone; a regular convex shell with an interior
center has a unique radial root; the demonstrated A-star GL638 map has an
affine pullback. Conditional maps recover the complement coordinates before
the dependent active coordinates.

For general intersections, choose geometric coordinates
(S_1(K),...,S_m(K),K_free). Evaluating the inverse geometric coordinates at a
given K is direct. Constructing K from prescribed residuals and free coordinates
requires the implicit solve. Full rank guarantees only local invertibility:
the forward branch and its domain must be controlled by a chart. A successful
root solve alone does not prove global uniqueness. Explicit patches, branch
identities and support predicates belong to the resolved map.

Exact inverses of every foreign channel are nevertheless not necessary for
unbiased multichanneling. Let each chart a be a one-to-one map K=T_a(u) with
true forward volume Jacobian J_a. Choose nonnegative scores rho_a(K), zero
outside that chart's image, with a positive sum wherever F requires coverage:

\[
w_a(K)=\frac{\rho_a(K)}{\sum_b\rho_b(K)},\qquad
I=\sum_a\int du\,F(T_a(u))w_a(T_a(u))J_a(u).
\]

Sampling a with probability pi_a and its cube coordinate with adaptive PDF g_a
gives the unbiased weight

\[
\frac{F(T_a(u))w_a(T_a(u))J_a(u)}{\pi_a g_a(u)}.
\]

The existing outer integration weights must be applied exactly once. The same
partition applies to the complete cancellation-bearing physical integrand,
including its requested cut sum, outside the LU residue differentiation.

The scores can be exact unadapted map densities or suitable geometric proxies;
they need not be normalized. They must be functions of the same raw point K,
independent of which channel generated it. Evaluating foreign Jacobians at the
selected channel's cube coordinates is not generally such a function and does
not establish this partition identity.

This is distinct from the mixture-density estimator F/(sum_a pi_a q_a), whose
denominator must contain the actual proposal densities. Substituting proxies
into that denominator would generally bias the result. Nor can a proxy replace
the selected map's true forward Jacobian. The fixed-partition organization is
consistent with the alternative discussed by Ohl.[^ohl]

Unknown multiplicity is a separate obstacle. If one map covers a point more
than once, its forward-Jacobian integral counts every preimage. Split controlled
branches into chart components and partition over those components, or supply
the exact multiplicity/pushforward correction. Arbitrary scores do not repair
an uncontrolled many-to-one map. Support predicates remain necessary even
when foreign inverse coordinates are not computed.

Proxy partitions preserve unbiasedness, but bounded weights require appropriate
singular behavior. A background channel assigned a finite fraction of a 1/R
singularity still generates unbounded weights, even if a focused channel is
also present. A sufficient local condition is

\[
c\,q_a^0(K)\le\rho_a(K)\le C\,q_a^0(K),\qquad 0<c\le C<\infty,
\]

for every relevant chart, where q_a^0 is its unadapted pushforward density.
Then the proxy-partition contribution differs from the exact-density partition
by at most C/c. Adaptive and channel probabilities also need the appropriate
positive bounds. Exact score normalization is unnecessary, but singular powers,
support and relative behavior matter.

### Proxy construction, asymptotic matching and positivity

Every proxy must capture the features its sampling density was designed for:
surface normals, angular/soft/collinear enhancements, correlated shrinking
widths, complement and tail factors, and the correct prepared or projected
frame. Conditional normalization factors can themselves carry singular powers
and must not be discarded without justification. Radial marginals and densities
with respect to Cartesian volume must be distinguished. A factor may be
approximated by a constant only where it is controlled above and below by
positive finite bounds.

Matching the actual powers is the default. Independently making a proxy more
singular than its own proposal is not sufficient for good weights: it can
allocate a singular region to a channel whose actual sampling is too weak.
For example, on (0,1), take normalized densities

\[
F=q_b=\tfrac14x^{-3/4},\qquad q_a=\tfrac34x^{-1/4}.
\]

The exact-density partition has F/(q_a+q_b)<=1. But positive proxies
rho_a=(3/4)x^(-7/8), rho_b=q_b satisfy the condition that each power is equal
or stronger than its density, while assigning w_a tending to one at zero.
Channel a then has weight asymptotic to (1/3)x^(-1/2), and its second moment
diverges as the integral of x^(-5/4). Stronger proxies are allowed only after
checking the complete partition, or strengthening the corresponding proposal
as well. One common positive factor h(K) multiplying every proxy cancels; the
comparability condition above may therefore be imposed on rho_a/h instead.

Scores should be nonnegative and strictly positive on the interiors of their
declared supports. Domain coverage must ensure that their sum is positive at
every ordinary point being integrated. A selected full-support component can
provide this guarantee. A local chart must still have zero score outside its
image; adding an unconditional floor to that chart would invalidate the
partitioned-integral identity. Positive scales, widths, coefficients, squared
norms and explicitly gated supports provide a controlled construction, rather
than relying on a general positivity proof for arbitrary expressions.

The resolved channel compiler should construct map kernels and proxy kernels
from one density-factor description, so changing a channel's profile or
conditional branch also changes its score consistently. Scalar score programs
reuse the same eager evaluator owner; they do not need dual components unless
a derivative is required to construct a geometric input. Shared energy and
projection expressions can be evaluated jointly. Runtime kinematics and branch
data are explicit inputs to the programs built during warmup.

For numerical normalization, compile log-scores where useful and use a scaled
sum/log-sum-exp reduction across the active supports. This preserves the
intended powers while avoiding overflow and underflow. A zero denominator at
an ordinary point is a coverage or evaluation error. Several infinite scores
at an exact singular endpoint need explicit limiting/branch treatment;
log-sum-exp alone cannot resolve infinity minus infinity. Score capping or
epsilon smoothing must not silently flatten a required enhancement.

Acceptance tests must verify positive denominators and partitions summing to
one, including absent-surface branches, extreme scales and competing singular
limits. Separate scaling tests check that proxies preserve the intended
bounded-weight or finite-variance behavior; integrating a reference density to
one only tests unbiasedness and cannot establish those stronger properties.

The implementation plan is therefore one chart/partition mechanism with exact
selected-map Jacobians, exact inverse densities where practical, and an
explicit geometric-proxy score option when foreign inversion is costly. The
proxy mode does not relax branch correctness. Validation includes nonconstant
analytic moments, forward/inverse and Jacobian checks, overlap multiplicities,
nonuniform discrete/adaptive probabilities, and singular-weight scaling.

For a simple nonconstant check, take T_0(u)=u, T_1(u)=u^2 on [0,1],
F(x)=1+x^2, and unnormalized proxy scores rho_0=1+x, rho_1=2-x. Using the
true forward Jacobians 1 and 2u with their proxy partition gives 4/3. A
40,000-point midpoint quadrature in each chart gives 1.3333333332899304,
an absolute error of 4.34e-11. This verifies the elementary partition example,
not the proposed multiloop map implementation.

## Acceptance harness using actual saved states

The acceptance harness should load an actual process/integrand state and use
its real graph-dependent map definitions, channel selection, routing, cut
preparation, centers and projected-target contexts. It runs the production
integration machinery with a temporary reference-function overlay. The saved
production state remains unchanged; test results use a separate workspace.

The substitution occurs at the raw physical-integrand boundary. In schematic
order:

1. Select the real graph/orientation/channel and construct the actual map.
2. Obtain the raw point K, including all geometry needed by its map.
3. Substitute a known normalized function for the physical scalar F(K).
4. Retain the true selected-map Jacobian, common channel partition, discrete
   probabilities and adaptive-grid corrections.

Substituting the entire final event weight would bypass the quantities under
test. Inserting a Gaussian numerator inside each physical cut would instead
retain inappropriate LU/residue factors and count the target repeatedly. The
reference measure is plain d^(3L)K: physical flux, symmetry, 2pi, observable and
residue factors are not part of it. Physical event selectors and acceptance
cuts are bypassed, while geometric chart-support predicates and fallback
decisions remain active. Actual geometry needed for the maps is retained even
when the corresponding physical evaluator is bypassed.

The existing UnitVolumeIntegrand replaces the process and uses only the global
parameterization, so it cannot provide this graph-dependent acceptance test.
The correct owners to extend are the process evaluation boundary and the
existing integrate request/warmup path. The outer integrator already owns
parameterization-Jacobian and adaptive/discrete weights.[^acceptance]

A possible request interface, still proposed, is

```text
integrate -p PROCESS -i INTEGRAND --sampling-test gaussian
```

The Python API should expose the same request. A structured test profile can
specify scale, shift, covariance, precision budget and seed without modifying
the loaded integrand's production settings. The known target is supplied
automatically. Standard execution options continue to select resources and
sampling settings.

### Reference functions and discrete normalization

In the fixed canonical raw/master LMB, with D=3L, use

\[
G_{\mu,\Sigma}(K)=
\frac{\exp[-\tfrac12(K-\mu)^T\Sigma^{-1}(K-\mu)]}
{(2\pi)^{D/2}\sqrt{\det\Sigma}},
\qquad \int_{\mathbb R^D}G\,d^D K=1.
\]

The default covariance is sigma^2 times the identity, with sigma comparable to
the positive characteristic energy scale, normally sqrt(s). An explicit
positive scale is available for amplitudes without a suitable physical sqrt(s).
The reference is evaluated independently in one common raw frame, never in a
different native frame for each channel.

When testing several graphs or orientations together, use known coefficients
b_(g,o)>0 with sum b_(g,o)=1 and target
F_ref(g,o,K)=b_(g,o) G_(D_g)(K). These coefficients are independent of adaptive
selection probabilities. Each label can also be tested separately. The
physical cut and CT multiplicities do not introduce extra copies of G.
Individual partitioned channel contributions need not integrate to one;
their total does.

Include centered Gaussians at several widths, shifted anisotropic Gaussians
with cross-loop covariance, and known moments. The first moment is mu and the
centered second-moment matrix is Sigma. These probes supplement normalization
and target the shared-frame transformations. A normalized Student-t reference
can additionally exercise tails and support beyond focused windows. Its
statistical behavior must be assessed with the actual proposal; full support
alone does not ensure a finite or manageable variance.

Existing rotation-stability operations must use a consistently transformed
reference mean/covariance or evaluate the reference in the canonical unrotated
frame. Otherwise anisotropic probes would trigger artificial instability.

### Systematic coverage and pass criteria

The harness should accept actual GL297, GL638 and other saved states, with both
amplitude and cross-section cases. Its case descriptions cover complete maps
alone, valid mixtures of partial-support patches, summed/Monte-Carlo channel
selection, exact/proxy partitions, automatic/manual catalogues, unequal channel
probabilities, adaptive grids, and direct/projected targets. Geometric tests
independently verify roots, selected target identities, finite-difference versus
dual Jacobians, inverse coordinates, support and branch multiplicities.

Record visits to focused and fallback branches, present/absent/pinched fibers,
unexpected solver failures, inverse residuals and weight extremes. A successful
run that never uses its focused branch does not validate that branch. Actual
states should be complemented with deliberately constructed boundary cases.

For a controlled unit fixture, take

\[
E_y(x)=2\sqrt{|x|^2+m^2}
-[2m+\Delta\tanh(y/R)],\qquad 0<\Delta<2m.
\]

This synthetic conditional residual has a regular shell for y>0, no shell for
y<0, and a pinch at zero. With a centered product Gaussian reference, the
total integral is one and each signed-y contribution is one half. Two
independent conditioning variables reproduce the four left/right existence
combinations, each with reference integral one quarter. This fixture tests
conditional transport directly; actual-graph state tests then validate the
production context and target resolution.

Acceptance requires both statistical agreement with the known result and an
uncertainty below a declared precision budget. Report the estimate, standard
error and absolute pull |estimate-target|/standard_error. Fixed seeds make
failures reproducible; independent validation seeds check adaptation artifacts.
Insufficient precision is inconclusive, not a pass. Structural failures remain
failures regardless of the central value. Statistical thresholds must account
for the number of cases; a unit-normalization result alone is not proof of
injectivity, correct target geometry or bounded physical-integrand weights.

## Sources

[^soper]: Davison E. Soper, [Choosing integration points for QCD calculations by numerical integration](https://arxiv.org/pdf/hep-ph/0103262), 2001, Sections V–IX.
[^toml]: [TOML v1.0.0 specification](https://toml.io/en/v1.0.0#inline-table), inline tables and keys.
[^symbolica]: [Symbolica pattern-matching documentation](https://symbolica.io/docs/pattern_matching.html), structural matching and ordering of function arguments.
[^parser]: GammaLoop [threshold expression parser](../../../crates/gammalooprs/src/integrands/process/threshold_multiplier.rs), [edge-set validation](../../../crates/gammalooprs/src/integrands/process/threshold_multiplier.rs), and [orientation matching](../../../crates/gammalooprs/src/settings/global.rs).
[^radial]: GammaLoop [three-dimensional parameterization](../../../crates/gammalooprs/src/utils/mod.rs).
[^lmb]: GammaLoop [amplitude LMB selection](../../../crates/gammalooprs/src/graph/mod.rs).
[^dual]: GammaLoop [GenericEvaluator](../../../crates/gammalooprs/src/integrands/process/evaluators.rs), [dual construction](../../../crates/gammalooprs/src/integrands/process/evaluators.rs), and pinned Numerica [first-derivative dual shape](/common/dev/symbolica_gl638/lib/numerica/src/domains/dual.rs:1406).
[^ohl]: Thorsten Ohl, [Vegas revisited: Adaptive Monte Carlo integration beyond factorization](https://arxiv.org/pdf/hep-ph/9806432), 1998, Section 4.3.
[^acceptance]: GammaLoop [existing UnitVolumeIntegrand](../../../crates/gammalooprs/src/integrands/mod.rs), [integrate warmup](../../../crates/gammaloop-api/src/commands/integrate.rs), and [outer sampling-weight application](../../../crates/gammalooprs/src/integrate/mod.rs).
