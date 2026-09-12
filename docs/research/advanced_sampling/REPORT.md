# E-surface-adapted sampling for GammaLoop

E-surface-adapted channels are mathematically sound and can be integrated into
GammaLoop's existing multichannel machinery. Their natural building blocks are
resolved energy surfaces, chosen loop-coordinate subspaces, conditional radial
maps, and normalized one-dimensional density profiles. The same framework can
include ordinary LMB channels and more advanced intersection coordinates.

Three qualifications determine a correct design. A surface is spherical after a
nonlinear change of variables only in an appropriate nondegenerate active space.
A Jacobian exactly proportional to distance from an entire non-pinched threshold
would generate a nonnormalizable density; fractional powers are needed there.
Finally, separate subspaces do not by themselves imply independent maps: their
external data can couple them. These restrictions do not defeat the proposal.

For GL638, the idea supplies concrete candidates for the measured hard H/Z
corner, including an unexpectedly simple exact channel for its native A-star
image. One fractional-power shell can make the local second moment finite.
Two fractional-power factors in a joint normal chart can also compensate the
leading large-weight growth. Neither statement establishes global convergence
of the complete graph without further angular, hierarchical and numerical tests.

| Question | Assessment |
|---|---|
| Can an E-surface be made spherical? | Yes for a bounded convex fiber with nonempty interior, after removing spectator directions; generally through a nonlinear radial map. |
| Can the forward Jacobian vanish exactly linearly across the whole shell? | Not as an ordinary normalized finite-probability map. A power strictly below one is allowed. |
| Can several surface maps be composed? | Yes with independent blocks, a valid conditional order, or a joint chart. Distinct edge lists or subspace names alone are insufficient. |
| Can ordinary and new channels coexist? | Yes. Each must have a defined inverse density in a common integration frame, and channel selection must retain the physical cut sum. |
| Does this help GL638? | There are explicit locally sufficient H, A-star H and joint H/Z constructions. Performance and global coverage remain to be demonstrated. |

## 1. Geometry: what becomes a sphere

Fix complement coordinates and external data, denoted by y. Let the active
Cartesian loop coordinates be x in R^d, with d=3l for l active loops. A standard
positive-energy threshold has the form

\[
 E_y(x)=\sum_e\sqrt{\lVert A_ex+b_e(y)\rVert^2+m_e^2}-Q(y).
\]

This is convex. If the stacked active routing map has no kernel, the positive
energy sum grows without bound as the active coordinates grow. If its minimum
is strictly negative, the region E_y<0 has a nonempty bounded interior.
Capatti's thesis develops precisely this threshold geometry and its radial
parameterization in Section 5.1.2, especially equations (5.51)–(5.69).[^zeno]

Choose any interior point c(y). For every unit direction n there is one positive
root r*(n,y) satisfying E_y(c+r*n)=0. Define

\[
 x=c(y)+\rho\,r_*(n,y)n,\qquad \rho\ge0.
\]

The threshold is now **rho=1**. Its enclosed region is the unit ball in these
new coordinates. The construction is generally nonlinear; it does not imply
that every multi-loop energy equation is a quadratic ellipsoid or that a linear
LMB transformation alone will make it spherical.

The exact active volume element is

\[
 d^dx=r_*^d\rho^{d-1}\,d\rho\,d\Omega_{d-1}.
\]

Angular derivatives of r* do not introduce an extra determinant factor: those
terms point in the radial direction, already present as the radial column.
Angular unit-cube coordinates still contribute their ordinary sphere Jacobian.
If c and r* depend only on the already sampled complement, those derivatives
also occupy off-diagonal blocks of the full Jacobian. A center that changes
with the active coordinates requires a different calculation.

At a regular root, E_y and r−r* are locally proportional with nonzero factor
n·gradient(E_y). They are not generally proportional far from the threshold.
Energy distance, radial distance and geometric normal distance should therefore
be distinguished in the API, including their units and width conventions.

The qualifications are substantive:

* A surface independent of some chosen variables is cylindrical in those
  directions. For example, GL638's A(p,s)=2E(s)−Q in its native two-loop CT
  space is R^3_p times a two-sphere, not a five-sphere. Its spherical sampling
  block should vary s and leave p as a complement.
* A pinched threshold has no strict interior; the above construction does not
  apply. Point and collinear limits require different charts or a declared
  fallback. A finite but zero enclosed volume is not sufficient.
* Massless zero-energy points can make the map nonsmooth. Convex radial
  parameterization remains useful almost everywhere, but smooth Jacobian and
  conditioning arguments must exclude or explicitly handle those loci.
* After nonlinear LU or CT-star projection, a target equation need not remain
  convex in the original integration variables. Its actual pullback must be
  examined; the original E-surface's convexity is insufficient.

A sampler center can differ from every subtraction center. It is a coordinate
choice for the unchanged integral. Matching a particular CT center becomes
necessary only when the intended target explicitly means that CT's star image.
Sampling groups must not alter threshold group IDs, common centers, maximal
overlaps, e-surface multichanneling or dual cancellation.

## 2. Density: the correction to a linearly vanishing Jacobian

Near a regular non-pinched shell, set delta=r−r*. If a map from a finite uniform
coordinate u obeyed

\[
 \left|\frac{dr}{du}\right|\propto|\delta|,
\]

then its conditional radial density would satisfy p_r proportional to
1/|delta|. Its integral diverges logarithmically on either side of the shell.
Equivalently, integrating du proportional to d(delta)/|delta| takes an infinite
coordinate interval to reach delta=0. The finite volume enclosed by the
surface does not remove this transverse divergence.

A useful normalized replacement is

\[
 p(\delta)=\frac{1-\beta}{a^{1-\beta}+b^{1-\beta}}
 |\delta|^{-\beta},\qquad -a<\delta<b,\quad0<\beta<1.
\]

The negative branch has probability a^(1−beta)/(a^(1−beta)+b^(1−beta));
conditional on it, sample delta=−a U^(1/(1−beta)). The positive branch is
analogous. In particular, **a quadratic map gives beta=1/2**, so its forward
radial Jacobian vanishes as sqrt(|delta|). This is a natural first profile.
An inner window must be clipped at r=0 and its branch probability recomputed;
equal branch probabilities for unequal intervals would not generate this stated
common-coefficient density. They can define a different valid profile when their
actual conditional densities are accounted for.

A compact focused channel can be added to ordinary full-support channels.
Alternatively, a single channel can contain a normalized focused branch and
an ordinary tail. The mixture and all inverse branches must be included in its
density. If a surface does not exist for a complement, a deterministic ordinary
conditional map is valid when both forward and inverse use that same decision.
Dropping the point, retrying complements until a surface exists, or returning
zero on a numerical failure would require additional normalization and can bias
the integral.

A finite-width profile 1/(|delta|+epsilon) is also possible. It helps resolve a
finite peak, but becomes smooth below epsilon and cannot by itself establish
an asymptotic variance cure. The analogous regulated shell focusing appears
explicitly in Soper's sampling construction.[^soper2001]

### One normal versus two normals

Suppose a regular intersection has scalar normal coordinates u,v, radius
R=sqrt(u^2+v^2), a bounded smooth measure factor, and f approximately C/R with
bounded angular and tangential coefficient C. Then

\[
 \int |f|\,R\,dR\,d\theta<\infty,
 \qquad
 \int |f|^2\,R\,dR\,d\theta\sim\int\frac{dR}{R}.
\]

The density power forbidden across a whole codimension-one shell is allowed
around a codimension-two intersection. Uniform R and angle give normal-plane
density 1/(2 pi Rmax R) on a disk. This compensates the leading 1/R growth.

There is also a particularly useful realization of the product-density idea:

\[
 q(u,v)\propto\frac1{\sqrt{|u|}\sqrt{|v|}}
 =\frac1{R\sqrt{|\cos\theta\sin\theta|}}.
\]

Each normal profile is separately normalizable. Together they yield bounded
leading weights proportional to sqrt(|cos(theta) sin(theta)|). This does not
require isotropic polar sampling, but it does require genuine joint normal
coordinates and their Jacobian. Simply adding two half-power channels is an
additive mixture, not this product, and retains weaker radial focusing.

| Local proposal for f approximately 1/R | Local second moment | Leading generic weight |
|---|---|---|
| Smooth ordinary proposal | Logarithmically divergent | R^-1 |
| One normal, q proportional to \|u\|^-1/2 | Finite | R^-1/2 |
| Two joint normals, q proportional to \|uv\|^-1/2 | Finite | Bounded |
| Joint polar proposal, q proportional to 1/R | Finite | Bounded |

These are local sufficiency results. Tangencies, vanishing coefficients in the
coordinate transformation, extra angular poles and unbounded complement tails
must be assessed separately.

## 3. Relation to Soper and adaptive multichannel integration

Soper's 1999 paper gives the general numerical setup and explains why sampling
must follow the residual singular geometry after cancellation. Its Section VI
discusses normalized mixtures and correlated soft/angular enhancement.[^soper1999]

The closest direct predecessor is **Soper, “Choosing integration points for QCD
calculations by numerical integration” (2001)**. It constructs explicit
ellipsoidal coordinates with the target at A+=S+, a normalized regulated shell
profile in equations (39)–(42), and a correlated enhancement at the special
soft point in equations (49)–(54). The reciprocal whole-shell factor has a
nonzero width for a non-pinched ellipsoid. Genuine singular focusing correlates
the shell and tangential distances. These are massless one-loop examples, not
a ready-made solution for every massive multi-loop or projected surface.[^soper2001]

Kleiss–Pittau provides the normalized mixture and channel-probability optimization
contract. Ohl extends adaptive grids across different coordinate charts: exact
adaptive mixtures need inverse maps and foreign-grid density evaluations.
His Section 4.3 also describes the cheaper fixed-map partition with independent
grid adaptation, which is a suitable first integration route here. His normalized
spherical/cylindrical/Cartesian benchmark is a useful precedent for independent
validation.[^kleiss][^ohl]

The application to GammaLoop below is an implementation proposal and the
GL638 conclusions are derived from its current geometry and recorded evidence,
rather than claims made in those papers.

## 4. Composing surfaces and complements

The expression

```text
cut(3,4,7) x cut(2,11) x complement(4,7)
```

is a good compact description after its meanings have been resolved. Surface
edges are the edges entering an energy equation; they are not automatically
the active LMB edges. Repetition of an edge between a surface list and the
complement is not by itself an error. The checks must use momentum signatures
and cycle ranks.

The canonical definition needs a master graph, a complete parent LMB, each
active block, and the actual energy equation, including its external shift.
It must also identify the evaluation frame: raw integration point, another
cut's LU-prepared point, or a specified CT-star point. Defaults can resolve
unique cases; ambiguity must list the competing equations, cuts, sides,
orientations, variants and subspaces. The complementary LMB edges can usually
be inferred from the complete parent and active blocks.

There are three different composition operations:

| Operation | Valid geometry | Jacobian and validation |
|---|---|---|
| Product | Disjoint coordinate blocks whose maps depend only on their block and an already fixed common complement | Product of block Jacobians; verify that one block does not move another block's target data. |
| Conditional | A directed acyclic dependency order: complements first, then each block using only previously fixed variables | Triangular Jacobian; order and conditional normalization are part of the definition. |
| Intersection | Several independent scalar surface normals within a common active block | Full-rank normal map, tangent coordinates, inverse branches and full Jacobian. |

“Different subspaces” is weaker than the product condition. A threshold in one
block can depend on momenta from another through its effective external data.
Mutually dependent radial maps require a joint construction or a clean
unsupported-dependency diagnostic. Repeatedly mapping one surface and then
another while multiplying guessed Jacobians is not sufficient.

A target named after a physical cut must still sample the complete requested
physical integrand. It must not cause evaluation of only that cut. Soper's
common-momentum cut sum is essential to the cancellation strategy.[^soper1998]
The same rule applies to all orientations requested by the runtime settings.

## 5. Proposed user API

All fields and commands introduced in this section are proposals. The current
repository does not implement them. The preferred design is one structured
runtime schema, shared by TOML, Python, validation, display and persistence.
Small readable shorthand should compile to that schema, rather than define
a separate sampler language with different semantics.

### A single useful GL638 channel

```toml
[sampling]
graphs = "monte_carlo"
orientations = "summed"
channel_mode = "monte_carlo"
channel_partition = "map_density"
lmb_basis_ids = { GL638 = [1, 2, 13, 14, 21, 22, 84] }
coordinate_system = "spherical"
mapping = "linear"
b = 0.3

[sampling.channels.GL638.H_cut1]
parent_lmb = [3, 6, 7, 10]
map = "radial"
lmb = [3]
target = { kind = "energy_surface", edges = [2, 4, 12], at = "lu_cut", host_cut = [2, 6, 10] }
profile = { type = "threshold_power", power = 2.0, distance = "radial", width_gev = 10.0, support = "window" }
```

This appends a named channel to the ordinary LMB channels. The active direction
is p; the remaining parent edges `[6,7,10]` define the complement. The host
cut's LU rescaling is solved from that complement before constructing the H
fiber. The 10-GeV width is an explicit trial parameter, not an established
optimal value. The sampler can choose its own interior center; H has an
analytic one. `distance="radial"` measures distance in the target's LU-prepared
frame; a CT-star target instead measures it in that star frame. The raw-coordinate
width follows from the corresponding transformation. GeV alone does not identify
which of these distances is intended.

The explicit first support contract is a normalized finite window: inner extent
min(width,r*) and outer extent width, with the sign probabilities from Section 2.
The focused channel has zero density outside that fiber's window. Retained
ordinary LMB channels supply global coverage and must keep positive sampling
probability. An absent or pinched fiber uses the deterministic ordinary conditional
branch described earlier. A full-support internal tail is a later alternative
requiring its own law and branch probability; it is not implicit in this example.

`channel_mode` generalizes the current `lmb_multichanneling`/`lmb_channels`
controls, and `channel_partition` generalizes weighting. The implementation
should migrate these concepts through their existing owners rather than retain
two competing controls for the same channel axis. The ordinary LMB list remains
a concise way to request existing channels. Existing cards and saved workspace
signatures would need an explicit migration policy.

The canonical chart always contains a full parent LMB. Reusable card blocks can
avoid repeating it; canonical exported settings should print it explicitly.
The existing selected-edge plus signed-cycle definition supplies the active
subspace identity. Sampler channel names and CT `group_id` values are separate.

A surface selector is independent of the CT catalogue. In particular, GL638's
useful H surface is another physical cut evaluated on cut 1 and has no ordinary
cut-1 threshold CT entry. A resolver limited to generated CT metadata would
exclude the most useful first target.

### Several blocks and an inferred complement

```toml
[sampling.channels.G.two_shells]
parent_lmb = [1, 2, 4, 7]
map = "product"
blocks = [
  { target = { kind = "energy_surface", edges = [3,4,7] }, lmb = [1] },
  { target = { kind = "energy_surface", edges = [2,11] }, lmb = [2] },
]
profile = { type = "threshold_power", power = 2.0, distance = "radial", width_gev = 10.0, support = "window" }
```

This illustrates the schema, not a verified assignment for a particular graph.
The resolver infers complement `[4,7]` and validates the proposed product.
An optional explicit complement acts as an assertion. The resolved display can
then render the original shorthand with active-edge annotations. `map="conditional"`
uses an explicit validated block order instead of claiming independence.

For H and Z, which share the same p block, the definition is different:

```toml
[sampling.channels.GL638.HZ_cut1]
parent_lmb = [3, 6, 7, 10]
map = "intersection"
lmb = [3]
targets = [
  { kind = "energy_surface", edges = [2,4,12], at = "lu_cut", host_cut = [2,6,10] },
  { kind = "energy_surface", edges = [3,10,13], at = "lu_cut", host_cut = [2,6,10] },
]
profile = { type = "product_threshold_power", powers = [2.0,2.0], distance = "energy", widths_gev = [10.0,10.0], support = "window" }
```

The active block has three scalar coordinates: two normals and one tangent,
not two three-vector radii. The displayed widths are energy residuals of the two
target equations, specified as requested upper bounds; the actual conditional
support must be certified feasible and reported. A polar-normal
profile can be an alternative on this same geometry. Coordinate construction
and density profile are deliberately separate concepts.

### Exact projected targets

An advanced target qualifier can select a particular CT projection:

```toml
target = { kind = "energy_surface", edges = [2,4,12], at = { ct_star = { host_cut = [2,6,10], threshold = [7,8], subspace_lmb = [3,7], centers = "all" } } }
```

Here the chart's active `[3]` block differs from the projection's native `[3,7]`
subspace. The shared parent resolves both. Any remaining variant ambiguity
requires an explicit variant selector. The geometry must prove that this
pullback admits the requested radial or intersection chart; a frame label alone
does not provide an inverse.

When several overlap centers contribute, represent them as normalized internal
branches with the sum of all inverse-branch densities. A complement-dependent
number of centers must not change the top-level adaptive grid's dimension or
channel count mid-run. Exact targeting should share the actual prepared center
context with the evaluator; independent SOCP solutions differing by numerical
noise can misalign the singular density at sufficiently small distances.

### Python and inspection

Python should accept the same dictionary shape through the existing transactional
settings owner. A proposed convenience is

```python
gl.set_integrand_settings(
    process="epem_a_tth", integrand="NNLO",
    settings={"sampling": {"channels": {"GL638": {"H_cut1": h_channel}}}},
)
```

The existing `gl.run(...)` method already permits the same settings through CLI
commands; a new fluent object hierarchy is not required. Optional builders can
later generate these dictionaries. Named processes, integrands and channels
should be first-class inputs and stable replay labels.

Proposed `display sampling --resolved` should print the master graph, target
equation and frame, full parent, signed cycles, complement, dependency order,
center rule, support and profile. `inspect sampling --channel NAME` should
show forward/inverse coordinates, root residual, full Jacobian, every channel
density and the partition sum. Definition errors must be transactional and
exhaustive. A malformed or unsupported channel must leave the integrand's
previous settings intact.

## 6. Integration into the current implementation

The source audit was performed at `451a36212` on
`ir_safe_threshold_subtraction`. The existing pieces are substantial, but the
extension is more than another `coordinate_system` enum value.[^code]

| Existing owner | Role in the extension |
|---|---|
| `SamplingSettingsParser`, runtime settings | Named channel schema, shared serialization, validation and defaults. |
| `LmbMultiChannelingSetup` | Generalize the existing resolved-channel owner to ordinary LMBs and surface charts. |
| `SubspaceData` | Reuse parent routing, selected edges, signed cycles and complement transformations. |
| `Esurface` and radial/existence routines | Reuse exact energy equations, shifts, values, radial derivatives and regularity classification. |
| Graph/channel selection and sample preparation | Construct graph-aware maps after the target graph/channel is known. |
| Existing continuous and discrete grids | Keep adaptation and the graph/orientation axes; generalize channel labels and replay identities. |
| LU and threshold geometry preparation | Supply correct native cut data and exact star contexts without changing the subtraction algebra. |

### The estimator boundary

Let T_a map a unit cube into the same raw graph integration variables K, with
Jacobian J_a. For a one-to-one chart with adaptive cube density g_a,

\[
 q_a(K)=g_a(T_a^{-1}K)/|J_a(T_a^{-1}K)|.
\]

Multiple inverse branches require their probability-weighted sum. In a full
adaptive mixture with probabilities pi_a, the event weight is
F(K)/(sum_a pi_a q_a(K)), with F the complete requested subtracted cut sum.
Foreign inverse maps, branch probabilities and adaptive-grid densities must
all be available to implement this formula.

A simpler first route fits the current engine. Let q_a^0 be each channel's
normalized density for an unadapted cube, including its branches and support.
Use the common partition

\[
 W_a(K)=\frac{q_a^0(K)}{\sum_bq_b^0(K)},\qquad
 I=\sum_a\mathbb{E}_{U\ \mathrm{uniform}}\left[
 \frac{F(T_a(U))W_a(T_a(U))}{q_a^0(T_a(U))}\right].
\]

For a bijective cube map, 1/q_a^0 equals its Jacobian. For overlapping internal
branches, it is the reciprocal of their summed pushforward density; substituting
the selected branch's Jacobian would count overlapping images incorrectly.
With the stated partition the latent integrand simplifies to F/(sum_b q_b^0).
The existing integrator applies its discrete and
continuous sampling corrections exactly once. This avoids foreign adaptive-grid
PDF lookup while retaining unbiasedness and useful focused partitions. Its
variance is not identical to the full adaptive balance heuristic. Ohl's fixed-map
alternative is the relevant precedent.[^ohl]

The current LU implementation evaluates the LMB partition at each cut's
rescaled momenta t*_c K. This is a legitimate per-cut representation partition,
but it is not automatically the mixture density of the raw proposal. A new
compact-support density cannot simply be substituted at that point. The new
map-density partition must use common raw coordinates and the actual supports,
and multiply the fully subtracted cut sum—or the same factor on each cut.
It must stay outside raised-residue derivatives and multiply generated-event
weights consistently. Otherwise one can double-count Jacobians, lose support,
or damage cancellation while still obtaining superficially reasonable numbers.

Map construction must also respect order. For the simple LU channels, generate
the cut-preserving complements, solve their host t*, construct the active
surface map, and then evaluate the full graph with all required cut kinematics.
Any shared threshold SOCP still requires all relevant cuts to have been
prepared. A proposal label never changes that physics requirement.

Runtime channel changes need no new CFF expressions when the graph geometry
is already available. They do change the adaptive coordinate meaning. Persist
a resolved map signature and reject incompatible workspace resume; use stable
names in maximum-weight records. Read-only state behavior must remain intact.

The cost is one selected forward chart plus inverse-density work for eligible
channels. General radial charts may solve r* again in each inverse; analytic
one-loop roots and complement-only center caches help. Exact star targeting
must reuse prepared geometry where possible. The overhead should be measured
against the full GL638 evaluation, not guessed from the number of new fields.

## 7. Concrete GL638 constructions

In prepared cut-1 coordinates, let p=q3, s=q7, t=q10 and a=q2. The cut is
(2,6,10), and its original-to-prepared scale tau is determined by

\[
 E_h(\tau a_0)+E(\tau(a_0+t_0))+E(\tau t_0)=Q.
\]

It depends on a0,t0, not on p0,s0. Here Q=1000 GeV, mt=173 GeV and mh=125 GeV.
Define B=E(a+t)+E(t). The relevant equations are

\[
 H=E(a+p)+E(p)-B,\quad Z=E(p)+E(t)+|p-t|-Q,
 \quad P=2E(p)-Q,\quad A=2E(s)-Q.
\]

H is physical cut 3 evaluated on cut 1. The WH multiplier uses
H^2/[H^2+(PZ/Q)^2] at the appropriate CT star. A and U have complementary
shared-one-loop and native-two-loop variants. The supplied deformation note
motivates this splitting, while the current research records identify the hard
H/Z exception and its propagated CT images.[^defo][^glrecords]

### Direct H shell

At fixed a,t, H=0 is an exact ellipsoid centered at b=−a/2, with

\[
 r_H(n)^2=\frac{B^2-|a|^2-4m_t^2}{4[1-(a\cdot n/B)^2]}.
\]

The numerator is positive except at the degenerate pair t=−a/2. This gives an
analytic, inexpensive radial map in active parent edge `[3]`, complement
`[6,7,10]`, parent `[3,6,7,10]`. In raw p0 coordinates the center and root are
divided by tau. Equivalently, mapping prepared p back to raw coordinates adds
tau^-3 to the Jacobian; its complement dependence is triangular.

This is a genuinely useful first channel. A half-power H density gives finite
local variance for the measured H/Z corner under the regularity assumptions.
It does not give bounded leading weights or automatically target every CT image.

### Native A-star H shell

At fixed a0,t0,s0 and the actual native common center c=(c_p,c_s), the A-star
projection scale ell is fixed independently of p. Therefore

\[
 p_*=c_p+(p-c_p)/\ell.
\]

H(p*)=0 is another exact ellipsoid in p, with center and radial root

\[
 b_{\mathrm{eff}}=c_p+\ell(-a/2-c_p),\qquad
 r_{\mathrm{eff}}(n)=\ell r_H(n).
\]

Divide these by tau for the raw-coordinate chart. The p* to raw-p Jacobian is
(ell/tau)^3. Thus this particular star target can reuse the one-block radial
sampler; it does not require a general six-dimensional star-first engine.
The conditional independence and actual center branch must be certified and
shared with the evaluator. Generic U-star geometry is more coupled and does
not inherit this simplification.

New arithmetic checks using recorded GL638 complements and its actual numerical
center tested 45 affine-map cases. Coordinate errors were at most 1.43e−13 GeV,
and H-star shell residuals at most 2.28e−13 GeV. On those exact H-star shells,
base H ranged from approximately −139 to +9927 GeV. This illustrates why
focusing only on the unprojected H equation can miss the intended image.
These are geometry checks, not new integrand integrations.[^newchecks]

### Joint H/Z normals

H and Z occupy the same p block. At the recorded hard intersection,
|gradient(H) cross gradient(Z)|=1.73380575, confirming rank two there.
A local `(H,Z,tangent)` chart is therefore available. Two independent
three-dimensional surface maps would overcount dimensions.

An explicit conditional chart can be constructed by setting h=H, z=Z,
u=E(p), S=B+h and D=Q+z−E(t). Squaring with the original energy signs retained
gives two linear equations for p at fixed u:

\[
 2a\cdot p=S^2-2Su-|a|^2,\qquad
 2t\cdot p=2Du-D^2-m_t^2+|t|^2.
\]

For a cross t nonzero, write the in-plane solution as d+e u and the remaining
component as w n0. The mass shell gives
w^2=u^2−mt^2−|d+eu|^2. Parameterizing the resulting ellipse by a full angle
phi includes both inverse branches. With kappa=|e|^2−1, its Jacobian is

\[
 \left|\frac{d^3p}{dh\,dz\,d\phi}\right|
 =\frac{u(S-u)(D-u)}{|a\times t|\sqrt{\kappa}}.
\]

The feasible h,z domain, positive energy signs, kappa and ellipse discriminant
must be checked. Restricted phi arcs require their exact total length in the
density. A certified compact normal patch plus ordinary coverage is sufficient;
there is no need to solve the entire global feasible domain at once.

Choose a fixed positive dimensionless normal metric alpha; for the studied
patch, alpha=|P0|/Q with P0=2E(p_intersection)−Q separated from zero is convenient.
Either uniform transverse R with h=R cos(theta), z=R sin(theta)/alpha, or
product half-power h,z profiles can compensate the leading 1/R weight on a
regular patch. The same chart can be applied to p* and pulled back affinely
for the A-star target. A previous independent determinant check supports the
conditional p-chart formula to 3.8e−9 at a regular point; global support remains
an implementation obligation.[^chart]

### What targeting Cutkosky cuts does and does not buy

In this e+e− CM setup, scaling all raw loop momenta by lambda changes tau to
tau/lambda and leaves the prepared physical point unchanged. Four new checks
with lambda=0.07,0.3,2,17 reproduced that point within 1.14e−13 GeV.
Targeting a cut's own equation after its own LU solve is identically zero.
Likewise, A evaluated at its own A-star is always zero.

Raw cut-shell focusing may match the auxiliary LU radial localization, but it
does not by itself change the sampled physical event shape or solve the H/Z
corner. The useful targets above vary the p shape at fixed cut-preserving
complements. The API must identify both the surface and the frame/block in
which distance to it actually varies.

## 8. Validation: beyond a unit-volume integral

A unit-volume test is valuable when it refers to a finite physical region with
an independently known volume. The whole Cartesian momentum space has infinite
volume. Also, checking q(T(u))J(u)=1 is potentially tautological when q and J
share the same mistake—as the recently corrected hyperspherical Jacobian
demonstrated.

An independent standalone calculation tested a shifted target ellipsoid with
axes (1,2,3), its radial sphere map, and a quadratic threshold profile. The
integration domain was rho<2, with semiaxes (2,4,6) and volume 64 pi. Numerical
Cartesian differentiation agreed with the analytic Jacobian within 1.07e−10.
Independent quadrature reproduced normalized volume
0.9999999999999966 and central second moments (0.8,3.2,7.2) within 1e−13.
Forward/inverse coordinate errors were at most 1.11e−16. No GammaLoop sampler
was implemented or executed for these checks.[^newchecks]

For the GL638 H ellipsoid there is also an exact volume oracle. Let
a_perp=sqrt(B^2−|a|^2−4mt^2)/2 and
a_parallel=a_perp/sqrt(1−|a|^2/B^2). Its volume is
4 pi a_perp^2 a_parallel/3. Raw H divides that volume by tau^3; the A-star
pullback adds ell^3. These can test the actual future sampler directly.

A transparent singular benchmark is f=1/sqrt(u^2+v^2) on `[-1,1]^2`, whose
integral is 8 asinh(1). The normalized product density
q=1/(16 sqrt(|uv|)) gives f/q at most 8 sqrt(2). On a unit disk,
q=1/(2 pi R) gives constant weight 2 pi for f=1/R. These examples distinguish
finite first moment, finite variance and bounded weights without physics
cancellation obscuring the result.

Independent exact moment calculations for the disk also make the distinction
quantitative. Relative variance means Var(weight)/I^2:

| Proposal and target domain | Relative variance |
|---|---:|
| Uniform disk, f set to zero below R=1e−4 | 3.60609 |
| Uniform disk, f set to zero below R=1e−12 | 12.81551 |
| One H half-power, normalized on the full disk | 0.697653 |
| Product half-powers normalized on the square, f supported on the disk | 0.373450 |
| Uniform transverse radius on the disk | 0 |

The ordinary proposal's second moment diverges as the artificial cutoff is
removed. The focused cases need no cutoff. Their domains and normalizations
are stated explicitly; this is an analytic benchmark, not a GL638 efficiency
measurement. Complete values are in the accompanying CSV.[^newchecks]

The implementation validation should then cover:

1. Resolver errors, signed-cycle ranks, inferred complements, dependent blocks,
   ambiguous targets and transactional CLI/Python parity.
2. Analytic roots/volumes, genuinely nonquadratic multi-loop fibers, independent
   Jacobians, inverse maps, all branches, units and precision escalation.
3. Empty/pinched fibers, root-origin boundaries, certified supports, conditional
   fallback and tails. Never infer a support certificate from a few probes.
4. Pointwise partition sums, equivalent/duplicated channels, nonuniform discrete
   probabilities, frozen and trained grids, and summed versus sampled channels.
5. The actual outer integration-weight path, events/observables, raised cuts,
   preserved CT group geometry and the full physical cut sum.
6. GL638 H/Z angular and tangential patches, hierarchical approaches, A/U-star
   images, overlap-center changes, ordinary soft limits, UV tails, old/new maxima,
   and independent seeds at equal sample and wall-time budgets.

All final GL638 tests must include its 936 orientations and complete UV setup.
Arb controls remain necessary: importance sampling deliberately visits difficult
regions more often and cannot define the exact coincident A/P residue expression.

## 9. Recommended implementation sequence and decision

First generalize the existing channel owner and raw-frame map-density partition,
then implement one regular radial fiber with normalized power focusing. The
direct GL638 H channel and its certified affine A-star image are concrete first
targets. Add independent products and conditional compositions through the same
map contract. Then add a rank-two intersection chart with reusable scalar normal
profiles, and handle the more general U-star and nonlinear pullbacks explicitly.

This is a contained architectural extension with reusable geometry; it is not a
trivial enum change and does not require replacing the Monte Carlo engine.
Automatic discovery of optimal metadata or channels is unnecessary for the
first implementation. Clear manual definitions, canonical resolution, inspection
and independent normalization tests provide a solid starting API.

The proposal therefore succeeds in its central objective: thresholds can define
useful sampling coordinates, ordinary complements can complete them, and such
channels can target GL638's remaining integrable corner. The exact-linear
whole-shell Jacobian must be replaced by a normalizable profile, composition
must be validated geometrically, and projected targets must name their actual
frame. With those conditions, both finite-local-variance and bounded-leading-weight
constructions are available for the studied corner. Complete GL638 convergence
remains a numerical and global-coverage question to test after implementation.

## Sources and supporting records

[^zeno]: Zeno Capatti, *Singularities of Feynman diagrams and their local cancellation in collider cross-sections* (ETH doctoral thesis, 2023), Chapter 3 for graph/LMB context and Section 5.1.2, printed pp. 107–111 for threshold geometry and radial parameterization. Supplied local file: [thesis_zeno.pdf](/common/dev/gammaloop_ir_safe_thresholds/TMP_TO_IGNORE/thesis_zeno.pdf). Publication identifier: 10.3929/ethz-b-000647466.
[^soper1999]: Davison E. Soper, [*Techniques for QCD calculations by numerical integration*](https://arxiv.org/pdf/hep-ph/9910292) (1999), Section VI, especially equations (71)–(77) and (89).
[^soper2001]: Davison E. Soper, [*Choosing integration points for QCD calculations by numerical integration*](https://arxiv.org/pdf/hep-ph/0103262) (2001), Sections IV and VI–VIII; equations (25)–(42), (49)–(54), and (62)–(65).
[^soper1998]: Davison E. Soper, [*QCD calculations by numerical integration*](https://arxiv.org/pdf/hep-ph/9804454) (1998), common-loop-momentum sum over cuts, especially equation (9).
[^kleiss]: Ronald Kleiss and Roberto Pittau, [*Weight optimization in multichannel Monte Carlo*](https://arxiv.org/pdf/hep-ph/9405257) (1994), equations (1)–(7).
[^ohl]: Thorsten Ohl, [*Vegas Revisited: Adaptive Monte Carlo Integration Beyond Factorization*](https://arxiv.org/pdf/hep-ph/9806432) (1998), equations (9), (15), (25)–(29), Sections 4.1–4.3.
[^defo]: Supplied [ttH_defo.pdf](/common/dev/gammaloop_ir_safe_thresholds/TMP_TO_IGNORE/ttH_defo.pdf), complementary one-/two-loop deformation-sector construction, especially the concluding example. This is contextual unpublished material, not a demonstrated global sampling prescription.
[^glrecords]: Current repository [GL638 validation](../gl638/README.md) and [runtime follow-up](../gl638/runtime-followup.md), including the frozen full-state H/Z and maximum-weight replays. Numerical provenance is recorded there; these data precede the proposed sampler.
[^code]: Current source: [runtime sampling schema](../../../crates/gammalooprs/src/settings/runtime.rs), [channel setup](../../../crates/gammalooprs/src/integrands/process/mod.rs), [LU partition boundary](../../../crates/gammalooprs/src/integrands/process/cross_section/mod.rs), [E-surface owner](../../../crates/gammalooprs/src/cff/esurface.rs), [transactional settings](../../../crates/gammaloop-api/src/commands/set.rs). Further line-level analysis: [API study](api.md).
[^chart]: Exact prior [H/Z and star-chart derivation](TWO_NORMAL_PROPOSAL.md), and [GL638 geometry assessment](gl638.md). The derivation and old finite-difference checks are distinct from a production implementation.
[^newchecks]: New standalone records: [ellipsoid geometry, volume and moments](geometry_checks.json), [GL638 rank, LU-scale and A-star checks](gl638_geometry_checks.json), [exact normal-plane moment table](normal_plane_moments.csv). These calculations change no GammaLoop code or configuration and include no new full-state integration. Independent primary-literature analysis: [literature study](literature.md).
