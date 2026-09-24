= Proposal: E-surface sampling channels in GammaLoop
<proposal-e-surface-sampling-channels-in-gammaloop>
Research only. No repository source, tests, configuration, or state was changed. This proposal is based on the current `ir_safe_threshold_subtraction` implementation and cross-checks with the independent geometry and literature agents. All new API names below are proposed, not currently implemented.

== Main conclusion
<main-conclusion>
A useful first implementation is feasible without replacing GammaLoop\'s integration engine. It is more than adding a `coordinate_system` variant: existing channels are LMB indices, their maps are graph-independent, and LU currently evaluates its channel partition at each cut-rescaled point. A new graph-aware map must own forward and inverse transformations, the physical measure, and a support-aware channel partition in one common integration coordinate system.

The user\'s compact expression is a good #emph[display language];, but `cut(3,4,7) x cut(2,11) x complement(4,7)` is not a complete unambiguous input in general. Edge lists do not specify the energy shift, host cut, routing, frame, or whether the two radial constructions remain independent. Named, validated surface selectors and a mandatory complete parent LMB supply these missing meanings; a short expression can then reference those names.

A sampling center does #strong[not] need to equal a CT center. Sampling channels partition the complete integrand and never change CT subspaces, common centers, groups, weights, or dual cancellation. A target evaluated at a particular CT star #emph[does] depend on that CT\'s projection definition and therefore must identify the CT variant and common-center branch.

== Actual implementation owners and consequences
<actual-implementation-owners-and-consequences>
#figure(
  align(center)[#table(
    columns: 3,
    align: (auto,auto,auto,),
    table.header([Existing owner], [What is available now], [Required extension],),
    table.hline(),
    [#link("../../../crates/gammalooprs/src/settings/runtime.rs")[`settings/runtime.rs:1095`];, `SamplingSettingsParser`], [Flat runtime controls for graph/orientation/channel summation, LMB subset IDs, coordinates, radial mapping and weighting], [Add optional named graph-scoped channel specifications here; preserve one serialization/validation route for TOML and Python],
    [#link("../../../crates/gammalooprs/src/integrands/process/mod.rs")[`integrands/process/mod.rs:1919`];, `LmbMultiChannelingSetup`], [Master graph, all parent bases, effective LMB channels, momentum routing and channel partitions], [Generalize the existing owner to resolved channels: ordinary LMB or graph-aware surface chart; do not invent a second unrelated integration framework],
    [#link("../../../crates/gammalooprs/src/integrands/process/gammaloop_sample.rs")[`gammaloop_sample.rs:650`];, `default_parametrize`], [Calls `global_parameterize` before receiving target graph geometry], [Construct surface maps after graph/channel selection, while keeping the ordinary map implementation reusable],
    [#link("../../../crates/gammalooprs/src/momentum/sample.rs")[`momentum/sample.rs:120`];, `SubspaceData`], [Complete parent routing; selected edge plus signed fundamental-cycle signatures; complement projection], [Reuse for exact active/complement coordinates; do not redefine subspace identity using edge names alone],
    [#link("../../../crates/gammalooprs/src/cff/esurface.rs")[`cff/esurface.rs:49`];, `Esurface`], [Energy edges #strong[and external shift];, existence classification, radial value and derivative], [Resolve geometry independently of whether a threshold CT exists; reuse root/value machinery],
    [#link("../../../crates/gammalooprs/src/cff/esurface.rs")[`cff/esurface.rs:632`];], [Active radial rescaling with complement held fixed], [Natural local map primitive for a regular surface fiber],
    [#link("../../../crates/gammalooprs/src/subtraction/lu_counterterm.rs")[`subtraction/lu_counterterm.rs:1665`];], [Prepares shared overlaps from already prepared per-cut kinematic points], [Exact star targets can reuse this physics geometry, but building a global invertible sampling pullback requires additional work],
    [#link("../../../crates/gammalooprs/src/integrands/process/mod.rs")[`process/mod.rs:3344`];], [Existing nested graph/orientation/channel discrete grids, per-channel continuous grids], [Keep graph/orientation axes; generalize the final axis from LMB channel to sampling channel],
    [#link("../../../crates/gammalooprs/src/integrate/mod.rs")[`integrate/mod.rs:940`];], [Channel labels currently assume basis edge lists], [Persist/display stable names and resolved chart identities, not fabricated LMB IDs],
    [#link("../../../crates/gammaloop-api/src/commands/set.rs")[`gammaloop-api/src/commands/set.rs:208`];], [Transactional runtime setting updates], [Use this validation boundary; a bad channel must leave every target integrand unchanged],
    [#link("../../../crates/gammaloop-api/src/python.rs")[`gammaloop-api/src/python.rs:3475`];], [Python `run(command)` invokes the same CLI parser/session], [A first API already obtains Python parity using the existing settings-string command; a typed convenience method can serialize the same schema later],
  )]
  , kind: table
  )

The key current LU distinction is explicit at #link("../../../crates/gammalooprs/src/integrands/process/cross_section/mod.rs")[`cross_section/mod.rs:1963`];: channel weights partition the fully subtracted LU-cut integrand after raised-residue derivatives. The inverse-Jacobian partition is evaluated at `tstar * raw_loop_momenta` (lines 1987--1995), not necessarily at the integration point from which the proposal was drawn. It is a legitimate per-cut partition, but it is #strong[not automatically the balance heuristic of the raw integration-variable proposal];. A compact-support new chart would be unsafe if its density were simply inserted at this rescaled argument.

== Geometry contract to expose
<geometry-contract-to-expose>
For a resolved positive-energy surface, fixed complement `y`, and active dimension `d`, select an interior point `c(y)` and solve one positive radius `rstar(n,y)` for every direction. Then

`k_active = c(y) + rho * rstar(n,y) * n`

sends the surface to `rho=1`. Its active volume is

`rstar(n,y)^d * rho^(d-1) d rho d Omega`.

Angular derivatives of `rstar` cancel from the determinant; derivatives of the center and radius with respect to already sampled complement coordinates lie in off-diagonal blocks. This simplification requires that the center be a function of the fixed complement, not an unaccounted function of the sampled active coordinates.

The surface must exist, have a nonempty interior, and bound a finite region in the #strong[chosen active space];. Directions annihilated by every surface momentum make a cylinder and belong in the complement. Pinched surfaces and nonsmooth massless loci require classification rather than a promise that every E-surface is globally a smooth sphere.

Near the threshold choose `rho-1 = sign(z) * C * |z|^nu`, with `nu>1`. This produces a normal proposal singularity `q proportional to |rho-1|^(-beta)`, `beta=1-1/nu`, hence `0<beta<1`. `power=2` means `beta=1/2`. The Jacobian vanishes like distance #strong[to a power strictly smaller than one];. A full-shell density proportional to `1/|rho-1|` is not normalizable; the literal linear vanishing of the Jacobian cannot be requested as an ordinary finite-probability shell channel.

Use a normalized full-support radial profile, with an ordinary decaying tail outside a finite threshold window. A piecewise inside/window/outside map needs explicit branch masses and inverse densities. If a window is specified in GeV, its inner half is clipped to the available `rstar` without silently extending to negative radii. A relative width may be offered separately, but do not make dimensional conventions implicit.

For LU, the target `E(tstar(k) k)` is not automatically a positive-energy convex surface in raw `k`. The elementary construction applies directly when the active coordinates leave all host-cut momenta fixed, so that `tstar` depends only on the complement. Otherwise the host root changes while constructing the radial map and its derivatives belong in the pullback Jacobian. A generic CT-star pullback may be still more complicated. These are explicit chart capabilities, not silently accepted enum combinations.

== Recommended compact runtime schema
<recommended-compact-runtime-schema>
Keep the current ordinary LMB list and coordinate controls. Add optional named channels, whose names are distinct from threshold `group_id` values. The canonical runtime names are `sampling_multichanneling`, `sampling_channels`, and `sampling_channel_weight`; the old `lmb_*` spellings remain input aliases for existing cards. This naming migration does not require a second sampler architecture. Example for the useful GL638 H target:

```toml
[sampling]
# Existing ordinary channels remain the coverage component.
graphs = "monte_carlo"
orientations = "summed"
sampling_multichanneling = true
sampling_channels = "monte_carlo"
lmb_basis_ids = { GL638 = [1, 2, 13, 14, 21, 22, 84] }
coordinate_system = "spherical"

# PROPOSED fields follow.
sampling_channel_weight = "map_density"

[sampling.channels.GL638.H_shell]
parent_lmb = [3, 6, 7, 10]
subspace_lmb = [3]
map = "radial"
lmb = [3]
target = { kind = "energy_surface", edges = [2, 4, 12], at = "lu_cut", host_cut = [2, 6, 10] }
profile = { type = "threshold_power", power = 2.0, width_gev = 10.0 }
```

Defaults should do useful work without hiding physics:

- The graph table resolves against the graph group\'s master topology. For the first implementation, allow single-master cross sections and ordinary existing amplitude groups; veto unsupported nonmaster geometry rather than imposing a new cross-section graph-group design.
- `parent_lmb` is mandatory for an explicit chart. Active `lmb` edges are a subset. The remaining parent edges are the ordinary complement, here `[6,7,10]`; no need to repeat them. An explicit `complement` is an optional assertion and must agree exactly with the resolved direct-sum decomposition.
- `center` defaults to a deterministic interior center of this surface fiber. This may differ from CT centers without changing the integral. A center defined from complement-only kinematics makes forward and inverse maps consistent.
- `profile.power` defaults to 2; a documented threshold window and ordinary tail ensure normalization and whole-space support. Keep width explicit in the first user card rather than inventing a physics-independent optimal value.
- Graph/host/energy-shift ambiguity is a clean error with all candidate equations, cuts, sides and orientations. `side` is only required where necessary to distinguish threshold instances. Do not require a CT side for a generic physical cut surface.
- If the surface is absent or pinched for a sampled complement, use the documented deterministic ordinary active map for that fiber. Decide the branch from the complement #strong[before consuming active target coordinates];. The inverse density recomputes the identical decision. Numerical failure to determine a supposedly regular map is an error/rescue request, not silently a zero-weight point.
- Main named channels have a static catalog. A runtime-varying set of overlap centers is an internal normalized mixture, not a changing top-level discrete-grid dimension: record the selected branch and sum all of its inverse densities. A missing surface/center set uses the declared ordinary branch consistently.
- The new map settings are runtime settings. Changing a channel shape should not regenerate CFF evaluators; it does invalidate the associated learned grid or require a new workspace.
- A useful first UI keeps adaptive channel probabilities managed by the existing discrete grid. Add named initial probabilities/floors only through that owner\'s API, not a separate set of weights with indistinguishable meanings.

`H` must be resolvable even if it is not an existing CT metadata entry. For GL638, `[2,4,12]` is the physical cut-3 energy surface evaluated after host cut 1. The existing generation option can deliberately omit thresholds that are physical cuts; a CT-only registry would therefore exclude precisely this useful sampling target.

=== Multiple surfaces
<multiple-surfaces>
The full canonical form uses named blocks. A compact display can then reproduce the user\'s language:

```toml
[sampling.channels.G.two_shells]
parent_lmb = [1, 2, 4, 7]
map = "product"
blocks = [
  { target = { kind = "energy_surface", edges = [3,4,7] }, lmb = [1] },
  { target = { kind = "energy_surface", edges = [2,11] }, lmb = [2] },
]
# complement [4,7] is inferred from the full parent LMB.
profile = { type = "threshold_power", power = 2.0, width_gev = 10.0 }
```

This is a schema illustration, not a claim that these invented edge assignments describe a valid graph. The resolver must verify the actual graph. Display:

`two_shells: cut(3,4,7)[lmb=1] x cut(2,11)[lmb=2] x complement(4,7)`.

Support three deliberately different constructions:

+ #strong[`product`];: active blocks form a direct sum; each target/center/radius is independent of the other active blocks at fixed common complement. This is the genuinely commuting product. Merely having different selected edges is insufficient: edge momenta in the two surfaces may depend on both blocks.
+ #strong[`conditional`];: an explicit ordered composition of invertible conditional maps, with full inverse applied in reverse order and product of conditional determinants. A dependency DAG gives a triangular sampling schedule. If a later map changes a surface targeted earlier, the overall map can remain valid but simultaneous final alignment is not established; report that distinction. A cyclic claim of simultaneous independent roots requires a joint solution, not a guessed Jacobian product.
+ #strong[`intersection`];: multiple surface equations in the same active block. Use a local chart of their independent normals plus tangential coordinates, requiring full normal rank and branches/patches. Separate the coordinate geometry from the normal density: either a normal-plane radial profile, or a product of one-dimensional threshold-power profiles in the exact normal coordinates. For two regular normals, power 2 in each yields `q proportional to 1/sqrt(|H Z|)`, which scales as `1/R` along generic rays and also bounds the leading `f proportional to 1/R` corner weights under bounded angular coefficients. This is the extension relevant to simultaneous H/Z focusing; it is not an ordinary product of two copies of the same three-vector coordinates.

A string grammar is optional sugar after the structured schema is stable. Resolve named targets first, then accept `H_shell x Z_shell x complement(...)` when its meaning is unique. Never make a handwritten mini-language the only route to specifying host cuts, shifts, signs, orientations, frames or profiles.

=== Base, cut and star target selectors
<base-cut-and-star-target-selectors>
Use one target family with a clear `at` qualifier:

- `at="integration"`: evaluate the specified graph energy equation at raw integration momenta.
- `at="lu_cut", host_cut=[...]`: evaluate it at that host cut\'s prepared `tstar` point. For an elementary radial chart require host-cut invariance under the chosen active displacement, otherwise select a supported implicit chart.
- `at={ ct_star = { host_cut=[...], threshold=[...], variant="...", center="..." } }`: evaluate the target equation at the named CT projection. This is an advanced capability whose selector must resolve the variant subspace and overlap-center branch. An `at` selector alone does not prove existence of a global inverse chart. Some useful star pullbacks are affine in a suitable active block and belong in the first implementation once that property is certified; generic nonlinear star pullbacks need not block them.

Internal surface identity contains the signed energy equation/external shift and topology identity. Internal subspace identity uses the existing selected-edge/signed-cycle representation. A stored selector resolves symbolic/named user input into these explicit identities and prints them for inspection. A raw ordinal cut or CT ID is acceptable as a convenience for one loaded state, but a named/edge-based canonical export is preferable for durable cards.

For exact star targets, forward/inverse/evaluator paths should share the same prepared-cut and group-center context, or certify that their complement-only computation is identical. Independently rerunning the SOCP can move a center by numerical noise: the resulting sampler may still normalize, but its density need not diverge on the actual runtime singular set below that displacement. Reusing prepared geometry also avoids repeated SOCP work.

For exact star targets, do not reuse a sample-dependent CT center as though it were fixed geometry throughout a new sampling map. Either parameterize with the actual center dependence and its Jacobian, or use a coordinate construction in which that center depends only on earlier coordinates. Which option is available is a property of the target, to be displayed/validated.

== CLI and Python usability
<cli-and-python-usability>
The existing settings entrypoint already supports a first implementation without a parallel public configuration system:

```python
# Existing Python method; new TOML fields would enter the same Rust validator.
gl.run("set process -p epem_a_tth -i NNLO string '" + sampling_toml + "'")
settings = gl.get_integrand_settings(process_id=0, integrand_name="NNLO")
```

The current `get_integrand_settings` accepts an optional integer process ID (`python.rs:3338`); named selection already works through `run(command)`. A convenience API should accept process names consistently, rather than pretending the existing getter does so. A better proposed convenience would be

```python
gl.set_integrand_settings(process="epem_a_tth", integrand="NNLO", settings={
    "sampling": {"channels": {"GL638": {"H_shell": h_shell_spec}}}
})
```

with a ordinary dictionary matching TOML, returning the resolved settings through the existing settings wrapper. Add this only by extending the existing transactional settings command owner. There is little reason to require a separate fluent class hierarchy in v1. Later small named builders could produce these same dictionaries.

Proposed inspection commands should expose rather than conceal geometry:

- `display sampling -p ... -i ... --resolved`: channel name, master graph, full parent LMB, signed active cycles, complement, exact energy equation, target frame, center rule, profile, dependencies and support.
- `inspect sampling ... --channel H_shell --point ...`: forward/inverse roundtrip, `rstar`, signed distance, full Jacobian, per-channel density and partition sum.
- `validate sampling ...`: static validation plus an optional pointwise geometry/normalization audit, with failures attached to resolved cut/surface names.

These are proposed commands. Existing display/inspect infrastructure and labels should own them; do not add a standalone sampler process with separate state-loading rules.

== The exact estimator to implement and test
<the-exact-estimator-to-implement-and-test>
Let `T_a` map unit-cube coordinates into the #emph[same raw graph integration frame];, with Jacobian `J_a`, inverse branch set `B_a(k)`, and optional adaptive cube density `g_a`. Its physical density is

`q_a(k) = sum_{u in B_a(k)} g_a(u) / |J_a(u)|`.

For a normalized mixture `q_mix=sum_a alpha_a q_a`, a point is weighted by `F(k)/q_mix(k)`, with `F` the complete physical-cut/orientation integrand sum selected by the existing runtime mode. Channel labels only select the proposal. They must not select which Cutkosky cuts are evaluated. This is precisely the distinction emphasized by Soper\'s 2001 scattering-singularity sampling paper, discussed in the literature companion.

The current integrator already attaches an outer `sample.get_weight()` for discrete and continuous grids (`integrate/mod.rs:2119`, `3199`). Do not multiply by a mixture denominator and then retain incompatible outer probabilities. A minimal compatible alternative uses a partition of unity built from #emph[unadapted] map densities, `w_a=q_a^0/sum_b q_b^0`, and integrates `F(T_a(u))*J_a(u)*w_a(T_a(u))` in each existing grid. This is unbiased with the current outer adaptive/discrete corrections and requires no foreign adaptive-grid PDF lookup. It is generally less variance-optimal than the full adaptive balance heuristic. State which estimator is used; both are valid, but their formulas must not be conflated.

For v1, extend the existing inverse-Jacobian partition path into an explicitly named `map_density` partition evaluated at raw integration coordinates. It should use each channel\'s actual inverse map and support. Retain the existing grid corrections. If the full adaptive balance heuristic is later wanted, first add/verify foreign-grid density evaluation and count discrete probabilities exactly once.

The sampling partition remains outside raised LU derivative operations and outside the threshold CT multiplier algebra. Apply one raw-frame partition to each fully subtracted physical cut or after their full sum, consistently including generated-event weights. Do not differentiate it as physical numerator data, and do not let it break the common grouping/center treatment that guarantees dual cancellation.

The extra cost is roughly one selected forward map plus inverse density evaluation for every eligible channel. In a radial chart the inverse recovers direction and complement directly but may still need `rstar` solving. Cache shared fibers/centers and exploit analytical one-loop roots. Benchmark this geometric overhead against the expensive full GL638 evaluator before judging the feature by channel count alone.

== Validation plan and scope
<validation-plan-and-scope>
A unit-volume integration is valuable but not sufficient. Existing unit-volume tests are at #link("../../../crates/gammalooprs/src/tests.rs")[`gammalooprs/src/tests.rs:317`];. It is easy for matching wrong forward/inverse Jacobians to pass a roundtrip, and it is easy for a finite run to miss a branch/support error.

Required checks, ordered by boundary:

+ #strong[Static resolver];: unknown graph/edge, invalid full LMB, duplicate active directions, uncovered complement, wrong host/side/orientation, ambiguous energy shifts, incompatible product dependencies, CT-star variant ambiguity, nonmaster restriction. Each error lists all resolved conflicting objects. No state mutation after validation failure.
+ #strong[One-loop analytic oracle];: massive two-body ellipsoid/sphere with known center, root and volume; empty and pinched fibers; both sides of the shell; near-origin and far-tail coverage. The GL638 H target has an especially useful analytic one-loop oracle.
+ #strong[Multiloop numerical Jacobian];: independent finite-difference determinant at regular points for a genuinely nonquadratic positive-energy surface, including complement-dependent centers. This must test the actual map, not a duplicate formula.
+ #strong[Forward/inverse and PDF];: map/inverse roundtrips, density from every branch, normalization and CDF tests, dimensions/units, arbitrary-precision consistency. Test piecewise branch transitions and deterministic missing-surface fallback explicitly.
+ #strong[Known integrals];: normalized finite-volume domains and smooth compact/noncompact test functions in 3, 6, 9, 12 dimensions; massive and shifted geometry. Use enough independent analytic domains so a shared Jacobian error cannot cancel.
+ #strong[Multichannel identity];: `sum_a w_a(k)=1` pointwise, including outside special supports; adding/removing/duplicating equivalent channels leaves the answer invariant; Monte Carlo versus summed channels; nonuniform discrete probabilities; trained versus frozen continuous grids. Test the estimator\'s actual current outer `sample.get_weight()` path.
+ #strong[LU/CT physics];: all cuts still evaluated under a sampling-channel label, one partition outside derivatives, unchanged CT grouping, native `tstar` preparation before target/star geometry, events and summed observables consistent. Compare the same physical points before and after channel introduction.
+ #strong[GL638 effectiveness];: first direct H shell in `[3]` with host cut 1 and fixed complement; measure `f/q`, first and second local moments on H/Z rays, multiple seeds at fixed evaluation/time budgets, absolute-integral/max-weight behavior, and replay all old and new extremes in Arb. Then evaluate own-A-star channels separately; the direct H shell does not prove control of those extrema.
+ #strong[State/workspace];: read-only state loading, channel definitions resolved identically in CLI/Python, stable names in extrema/replay records; changed map geometry or channel order cannot silently reuse an incompatible learned grid.

Phases should be correctness-driven: (I) ordinary LMB plus one regular energy-surface fiber, including certified affine CT-star pullbacks, shared runtime schema/inspection, exact unadapted map-density partition and tests; (II) validated independent products and ordered conditional compositions; (III) local same-subspace intersection charts and generic nonlinear CT-star pullbacks. Both the GL638 direct H channel and its native-A-star affine pullback are meaningful phase-I experiments. A claim that phase I cures every GL638 maximum would exceed the present evidence.

== GL638 star-target refinement from the independent geometry audit
<gl638-star-target-refinement-from-the-independent-geometry-audit>
The current hard native A-star target has a useful simplification. Work in parent `[3,6,7,10]`, and fix the complementary raw coordinates `(a0,t0,s0)`. The host root `tau`, the relevant native A-group center `c=(c_p,c_s)`, and the A-projection scale `ell` are then independent of the active `p0`. In prepared coordinates,

`p_star = c_p + (p-c_p)/ell`.

Consequently `H(p_star)=0` is an ellipsoid in the same three-dimensional `p` block. If the original H ellipsoid has center `-a/2` and radial root `rH(n)`, its affine pullback has center

`c_p + ell*(-a/2-c_p)`

and radial root `ell*rH(n)`. Transform to raw coordinates by dividing center and root by `tau`; include the corresponding scale in the full raw-coordinate Jacobian. This gives a concrete `at=ct_star` one-surface channel using the same radial owner, with no requirement to implement a general six-dimensional cone sampler first.

The selector still has to identify A=(7,8), its native two-loop variant `[3,7]`, and its actual overlap-center branch. If several centers contribute, all relevant pullbacks need representation. The affine independence must be certified against the actual runtime center dependence; it cannot be inferred from the A energy-edge list alone. The U-star projection has a different dependency and is not covered by this argument.

The same-subspace H/Z intersection chart can also be more powerful than either one-shell map. Exact local `(H,Z,phi)` coordinates in the active p block, with a product of power-2 normal maps, realize the user\'s product-density intuition despite H and Z occupying the same loop subspace. The coordinate inversion and rank/branch checks belong to `map="intersection"`; the density can remain a simple product of one-dimensional profiles. This is a useful separation of API concepts.

== Primary literature connections supplied by the literature audit
<primary-literature-connections-supplied-by-the-literature-audit>
- #link("https://arxiv.org/abs/hep-ph/9910292")[Soper, #emph[Techniques for QCD calculations by numerical integration] (1999)] motivates surface-aware coordinate and density choices.
- #link("https://arxiv.org/abs/hep-ph/0103262")[Soper, #emph[Choosing integration points for QCD calculations by numerical integration] (2001)];, especially the separate cut/scattering-surface channel indices and discussion around Eqs. 17--18: the sampled channel\'s cut label does not restrict the sum over physical cuts in the evaluated integrand; every other channel density must be evaluable at the same point.
- #link("https://arxiv.org/abs/hep-ph/9806432")[Ohl, #emph[Vegas revisited: Adaptive Monte Carlo integration beyond factorization] (1998)];, Eqs. 9/15 versus Sec. 4.3 Eq. 27: exact adaptive mixtures need foreign inverse maps/Jacobians/grid densities; fixed analytic partitions with independent adaptive integrations are a legitimate cheaper alternative. Eqs. 25--26 give an analytic unit-integral benchmark combining spherical, cylindrical and Cartesian structures, useful beyond a single volume check.
