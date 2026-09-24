#import "../../shared.typ": source-link

#let sampling = [
= Sampling channels and maps

Sampling channels map unit-cube coordinates to the graph's unbounded spatial loop momenta.
They can use ordinary loop-momentum bases (LMBs), concentrate near an energy surface, or combine
a physical Cutkosky cut with thresholds in its remaining loop variables. They change the
integration measure and its partition, not the integrand's physical cut sum or its threshold
counterterms. Amplitudes and cross sections use the same channel catalogue and integration engine.
Channel names and IDs are independent of threshold `group_id` values. The
#link("guides/threshold-subtraction/")[threshold metadata guide] describes the common solve
spaces and centers needed for cancellation; sampling is an additional change of variables.

The settings below are *runtime overlays*: apply them to a loaded integrand with
`set process -p PROCESS -i INTEGRAND file sampling.toml`. In a run card, put the same settings
under `[default_runtime_settings.sampling]`, with the corresponding prefix on its subtables,
or apply the overlay in a command block. Changed channel definitions invalidate warm-up caches.
Start a new integration workspace when changing channels or their weighting so an old grid is
not confused with the new channel order.

== Select channels

This baseline uses optimized ordinary LMB channels and samples one channel per point:

// docs-example: syntax
```toml
[sampling]
graphs = "monte_carlo"
orientations = "summed"
sampling_multichanneling = true
sampling_channels = "monte_carlo"
default_channel_selection = ["auto:optimized_lmb"]
sampling_channel_weight = "ose"
alpha = 3.0
coordinate_system = "spherical"
mapping = "linear"
b = 1.0
```

`graphs`, `orientations`, and `sampling_channels` are separate choices, each accepting `summed`
or `monte_carlo`. With `sampling_channels = "summed"`, every selected channel maps the same
outer cube point and their contributions are summed before statistics are accumulated.
With `monte_carlo`, the channel is another adaptive discrete coordinate. Keep
`sampling_multichanneling = true` for either mode: setting it to false disables the channel
partition rather than merely disabling advanced maps. By default multichanneling is off and
all three sum modes are `summed`.

#table(
  columns: (auto, 1fr),
  table.header([*Selector*], [*Expansion*]),
  [`auto:lmb`], [All admissible generated LMB channels.],
  [`auto:optimized_lmb`], [The generated optimized subset, augmented for elementary massless-edge
    coverage when no explicit basis restriction is supplied.],
  [`auto:surfaces`], [Currently the same optimized-LMB coverage. Production does not yet discover
    automatic surface candidates; add named surface channels explicitly.],
  [`lmb(6,12,13,14)`], [One generated basis with exactly these ordered edges. No named definition
    is required.],
  [`cut1_joint_HZ`], [A user-defined channel in this graph's `channel_definitions` table.],
)

A graph-specific selection replaces the default list. Entries *within* a list are additive:
automatic selectors, explicit generated LMBs and names can be combined. Generated entries are
deduplicated in first-selection order; selecting an automatic preset twice does not add copies.
Named definitions keep their identity, even if two happen to describe equal maps.

// docs-example: syntax
```toml
[sampling]
default_channel_selection = ["auto:optimized_lmb"]
channel_selection = { GL638 = ["auto:surfaces", "cut1_joint_HZ"] }
```

`GL638` is the actual master graph name; `cut1_joint_HZ` is any name chosen for a definition.
TOML inline tables use `=`, not JSON's `:`. Bare keys such as `GL638` are valid; quote names
containing spaces or other special characters. If neither a graph override nor a nonempty
default is provided, catalogue resolution uses `auto:optimized_lmb`. An explicitly empty
graph selection stays empty and is not a request for automatic fallback.

`lmb_basis_ids = { GL638 = [1, 2] }` restricts admissible generated bases by their displayed
basis IDs, in the requested order. The restriction also applies to explicit `lmb(...)`
selectors. Missing, excluded or ambiguous ordered matches are errors; edge order is never
silently changed. Basis IDs, physical cut IDs and zero-based *sampling channel IDs* are
different index spaces. Use edge selectors or named definitions for portable steering and
inspect the resolved catalogue before using a channel ID.

To use only advanced maps, select only their names. They still must cover the integration
domain: a compact joint map requires an explicitly selected full-support sibling, which may
be a full-volume physical-cut channel rather than an LMB channel. Elementary massless-edge
coverage is useful even when no regular E-surface exists; it does not by itself prove coverage
of every simultaneous soft or collinear limit.

== Define a channel

A named definition has one Symbolica expression in `around`. Its structural validation and
graph routing determine the map; this is not a regular-expression mini-language.

#table(
  columns: (auto, 1fr),
  table.header([*Field*], [*Meaning*]),
  [`around`], [Map expression. Its `surface` or `cut` edges identify energy equations, not the
    independent variables used to solve them.],
  [`parent_lmb`], [Complete ordered parent basis defining the routing of the independent variables.],
  [`subspace_lmb`], [Active parent edges for a simple surface map. For compositions,
    `block(lmb(...), ...)` assigns the individual active blocks.],
  [`on_cut`], [Optional list of physical cut IDs used to check an explicitly identified host.
    It does not select the cuts evaluated by the integrand.],
  [`channel_weight`], [Optional override of the global channel partition score.],
  [`radial_profile`], [Optional `lu_h` proposal for the unique physical-cut block.],
  [`singularity_proxy`], [Explicit Symbolica score expression, required when using
    `singularity_proxy` weighting.],
)

Surface definitions require a complete valid parent and active edges inside it. The remaining
coordinates form the complement. The full parent matters: an edge label alone does not
identify the cycle generated by varying that momentum. A requested parent ordering can be
resolved from a generated basis without changing the generated catalogue's IDs.

#table(
  columns: (auto, 1fr),
  table.header([*Expression*], [*Meaning*]),
  [`surface(a,b,...)`], [A routed positive-energy surface, including its external energy shift.],
  [`phase_space(cut(a,b,...))`], [A physical-cut radial map retaining the auxiliary LU scale.],
  [`lmb(a,b,...)`], [An ordinary complete LMB map when used as a named channel's expression.],
  [`complement(a,b,...)`], [Ordinary coordinates for these remaining parent edges.],
  [`block(lmb(a,...), MAP)`], [Assign these active coordinates to `MAP`. The descriptor itself
    consumes no extra coordinates.],
  [`then(A,B,...)`], [Ordered conditional blocks; later geometry may depend on earlier outputs.],
  [`product(A,B,...)`], [Independent blocks with a direct-sum dependency check.],
  [`at_cut(cut(a,...), MAP)`], [Evaluate the target in the kinematics of this physical host.],
  [`left(MAP)`, `right(MAP)`], [Resolve a target on the specified side of its host.],
  [`intersect(surface(...),surface(...))`], [One joint block for two supported, distinct energy
    equations. It does not allocate two independent copies of the same loop variables.],
)

The parser also recognizes `soft(edge)` and `collinear(a,b)`, but numerical compilation of
these dedicated primitives is not implemented. Use ordinary LMB coverage for soft sampling.
General nested compositions, arbitrary multi-surface intersections and CT-projected
star-target maps are likewise not implemented. A parsed expression is not a promise that
every topology admits its geometry; unsupported or ambiguous requests fail during resolution
or warm-up.

== Amplitude example: a massive two-loop kite

The #source-link("tests/resources/graphs/massive_kite.dot", label: "massive kite fixture")
has parent `[4,6]`. With unit internal masses and incoming energy 5 in its rest frame, the
following is the tested surface-plus-LMB configuration:

// docs-example: syntax
```toml
[sampling]
graphs = "summed"
orientations = "summed"
sampling_multichanneling = true
sampling_channels = "summed"
sampling_channel_weight = "map_density"
power = 2.0
channel_selection = { massive_kite = ["C", "D", "ordinary"] }

[sampling.channel_definitions.massive_kite.C]
around = "surface(2,4,6)"
parent_lmb = [4,6]
subspace_lmb = [4,6]

[sampling.channel_definitions.massive_kite.D]
around = "surface(3,5,6)"
parent_lmb = [4,6]
subspace_lmb = [4,6]

[sampling.channel_definitions.massive_kite.ordinary]
around = "lmb(4,6)"
parent_lmb = [4,6]
```

These three-energy surfaces genuinely depend on both loop vectors. In contrast, a surface
with only one independent varying energy direction must use a proper three-dimensional
subspace; requesting it as a bounded six-dimensional surface is invalid. The simple
`subspace_lmb` form samples its complement first, then the conditional surface block.
There is no Cutkosky host, LU rescaling or LU h-function for an amplitude.

For an existing regular surface and fixed complement, the map chooses an interior center
and solves the positive radial root in each direction. `power > 1` concentrates on both
sides of that root. Locally its density scales as
`q ~ |r-r_star|^(1/power-1)`; `power = 2` gives an inverse square-root enhancement.
`b * e_cm` sets the radial scale. The proposal retains normalized tails. A density proportional
to `1/|r-r_star|` across an entire smooth shell would not be normalizable and is not supplied
by taking a finite value of `power`.
A sampling chart's interior center need not be the subtraction group's center. This freedom
does not extend to a hypothetical CT-star target: focusing the projected singularity would
require the actual counterterm projection, which the current map language does not implement.

Some conditional surfaces exist only for part of the complement phase space. Their existence
is decided from the already sampled complement. Certified absent or pinched fibers use a
normalized ordinary conditional fallback; forward and inverse use the same decision. An
ambiguous numerical solve is an error, not permission to discard a point or reinterpret a
regular surface as absent. This keeps the channel list and discrete grid dimensions fixed.

== Cross-section example: a cut and a joint threshold block

This GL638 overlay combines optimized LMB coverage with one cut channel and a composed
Cut 1/HZ channel. The short labels H and Z are mnemonic names, not particle identities:
here H denotes `surface(2,4,12)` and Z denotes `surface(3,10,13)`.

// docs-example: syntax
```toml
[sampling]
graphs = "monte_carlo"
orientations = "summed"
sampling_multichanneling = true
sampling_channels = "monte_carlo"
sampling_channel_weight = "map_density"
power = 2.0
b = 0.3
channel_selection = { GL638 = ["auto:optimized_lmb", "cut1", "cut1_joint_HZ"] }

[sampling.channel_definitions.GL638.cut1]
around = "phase_space(cut(2,6,10))"
parent_lmb = [3,4,7,10]
subspace_lmb = [3,4,7,10]
on_cut = [1]
radial_profile = "lu_h"

[sampling.channel_definitions.GL638.cut1_joint_HZ]
around = "then(block(lmb(6,10),phase_space(cut(2,6,10))),complement(7),block(lmb(3),at_cut(cut(2,6,10),intersect(surface(2,4,12),surface(3,10,13)))))"
parent_lmb = [3,6,7,10]
on_cut = [1]
radial_profile = { kind = "lu_h", approximation = "log_logistic", broad_fraction = 0.02 }
```

The composed channel first prepares Cut 1 on variables `(6,10)`, samples complement `7`, then
focuses both thresholds in variable `3`. The active displacement leaves the host cut fixed;
the physical-to-raw affine transformation and its determinant are part of the complete map.
Neither its name nor `on_cut = [1]` removes the other physical cuts from the evaluation.
The full #source-link("examples/cli/epem_a_ttxh/NNLO/gl638_advanced_sampling.toml",
label: "GL638 run card") supplies the other cut channels, a right-side threshold channel,
and the generation and physics settings needed for a complete calculation.

For example, Cut 3 and its right threshold use
`then(block(lmb(4,12),phase_space(cut(2,4,12))),complement(7),block(lmb(5),at_cut(cut(2,4,12),surface(5,10))))`
with parent `[4,5,7,12]`. Side qualifiers can disambiguate multiple eligible associations;
they do not replace checking the actual energy equation and its shift. Compositions may
prepare left and right targets after one host when the routed dependency checks permit it.
`product` is suitable only when the blocks are independent: distinct active edge lists
alone do not establish independence.

=== Matching the Cutkosky h-function

`radial_profile = "lu_h"` is shorthand for the detailed table above. Optional positive
`scale` and `shape` override the automatic fit. The default `broad_fraction = 0.02` retains
a positive ordinary-radius component. The implementation fits a log-logistic proposal to
the runtime `[h_function]`; it does *not* invert the exact physical h-function.

For a simple multiplicative-flow cut, let the physical radial root be R and the raw radius
be r, so `t = R/r`. At fixed angular shape the physical radial measure is proportional to
`h(t) dt`. Sampling t approximately according to h therefore reduces this auxiliary
variance. The focused CDF is `t^a/(t^a+s^a)` for positive shape a and scale s. A normalized
broad component is added, and its complete CDF is inverted with a safeguarded solve in
`log(t)`. The Jacobian uses the derivative of this *actual proposal*, including its mixture,
and inverse density is evaluated at the supplied point.

Changing this proposal does not alter the physical h-function, its raised-cut derivatives,
or the symmetric localization multiplying local threshold counterterms. Matching h alone
does not flatten angular structure, soft limits or threshold intersections. In particular,
raised-cut derivatives may require stronger tail optimization than this fit provides.

=== What a joint map improves

The implemented joint primitive supports two energy sums sharing one routed on-shell energy
in the same three-dimensional active block. It retains the original unsquared equations,
their masses and shifts, and both branches of the tangential circle. A certified compact
disk in the two normal coordinates gives density `q ~ 1/R`, where R is their normal-plane
radius. Unsupported energy classes, rank loss and tangencies are not solved by this primitive.

A locally integrable corner `f ~ 1/R` in two normal dimensions has a logarithmically divergent
second moment under a smooth nonzero proposal: its first moment contains `dR`, but its
second contains `dR/R`. A regular joint `q ~ 1/R` removes that leading weight growth.
A one-normal inverse square-root map can improve local variance while leaving generic-ray
weights unbounded. This is why a direct H map and a joint H/Z map have different purposes.

The GL638 local scans found flat joint-map weights on the investigated regular H/Z approach
directions. This is not a global bounded-weight or variance guarantee: other thresholds,
thresholds approached only after a counterterm projection, soft endpoints, tangencies and
shrinking chart support need separate analysis. Sampling a physical threshold does not
automatically focus its displaced counterterm-star image.

== Partition scores and support

For channel c, the integrand contribution contains its map determinant `J_c` times a
partition `w_c = rho_c / sum_b rho_b`, evaluated in the *same raw master momentum frame*.
The integrator separately supplies adaptive cube and discrete-selection probabilities.
This is a partition of unity using static map scores, not a balance heuristic containing
the foreign channels' adaptive grid densities.

#table(
  columns: (auto, 1fr),
  table.header([*Weight setting*], [*Score*]),
  [`map_density`], [Actual unadapted inverse-map density, including branch and support information.],
  [`inverse_jacobian`], [Another spelling for the same exact-density strategy.],
  [`ose`], [On-shell-energy score for a complete ordinary LMB map.],
  [`singularity_proxy`], [The supplied positive Symbolica expression in raw components
    `x0`, `x1`, ...; it is not inferred from the map.],
)

For L loops and selected LMB edges B, the OSE score is
`rho = e_cm^(-3L) product_(e in B) (e_cm/E_e)^alpha`. Energies include masses and full
affine routing at the unrescaled point. `alpha` is finite and nonnegative, defaults to 3,
and zero gives a constant score. The energy-dimension factor matters when mixing OSE
and exact map densities. OSE changes the partition, not the map's Jacobian or its radial law.

For a mixed catalogue, keep the global `map_density` default and override only an ordinary
soft-coverage channel:

// docs-example: syntax
```toml
[sampling.channel_definitions.GL638.soft_6_12]
around = "lmb(6,12,13,14)"
parent_lmb = [3,4,7,10]
channel_weight = "ose"
```

Add `soft_6_12` to that graph's selection to activate it. A global `ose` choice requires
advanced channels to override it explicitly; OSE is not supported on a surface or composed
map. A user proxy must be finite and positive wherever its channel has support, and should
retain the singular features the map was designed to sample. Actual support still gates
proxy scores: a proxy cannot assign weight to a point outside its map. A selected full-support
channel keeps the denominator positive; numerical uncertainty is not treated as zero support.
Stronger OSE or proxy singularities do not automatically improve variance, since they may
assign an overlapping region to a channel that samples its other directions poorly.

== Inspect and export the effective catalogue

// docs-example: syntax
```text
display integrands -p epem_a_tth -i NNLO --graph GL638 --show_sampling pretty
display integrands -p epem_a_tth -i NNLO --graph GL297 GL638 --show_sampling toml
```

The pretty view shows ordered selector expansion, canonical channel IDs, generated basis
and physical cut IDs, complete expressions, resolved blocks and complements, profiles and
scores. It states whether channels are sampled, summed or inactive. These are static recipes;
point-dependent roots, centers and existence decisions are not fabricated for display.

The TOML view serializes the actual settings types, replacing generated LMB expansions with
explicit `lmb(...)` selectors and retaining named definitions. Its clean stdout can be
redirected to a file and applied through `set process ... file`. Multiple graphs produce one
valid document; separate exports can be applied successively without resetting unrelated
graphs. Do not concatenate complete TOML documents. The graph-specific entries belong in
the existing run-card sampling tables, with their `default_runtime_settings` prefix.

== Validate a map on a known integral

Use the normal integration workflow with a Gaussian reference target before interpreting
a variance comparison on an expensive physical integrand:

// docs-example: syntax
```text
integrate -p epem_a_tth -i NNLO -c 4 -w reference_workspace --reference-gaussian '{"width":300}' --show-phase both
```

This substitutes a normalized Gaussian in all D spatial raw-momentum components while
retaining the real graph routing, selected maps, channel partition, adaptive grid and
stability machinery. Re reports its normalization; Im reports the raw `|K|^2` moment
divided by `|center|^2 + D*width^2`. *Both targets are one.* Im here is a second acceptance
observable, not a physical imaginary amplitude. The optional `center` is a D-element list;
omitting it means zero. `width` is positive and in the momentum units of the calculation.

Disable active physical selectors and observables for this substitution. All loaded graphs
in one integrand slot must have the same nonzero loop count. Use a separate workspace and
repeat the same reference descriptor when resuming; it cannot be silently resumed as a
physics integration. Through Python, use the same command with `api.run(...)`.

Check normalization *and* the nonconstant moment, their reported errors, invalid-point
counts and stability behavior. Repeat with enough statistics and suitable width; agreement
inside a very large error bar is weak evidence. A normalization test alone does not certify
individual branches, support or every singular limit: forward/inverse, independent
Jacobian and physical-limit tests remain necessary.

For physical runs, the absolute monitor takes componentwise absolute values after summing
physical cuts, counterterms, group members and any summed orientations at each mapped point.
When channels are explicitly summed, their positive absolute contributions are added before
the outer-sample statistics update; the code does not take the absolute value of the final
signed sum over different mapped points. If orientations are sampled, the monitor concerns
that sampled-orientation domain. Report these choices when comparing signed and absolute
errors, sample counts or maximum-weight influence.
]
