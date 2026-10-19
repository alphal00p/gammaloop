#import "../../shared.typ": callout, source-link

#let threshold-subtraction = [
= Threshold subtraction metadata

Threshold metadata tells GammaLoop which loop variables to vary when constructing each
threshold counterterm and how to combine alternative counterterms for one threshold.
It is generation input attached to a DOT graph. Changing it requires regenerating the
integrand; changing runtime sampling channels does not change the subtraction prescription.

The interface applies to amplitudes and cross sections. GammaLoop discovers eligible
thresholds, resolves the supplied directives, and checks their geometry. It does not
automatically derive an IR-safe splitting or prove that a user multiplier preserves the
integral. The examples below explain the implemented interface and the checks needed when
designing such a prescription.

== Start from the generated metadata

Generate the desired process and orientation selection, then inspect the effective metadata:

// docs-example: syntax
```text
display integrands -p epem_a_tth -i NNLO --graph GL638 --show_threshold_subtraction pretty
display integrands -p epem_a_tth -i NNLO --graph GL638 --show_threshold_subtraction toml
```

The pretty report shows physical cuts, sides, threshold edge sets, named variants,
generation status, selected subspaces, full parent bases, solve groups, multipliers and
shared definitions. This is persisted generation information; changing the current beam
energy does not retrospectively change the exported declarations.

TOML mode emits a complete `threshold_counterterms = "...";` DOT attribute for exactly
one master graph. It expands eligible defaults and fills resolved fields on explicit
variants, preserving skipped and disabled declarations, function-map scope, and the absence
of automatic group IDs or identity multipliers. To capture it from a saved state:

// docs-example: syntax
```sh
./gammaloop --read-only-state -s ./saved_state \
  run -c "display integrands -p epem_a_tth -i NNLO --graph GL638 --show_threshold_subtraction toml; quit" \
  > threshold-attribute.dot
```

Replace the graph's existing attribute with that fragment; it is not a complete graph file.
The embedded TOML's double quotes are escaped as `\"` inside the DOT string. Linnet removes
that DOT quoting before the TOML parser runs. Prefer this exporter or `save dot` to manual
escaping. Numerical edge arrays are emitted on one line.

Threshold generation is controlled by `global.generation.threshold_subtraction.enable_thresholds`.
If generation was disabled, the pretty report says so and TOML export fails clearly: a
DOT threshold attribute cannot encode that global toggle. Enabled generation with no
eligible thresholds exports a valid empty specification. Runtime
`subtraction.disable_threshold_subtraction` is a separate diagnostic switch.

== Schema and defaults

The attribute contains versioned TOML with the hierarchy
`cuts → thresholds → counterterms`. A minimal partial GL297 declaration is:

// docs-example: syntax
```text
threshold_counterterms = "
schema_version = 1

[[cuts]]
edges = [2, 11, 14]
[[cuts.thresholds]]
edges = [5, 7]
[[cuts.thresholds.counterterms]]
name = \"left_1l\"
subspace = [5]
parent_lmb = [3, 5, 11, 14]
";
```

This is one association, not the complete GL297 prescription; use the maintained DOT
or its expanded export when reproducing the full calculation.

All IDs are graph edge IDs, not half-edge indices, E-surface IDs, cut ordinals or integration
channel indices. The declared cut and threshold edge lists are sets: their order is
canonicalized and duplicates are errors. `subspace` and `parent_lmb` are ordered lists of
defining edges; their order is retained. The left/right side is inferred from topology and
must not be supplied as another field. For an amplitude use `cuts.edges = []`.

#table(
  columns: (auto, 1fr),
  table.header([*Field*], [*Meaning*]),
  [`schema_version`], [Required; currently `1`. Unknown fields and versions are rejected.],
  [`function_map`], [Optional graph-wide map of scalar definitions used by multipliers.],
  [`cuts`], [Physical cuts; each has required `edges` and optional `thresholds`.],
  [`thresholds`], [Each has required nonempty `edges` and optional `counterterms`.],
  [`name`], [Variant name. A single unnamed variant becomes `default`; multiple variants
    require distinct nonblank names within that threshold.],
  [`parent_lmb`], [Required on every explicitly supplied variant, including disabled ones.
    It must specify a complete valid parent loop-momentum basis.],
  [`subspace`], [Optional nonempty list of active defining edges from that parent.
    Omission requests the maximal inferred cut-side subspace in the supplied parent.],
  [`group_id`], [Optional nonnegative refinement label; explicit IDs must start at zero
    and be contiguous in their graph/group namespace.],
  [`disable`], [Defaults to `false`; disables this variant when `true`.],
  [`multiplier`], [Optional scalar weight. Absence is identity and allocates no evaluator.],
)

With no metadata, GammaLoop uses the maximal available loop subspace on each side of
each cut. All thresholds on one side of that cut consequently share an automatic solve
group. Missing cuts, missing thresholds and an empty or omitted `counterterms` list retain
those defaults. An explicit nonempty variant list replaces the default for that association;
it does not append to it. If every listed variant is disabled, no default is restored.

A topology-valid declaration filtered out by the selected orientation remains dormant and
is retained in the export. This differs from a malformed edge set or invalid parent, which
is an error. A selected-orientation run cannot establish coverage of every orientation.

== What identifies a common solve space

A defining edge alone does not identify the variation: its fundamental cycle depends on
the parent LMB. GammaLoop derives a canonical signature from every selected edge and the
complete signed cycle generated by varying its momentum with the other parent momenta
fixed. Each cycle is oriented positively along its defining edge and sorted by edge ID;
the list of defining-edge/cycle pairs is also sorted.

For example, the pretty report's
`e7 → +e7 +e8 +e14` means that an increment of the selected edge-7 momentum induces
the same increment on edges 7, 8 and 14, with all omitted edges unchanged. Minus signs
would mean the opposite increment. These are coefficients of a vector-momentum variation,
not energies or an additional threshold equation.

Two different parents can therefore define the same solve space. They share a group when
their selected defining edges and normalized signed cycles coincide. Equal parent names
are unnecessary; merely spanning the same abstract vector space with different defining
edges is insufficient. Parent order and changes of the fixed complement do not enter the
group key. Each association still retains its own native parent and sampled complement.

An explicit `group_id` can only subdivide these automatic groups:

- Equal signatures with the same explicit ID share a group.
- Equal signatures with different IDs, including an absent ID, are separate groups.
- Different signatures never share a group. Reusing an explicit ID across them is an error,
  with the conflicting cuts, thresholds, variants, parents and signed cycles reported.

Use such subdivision only when the physics permits independent subtraction. Constituents
of one degenerate higher-power residue must retain compatible variant names, group IDs,
solve cycles, disabling and multiplier directives. GammaLoop rejects incompatible
declarations instead of splitting that residue across independently treated groups.

== Cross-cut kinematics and cancellation

Compatible thresholds from different Cutkosky cuts of one graph participate in the same
group. *Every participating cut's kinematics, including its own LU rescaling root, must be
prepared before constructing that group's SOCP.* The common problem uses those distinct
cut-specific external data. Reusing one cut's rescaling for another would change the
surfaces being solved.

Each group independently owns surface classification, maximal-overlap discovery, SOCP
center construction, radial projections and multichanneling over its overlap structure.
More than one maximal overlap can require more than one center within a group. A group's
center is expressed in the common selected-edge coordinates and transported to each
association's native frame while retaining that association's complement. It is not a
common full-graph momentum point across different cuts. The radial measure uses the number
of active loops of that variant. Changing the parent chart does not authorize dropping
the existing projection or coordinate Jacobians.

This threshold multichanneling is distinct from the integration sampling channels in the
#link("guides/sampling/")[sampling guide]. A sampling-map change cannot repair an
incorrect threshold-group assignment.

There are two competing requirements when choosing subspaces:

- Cuts approaching a common pinched IR configuration need compatible threshold
  counterterms in a common solve space so their leading soft terms cancel.
- Intersecting thresholds within one amplitude need common coordinates and centers for
  dual cancellation. They are not repaired by adding same-amplitude iterated threshold
  counterterms.

Split variants and weights allow one threshold to participate in different regions with
different subspaces. GammaLoop also constructs the required left/right Cartesian-product
terms for cross sections; users do not enumerate those pairs in DOT metadata.
For such pairs each selected cycle must leave the other variant's defining edges fixed;
incompatible pairs are rejected. This linear compatibility is necessary but does not prove
that arbitrary nonlinear projections preserve both threshold equations. Coupled choices
still require a physics check of their merged roots.

For amplitude graph groups, supplied member metadata is interpreted using the master
graph's edge IDs and topology. The user must align those IDs meaningfully; GammaLoop does
not infer graph correspondences. Native member masses and residue evaluators remain
member-specific. Implicit amplitude defaults retain native full-space treatment. For
cross-section graph groups, explicit threshold metadata on non-master graphs is currently
rejected. This does not restrict sharing between cuts of the same master graph.

== Multiplier expressions and reusable functions

Each multiplier contains `expression`, optional `function_map`, `symmetrize = false` and
`opaque_derivatives = true`. Only the latter two values are implemented:
`symmetrize = true` is not an automatic cure, and `opaque_derivatives = false` is rejected.
Multipliers are applied after the raised-residue derivatives have been formed, as constants
with respect to those derivatives. They must evaluate to a finite real scalar in the active
precision. Tensor expressions must contract to a scalar without free indices.

#table(
  columns: (auto, 1fr),
  table.header([*Expression*], [*Value*]),
  [`E(edge)`, `E(effective,edge)`, `E(star,edge)`],
    [Positive on-shell energy of a graph edge in the selected view.],
  [`Q3(edge,cind(i))`, `Q3(star,edge,cind(i))`],
    [Spatial momentum component, with `i = 1, 2, 3`; explicit `effective` also works.
    The Lorentz-indexed tensor has zero temporal component.],
  [`eta(eset(a,b,...))`, `eta(star,eset(a,b,...))`],
    [Cut-local E-surface residual; explicit `effective` also works. The symmetric
    `eset` wrapper is mandatory. Ambiguous external shifts are errors.],
  [`P(position,cind(i))`], [External four-momentum component in the existing external order.],
  [`UFO::MT`], [Example of a model parameter, using its registered model symbol.],
)

`Q` and `K` are not multiplier primitives; use `Q3` for edge spatial momenta. A Minkowski
contraction of two `Q3` tensors is minus the Euclidean spatial dot product. The default
view is `effective`, not `star`. `eta` only addresses an unambiguous surface in the owning
cut's catalogue; a foreign-cut residual may need an equivalent energy expression.

Graph-wide definitions and per-multiplier definitions are merged before constructing
the Symbolica evaluator. A local entry overrides the identical shared key. Use the same
formal signature for an override: multiple signatures for one function head are rejected.
Bare scalar names and zero-argument functions are supported, as are functions with distinct
scalar formal arguments and nested calls. Edge IDs and frame tags are static inputs to the
kinematic primitives, not generic function formals. Reserved primitives cannot be redefined;
recursive reachable functions and undefined scalar symbols are errors.

Here is an amplitude syntax example; edge IDs and the parent must match the actual graph:

// docs-example: syntax
```toml
schema_version = 1

[function_map]
"scale()" = "P(0,cind(0))+P(1,cind(0))"
"damping(x)" = "1/(1+(x/scale())^2)"

[[cuts]]
edges = []
[[cuts.thresholds]]
edges = [3, 4]
[[cuts.thresholds.counterterms]]
name = "weighted"
subspace = [3]
parent_lmb = [3]
[cuts.thresholds.counterterms.multiplier]
expression = "damping(eta(effective,eset(3,4)))"
opaque_derivatives = true
[cuts.thresholds.counterterms.multiplier.function_map]
"scale()" = "UFO::MT"
```

This illustrates binding and local scope, not a complete IR-safe partition: another variant
or an analytic argument is needed before changing a physical subtraction weight.

Graph-declared scalar `params` can also appear in multipliers. For example, declare
`params = "threshold_weight_bias";` in DOT and use that symbol in the expression. After
generation its value is supplied in declaration order through the ordinary runtime field:

// docs-example: syntax
```toml
[general]
additional_param_values = [1.0]
```

The number of supplied values must match the declared parameters. Adding a parameter or
changing an expression requires regeneration; varying an already bound parameter does not.

The built-in, overridable scalar helper
`dE(ma,mb,Ea,Eb,ax,ay,az,bx,by,bz)` is the rationalized on-shell energy difference

// docs-example: syntax
```text
((ma-mb)*(ma+mb)+(ax-bx)*(ax+bx)+(ay-by)*(ay+by)+(az-bz)*(az+bz))/(Ea+Eb)
```

For consistent on-shell inputs it equals `Ea-Eb`, while avoiding subtraction of nearly
equal rounded square roots. Its denominator must still be nonzero. A compact edge-specific
alias can capture the static geometry:

// docs-example: syntax
```toml
[function_map]
"de46()" = '''
dE(UFO::MT,UFO::MT,E(star,4),E(star,6),
   Q3(star,4,cind(1)),Q3(star,4,cind(2)),Q3(star,4,cind(3)),
   Q3(star,6,cind(1)),Q3(star,6,cind(2)),Q3(star,6,cind(3)))
'''
```

=== Which point does a multiplier see?

Let `base` be the cut-rescaled sample and let `left` and `right` denote the selected
threshold projections. `L` is a local counterterm, `I` an integrated contribution, and
`O` the unmodified side. Every variant evaluates its own projection; the two variants of
one threshold need not have equal star points.

#table(
  columns: (auto, 1fr, 1fr),
  table.header([*Component*], [*effective*], [*star*]),
  [`L × O`], [base], [left],
  [`I × O`], [left], [left],
  [`O × L`], [base], [right],
  [`O × I`], [right], [right],
  [`L × L`], [base], [merged left/right],
  [`I × L`], [left], [merged left/right],
  [`L × I`], [right], [merged left/right],
  [`I × I`], [merged left/right], [merged left/right],
)

Both multipliers in a pair see the same complete component-specific context. A missing
counterterm side leaves its spectator parent coordinates unchanged; cross-dependent
momenta are nevertheless rebuilt from the complete context. The coefficient of a pair is
the product of its two multiplier values. Exact zero can skip a completed component but
does not remove its threshold from overlap discovery or change a group's normalization.
All participating factors must be finite even if another factor is zero.

In particular, `eta(star,eset(...))` for the counterterm's own threshold vanishes at its
root. A partition built from that quantity alone may select a variant identically or
produce `0/0`, rather than distinguish the desired physical regions.

== Worked GL638 splitting

The maintained #source-link("examples/cli/epem_a_ttxh/NNLO/graphs/GL638.dot", label: "GL638 DOT")
uses explicit variants on cut `(2,6,10)` with full parent `[3,6,7,10]`.
Other cut sides retain maximal one-loop defaults. The selected one-loop space `[7]`
supplies compatible partners across cuts; the two-loop space `[3,7]` retains the
same-amplitude dual cancellation needed away from the relevant soft region.

In the center-of-mass frame, at each variant's own cut-preserving star, write
`Q = P0^0 + P1^0`, and label these residuals:

// docs-example: syntax
```text
H  = eta(2,4,12)       P  = eta(3,12)        Z = eta(3,10,13)
G0 = eta(2,6,12,13)    G4 = eta(2,4,10,13)
WH = H^2 / (H^2 + (P*Z/Q)^2)
WF = G0^2*G4^2 / (G0^2*G4^2 + P^2*Z^2)
```

This shorthand describes the host-shell physics; the actual definitions use `E`, `Q3`
and `dE` because foreign surfaces need not belong to the local catalogue. Off the host
shell, `wH`, `g0`, and `g4` are the corresponding foreign residual *minus* the host
`eta(2,6,10)` residual, and `wP = 2*E(star,3)-wQ()`. Replacing these definitions with
foreign `eta` calls would change their off-shell meaning.

#table(
  columns: (1fr, auto, auto),
  table.header([*Threshold edges on cut `(2,6,10)`*], [*`[7]` weight*], [*`[3,7]` weight*]),
  [`(7,8)` and `(8,12,14)`], [`WHc()`], [`WH()`],
  [`(8,10,13,14)`], [`1-WF()`], [`WF()`],
  [`(3,7,14)`, `(3,12)`, `(3,10,13)`], [absent], [1],
)

For instance the first threshold is declared in the embedded TOML as follows, using the
shared definitions already present in that DOT:

// docs-example: syntax
```toml
[[cuts]]
edges = [2, 6, 10]
[[cuts.thresholds]]
edges = [7, 8]
[[cuts.thresholds.counterterms]]
name = "shared_1l"
subspace = [7]
parent_lmb = [3, 6, 7, 10]
[cuts.thresholds.counterterms.multiplier]
expression = "WHc()"
[[cuts.thresholds.counterterms]]
name = "native_2l"
subspace = [3, 7]
parent_lmb = [3, 6, 7, 10]
[cuts.thresholds.counterterms.multiplier]
expression = "WH()"
```

No explicit `group_id` is needed. In the DOT, `WH` and `WHc` use complementary lazy
common-zero values zero and one. That assigns a value at their common zero; it is neither
a denominator regulator nor a proof of continuity. Algebraic complementarity at one point
also does not establish a pointwise partition between variants evaluated at different
stars. The shared soft and intersecting-threshold limits must be checked together.

A possible continuous tuning of a complementary pair is `w/(w+c*(1-w))` and
`c*(1-w)/(w+c*(1-w))`, with positive dimensionless `c`. This preserves the endpoints
and their vanishing orders for `0 <= w <= 1`; apply it consistently to partners whose
cancellation requires common weights. It does not remove a singular common-zero limit or
establish an optimal variance. Changing a valid subtraction prescription may change the
integral of the absolute value while preserving the signed integral. A sampling-only change
should preserve both for the same physical integrand and orientation-summation convention.

The top and Higgs masses separate some otherwise apparently competing soft/threshold
regimes. This motivates the Higgs-aware partition; it is not a topology-independent
recipe. The maintained #source-link("examples/cli/epem_a_ttxh/NNLO/graphs/GL297.dot", label: "GL297 DOT")
provides the simpler example with explicit one-loop associations across multiple cuts.

== Runtime localization and validation boundaries

The local threshold envelope uses the existing
`subtraction.local_ct_settings.uv_localisation` controls, despite that historical name:

// docs-example: syntax
```toml
[subtraction.local_ct_settings.uv_localisation]
gaussian_width = 1.0
sliver_width = 10.0
dynamic_width = false
smooth_sliver = false
force_uv_dampers_to_one = false
```

For signed radial displacement `delta = r-r_star`, the envelope is a Gaussian times
compact support `|delta| <= sliver_width*S`, with Gaussian width `gaussian_width*S`.
Here `S = E_cm` unless `dynamic_width = true`, when `S = r_star`. These settings are
dimensionless. `smooth_sliver = true` adds an even bump that is flat at the support
boundary and remains one at the pole.

The even envelope and both radial mirror terms preserve the local counterterm's zero
principal-value integral. An arbitrary asymmetric extra radial weight does not inherit
that property. Nor does evenness alone guarantee cross-cut IR cancellation: hard support
boundaries can separate nearby projected cut contributions. Inspect boundary approaches
as well as the interior; a smooth sliver removes the abrupt boundary but is not a general
variance optimization. Do not infer higher-power moving-boundary correctness from a
simple-cut check.

The integrated counterterm's normalized radial profile is
`h(r/r_star)/r_star`. Its `subtraction.integrated_ct_settings.range` and nested
`h_function_settings` are distinct from both the local envelope and the LU Cutkosky
h-function. The implemented integrated range is `type = "infinite"`; the compact variant
is not implemented. Preserve normalization when changing profile parameters.

Validation of a new prescription should compare CTs off, default CTs and the proposed
metadata at identical physical coordinates. Examine each cut and their sum, each split
variant, relevant left/right products, soft limits, intersecting-threshold limits and
localization boundaries. Include normal and higher-precision replays. A generic rotation
can expose cancellations missed by a special axis permutation; a stable complex norm does
not by itself bound the relative error in a tiny real component.

For a single soft three-momentum radius `lambda`, include `lambda^2 d lambda`; for two
simultaneously soft three-momenta at fixed ratio include `lambda^5 d lambda`. A fitted ray
power without that measure cannot decide integrability. Rays also do not certify angular,
hierarchical or tangential limits. All-orientation generation verifies construction, not
numerical cancellation; evaluate the complete orientation sum as a separate acceptance gate.

Two durable GL638 cautions remain relevant when interpreting such tests. Coincident radial
roots can expose a removable singularity as separate undefined residues; increasing
precision rescues nearby cancellation but does not define the exact common-root limit.
Separately, a hard two-normal corner with local behavior `f ~ 1/R` has a finite first
moment under `R dR`, but a smooth nonzero sampling density gives a logarithmic second
moment if the leading coefficient persists on an angular/tangential neighborhood. Correct
subtraction and effective sampling are therefore separate requirements. Surface-focused
sampling can address that corner; it does not prove bounded weights everywhere or replace
the threshold-prescription checks above.
]
