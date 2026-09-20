#import "../../shared.typ": callout, boundary, product-link, source-link

#let dot-input = [
= GammaLoop DOT input

GammaLoop imports Feynman graphs from Graphviz `digraph` files. The generic DOT syntax and parser
behavior are covered in Linnet's #product-link("linnet", page: "guides/dot-input/", label: "DOT input and graph indices");
this page records GammaLoop's physics contract.

== Explore a DOT diagram

Edit an amplitude or cross-section diagram using the same rendering pipeline as `save dot`
followed by `just draw`.

The notebook loads automatically when it comes into view. The first visit downloads browser
Python and its dependencies; rendering then runs locally in your browser.

#context if target() == "html" {
  html.elem("div", attrs: (
    class: "live-notebook",
    "data-notebook": "physics_render_settings",
    "aria-label": "Physics DOT rendering notebook",
  ))[
    #html.elem("p", attrs: (class: "live-notebook-fallback"))[
      This interactive notebook requires JavaScript and this documentation build's notebook
      assets. The #source-link("crates/linnet-py/examples/physics_render_settings.py", label: "physics rendering notebook source")
      is also available to run locally.
    ]
  ]
} else {
  [The online guide includes an interactive notebook. Open the
  #source-link("crates/linnet-py/examples/physics_render_settings.py", label: "physics rendering notebook source")
  to run the same example locally.]
}

Try changing an internal edge's `particle`: `a` draws a photon wave, `g` a gluon coil, `H` a
dashed scalar, and `t` a fermion with its model label. The same generated particle map and
Typst callbacks control these styles in exported GammaLoop drawings.

The initial settings match ordinary `just draw`: 100 steps per epoch, 30 epochs, seed 42, and momentum arrows
and labels disabled. The preview draws edited DOT without importing it into a GammaLoop
calculation, so it does not validate interactions or numerators against the model.

The amplitude assigns half-edge ports `5:12` and `4:15` to put the outgoing leg at vertex `5`
above the leg at vertex `4`. This external order allows the default layout to draw the internal
cycle without crossings. The ports are global half-edge indices; the remaining indices are
inferred. Reordering the DOT statements alone would not establish this order.

== Graph shape and half-edges

// docs-example: syntax
```text
digraph box {
  ext [style=invis];
  ext -> v0:0 [particle="e-"];
  v0 -> v1 [id=1, particle="a"];
  v1 -> ext [particle="e+"];
}
```

Use `digraph NAME { ... }`. DOT arrows define the underlying edge orientation. An invisible node
represents an external dangling half-edge: invisible-node to vertex is incoming, and vertex to
invisible-node is outgoing. An edge between two invisible nodes is invalid.

Use one shared `ext [style=invis]` node for all external legs. Each connection still creates a
separate dangling half-edge. Keep the shared node otherwise attribute-free: attributes on a
dangling node are merged into every external edge, so put `particle`, `num`, and other per-leg data
on the edge.

Numeric ports (`vertex:0`) fix half-edge indices; an edge `id=1` fixes its edge index. These are
separate index spaces. Assign only indices needed by momentum or numerator references; omitted
indices are allocated automatically. Either endpoint of an internal edge can be indexed
independently. Supplied indices must be unique and in range across the entire graph, including
external legs. A numeric port is not a slot numbered locally within one vertex. Omitted indices
follow parser order, which need not be DOT statement order; see
#product-link("linnet", page: "guides/dot-input/", label: "DOT input and graph indices") for the allocation rules.

#table(
  columns: (auto, 1fr, 1fr),
  table.header([*Input*], [*Drawing*], [*GammaLoop calculation*]),
  [Edge `id=n`], [Names the momentum label $q_n$.],
    [Identifies an edge in momentum expressions such as `Q(n, ...)`.],
  [Numeric port `vertex:h`], [Sets initial amplitude row order and identifies half-edge debug labels.],
    [Orders amplitude external momentum and helicity entries; also supplies tensor-index identities.],
  [Edge `is_cut=k`], [Orders cross-section rows numerically and pairs legs at the same vertical position.],
    [Marks an initial-state cut and orders its external momentum entries; equal tags sew split legs.],
  [Edge `lmb_id=k`], [Does not rename the edge's $q$ label.],
    [Selects this edge as the carrier of loop momentum `K(k, ...)`.],
)

DOT arrow direction determines incoming versus outgoing, independently of these numbers and of
`dir`, which controls the displayed or particle orientation. For an amplitude, the runtime
external momentum and helicity lists follow increasing dangling half-edge index. Incoming
momenta enter the conservation equation with a plus sign and outgoing momenta with a minus sign.
Automatic amplitude placement uses this same half-edge order, separately for the incoming left
column and outgoing right column. Gaps in the indices do not create empty rows. These are initial
vertical positions that layout may adjust; explicit positions still take precedence. Moving a
port index to another leg changes its kinematic-list position and its initial drawing order.
See #product-link("gammaloop", page: "guides/conventions/", label: "Kinematics, normalization, and weights")
for the four-vector and sign conventions.

For example, the edge and half-edge indices deliberately differ here:

// docs-example: syntax
```text
digraph amplitude_order {
  ext [style=invis];
  ext -> v:3 [id=0, particle="H"];
  v:2 -> ext [id=1, particle="H"];
  v:1 -> ext [id=2, particle="H"];
  ext -> v:0 [id=3, particle="H"];
}
```

The incoming half-edges are `0, 3`, so the left column starts with $q_3$ above $q_0$.
The outgoing half-edges are `1, 2`, so the right column starts with $q_2$ above $q_1$.
The runtime external list follows half-edges `0, 1, 2, 3`, belonging to edges `3, 2, 1, 0`.
This ordering does not renumber either collection: edge IDs continue to determine the $q$ labels.

== Cross-section cuts and drawing order

A cross-section input may describe an initial-state cut as an internal edge carrying `is_cut`,
or as two dangling legs with the same `is_cut` value. In the split form, one leg must be incoming
and the other outgoing. They still have distinct edge IDs and distinct half-edge IDs. For example,
these two statements in a larger graph identify one pair:

// docs-example: syntax
```text
ext -> left:0 [id=0, is_cut=0, particle="e-"];
right:7 -> ext [id=5, is_cut=0, particle="e-"];
```

Here `ext` must be declared invisible, and the enclosing graph must have enough edges and
half-edges for the explicit indices. Sharing `id=0` between the two statements would be an error;
sharing `is_cut=0` is what associates them. Use one tag per initial-state pair, with consistent
particle data on its two legs. A cut-tagged dangling leg without its opposite-flow partner is
invalid for a GammaLoop calculation.

On import, GammaLoop sews matching legs into one edge, then sorts the cut tags numerically and
moves the resulting cut edges and their incoming half-edges to leading indices `0, 1, ...`.
These cut entries determine the cross-section's external momentum order. Sewing and this
normalization can change other edge and half-edge indices too. Exports conventionally use the
incoming cut half-edge index as the tag. On import, however, tags are non-negative integer matching and ordering
values: they need not be contiguous or name an existing half-edge. The index bounds described
above apply to `id` and numeric ports, not to these tags.

The notebook and `just draw` render DOT directly, without this physics import or sewing.
Cross-section placement follows the same numeric `is_cut` order as the calculation: incoming
legs occupy rows from top to bottom, and outgoing legs find their row through the matching tag.
For example, tag `2` precedes tag `10`; gaps do not create empty rows. Incoming legs go on the left
and outgoing legs on the right. These cross-section rows are pinned vertically, while amplitude
rows use the adjustable starting positions described above. All dangling depths are pinned to
zero, and explicit horizontal and vertical placements retain precedence. An outgoing cut tag
without an incoming counterpart can still be drawn, but receives no automatic row placement.

Generated momentum labels use the current drawing edge's `eid`, including on split cut legs: the
example's labels are $q_0$ and $q_5$. With edge IDs fixed, changing `is_cut`, a half-edge port, or
`lmb_id` does not rename these labels. After importing and saving a graph, the labels instead
reflect the edge IDs in that newly exported DOT. Do not infer kinematic-list positions from the drawing's $q$ subscripts.

== Physics attributes

#table(
  columns: (auto, auto, 1fr),
  table.header([*Location*], [*Fields*], [*Meaning*]),
  [Graph], [`num`], [Global numerator factor, default `1`.],
  [Graph], [`overall_factor`, `projector`, `params`, `group_id`, `is_group_master`], [Advanced imported-state and grouping metadata.],
  [Node], [`int_id`], [UFO interaction identifier when incident particles do not uniquely select a rule.],
  [Node], [`num`, `dod`], [Vertex numerator override and degree-of-divergence override.],
  [Edge], [`particle`, `pdg`], [Model particle name or signed UFO PDG code.],
  [Edge], [`dir`], [Explicit orientation; otherwise fermion orientation is derived from the signed particle.],
  [Edge], [`name`, `mass`], [Optional stable name and Symbolica mass override.],
  [Edge], [`lmb_id`], [Loop-momentum slot carried by this edge.],
  [Edge], [`num`, `is_cut`, `is_dummy`, `momtrop_edge_power`, `vakint_edge_power`], [Numerator and evaluator controls.],
)

For fermions, signed particles and edge orientation must give a consistent fermion flow around
closed loops and along open fermion lines. External momentum labels are not a current GammaLoop
physics attribute; older examples using `mom` should use edge `particle`/`pdg` and half-edge ports.
The active model supplies propagator tensors; the legacy `propagator` field is not a current
GammaLoop input field. `multiplicity_factor` is likewise a legacy field; current grouping uses
`group_id`, `is_group_master`, and `overall_factor`.

== Numerators and loop bases

A graph, vertex, or edge `num` is parsed as a Symbolica expression. Quote expressions containing
punctuation or spaces. Keep vertex and edge numerators attached to their owning statement when
editing an exported graph. In an edge numerator, `eid`, `source`, and `sink` refer to that edge's
index and its existing endpoint indices; their function forms, such as `source(0)`, construct
corresponding tensor-index identities. This local shorthand is resolved for edge numerators,
not for an arbitrary graph or vertex expression.

By contrast, generated `hedge(h)` or `hedge(h, k)` expressions contain explicit tensor-index
identities: `h` is not a local port number within the statement. Preserve the exported topology,
indices, and numerator expressions together. Sewing can move graph records without renaming
those abstract tensor identities, so changing only the printed ports or only numbers in `num`
is not a safe way to renumber an exported calculation.

Edge momenta `Q(e, ...)`, loop momenta `K(k, ...)`, and external momenta `P(p, ...)` use different
indices. Setting `lmb_id=k` chooses the internal edge carrying `K(k, ...)`; it does not set the
edge ID or an external momentum slot. The remaining edge momenta are expressed through the
chosen basis and external flows. Exports show these expressions in `lmb_rep`; they are generated
routing information, not a replacement for the `lmb_id` input.

Load a DOT file after importing its model:
// docs-example: syntax

```text
import model assets/models/json/sm/sm.json
import graphs ./my_graphs.dot
```

Save loaded graphs back to DOT with `save dot ./my_graphs_exported`. GammaLoop fills omitted names,
indices, factors, and loop-momentum basis entries during import/export. The current spelling is
`lmb_id`; older wiki examples may use `lmb_index`.

#callout("Keep external data on edges", [
  A shared invisible node is deliberately compact. Put each external particle, numerator, and
  cut flag on its own connecting edge so those attributes do not get copied to every leg.
])

For PDF figures, run `just draw` from the directory written by `save dot`. See
#product-link("linnet", label: "Clinnet DOT rendering", page: "guides/clinnet/") for renderer
configuration, or #product-link("linnet", label: "Rendering graphs from Python", page: "guides/python-rendering/")
to build graphs with arbitrary Python data and custom drawing selectors.
]
