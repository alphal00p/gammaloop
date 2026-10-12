#import "../../shared.typ": callout, boundary, product-link, source-link

#let dot-input = [
= GammaLoop DOT input

GammaLoop imports Feynman graphs from Graphviz `digraph` files. The generic DOT syntax and parser
behavior are covered in Linnet's #product-link("linnet", page: "guides/dot-input/", label: "DOT input and graph indices");
this page records GammaLoop's physics contract.

== How GammaLoop reads DOT

`import graphs PATH` reads every `digraph` of a DOT file, or of every `.dot` file in a directory
in file-name order; `import graphs --inline-dot '<DOT>'` reads a string instead. GammaLoop reads
each graph with FeynKit's DOT importer, the one behind `FeynmanDiagram.from_dot`, against the
active model. It accepts two dialects: the compact physics DOT this page describes, and the
annotated DOT written by FeynKit's `FeynmanDiagram.to_dot()`. The annotated dialect is
recognized by graph attributes such as `model_fingerprint`; the fingerprint must match the active
model, and the diagram is restored exactly, including its loop momentum basis and physical cuts.

FeynKit matches particles and interactions, sews cross-section initial states, finds physical
cuts and routes momenta. It then applies the model's Feynman rules to every factor the input
omits: each vertex takes its interaction's vertex factor, each internal edge its propagator
numerator, and the graph the external-state projector built from the polarization vectors and
spinors of its external particles. A vertex or internal-edge `num`, or a graph `projector`,
replaces the corresponding factor; it must use the tensor indices of the model's rules, described
under Numerators and loop bases. GammaLoop refuses a graph whose numerators, projector and
factors together leave a tensor index uncontracted, apart from the colour indices of external
legs. It finally applies the attributes FeynKit does not
interpret: `group_id`, `is_group_master` and `params` on the graph, `dod` on vertices, and
`name`, `mass`, `dod`, `momtrop_edge_power` and `vakint_edge_power` on edges. Any other attribute
of compact DOT, apart from Graphviz presentation attributes, produces a warning that suggests the
canonical spelling of a likely typo or older alias.

== Explore a DOT diagram

Edit an amplitude or cross-section diagram with FeynKit. The notebook parses compact model-aware DOT
into a `FeynmanDiagram` and uses its interactive renderer for the preview.

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

The editor accepts compact physics DOT and the annotated format emitted by
`FeynmanDiagram.to_dot()`. FeynKit resolves `particle` or `pdg` against the bundled Standard
Model and infers interaction slots. Supply `int_id` if the incident particles match several
interactions. Invalid particles, interactions, cut pairs, or loop bases produce a parsing
error in the notebook.

The amplitude retains the photon/top/gluon example and its explicit external half-edge order.
The cross-section describes electron–positron annihilation into muons. Matching `is_cut` tags
pair its initial-state legs, while `graph [final_state="mu-,mu+"]` requests physical cuts with
those final-state particles. FeynKit sews the initial-state legs, finds matching separating
cuts, and routes momenta through the resulting diagram. The preview opens initial-state
connections into matching left/right legs by default. Turn off *Split initial state* to
draw the sewn graph; the parsed diagram and physical final-state cut remain unchanged. The
preview shows the parsed diagram before the model's Feynman rules fill in omitted numerators and
the projector, and without GammaLoop's evaluator attributes; `import graphs` adds both.

The initial layout uses 100 steps per epoch, 30 epochs, and seed 42, with momentum arrows and
labels disabled. Enable momentum labels to display the parsed diagram's stored loop momentum
basis. The preview and edge table both read the parsed FeynKit diagram. Hover to inspect
physics information, and Shift-click to select vertices, edges, or half-edges near an internal
edge's ends. The SVG has a transparent background and follows the documentation theme.

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
separate index spaces, and an amplitude keeps both; sewing a cross section renumbers them, as
described below. Assign only indices needed by momentum or numerator references; omitted
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
  [Numeric port `vertex:h`], [Sets initial row order in both drawing modes and identifies half-edge debug labels.],
    [Orders amplitude external momentum and helicity entries; also supplies tensor-index identities.],
  [Edge `is_cut=k`], [Groups matching cross-section legs at a shared, movable vertical position.],
    [Marks a cross-section initial state and orders its external momentum entries; equal tags sew split legs.],
  [Edge `lmb_id=k`], [Does not rename the edge's $q$ label.],
    [Selects this edge as the carrier of loop momentum `K(k, ...)`.],
)

DOT arrow direction determines incoming versus outgoing, independently of these numbers and of
`dir`. The particle fixes the edge orientation: an edge's own `dir` must agree with it, while a
default `edge [dir=...]` only styles the drawing. For an amplitude, the runtime
external momentum and helicity lists follow increasing dangling half-edge index. Incoming
momenta enter the conservation equation with a plus sign and outgoing momenta with a minus sign.
Automatic placement in both modes uses this half-edge order, separately for the incoming left
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

A cross-section input marks each initial state with an `is_cut` tag, either on one internal edge
or on two dangling legs with the same tag. Every external leg of a cross section needs a tag. In
the split form, one leg must be incoming and the other outgoing, both of the same particle. They
still have distinct edge IDs and distinct half-edge IDs. A tagged internal edge is split at
import: its target receives the incoming leg and its source emits the outgoing one. For example,
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
invalid.

On import, FeynKit sews matching legs into one edge, then sorts the cut tags numerically and
moves the resulting cut edges and their incoming half-edges to leading indices `0, 1, ...`.
These cut entries determine the cross-section's external momentum order. The other edges follow
in DOT order, and the other half-edges are renumbered too. Indices in `num`, `projector` and
`overall_factor` expressions refer to the DOT input; FeynKit transports them to the sewn
indices. Exports conventionally use the incoming cut half-edge index as the tag. On import, however, tags are non-negative integer matching and ordering
values: they need not be contiguous or name an existing half-edge. The index bounds described
above apply to `id` and numeric ports, not to these tags.

The graph attribute `final_state` lists the particles of the physical cuts, comma-separated, for
example `graph [final_state="mu-,mu+"]`. FeynKit keeps the cuts separating the incoming from the
outgoing initial-state legs whose cut edges carry exactly these particles, and refuses a graph
without any. Without `final_state`, every separating cut is physical. `import graphs
--process-spec '<process>'` then keeps only the cuts matching that process definition, for
example `--process-spec 'e+ e- > mu+ mu- [{{1}}]'`.

`just draw` renders compact DOT directly, without this physics import or sewing. The notebook
instead imports compact DOT with FeynKit and selects physical cuts using `final_state`. Both renderers give incoming legs one shared, movable X coordinate on the left and outgoing legs
another on the right. Initial Y positions follow half-edge order separately on each side.
Cross sections additionally group Y coordinates by numeric `is_cut` identity: matching tags
share a freely moving row, seeded from the mean of their starting positions. The tag's numeric
order does not determine drawing order. The existing dangling-centroid repulsion spreads
external endpoints in both modes; its default strength is `gamma-dangling-centroid=1.25`.
Neither mode automatically pins X or Y. All dangling depths are pinned to zero, and explicit
horizontal and vertical placements retain precedence. A tag without a counterpart still gets
its side's X group and its own movable Y group; an untagged leg keeps an independent Y coordinate.

In `just draw`, generated momentum labels use the current drawing edge's `eid`, including on split cut legs: the
example's labels are $q_0$ and $q_5$. With edge IDs fixed, changing `is_cut`, a half-edge port, or
`lmb_id` does not rename these labels. After importing and saving a graph, the labels instead
reflect the edge IDs in that newly exported DOT. Do not infer kinematic-list positions from the drawing's $q$ subscripts.

== Physics attributes

#table(
  columns: (auto, auto, 1fr),
  table.header([*Location*], [*Fields*], [*Meaning*]),
  [Graph], [`num`], [Global numerator factor, default `1`.],
  [Graph], [`projector`], [Replaces the external-state projector built from the model.],
  [Graph], [`overall_factor`], [Graph weight, default `1`.],
  [Graph], [`final_state`], [Particles of a cross section's physical cuts; default: every separating cut.],
  [Graph], [`params`, `group_id`, `is_group_master`], [Additional parameters (separated by `;`) and grouping metadata.],
  [Node], [`int_id`], [UFO interaction name; required when the incident particles match several rules. It must name a rule of the incident particles, also when the vertex has a `num`.],
  [Node], [`num`], [Replaces the interaction's vertex factor; required for a vertex matching no rule.],
  [Node], [`dod`], [Degree-of-divergence override.],
  [Edge], [`particle`, `pdg`], [Model particle name or signed UFO PDG code; every edge needs one.],
  [Edge], [`dir`], [Must agree with the orientation that the particle fixes.],
  [Edge], [`name`, `mass`], [Optional stable name, and a Symbolica mass replacing the particle's in the propagator denominator; a numerator the model fills in keeps the model's mass, so give `num` to change it there too.],
  [Edge], [`lmb_id`], [Loop-momentum slot carried by this internal edge.],
  [Edge], [`num`], [Replaces an internal edge's propagator numerator; ignored on external legs.],
  [Edge], [`is_cut`], [Cross-section initial-state tag.],
  [Edge], [`is_dummy`], [Dummy edge, outside interaction matching and the external legs.],
  [Edge], [`dod`, `momtrop_edge_power`, `vakint_edge_power`], [Evaluator controls.],
)

For fermions, signed particles and edge orientation must give a consistent fermion flow around
closed loops and along open fermion lines. A dummy edge takes the vertex slots after the
interaction's legs. External momentum labels are not a current GammaLoop
physics attribute; older examples using `mom` should use edge `particle`/`pdg` and half-edge ports.
The active model supplies propagator tensors; the legacy `propagator` field is not a current
GammaLoop input field. `multiplicity_factor` is likewise a legacy field; current grouping uses
`group_id`, `is_group_master`, and `overall_factor`. GammaLoop no longer reads the `particle_id`,
`vertex_rule_id`, `initial_state_connection` and half-edge `ufo_order` fields of its earlier DOT
exports: name particles with `particle` or `pdg` and interactions with `int_id`, and let FeynKit
order the vertex slots.

== Numerators and loop bases

A graph, vertex, or edge `num` is parsed as a Symbolica expression, as are `projector` and
`overall_factor`. Quote expressions containing punctuation or spaces. Unqualified symbols resolve
in the `gammalooprs` namespace, where `Q`, `K`, `P`, `hedge`, `edge` and the polarization heads
live; qualify the others, as in `spenso::mink(4, hedge(3, 1))`. Keep vertex and edge numerators
attached to their owning statement when editing an exported graph. In an edge numerator, `eid`,
`source`, and `sink` refer to that edge's index and its endpoint half-edge indices in the
imported graph, after any sewing; their function forms construct the corresponding tensor-index
identities: `source(1)` is `hedge(h, 1)` for the source half-edge `h`, and `source(0)` is
`hedge(h)`. This local shorthand is resolved for edge numerators, not for an arbitrary graph or
vertex expression.

An explicit factor contracts with the factors the model supplies only through the same tensor
indices. The model's rules give every half-edge its own indices: the k-th index of a Lorentz,
spinor or colour representation on half-edge `h` is `hedge(h, k)`, counted from one, so a gluon
leg on half-edge `4` carries `spenso::mink(4, hedge(4, 1))` and `spenso::coad(8, hedge(4, 1))`.
The external-state projector uses the same identities, for example
`u(0, spenso::bis(4, hedge(0, 1)))` for an incoming fermion on edge `0` and half-edge `0`.
`hedge(h)` is a different index, which older GammaLoop graphs used for the first one: a
`projector` or `num` written with it leaves the model's indices open, and import stops with an
error naming the uncontracted indices. Graphs whose factors are all explicit only need to be
consistent among themselves. Dummy indices local to one factor use `vertex(v, k)` or
`edge(e, k)`. The model's projector holds no colour wave functions, so, as in a generated
amplitude, the colour indices of external legs may stay open; evaluating such an amplitude
needs a colour sum or average. For an incoming quark on half-edge `0` and an incoming antiquark
on half-edge `1`, the graph attribute
`num="spenso::g(spenso::dind(spenso::cof(3, hedge(0, 1))), spenso::cof(3, hedge(1, 1)))/3"`
averages it. Cross sections contract colour by sewing.

By contrast, generated `hedge(h)` or `hedge(h, k)` expressions contain explicit tensor-index
identities: `h` is not a local port number within the statement. Preserve the exported topology,
indices, and numerator expressions together. Sewing can move graph records without renaming
those abstract tensor identities, so changing only the printed ports or only numbers in `num`
is not a safe way to renumber an exported calculation.

Edge momenta `Q(e, ...)`, loop momenta `K(k, ...)`, and external momenta `P(p, ...)` use different
indices. Setting `lmb_id=k` chooses the internal edge carrying `K(k, ...)`; it does not set the
edge ID or an external momentum slot. The slots must run from `0` without gaps and complete a
loop momentum basis. Without any `lmb_id`, FeynKit routes momenta through its canonical basis for
the physical cuts it keeps, which keeps beam momenta off cut edges where possible. This is
generation's basis when those cuts are the generated ones. A cross section imported without
`final_state` keeps every separating cut while routing, and `--process-spec` selects among them
only afterwards, so its basis can differ from the generated one; `lmb_id`, which `save dot`
writes, restores the exact basis. The remaining edge
momenta are expressed through the chosen basis and external flows. Exports show these
expressions in `lmb_rep`; they are generated routing information that import ignores, not a
replacement for the `lmb_id` input.

A cross section rescales its loop momenta about the origin of this basis, where each cut edge
carries only its external momentum. The cut energies there must add up to less than the incoming
energy for every cut: a basis whose cut edges carry part of the beam momenta, such as `K + P1`,
can violate this, and GammaLoop then refuses the integrand before evaluating it. A kinematically
closed cut, or incoming momenta outside their centre-of-mass frame, can violate it whatever the
basis; the refusal names which cause applies.

Load a DOT file after importing its model:
// docs-example: syntax

```text
import model assets/models/json/sm/sm.json
import graphs ./my_graphs.dot
```

Save loaded graphs back to DOT with `save dot ./my_graphs_exported`. The export is compact DOT
with every factor explicit, which `import graphs` reads again: particle names, `int_id`
interaction names, vertex and edge `num`, `projector`, `overall_factor`, `is_cut` tags, `lmb_id`,
`params`, and the vertex and edge attributes GammaLoop applies. It does not record `group_id`,
`is_group_master` or `final_state`: a re-imported cross section keeps every separating cut until
`--process-spec` or its process definition selects the physical ones. The current spelling is
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
