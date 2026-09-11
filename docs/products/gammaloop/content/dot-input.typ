#import "../../shared.typ": callout, boundary, product-link, source-link

#let dot-input = [
= GammaLoop DOT input

GammaLoop imports Feynman graphs from Graphviz `digraph` files. The generic DOT syntax and parser
behavior are covered in Linnet's #product-link("linnet", page: "reference/rust/linnet/parser/index.html", label: "Linnet DOT parser reference");
this page records GammaLoop's physics contract.

== Explore a DOT diagram

The notebook parses DOT with `linnet-py` and renders it through Linnest and Typst. It loads
when it comes into view, downloading browser Python and its dependencies on the first visit.
All subsequent Python execution and rendering happen in your browser.

Edit the DOT, choose *Force* or *Stable layered*, or toggle the Feynman styling and annotations.
The diagram updates after an edit; the force controls update when a slider is released. Below
the drawing, inspect the generated Typst and a table of the parsed edges.

The notebook supplies its own Python `DotCodec` and `RenderConfig`, including its `pᵢ`
annotations. GammaLoop's exported `just draw` bundle uses its DOT-only Typst figure template and
generated particle styles. The notebook does not select that template or validate particles and
interactions against the active GammaLoop model.

#context if target() == "html" {
  html.elem("div", attrs: (
    class: "live-notebook",
    "data-linnet-notebook": "physics_render_settings",
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

Try changing an internal edge's `particle`: the notebook draws `a` as a photon wave, `g` as a
gluon coil, `H` as a dashed scalar, and `t` as a fermion. Its codec translates grouped external
coordinates such as `pin="x:@-left"` and `pin="x:@+right"` into typed placements. Turning off
Feynman styling retains those placements while showing generic structural labels. Half-edge
IDs work in either style; the momentum and `pᵢ` / `nᵢ` controls apply to its Feynman view.

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
independently. Supplied indices must be unique and in range.

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
  [Edge], [`lmb_id`], [Edge index in a chosen loop-momentum basis.],
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
editing an exported graph; generated exports may refer to local `hedge(...)` placeholders.

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
