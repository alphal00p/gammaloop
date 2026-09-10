#import "../../shared.typ": callout, boundary, product-link

#let dot-input = [
= GammaLoop DOT input

GammaLoop imports Feynman graphs from Graphviz `digraph` files. The generic DOT syntax and parser
behavior are covered in Linnet's #product-link("linnet", page: "reference/rust/linnet/parser/index.html", label: "Linnet DOT parser reference");
this page records GammaLoop's physics contract.

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
]
