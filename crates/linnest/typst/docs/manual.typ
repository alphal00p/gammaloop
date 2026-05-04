#import "@preview/tidy:0.4.3"
#import "../src/lib.typ": draw, graph, layout, subgraph

#set document(title: "Linnest Typst API")
#set page(margin: 22mm)
#set text(size: 10pt)

= Linnest Typst API

Linnest exposes the `linnest.wasm` graph layout plugin through a small Typst
API. Its drawing layer bundles `kurvst.wasm`; the standalone Kurvst package
manual documents the full curve and path-pattern API.
The public surface is intentionally narrow:

- `graph` for construction, parsing, inspection, joins, and graph algorithms.
- `subgraph` for subgraph object construction and inspection.
- `layout` for the separate layout pass.
- `draw` for rendering a laid-out graph object with CeTZ.

== Minimal Builder Example

```typ
#import "../src/lib.typ": draw, graph, layout, subgraph

#let b = graph.builder(name: "demo")
#let (node: a, builder: b) = graph.node(b, name: "a")
#let (node: c, builder: b) = graph.node(b, name: "c")
#let b = graph.edge(
  b,
  source: (node: a, compass: "e"),
  sink: (node: c, compass: "w"),
  statements: (
    color: "0055ff",
    source-color: "d72638",
    sink-color: "1b7f4c",
    label: "a-c",
  ),
)
#let g = graph.finish(b)
#let g = layout(g)
#let east = subgraph.compass(g, "e")
#let edges = graph.edges(g, subgraph: east)
#let dot = graph.dot(g)
#let edge-label(edge) = text(fill: rgb("#" + edge.color))[#edge.label]
#let source-style(edge) = (stroke: rgb("#" + edge.source-color) + 0.5pt)
#let sink-style(edge) = (stroke: rgb("#" + edge.sink-color) + 0.5pt)
#draw(g, subgraph: east, edge-label: edge-label, source-style: source-style, sink-style: sink-style)
```

#import "../src/lib.typ": draw, graph, layout, subgraph

#let b = graph.builder(name: "demo")
#let (node: a, builder: b) = graph.node(b, name: "a")
#let (node: c, builder: b) = graph.node(b, name: "c")
#let b = graph.edge(
  b,
  source: (node: a, compass: "e"),
  sink: (node: c, compass: "w"),
  statements: (
    color: "0055ff",
    source-color: "d72638",
    sink-color: "1b7f4c",
    label: "a-c",
  ),
)
#let g = graph.finish(b)
#let g = layout(g)
#let east = subgraph.compass(g, "e")
#let edges = graph.edges(g, subgraph: east)
#let dot = graph.dot(g)
#let edge-label(edge) = text(fill: rgb("#" + edge.color))[#edge.label]
#let source-style(edge) = (stroke: rgb("#" + edge.source-color) + 0.5pt)
#let sink-style(edge) = (stroke: rgb("#" + edge.sink-color) + 0.5pt)
#draw(g, subgraph: east, edge-label: edge-label, source-style: source-style, sink-style: sink-style)

== Graph Objects

Graph, builder, and subgraph objects are opaque zero-copy values. Build or parse
graph objects with `graph`, transform graph objects with `layout`, and pass
objects back to `graph` or `subgraph` for inspection.

- `graph.parse(input)` parses one or more DOT digraphs and returns an array of
  graph objects.
- `graph.build(..)` constructs one graph object from named arguments. This is
  convenience sugar over the builder functions.
- `graph.builder(..)` starts a builder object.
- `graph.node(builder, ..)` returns `(node: index, builder: builder)`.
- `graph.edge(builder, ..)` returns the updated builder.
- `graph.finish(builder)` turns a builder into a graph.
- `layout(graph, ..)` runs layout as an explicit
  second step. Its settings are named parameters so calls stay descriptive and
  Tidy can document each field.
- `draw(graph, ..)` draws a laid-out graph object with CeTZ.
- `graph.dot(graph)` returns a DOT string for inspection or export.

== Builder Spec

The builder API is the preferred construction API when nodes need to be reused.
Use destructuring to keep the builder value moving:

```typ
#let (node: a, builder: b) = graph.node(b, name: "a")
```

`graph.edge` accepts `source` and `sink` endpoint dictionaries. An endpoint can
contain `node`, `statement`, `id`, `port_label`, `compass`, and `in_subgraph`.
Set `source: none` or `sink: none` to create an external half edge.

String-valued DOT statements may contain `{name}` placeholders. `graph.build`
and the builder object expand each placeholder while constructing the graph
object. Use `{{` and `}}` for literal braces. Unknown placeholders are left
unchanged.

This is useful for graph-level `edge-statements`: the default edge statement is
merged into each edge, then expanded against that edge's complete statement
dictionary. The following default metadata records a display label derived from
the per-edge `label` statement:

```typ
#let b = graph.builder(
  edge-statements: (
    color: "000000",
    display_label: "{label}",
  ),
)
#let (node: a, builder: b) = graph.node(b, name: "a")
#let (node: c, builder: b) = graph.node(b, name: "c")
#let b = graph.edge(
  b,
  source: (node: a),
  sink: (node: c),
  statements: (
    color: "0055ff",
    label: "a-c",
  ),
)
```

Draw styling is Typst-native. Pass dictionaries or callbacks to `draw`; edge
callbacks receive the merged `scope`, edge statements, endpoint records, and the
edge index:

```typ
#let edge-label(edge) = text(fill: rgb("#" + edge.color))[#edge.display_label]
#let source-style(edge) = (stroke: red + 0.5pt)
#let sink-style(edge) = (stroke: blue + 0.5pt)
#draw(layout(graph.finish(b)), edge-label: edge-label, source-style: source-style, sink-style: sink-style)
```

`graph.build` accepts the same data in one dictionary:

```typ
#let g = graph.build(
  name: "demo",
  statements: (full_num: "x + y"),
  node-statements: (shape: "circle"),
  edge-statements: (
    color: "000000",
    display_label: "{label}",
  ),
  nodes: ((name: "a"), (name: "b")),
  edges: ((
    source: (node: 0, compass: "e"),
    sink: (node: 1, compass: "w"),
    statements: (color: "0055ff", label: "ab"),
  ),),
)
```

== Graph Queries

`graph.info(g)` returns graph metadata. `graph.nodes(g)` returns node records,
and `graph.edges(g)` returns edge records. Pass `subgraph: sg` to filter nodes
or edges by an subgraph object.

`graph.join(left, right, key: "statement")` joins matching dangling half edges.
The key is read from half-edge data and can be `"statement"`, `"port_label"`,
`"compass"`, or `"id"`.

`graph.cycles(g)` returns subgraph objects for a cycle basis.
`graph.forests(g)` returns subgraph objects for spanning forests.

== Layout Model

`layout` starts from a traversal-tree placement, then optimizes the positions of
graph nodes and edge control points. The initial tree spacing is

$ L = lambda sqrt((W H) / max(n, 1)) $,

with horizontal spacing $tau_x L$ and vertical spacing $tau_y L$. Here
$lambda$ is `length_scale`, $W$ is `viewport_w`, $H$ is `viewport_h`, $tau_x$
is `tree_dx`, and $tau_y$ is `tree_dy`. These fields set the geometry scale for
both layout modes.

The shared spring/charge model uses the following coefficients:

$ c_("vv") = beta L^2 $
$ c_("ev") = beta gamma_("ev") L^2 $
$ c_("ee") = beta gamma_("ee") L^2 $
$ c_("center") = beta g_("center") L^2 $
$ c_("dangling") = beta gamma_("dangling") L^2 $

The Typst parameter names are `beta` for $beta$, `gamma_ev` for
$gamma_("ev")$, `gamma_ee` for $gamma_("ee")$, `g_center` for $g_("center")$,
and `gamma_dangling` for $gamma_("dangling")$. The spring stiffness $k$ is
`k_spring`, and the softening constant $epsilon$ is `eps`.

In `layout_algo: "anneal"`, linnest minimizes an energy:

$ E =
sum_(i < j) 1/2 c_("vv") / (d(v_i, v_j) + epsilon)
+ sum_(i, e) c_("ev") / (d(v_i, e) + epsilon)
+ sum_((v, e) " incident") 1/2 k (ell_e - d(v, e))^2
\ + sum_("local edge pairs") 1/2 c_("ee") / (d(e_i, e_j) + epsilon)
+ sum_("dangling pairs") 1/2 c_("dangling") / (d(e_i, e_j) + epsilon)
+ sum_(i) 1/2 c_("center") / (d(v_i, 0) + epsilon)
+ p_("cross") N_("cross") $.

Here $p_("cross")$ is `crossing_penalty` and $N_("cross")$ is the number of
detected edge crossings. `temp`, `step`, `seed`, `steps`, `epochs`, `cool`,
`accept_floor`, `step_shrink`, and `incremental_energy` belong to this
simulated annealing mode. `crossing_penalty` is also anneal-only; force mode
does not currently add a crossing force.

In `layout_algo: "force"`, linnest applies the direct forces corresponding to
the same vertex-vertex, edge-vertex, incidence spring, local edge-edge,
dangling-edge, and center terms. `step` is the integration step, `delta` clamps
per-step movement, `steps` and `epochs` set the iteration budget, `cool` shrinks
the step after each epoch, and `early_tol` stops when movement is small.
`z_spring` and `z_spring_growth` are force-only helpers: the integrator gives
points temporary z coordinates to break overlaps and pulls them back toward the
2D plane.

`directional_force` is applied in both modes as an extra bias derived from
pin/port direction constraints.

After either graph layout mode, labels are relaxed separately. If $L_l$ is the
label target distance and $q_l$ is the label repulsion strength, then
$L_l = alpha_l L$ and $q_l = beta_l L^2$. The Typst names are
`label_length_scale` for $alpha_l$ and `label_charge` for $beta_l$.
`label_spring` is the spring constant pulling each label toward its target.
`label_steps`, `label_step`, `label_early_tol`, and `label_max_delta_scale`
control the label relaxation iteration.

== Subgraphs

Subgraph objects are opaque zero-copy values.

- `subgraph.label(g, label)` constructs a subgraph from a base62 label.
- `subgraph.bits(g, bits)` constructs a subgraph from a boolean hedge array.
- `subgraph.compass(g, compass)` selects half edges with a DOT compass point.
- `subgraph.to-label(sg)` returns the base62 label.
- `subgraph.hedges(sg)` returns included hedge indices.
- `subgraph.contains(sg, hedge)` tests hedge membership.

== Generated Reference

#let tidy-style = dictionary(tidy.styles.default)
#let _ = tidy-style.insert("show-example", tidy-style.show-example.with(scale-preview: 100%))

#let docs = tidy.parse-module(
  read("../src/lib.typ"),
  scope: (draw: draw, graph: graph, layout: layout, subgraph: subgraph),
)
#tidy.show-module(docs, style: tidy-style)

#let graph-docs = tidy.parse-module(
  read("../src/graph.typ"),
  name: "graph",
  scope: (draw: draw, graph: graph, layout: layout, subgraph: subgraph),
)
#tidy.show-module(graph-docs, style: tidy-style)

#let draw-docs = tidy.parse-module(
  read("../src/draw.typ"),
  name: "draw",
  scope: (draw: draw, graph: graph, layout: layout, subgraph: subgraph),
)
#tidy.show-module(draw-docs, style: tidy-style)

#let subgraph-docs = tidy.parse-module(
  read("../src/subgraph.typ"),
  name: "subgraph",
  scope: (draw: draw, graph: graph, layout: layout, subgraph: subgraph),
)
#tidy.show-module(subgraph-docs, style: tidy-style)
