#import "@preview/tidy:0.4.3"
#import "../src/lib.typ": draw, graph, layout, physics, subgraph

#set document(title: "Linnest Typst API")
#set page(margin: 22mm)
#set text(size: 10pt)

= Linnest Typst APIs

Linnest exposes the `linnest.wasm` graph layout plugin through a small Typst
API. Its drawing layer bundles `kurvst.wasm`; the standalone Kurvst package
manual documents the full curve and path-pattern API.
The public surface is intentionally narrow:

- `graph` for construction, parsing, inspection, joins, and graph algorithms.
- `subgraph` for subgraph object construction and inspection.
- `layout` for the separate layout pass.
- `draw` for rendering a laid-out graph object with CeTZ.
- `physics` for reusable particle-line edge styles.

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

`graph.edge` accepts `source` and `sink` half-edge dictionaries. A half edge can
contain `node`, `statement`, `id`, `port-label`, `compass`, and `in-subgraph`.
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
    display-label: "{label}",
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

== Layout Pins

`graph.pin` formats the `pin` statement consumed by layout. Numeric `x` and
`y` values fix coordinates for nodes or edge control points:

```typ
#let (node: a, builder: b) = graph.node(b, name: "a", pin: graph.pin(x: -2, y: 0))
#let (node: c, builder: b) = graph.node(b, name: "c", pin: graph.pin(x: 2, y: 0))
#let b = graph.edge(b, source: (node: a), sink: (node: c), pin: graph.pin(x: 0, y: 1.2))
```

`graph.pin-group` links one coordinate across several nodes or edge control
points. A `side` of `"+"` keeps the coordinate positive, and `"-"` keeps it
negative. GammaLoop external-edge columns use this to keep incoming and outgoing
external legs on opposite sides while pairing rows by a shared `y` group:

```typ
#let b = graph.edge(
  b,
  source: (node: right-ext),
  sink: (node: center),
  pin: graph.pin(x: graph.pin-group("right", side: "+"), y: graph.pin-group("edgee0")),
)
#let b = graph.edge(
  b,
  source: (node: center),
  sink: (node: left-ext),
  pin: graph.pin(x: graph.pin-group("left", side: "-"), y: graph.pin-group("edgee0")),
)
```

Draw styling is Typst-native. Pass dictionaries or callbacks to `draw`; edge
callbacks receive the merged `scope`, edge statements, source/sink half-edge
records, and the edge index:

```typ
#let edge-label(edge) = text(fill: rgb("#" + edge.color))[#edge.display-label]
#let source-style(edge) = (stroke: red + 0.5pt)
#let sink-style(edge) = (stroke: blue + 0.5pt)
#draw(layout(graph.finish(b)), edge-label: edge-label, source-style: source-style, sink-style: sink-style)
```

== Physics Edge Styles

`physics.style(..)` returns `source-style`, `sink-style`, and `edge-label`
callbacks for `draw`. The default source/sink strokes use darker/lighter halves
to encode the graph source/sink split. Set `momentum-arrows: true` to draw
centered parallel arrow decorations while the main edge remains connected to
the nodes. The arrows flow from source to sink independently of the DOT
`dir`/`orientation` value; there is one momentum marker per edge, even though
the base edge is internally styled as source and sink halves. On paired edges
the marker lives on the sink half; on dangling edges it lives on the existing
source or sink half-edge. The decoration defaults to a plain black `0.55pt`
stroke and a scaled CeTZ `"barbed"` arrowhead, so it does not inherit the base
particle stroke. The arrow length is capped by both `momentum-arrow-length` and
`momentum-arrow-ratio`, so the shorter one wins. When the layout provides an
edge-label position, the arrow offset is signed so the momentum marker is drawn
on the same side of the edge as that label. Use `momentum-arrow-offset` to set
the normal displacement. `momentum-arrow-mark: auto` uses the default barbed
marker, `momentum-arrow-mark: none` draws only the momentum line, and any other
value overrides the CeTZ mark style, for example `"stealth"` or `(end: "barbed",
scale: 1.6)`. Only the `end` mark is kept so source/sink splitting cannot
create a second momentum head.

Optional labels can be built from edge metadata with `show-edge-index`,
`show-half-edge-index`, `show-particle`, and, for explicit momentum fields,
`show-momentum`.
`show-half-edge-index` only emits a value for dangling half edges.

```typ
#import "../src/lib.typ": draw, graph, layout, physics

#let b = graph.builder(name: "physics")
#let (node: a, builder: b) = graph.node(b, name: "a")
#let (node: c, builder: b) = graph.node(b, name: "c")
#let b = graph.edge(
  b,
  source: (node: a),
  sink: (node: c),
  statements: (particle: "g", id: "7"),
)
#let callbacks = physics.style(
  momentum-arrows: true,
  show-edge-index: true,
  show-particle: true,
)
#draw(
  layout(graph.finish(b)),
  source-style: callbacks.source-style,
  sink-style: callbacks.sink-style,
  edge-label: callbacks.edge-label,
)
```

Set `edge-offset` on `draw`, or `offset` in a `source-style`/`sink-style`
dictionary, to draw a fitted parallel path. Style callbacks may also return an
array of dictionaries; the layers are drawn in order on the same graph edge, so
parallel strokes do not require duplicate graph edges.

```typ
#let base-style(edge) = (
  stroke: (paint: gray, thickness: 0.7pt, cap: "round"),
)
#let offset-style(edge) = (
  offset: 0.18,
  length: 1.4,
  ratio: 0.5,
  resolve-length: "min",
  stroke: (paint: rgb("#2f6f4e"), thickness: 1pt, cap: "round"),
)
#let source-style(edge) = (base-style(edge), offset-style(edge))
#let sink-style(edge) = (base-style(edge), offset-style(edge) + (mark: (end: ">")))
```

The parallel path is applied to the base edge geometry before patterns and other
decorations; node outsets then trim the shifted path, so shifted paths still
start and end outside fitted node circles. Add `edge-length` or `length` to
center-trim the shifted path to a fixed arc length, and add `edge-ratio` or
`ratio` to cap it by a fraction of the base edge length. `edge-resolve-length`
/ `resolve-length` decides how to combine both limits: `"min"`/`"shorter"` (default), `"max"`/`"longer"`,
`"length"`/`"fixed"`, `"ratio"`/`"relative"`, `"none"`/`"full"`, or a function
receiving `(base-length, length, ratio)`.
Set `offset-side: "label"` on an offset layer to choose the sign of `offset`
so the layer is drawn on the same side of the curve as the edge label.

The builder keeps `source-style-eval`, `sink-style-eval`, and `label-eval` as
ordinary kebab-case edge metadata for downstream renderers or DOT export. They
can be set as direct arguments on `graph.build`, `graph.builder`, or
`graph.edge`, instead of being manually nested under `edge-statements` or
per-edge `statements`:

```typ
#let b = graph.builder(
  name: "demo",
  source-style-eval: "(stroke: red + 0.5pt)",
  sink-style-eval: "(stroke: blue + 0.5pt)",
)
#let b = graph.edge(
  b,
  source: (node: a),
  sink: (node: c),
  label-eval: "(text(fill: rgb(\"#{color}\"))[{label}])",
  statements: (color: "0055ff", label: "a-c"),
)
```

`graph.build` accepts the same data in one dictionary:

```typ
#let g = graph.build(
  name: "demo",
  statements: (full_num: "x + y"),
  node-statements: (shape: "circle"),
  edge-statements: (
    color: "000000",
    display-label: "{label}",
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
The key is read from half-edge data and can be `"statement"`, `"port-label"`,
`"compass"`, or `"id"`.

`graph.cycles(g)` returns subgraph objects for a cycle basis.
`graph.forests(g)` returns subgraph objects for spanning forests.

== Layout Model

`layout` starts from a traversal-tree placement, then optimizes the positions of
graph nodes and edge control points. The initial tree spacing is

$ L = lambda sqrt((W H) / max(n, 1)) $,

with horizontal spacing $tau_x L$ and vertical spacing $tau_y L$. Here
$lambda$ is `length-scale`, $W$ is `viewport-w`, $H$ is `viewport-h`, $tau_x$
is `tree-dx`, and $tau_y$ is `tree-dy`. These fields set the geometry scale for
both layout modes.

The shared spring/charge model uses the following coefficients:

$ c_("vv") = beta L^2 $
$ c_("ev") = beta gamma_("ev") L^2 $
$ c_("ee") = beta gamma_("ee") L^2 $
$ c_("center") = beta g_("center") L^2 $
$ c_("dangling") = beta gamma_("dangling") L^2 $

The Typst parameter names are `beta` for $beta$, `gamma-ev` for
$gamma_("ev")$, `gamma-ee` for $gamma_("ee")$, `g-center` for $g_("center")$,
and `gamma-dangling` for $gamma_("dangling")$. The spring stiffness $k$ is
`k-spring`, and the softening constant $epsilon$ is `eps`.

In `layout-algo: "anneal"`, linnest minimizes an energy:

$ E =
sum_(i < j) 1/2 c_("vv") / (d(v_i, v_j) + epsilon)
+ sum_(i, e) c_("ev") / (d(v_i, e) + epsilon)
+ sum_((v, e) " incident") 1/2 k (ell_e - d(v, e))^2
\ + sum_("local edge pairs") 1/2 c_("ee") / (d(e_i, e_j) + epsilon)
+ sum_("dangling pairs") 1/2 c_("dangling") / (d(e_i, e_j) + epsilon)
+ sum_(i) 1/2 c_("center") / (d(v_i, 0) + epsilon)
+ p_("cross") N_("cross") $.

Here $p_("cross")$ is `crossing-penalty` and $N_("cross")$ is the number of
detected edge crossings. `temp`, `step`, `seed`, `steps`, `epochs`, `cool`,
`accept-floor`, `step-shrink`, and `incremental-energy` belong to this
simulated annealing mode. `crossing-penalty` is also anneal-only; force mode
does not currently add a crossing force.

In `layout-algo: "force"`, linnest applies the direct forces corresponding to
the same vertex-vertex, edge-vertex, incidence spring, local edge-edge,
dangling-edge, and center terms. `step` is the integration step, `delta` clamps
per-step movement, `steps` and `epochs` set the iteration budget, `cool` shrinks
the step after each epoch, and `early-tol` stops when movement is small.
`z-spring` and `z-spring-growth` are force-only helpers: the integrator gives
points temporary z coordinates to break overlaps and pulls them back toward the
2D plane.

`directional-force` is applied in both modes as an extra bias derived from
pin/port direction constraints.

After either graph layout mode, labels are relaxed separately. If $L_l$ is the
label target distance and $q_l$ is the label repulsion strength, then
$L_l = alpha_l L$ and $q_l = beta_l L^2$. The Typst names are
`label-length-scale` for $alpha_l$ and `label-charge` for $beta_l$.
`label-spring` is the spring constant pulling each label toward its target.
`label-steps`, `label-step`, `label-early-tol`, and `label-max-delta-scale`
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

#let draw-docs = tidy.parse-module(
  read("../src/draw.typ"),
  name: "draw",
  scope: (draw: draw, graph: graph, layout: layout, physics: physics, subgraph: subgraph),
)
#let docs = tidy.parse-module(
  read("../src/lib.typ"),
  scope: (draw: draw, graph: graph, layout: layout, physics: physics, subgraph: subgraph),
)
#let graph-docs = tidy.parse-module(
  read("../src/graph.typ"),
  name: "graph",
  scope: (draw: draw, graph: graph, layout: layout, physics: physics, subgraph: subgraph),
)

#let physics-docs = tidy.parse-module(
  read("../src/physics-edge-style.typ"),
  name: "physics",
  scope: (draw: draw, graph: graph, layout: layout, physics: physics, subgraph: subgraph),
)
#let subgraph-docs = tidy.parse-module(
  read("../src/subgraph.typ"),
  name: "subgraph",
  scope: (draw: draw, graph: graph, layout: layout, physics: physics, subgraph: subgraph),
)



#tidy.show-module(draw-docs, style: tidy-style)
#tidy.show-module(docs, style: tidy-style)
#tidy.show-module(graph-docs, style: tidy-style)
#tidy.show-module(physics-docs, style: tidy-style)
#tidy.show-module(subgraph-docs, style: tidy-style)
