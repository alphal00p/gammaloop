#import "@preview/tidy:0.4.3"
#import "../src/lib.typ": graph, layout, subgraph, config

#set document(title: "Linnest Typst API")
#set page(margin: 22mm)
#set text(size: 10pt)

= Linnest Typst API

Linnest exposes the `linnest.wasm` graph layout plugin through a small Typst API.
The public surface is intentionally narrow:

- `graph` for construction, parsing, inspection, joins, and graph algorithms.
- `subgraph` for archived subgraph construction and inspection.
- `layout` for the separate layout pass.
- `config` for default layout settings.

== Minimal Builder Example

```typ
#import "../src/lib.typ": graph, layout, subgraph

#let b = graph.builder(name: "demo")
#let (node: a, builder: b) = graph.node(b, name: "a")
#let (node: c, builder: b) = graph.node(b, name: "c")
#let b = graph.edge(b, source: (node: a), sink: (node: c))
#let g = graph.finish(b)
#let g = layout(g, seed: "2", steps: "5")
#let north = subgraph.compass(g, "n")
#let edges = graph.edges(g, subgraph: north)
#let dot = graph.dot(g)
```

== Graph Values

Linnest graph values are archived byte arrays. Treat them as opaque values: build
or parse them with `graph`, transform them with `layout`, and pass them back to
`graph` or `subgraph` for inspection.

- `graph.parse(input)` parses one or more DOT digraphs and returns an array of
  archived graphs.
- `graph.build(spec)` constructs one graph from a Typst dictionary. This is
  convenience sugar over the builder functions.
- `graph.builder(..)` starts an archived builder.
- `graph.node(builder, ..)` returns `(node: index, builder: builder)`.
- `graph.edge(builder, ..)` returns the updated builder.
- `graph.finish(builder)` turns a builder into a graph.
- `layout(graph, seed: "2", steps: "5", ..)` runs layout as an explicit
  second step. Its settings are named parameters so calls stay descriptive and
  Tidy can document each field.
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

`graph.build` accepts the same data in one dictionary:

```typ
#let g = graph.build(
  name: "demo",
  statements: (full_num: "x + y"),
  node-statements: (shape: "circle"),
  edge-statements: (color: "black"),
  nodes: ((name: "a"), (name: "b")),
  edges: ((
    source: (node: 0, compass: "e"),
    sink: (node: 1, compass: "w"),
    statements: (label: "ab"),
  ),),
)
```

== Graph Queries

`graph.info(g)` returns graph metadata. `graph.nodes(g)` returns node records,
and `graph.edges(g)` returns edge records. Pass `subgraph: sg` to filter nodes
or edges by an archived subgraph.

`graph.join(left, right, key: "statement")` joins matching dangling half edges.
The key is read from half-edge data and can be `"statement"`, `"port_label"`,
`"compass"`, or `"id"`.

`graph.cycles(g)` returns archived subgraphs for a cycle basis.
`graph.forests(g)` returns archived subgraphs for spanning forests.

== Subgraphs

Subgraph values are archived byte arrays.

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
  scope: (graph: graph, layout: layout, subgraph: subgraph, config: config),
)
#tidy.show-module(docs, style: tidy-style)

#let graph-docs = tidy.parse-module(
  read("../src/graph.typ"),
  name: "graph",
  scope: (graph: graph, subgraph: subgraph),
)
#tidy.show-module(graph-docs, style: tidy-style)

#let subgraph-docs = tidy.parse-module(
  read("../src/subgraph.typ"),
  name: "subgraph",
  scope: (graph: graph, subgraph: subgraph),
)
#tidy.show-module(subgraph-docs, style: tidy-style)
