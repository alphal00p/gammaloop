# linnest

Typst wrapper package for the `linnest.wasm` graph layout plugin.

## Usage

```typ
#import "src/lib.typ": draw, graph, layout, subgraph

#let a = 0
#let c = 1
#let g = graph.build(
  name: "demo",
  nodes: (graph.node(name: "a"), graph.node(name: "c")),
  edges: (graph.edge(source: (node: a), sink: (node: c), statements: (color: "0055ff", label: "a-c")),),
)
#let g = layout(g)
#let dot = graph.dot(g)
#let edge-label(edge) = text(fill: rgb("#" + edge.color))[#edge.label]
#draw(g, edge-label: edge-label)
```

See `docs/manual.typ` and the rendered `docs/manual.pdf` for the full API reference.
