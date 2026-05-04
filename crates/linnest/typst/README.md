# linnest

Typst wrapper package for the `linnest.wasm` graph layout plugin.

## Usage

```typ
#import "src/lib.typ": draw, graph, layout, subgraph

#let b = graph.builder(name: "demo")
#let (node: a, builder: b) = graph.node(b, name: "a")
#let (node: c, builder: b) = graph.node(b, name: "c")
#let b = graph.edge(b, source: (node: a), sink: (node: c), statements: (color: "0055ff", label: "a-c"))
#let g = layout(graph.finish(b))
#let dot = graph.dot(g)
#let edge-label(edge) = text(fill: rgb("#" + edge.color))[#edge.label]
#draw(g, edge-label: edge-label)
```

See `docs/manual.typ` and the rendered `docs/manual.pdf` for the full API reference.
