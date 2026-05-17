# linnest

Typst wrapper package for the `linnest.wasm` graph layout plugin.

## Usage

```typ
#import "src/lib.typ": draw, edge, graph, layout, node, sink, source, subgraph

#let g = graph.build({
  node(<a>, label: [A], statements: (fill-color: "cfe8ff"))
  node(<c>, label: [C], statements: (fill-color: "d6f5d6"))
  edge(
    source(<a>, compass: "e"),
    <a-c>,
    sink(<c>, compass: "w"),
    label: [a-c],
    statements: (color: "0055ff"),
  )
},
  name: "demo",
)
#let g = layout(g)
#let dot = graph.dot(g)
#let edge-label(edge) = text(fill: rgb("#" + edge.color))[#edge.label]
#draw(g, edge-label: edge-label)
```

See `docs/manual.typ` and the rendered `docs/manual.pdf` for the full API reference.
