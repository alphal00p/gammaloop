#import "../../shared.typ": callout

#let quickstart-typst = [
= Using Linnet from Typst

Linnest brings Linnet's graph model, layout algorithms, and physics-aware drawing styles into
Typst through a bundled WebAssembly plugin. This first document draws two nodes joined by one
internal edge.

#callout("Bundled source package", [
  Linnest `0.1.0` is not currently published on Typst Universe. Keep the repository's complete
  `crates/linnest/typst/` directory together: its Typst sources load the adjacent `linnest.wasm`
  plugin. Typst `0.15.0` or newer is required.
])

== Create the document

From a GammaLoop checkout, save this as `linnest-quickstart.typ` in the repository root:

// docs-example: syntax linnet-typst-quickstart
```typ
#import "crates/linnest/typst/src/lib.typ": draw, graph, layout

#let g = graph.build({
  graph.node(<left>, label: [left])
  graph.node(<right>, label: [right])
  graph.edge(
    graph.source(<left>),
    <propagator>,
    graph.sink(<right>),
    label: [$p$],
  )
})

#draw(layout(g))
```

Compile it from the same repository root:

// docs-example: syntax
```sh
typst compile --root . linnest-quickstart.typ linnest-quickstart.pdf
```

The PDF should contain two labeled nodes and one labeled edge. Construction, layout, and drawing
are separate stages: keep the graph value when you need queries or subgraphs before rendering it.

Continue with the #link("guides/linnest/")[Linnest guide] for DOT parsing, explicit styling, and
subgraph drawing, then use the #link("reference/typst/")[Typst API reference] for graph, layout,
drawing, physics, and subgraph functions.
]
