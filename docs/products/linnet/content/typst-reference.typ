#import "../../shared.typ": boundary, product-link

#let typst-reference = [
= Linnest Typst API

Linnest documents a small supported surface around Linnet's graph model. Start with the
#link("guides/linnest/")[workflow guide] if you want to build and draw a graph; use these
pages when you need an exact function, parameter, or default.

== Supported surface

- #link("reference/typst/graph/")[`graph`] is a module for construction, parsing,
  inspection, data transforms, joins, cycles, and forests.
- #link("reference/typst/layout/")[`layout` and `layout-sequence`] are top-level
  functions. They are documented together because both come from `layout.typ`.
- #link("reference/typst/drawing/")[`draw`, `edge-halves`, and
  `to-cetz-edge-halves`] are top-level drawing functions from `draw.typ`.
- #link("reference/typst/subgraph/")[`subgraph`] is a module for constructing and
  inspecting zero-copy subgraph values.
- #link("reference/typst/templates/")[Domain styles and templates] explains how
  application-owned callbacks compose the generic drawing surface.
- `curve` is a re-export of Kurvst, not a second Linnest API. Use the canonical
  #product-link("gammaloop", label: "Kurvst curve and path reference", path: "guides/kurvst/").

#boundary("Modules and functions import differently", [
  Import `graph` and `subgraph` as modules when you want qualified calls such as
  `graph.build`. Import `layout`, `layout-sequence`, and `draw` directly from `lib.typ`;
  they are functions rather than nested modules.
])

```typ
#import "crates/linnest/typst/src/lib.typ": draw, graph, layout, subgraph
#import graph: build, edge, node, sink, source
```

The same package also exports the `layouts` module for callers that need the underlying
layout namespace; its `layouts.layout` and `layouts.sequence` aliases are linked from the
layout page. Bindings exposed only as a consequence of Typst module loading, such as
serialization helpers or imported dependency modules, remain implementation details rather
than supported API. The focused reference uses the supported functions and modules above.
]
