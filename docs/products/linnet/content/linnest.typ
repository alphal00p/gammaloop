#import "../../shared.typ": boundary
#import "/crates/linnest/typst/docs/manual.typ": linnest-guide

#let linnest = [
= Build and draw with Linnest

Linnest is the Typst and WebAssembly interface to Linnet's graph model, DOT parser,
layout algorithms, graph queries, and generic drawing styles. This guide follows one graph
from construction through layout and drawing; the separate
#link("reference/typst/")[Typst API reference] carries the supported symbol inventory.

#boundary("Choose the interface deliberately", [
  Use Linnet's Rust or Python interfaces for graph construction and algorithms in an
  application. Use Linnest when the graph, layout, or final drawing is owned by a Typst
  document.
])

#linnest-guide

== Continue into the API

Use the #link("reference/typst/graph/")[graph module] for construction and queries,
the #link("reference/typst/layout/")[layout functions] for placement, and the
#link("reference/typst/drawing/")[drawing functions] for CeTZ output. Domain-template
composition and subgraph objects have their own focused reference pages.
]
