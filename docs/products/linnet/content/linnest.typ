#import "../../shared.typ": boundary
#import "/crates/linnest/typst/docs/manual.typ": manual as linnest-manual

#let linnest = [
= Linnest Typst API

Linnest is the Typst and WebAssembly interface to Linnet's graph model, DOT parser,
layout algorithms, graph queries, and physics-aware drawing styles. This page embeds
the package's complete maintained manual so its examples and generated API reference
remain part of the same searchable Linnet documentation corpus.

#boundary("Choose the interface deliberately", [
  Use Linnet's Rust or Python interfaces for graph construction and algorithms in an
  application. Use Linnest when the graph, layout, or final drawing is owned by a Typst
  document.
])

#linnest-manual
]
