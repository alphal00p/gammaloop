#import "/crates/linnest/typst/docs/manual.typ": subgraph-concepts, subgraph-reference

#let typst-subgraph = [
= subgraph module

Use `subgraph` to construct and inspect opaque zero-copy selections of half-edges. Pass
those selections to graph queries, drawing, and layout functions that accept a `subgraph`
argument.

== Concepts

#subgraph-concepts

#subgraph-reference
]
