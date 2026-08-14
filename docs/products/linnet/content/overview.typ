#import "../../shared.typ": callout, boundary, product-link, source-link

#let overview = [
= Overview

Linnet is a graph library for tensor networks, Feynman diagrams, and other workflows where
the boundary of a subgraph matters. Its central representation is a half-edge graph. An
ordinary internal edge is a pair of half-edges, while an external edge has one dangling
half-edge. A node owns incident half-edges, and a subgraph is a selection of half-edges rather
than a detached copy of the graph.

#callout("Why half-edges", [
  Cutting an internal edge can turn its two halves into boundary edges without changing the
  degree of either endpoint. This makes operations such as extraction, excision, contraction,
  sewing, traversal, and cut analysis natural for physics graphs. Linnet algorithms therefore
  accept subgraph views even when the caller wants to operate on the full graph.
])

== A working model

A Linnet workflow has four parts:

- construct a graph with `HedgeGraphBuilder`, or parse a DOT graph;
- retain the typed identifiers returned for half-edges, nodes, and edges;
- describe the region of interest with one of the subgraph representations;
- run traversal or mutation operations against that subgraph and the graph's involution.

The involution is the pairing map for half-edges. A paired entry represents an internal edge;
an identity entry represents a dangling edge and records its flow. Edge orientation and flow
are related concepts, but they are not interchangeable. Use the generated API descriptions for
the exact variants and conversion rules in the selected release.

#boundary("Indexes belong to one graph", [
  `Hedge`, `NodeIndex`, and `EdgeIndex` are compact typed indexes, not durable object IDs.
  Graph edits can change storage and mappings. Keep the result returned by an operation, and do
  not apply an index or subgraph bitset to an unrelated graph just because its numeric values
  happen to fit.
])

== Construction and interchange

The Rust package can be added directly from crates.io:

```toml
[dependencies]
linnet = "0.17"
```

Programmatic construction uses `HedgeGraphBuilder`. DOT support provides a convenient
human-readable interchange and debugging format; it can carry global, node, edge, and
half-edge attributes. Treat DOT as a representation boundary, not as a substitute for the
typed graph invariants. Validate or handle parser errors before running algorithms.

The `drawing` feature adds graph drawing and layout support. Layout is separate from graph
identity: coordinates and rendering attributes may change without changing the topology.

== Ecosystem ownership

#boundary("One graph layer, several consumers", [
  #product-link("spenso", label: "Spenso") uses Linnet to represent tensor-network connectivity
  and owns tensor contraction semantics. #product-link("gammaloop", label: "GammaLoop") uses
  graph concepts in collider workflows and owns its process/state lifecycle. Those manuals link
  here for graph operations instead of redefining half-edge behavior.
])

The crate's maintained entry points are the
#source-link("crates/linnet/README.md", label: "Linnet README") and generated Rust catalog.
]
