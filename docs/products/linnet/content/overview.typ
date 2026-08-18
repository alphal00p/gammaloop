#import "../../shared.typ": callout, boundary, developer-link, product-link

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

== Choose a task

- To construct and validate a graph with a boundary, follow the
  #link("tutorial/")[first-graph tutorial] and the exact
  #link("reference/rust/supported/linnet/#supported-hedgegraphbuilder")[`HedgeGraphBuilder`
  reference].
- To enumerate cycles, cuts, or connected regions, use the
  #link("guides/algorithms/")[graph-algorithms guide] with the
  #link("reference/rust/supported/linnet/#supported-hedgegraph")[`HedgeGraph` reference].
- To render a tree of existing DOT files and rebuild only changed figures, use the
  #link("guides/clinnet/")[Clinnet command-line guide].
- To lay out or render a graph, continue to the #link("guides/linnest/")[Linnest Typst guide];
  layout coordinates are deliberately separate from graph identity.

== A working model

A Linnet workflow has four parts:

- construct a graph with `HedgeGraphBuilder`, or parse a DOT graph;
- retain the typed identifiers returned for half-edges, nodes, and edges;
- describe the region of interest with one of the subgraph representations;
- run traversal or mutation operations against that subgraph and the graph's involution.

The involution is the pairing map for half-edges. A paired entry represents an internal edge;
an identity entry represents a dangling edge and records its flow. Edge orientation and flow
are related concepts, but they are not interchangeable. Consult the API reference for the exact
variants and conversion rules in your Linnet release.

#boundary("Indexes belong to one graph", [
  `Hedge`, `NodeIndex`, and `EdgeIndex` are compact typed indexes, not durable object IDs.
  Graph edits can change storage and mappings. Keep the result returned by an operation, and do
  not apply an index or subgraph bitset to an unrelated graph just because its numeric values
  happen to fit.
])

== Construction and interchange

The Rust package can be added directly from crates.io:

// docs-example: syntax
```toml
[dependencies]
linnet = "0.17"
```

Programmatic construction uses `HedgeGraphBuilder`. DOT support provides a convenient
human-readable interchange and debugging format; it can carry global, node, edge, and
half-edge attributes. Treat DOT as a representation boundary, not as a substitute for the
typed graph invariants. Validate or handle parser errors before running algorithms.

The `drawing` feature adds layout support. Layout is separate from graph identity: coordinates
may change without changing the topology. Linnet itself does not provide a supported renderer;
use DOT, #link("guides/clinnet/")[Clinnet], or #link("guides/linnest/")[Linnest] at the
presentation boundary.

== Related crates

#boundary("One graph layer, several consumers", [
  #product-link("spenso", label: "Spenso") uses Linnet for tensor-network connectivity and adds
  tensor contraction semantics. #product-link("gammaloop", label: "GammaLoop") uses the graph
  model in collider workflows and adds process generation, integration, and persistent state.
])

The rendered guides and generated APIs are the canonical usage reference. Contributors can
follow the stores, invariants, mutation flow, and feature boundaries in the
#developer-link("linnet-architecture", "linnet-architecture.typ", "Linnet implementation architecture").
]
