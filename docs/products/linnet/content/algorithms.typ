#import "../../shared.typ": callout, boundary

#let algorithms = [
= Subgraphs, traversal, cuts, and drawing

Linnet algorithms operate on a graph together with an explicit subgraph view. This lets one
algorithm work on the full graph, a selected region, or the boundary created by cutting paired
half-edges. Choose the subgraph representation by operation: bit-backed sets are efficient for
repeated set algebra, while borrowed or mapped views avoid copying when identity must remain
attached to the owning graph.

== Traversal and connectivity

A traversal needs a starting half-edge or node, the active subgraph, and a policy describing
which adjacent half-edges may be crossed. Results retain typed indexes into the original graph.
Breadth-first traversal is useful for distance layers and connectivity; depth-first traversal
is useful for trees, back edges, and decompositions. A traversal tree is a derived view, not a
new graph, so mutations must be coordinated with the indexes it contains.

Connected components, bridges, articulation structure, and biconnected regions all depend on
the chosen subgraph. Run them on the full graph only when external half-edges and excluded nodes
should participate in the answer.

== Cycles and cuts

`cycle_basis` returns an independent basis and rank information. `all_cycles` expands that basis
and can grow exponentially; prefer a basis when assigning loop momenta or computing a loop
count. Symmetric differences combine cycle bitsets without inventing new graph identities.

Cut enumeration separates required left and right node sets and returns the left region, cut
content, and right region. A cut edge is represented through its half-edges, preserving which
side owns each endpoint. Validate direction/flow after enumeration when a physical Cutkosky or
tensor-network interpretation depends on orientation.

#callout("Complexity is part of the API choice", [
  Enumerating every cycle or every admissible cut is inherently combinatorial. Filter the
  subgraph and endpoint sets first, and use basis or connectivity operations when enumeration
  is not required by the caller.
])

== Enumerate separating cuts

This complete example builds a square with a diagonal and asks for every cut separating opposite
nodes. It is the smallest graph that makes the distinction between a path and a cut family
visible without external data or drawing support.

// docs-example: run
```rust
use linnet::half_edge::{HedgeGraph, builder::HedgeGraphBuilder};

fn main() {
    let mut builder = HedgeGraphBuilder::new();
    let a = builder.add_node(());
    let b = builder.add_node(());
    let c = builder.add_node(());
    let d = builder.add_node(());

    builder.add_edge(a, b, (), true);
    builder.add_edge(b, c, (), true);
    builder.add_edge(c, d, (), true);
    builder.add_edge(d, a, (), true);
    builder.add_edge(b, d, (), true);

    let graph: HedgeGraph<(), (), ()> = builder.build();
    graph.check().expect("the half-edge involution is valid");
    let cuts = graph.all_cuts_from_ids(&[a], &[c]);

    assert_eq!(cuts.len(), 4);
    println!("separating cuts: {}", cuts.len());
}
```

The expected invariant is four distinct left/cut/right partitions separating `a` from `c`.
Changing the diagonal changes that count, while changing only node or edge payloads does not.
The generated #link("reference/rust/supported/linnet/#supported-hedgegraph")[`HedgeGraph`
reference] and
#link("reference/rust/supported/linnet/#supported-hedgegraphbuilder")[builder reference] give the
exact supported type boundaries.

#callout("Interpret a cut mismatch structurally", [
  A failed `graph.check()` indicates a storage or involution invariant and must be resolved before
  enumeration. Zero cuts usually means the endpoint sets overlap or are disconnected from the
  selected region. An unexpected nonzero count means the constructed edges or endpoint sets do
  not describe the graph you intended; inspect those before adding post-enumeration filters.
])

== Mutation boundaries

Extraction copies a selected region, excision separates a graph along its boundary, sewing
joins compatible dangling half-edges, and contraction identifies structure. These operations
return mappings or replacement indexes where needed. Keep those results: a numeric `Hedge` or
`NodeIndex` from before a storage-changing operation is not automatically a valid index into the
resulting graph.

== DOT and drawing

DOT parsing and emission are interchange/debugging boundaries. Attributes can carry graph,
node, edge, and half-edge data, but the typed Rust invariants are reconstructed by the parser;
handle parser errors before running algorithms. The `drawing` feature adds layout and rendering.
Layout coordinates are presentation data and may vary without changing topology, subgraph
membership, or edge flow.

#boundary("Tensor contraction uses Spenso", [
  Linnet determines connectivity and graph transformations. Tensor compatibility, contraction
  cost, and execution belong to Spenso even when a Spenso network delegates its graph
  algorithms to Linnet.
])
]
