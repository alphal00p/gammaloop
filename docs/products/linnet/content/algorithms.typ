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

Connectivity and traversal results depend on the chosen subgraph. A traversal tree retains its
parent relation, so callers can inspect immediate children, strict ancestors, and the fundamental
cycle closed by a non-tree internal edge without reconstructing an adjacency model. Boundary
queries distinguish half-edges already selected at a split from every boundary half-edge incident
to a touched node; dangling half-edges are available as a separate structural selection.

Linnet has two independent notions of direction. `DirectionBasis::Underlying` follows the
source/sink roles stored by the involution. This is the default in Python and is the basis used by
GammaLoop momentum routing. `DirectionBasis::Superficial` additionally applies an edge's
`Default`, `Reversed`, or `Undirected` orientation; a superficially undirected edge contributes no
directed arc. Reachability, topological order, transitive closure, and transitive reduction take
this choice explicitly. Neither mode rewrites stored flow or orientation.

Directed algorithms ignore dangling half-edges and edges split by the active subgraph boundary,
because neither represents a complete arc inside that selection. Topological ties use node-index
order for the initial ready set and deterministic edge-discovery order thereafter. Linnet does not
currently expose a dedicated articulation-point or biconnected-component API; derive those from
traversal data only when your application owns that algorithm and its tests.

== Cycles and cuts

`cycle_basis` returns an independent cycle vector together with the spanning forest used to
construct it. Use `cyclotomatic_number` when the numerical cycle rank itself is the result.
`all_cycles` expands a basis and can grow exponentially; prefer a basis when assigning loop
momenta. Native callers supply a maximum accepted cycle rank. Python callers supply
`max_results`, which bounds the complete non-empty basis-combination search space rather than
only the circuits that survive filtering. Both guards are checked before powerset construction.
Symmetric differences combine cycle bitsets without inventing new graph identities.

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

// docs-example: run linnet-separating-cuts
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
Changing topology can change that count, while changing only node or edge payloads does not.
The native #link("reference/rust/linnet/half_edge/struct.HedgeGraph.html")[`HedgeGraph`
Rustdoc] and
#link("reference/rust/linnet/half_edge/builder/struct.HedgeGraphBuilder.html")[builder Rustdoc]
give the exact compiled type boundaries for this revision.

The Python binding delegates to that same native enumeration and returns named graph-bound
results rather than nested tuples:

```python
partitions = graph.all_cuts(["a"], ["c"])
for partition in partitions:
    source_side = partition.source_side
    boundary = partition.boundary
    cut_edges = partition.edges
    target_side = partition.target_side
```

The source and target sides are structural `Subgraph` views. `boundary.left` and
`boundary.right` select the two oriented half-edges of every cut edge. The result is eagerly
materialized and sorted by half-edge membership at the Python boundary; do not use all-cut
enumeration where one witness or a minimum-cardinality cut is the actual question.

`all_bonds(min_size=..., max_size=...)` is the related native enumeration of
inclusion-minimal cutsets whose two induced sides are connected. A bond size is its boundary
half-edge count; it is not a max-flow result. `find_bond` short-circuits after one witness in the
same inclusive size range. Both full enumerations can be combinatorial.

`tadpoles(externals)` identifies a validated terminal set structurally and returns the
native tadpole components in stable order. Terminals must resolve to distinct structural nodes;
empty or duplicate sets, out-of-range nodes, and zero-crown terminals are errors instead of panic
paths.

#callout("Interpret a cut mismatch structurally", [
  A failed `graph.check()` indicates a storage or involution invariant and must be resolved before
  enumeration. Overlapping, empty, isolated, or closed terminal groups are rejected. For validated
  groups, zero cuts means the native enumerator found no admissible separating partition, commonly
  because the groups occupy disconnected components. An unexpected nonzero count means the
  constructed edges or endpoint sets do not describe the graph you intended; inspect those before
  adding post-enumeration filters.
])

== Mutation boundaries

Insertion extends a graph with the same declarative node and edge specs accepted by the builder.
`split_edge` separates one paired edge into two dangling edges, while `connect` joins two exact
dangling half-edges. Extraction copies a selected region, excision separates a graph along its
boundary, sewing joins compatible dangling half-edges, and contraction identifies structure.
These operations return replacement views or indexes where needed. Keep only the topology-independent
payloads across an edit: graph-bound Python views are revision checked and become stale after every
successful topology change.

Transitive reduction is an owning transformation over a directed acyclic graph. The Python method
returns a new graph, preserves arbitrary payload identity and copied drawing configuration, and
leaves views into the source graph live. Transitive closure remains a Rust-only operation until the
Python API has an explicit factory for every new edge and half-edge value it synthesizes.

== DOT and drawing

DOT parsing and emission are interchange/debugging boundaries. Attributes can carry graph,
node, edge, and half-edge data, but the typed Rust invariants are reconstructed by the parser;
handle parser errors before running algorithms. The `drawing` feature enables layout; rendering
belongs to a DOT consumer, Clinnet, or Linnest. Layout coordinates are presentation data and may
vary without changing topology, subgraph membership, or edge flow.

#boundary("Tensor contraction uses Spenso", [
  Linnet determines connectivity and graph transformations. Tensor compatibility, contraction
  cost, and execution belong to Spenso even when a Spenso network delegates its graph
  algorithms to Linnet.
])
]
