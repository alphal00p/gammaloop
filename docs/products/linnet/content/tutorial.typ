#import "../../shared.typ": callout

#let tutorial = [
= Tutorial

In this tutorial, you will construct a half-edge graph in Rust, validate its invariants, find
its one independent cycle, and print DOT. The builder returns Linnet's typed node and half-edge
identifiers directly, so you can see where graph identity enters the workflow.

== Prerequisites

Create a Rust binary project and add the current Linnet release:

```sh
cargo new linnet-first-graph
cd linnet-first-graph
cargo add linnet@0.17
```

No drawing feature or external layout program is needed to construct the graph and emit DOT.

== Build a graph with a boundary

Replace `src/main.rs` with the following program:

```rust
use linnet::half_edge::{
    builder::HedgeGraphBuilder,
    involution::Flow,
    HedgeGraph,
};

fn main() {
    let mut builder = HedgeGraphBuilder::new();
    let a = builder.add_node("a");
    let b = builder.add_node("b");
    let c = builder.add_node("c");

    builder.add_edge(a, b, "ab", false);
    builder.add_edge(b, c, "bc", false);
    builder.add_edge(c, a, "ca", false);
    builder.add_external_edge(a, "incoming", false, Flow::Sink);
    builder.add_external_edge(c, "outgoing", false, Flow::Source);

    let graph: HedgeGraph<&str, &str> = builder.build();
    graph.check().expect("the half-edge involution is valid");

    let (basis, _) = graph.cycle_basis();
    assert_eq!(basis.len(), 1);
    println!("{}", graph.base_dot());
}
```

Run it with `cargo run`. Success means the assertion passes and stdout contains a DOT
`digraph` with nodes `a`, `b`, and `c`, three paired internal edges, and two dangling boundary
edges. The external edges do not add a cycle, so the triangle's cycle-basis rank remains one.

#callout("Keep the returned indexes with their graph", [
  `a`, `b`, and `c` are compact `NodeIndex` values owned by this graph construction. Do not
  serialize their numeric representation as a durable identity or reuse them against another
  graph. Mutating operations return the mappings needed to relate old and new storage.
])

== Make the first useful variation

Add a fourth internal edge between `a` and `b`, run the program again, and change the assertion
to `basis.len() == 2`. This demonstrates the key distinction: parallel paired half-edges add an
independent cycle, while dangling identity half-edges describe the graph boundary.

From here, use a subgraph view to run traversal, cuts, or contraction against only the region
you intend to modify. Enable Linnet's `drawing` feature only when you need layout rather than
plain DOT interchange.

== Troubleshooting and next steps

- A type-inference error around `build()` means Rust cannot choose node storage or data types;
  keep the explicit `HedgeGraph<&str, &str>` annotation from the example.
- A failed `check()` points to an involution or storage invariant. Validate immediately after
  parsing or construction, before applying graph algorithms.
- If a traversal includes an unexpected dangling edge, check the selected subgraph and the
  `Flow` assigned to each external edge.
- Continue with the subgraph-and-algorithms manual, then consult the Rust API reference for the
  exact return mappings of mutating operations and additional builder methods.
]
