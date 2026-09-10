#import "../../shared.typ": callout

#let tutorial = [
= Tutorial

In this tutorial, you will construct a half-edge graph in Rust, validate its invariants, find
its one independent cycle inside an explicitly selected subgraph, and print DOT. The builder
returns typed node indexes; Linnet's subgraph filters then select the half-edges on which an
algorithm operates.

For the shortest installation and first-run path, begin with
#link("quickstart/rust/")[Using Linnet from Rust]. Then return here to work with a boundary,
subgraphs, and cycle algorithms.

== Prerequisites

Create a Rust binary project and add the current Linnet release:

// docs-example: syntax
```sh
cargo new linnet-first-graph
cd linnet-first-graph
cargo add linnet@0.17
```

No drawing feature or external layout program is needed to construct the graph and emit DOT.

== Build a graph with a boundary

Replace `src/main.rs` with the following program:

// docs-example: run linnet-first-graph
```rust
use linnet::half_edge::{
    builder::HedgeGraphBuilder,
    involution::Flow,
    subgraph::{SubSetOps, SuBitGraph},
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

    let boundary: SuBitGraph = graph.external_filter();
    let internal = graph.full_filter().subtract(&boundary);
    let (basis, _) = graph.cycle_basis_of(&internal);
    assert_eq!(basis.len(), 1);
    println!("{}", graph.base_dot());
}
```

Run it with `cargo run`. Success means the assertion passes and stdout contains a DOT
`digraph` with nodes `a`, `b`, and `c`, three paired internal edges, and two dangling boundary
edges. The external edges do not add a cycle, so the triangle's cycle-basis rank remains one.

#callout("Verification scope and cost", [
  The docs harness compiles and runs this Rust program; it syntax-checks the setup commands
  without creating a project or using the network. Success requires `graph.check()` and the
  subgraph cycle-rank assertion to pass. The run is small; a clean external Cargo build is
  dependency-dominated and can take minutes, while the graph operation itself should finish in
  well under a second.
])

#callout("Keep the returned indexes with their graph", [
  `a`, `b`, and `c` are compact `NodeIndex` values owned by this graph construction. Do not
  serialize their numeric representation as a durable identity or reuse them against another
  graph. Mutating operations return the mappings needed to relate old and new storage.
])

== Make the first useful variation

Add a fourth internal edge between `a` and `b`, run the program again, and change the assertion
to `basis.len() == 2`. This demonstrates the key distinction: parallel paired half-edges add an
independent cycle, while dangling identity half-edges describe the graph boundary excluded by
`internal`.

From here, use a subgraph view to run traversal, cuts, or contraction against only the region
you intend to modify. Enable Linnet's `drawing` feature only when you need layout rather than
plain DOT interchange.

== Troubleshooting and next steps

- A type-inference error around `build()` means Rust cannot choose node storage or data types;
  keep the explicit `HedgeGraph<&str, &str>` annotation from the example.
- A failed `check()` points to an involution or storage invariant. Validate immediately after
  parsing or construction, before applying graph algorithms.
- If an algorithm includes an unexpected dangling edge, check that you passed `internal` rather
  than `graph.full_filter()`, then inspect the `Flow` assigned to each external edge.
- Continue with the subgraph-and-algorithms manual, then consult the Rust API reference for the
  exact return mappings of mutating operations and additional builder methods.
]
