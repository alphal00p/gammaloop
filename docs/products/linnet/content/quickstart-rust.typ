#import "../../shared.typ": callout

#let quickstart-rust = [
= Using Linnet from Rust

Use Linnet when a problem is naturally expressed as a half-edge graph with an explicit boundary.
This example builds the smallest such graph, validates its involution, and emits DOT that can be
inspected or handed to a renderer.

== Create a Rust project

// docs-example: syntax
```sh
cargo new linnet-quickstart
cd linnet-quickstart
cargo add linnet@0.17
```

Replace `src/main.rs` with:

// docs-example: run linnet-quickstart-boundary
```rust
use linnet::half_edge::{
    builder::HedgeGraphBuilder, involution::Flow, HedgeGraph,
};

fn main() {
    let mut builder = HedgeGraphBuilder::new();
    let source = builder.add_node("source");
    let sink = builder.add_node("sink");

    builder.add_edge(source, sink, "propagator", false);
    builder.add_external_edge(source, "incoming", false, Flow::Sink);
    builder.add_external_edge(sink, "outgoing", false, Flow::Source);

    let graph: HedgeGraph<&str, &str> = builder.build();
    graph.check().expect("the half-edge involution is valid");

    assert_eq!(graph.n_externals(), 2);
    println!("{}", graph.dot_label(&graph.full_filter()));
}
```

Run it with `cargo run`. The output is a DOT graph with one paired internal edge and two dangling
boundary edges. `graph.check()` is the important first invariant: run it immediately after
constructing or parsing a graph, before applying graph algorithms.

#callout("Half-edges make the boundary explicit", [
  An internal edge is a pair under the graph's involution. A dangling external edge maps to
  itself and carries a `Flow`. The `NodeIndex` values returned by the builder belong to this
  graph; do not reuse their numeric representation as a durable identifier.
])

Continue with the #link("tutorial/")[graph tutorial] to select a subgraph and compute a cycle
basis. Use the #link("guides/clinnet/")[Clinnet guide] to render existing DOT files, or the
#link("quickstart/typst/")[Typst guide] to construct and draw an editable graph with Linnest.

`linnet` is a library, so add it to a project with `cargo add`; `cargo install` and Cargo Binstall
apply to executables. The separately versioned Clinnet renderer can be installed with
`cargo install clinnet --locked`.
]
