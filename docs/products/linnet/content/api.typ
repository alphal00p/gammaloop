#import "../../shared.typ": boundary, catalog-contract, source-link

#let api = [
= Rust and Python API boundaries

#catalog-contract(
  rust-scope: "linnet",
  python-scope: "linnet_py",
)

== Rust package

The `linnet` package exposes these principal areas:

- `half_edge` owns `HedgeGraph`, builders, the involution and edge data, subgraph types,
  traversal trees, and graph algorithms;
- `parser` owns DOT parsing and DOT-backed graph data;
- `drawing` owns rendering data and layout helpers;
- `permutation`, `tree`, and `union_find` provide supporting algorithms and stores.

The generic graph types separate edge, vertex, and half-edge payloads from the node-storage
implementation. Generated reference pages own the full generic bounds and method signatures.
Guides should name a concrete type only when that choice matters to the workflow.

Feature gates affect the compiled reference surface:

- `serde` adds serialization traits;
- `bincode` adds binary encoding support;
- `drawing` adds the drawing stack;
- `symbolica` adds Symbolica interoperability.

Reference generation must record feature availability on each affected item. An all-features
catalog is useful for discovery, but it must not imply that every item is present in a default
build.

== Standalone Python distribution

#boundary("Distribution and import have different names", [
  Install the Maturin distribution `linnet-py`, then import `linnet_py`. It requires Python
  3.10 or newer. It is a standalone extension and is not a `symbolica.community` module.
  Its package version is independent of the Rust `linnet` version, so documentation snapshots
  must display both component versions.
])

The Python binding focuses on DOT-backed graphs and mirrors typed identifiers and subgraphs:

```python
from linnet_py import DotGraphBuilder

builder = DotGraphBuilder()
left = builder.add_node()
right = builder.add_node()
builder.add_edge(left, right)
graph = builder.build()

print(graph.dot())
```

`DotGraph` can also parse one graph from a string or file, or a set of graphs from a string.
The exported classes include `Hedge`, `NodeIndex`, `EdgeIndex`, `Flow`, `Orientation`, graph
attribute wrappers, `HedgePair`, `Subgraph`, `Cycle`, `OrientedCut`, `TraversalTree`,
`DotGraph`, and `DotGraphBuilder`.

Mutable statement wrappers expose mapping-like access to DOT attributes. Their lifetime and
write-through behavior are part of the generated Python descriptor; examples should avoid
presenting a copied mapping as though it were always attached to its graph.

Source starting points are #source-link("crates/linnet/src", label: "the Rust package") and
#source-link("crates/linnet-py/src/lib.rs", label: "the Python extension").
]
