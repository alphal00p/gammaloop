#import "../../shared.typ": boundary

#let api = [
= Rust and Python APIs

== Rust package

The `linnet` package is organized into these principal areas:

- `half_edge` provides `HedgeGraph`, builders, the involution and edge data, subgraph types,
  traversal trees, and graph algorithms;
- `parser` provides DOT parsing and DOT-backed graph data;
- `drawing` provides rendering data and layout helpers;
- `permutation`, `tree`, and `union_find` provide supporting algorithms and stores.

The generic graph types separate edge, vertex, and half-edge payloads from the node-storage
implementation. Consult the Rust API reference for complete generic bounds and method
signatures; choose a concrete storage type according to the workflow and mutation behavior you
need.

Cargo features enable optional capabilities:

- `serde` adds serialization traits;
- `bincode` adds binary encoding support;
- `drawing` adds the drawing stack;
- `symbolica` adds Symbolica interoperability.

An item shown in an all-features API reference may not be available in a default build. Check
its required feature and your `Cargo.toml` before using it.

== Standalone Python distribution

#boundary("Distribution and import have different names", [
  Install the Python distribution `linnet-py`, then import `linnet_py`. It requires Python 3.10
  or newer. It is a standalone extension and is not a `symbolica.community` module.
  Its package version is independent of the Rust `linnet` version. Record both versions when
  diagnosing compatibility between Rust and Python code.
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

Mutable statement wrappers expose mapping-like access to DOT attributes. Consult the Python API
reference for their lifetime and write-through behavior, and do not assume that a copied mapping
remains attached to its graph.
]
