#import "../../shared.typ": boundary

#let api = [
= Choose a Linnet interface

Use Rust when Linnet is part of an application or library, Python for standalone DOT-backed
graph workflows, and Linnest when a Typst document owns layout and drawing. Each interface
shares the half-edge model while exposing the conventions natural to its language.

== Typst package

Linnest combines Typst functions with a WebAssembly plugin. Follow the
#link("guides/linnest/")[Linnest workflow guide] for a first drawing, then use the focused
#link("reference/typst/")[Typst API reference] for its graph, layout, drawing, physics, and
subgraph surfaces. The curve helpers are Kurvst re-exports and stay documented in their
canonical GammaLoop guide.

== Rust package

The `linnet` package is organized into these principal areas:

- `half_edge` provides `HedgeGraph`, builders, the involution and edge data, subgraph types,
  traversal trees, and graph algorithms;
- `parser` provides DOT parsing and DOT-backed graph data;
- `half_edge::layout` provides layout helpers when the `drawing` feature is enabled;
- `permutation`, `tree`, and `union_find` provide supporting algorithms and stores.

The generic graph types separate edge, vertex, and half-edge payloads from the node-storage
implementation. Start with the #link("reference/rust/")[Rust orientation], then use the
revision-specific Rustdoc for complete generic bounds, implementations, and method signatures;
choose a concrete storage type according to the workflow and mutation behavior you need.

Cargo features enable optional capabilities:

- `serde` adds serialization traits;
- `bincode` adds binary encoding support;
- `drawing` enables Linnet's layout module and its optional geometry dependency;
- `symbolica` adds Symbolica interoperability.

An item shown in an all-features API reference may not be available in a default build. Check
its required feature and your `Cargo.toml` before using it.

#boundary("Rendering runs through Typst", [
  Linnet can compute layout coordinates, but the supported renderer is not a native Rust drawing
  backend. Clinnet's shared renderer invokes an external Typst 0.15 or newer executable over
  Linnest and Kurvst assets. Use #link("guides/clinnet/")[Clinnet] for command-line figure
  batches, or use #link("guides/linnest/")[Linnest] when a Typst document owns the final drawing.
])

== Standalone Python distribution

#boundary("Distribution and import have different names", [
  Install the Python distribution `linnet-py`, then import `linnet_py`. It requires Python 3.10
  or newer. It is a standalone extension and is not a `symbolica.community` module.
  Its package version is independent of the Rust `linnet` version. Record both versions when
  diagnosing compatibility between Rust and Python code.
])

The Python binding focuses on DOT-backed graphs and mirrors typed identifiers and subgraphs. Its
builder accepts native names and mappings directly; string values are accepted for the common
orientation and external-flow choices:

```python
from linnet_py import DotGraphBuilder

builder = DotGraphBuilder()
left = builder.add_node("left", {"label": "incoming vertex"})
right = builder.add_node("right", {"label": "outgoing vertex"})
builder.add_external_edge(left, {"label": "p0"}, flow="sink")
builder.add_edge(left, right, {"label": "p1"}, orientation="default")
graph = builder.build()

print(graph.dot())
```

`orientation` accepts `"default"`, `"reversed"`, or `"undirected"`, while `flow` accepts
`"source"` or `"sink"`; the corresponding `Orientation` and `Flow` objects remain available
when typed values are preferable.

`DotGraph.render(output, template=None, inputs=None, typst='typst') -> Path` writes a PDF, SVG,
or PNG according to the output suffix. `DotGraph.to_svg(template=None, inputs=None,
typst='typst') -> str` returns SVG text without requiring the caller to manage a temporary file.
In a notebook, displaying a graph calls the same SVG path automatically through `_repr_svg_`:

```python
output = graph.render("diagram.pdf", inputs={"title": "Example"})
svg = graph.to_svg()
graph  # the final expression in a notebook renders inline
```

These methods invoke the external `typst` command, which must be Typst 0.15 or newer; `typst`
can name another executable or an absolute path. The default render uses the generic Clinnet
figure template and the embedded Linnest and Kurvst package assets. Those assets describe generic
graph layout and drawing, not GammaLoop's model-specific particle decorations. To preserve photon,
gluon, fermion, and other generated styles, pass the GammaLoop bundle's `figure.typ` as `template`
and keep its generated `edge-style.typ`, layout files, and package tree together beneath the same
Typst project root.

`DotGraph` can also parse one graph from a string or file, or a set of graphs from a string.
The exported classes include `Hedge`, `NodeIndex`, `EdgeIndex`, `Flow`, `Orientation`, graph
attribute wrappers, `HedgePair`, `Subgraph`, `Cycle`, `OrientedCut`, `TraversalTree`,
`DotGraph`, and `DotGraphBuilder`.

Mutable statement wrappers expose mapping-like access to DOT attributes. Consult the structured
#link("reference/python/")[Python API] for their lifetime and write-through behavior, and do not
assume that a copied mapping remains attached to its graph.
]
