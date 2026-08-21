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

The Python binding separates arbitrary application data from drawing data. `Node.data`,
`Edge.data`, and `HalfEdge.data` may contain any Python object. Linnet retains those objects by
reference: it does not copy, compare, hash, pickle, stringify, or send them to Typst. The matching
`drawing` properties contain only validated, typed values that the renderer can stage safely.

Build a graph declaratively with node specifications, endpoint specifications, and edge
specifications:

```python
from linnet_py import Compass, MathSymbol, build, edge, node, sink, source

incoming_payload = object()
edge_payload = {"model id": 17}

incoming = node(
    "in",
    data=incoming_payload,
    label=MathSymbol("n", subscript=0),
)
outgoing = node("out", data={"kind": "final"})
graph = build(
    incoming,
    outgoing,
    edge(
        source(incoming, data={"side": "source"}, compass=Compass.E),
        "propagator",
        sink(outgoing, data={"side": "sink"}, compass=Compass.W),
        data=edge_payload,
        particle="e-",
        label=MathSymbol("p", subscript=0),
    ),
    name="example",
)
```

`graph.node("in")`, `graph.edge("propagator")`, and `graph.half_edge(0)` return live views.
Names and integer indices are accepted for node and edge lookup. Topology properties such as
`index`, `incidence`, `source`, `sink`, `pair`, and `flow` are read-only, while `data` and every
field on `drawing` are mutable. Drawing fields include labels and styles as well as placement,
routing, rank, size, bend, momentum, particle, cut, and half-edge anchor information. Put a
template-specific field in `drawing.extensions`; built-in field names are reserved.

`Graph.map(node=..., edge=..., source=..., sink=...)` evaluates the supplied callbacks against
live element views and returns a new graph. An omitted callback preserves the exact Python data
object. Drawing and rendering configuration are copied by value. Reordering nodes or edges and
reversing an edge preserve each data/drawing association, but invalidate every older live view;
using a stale view raises `ReferenceError`.

=== Typed Typst configuration

Plain Python strings are always escaped Typst strings, never source expressions. `None` maps to
Typst `none`, `AUTO` maps to `auto`, and `INHERIT` leaves the value from the preceding
configuration layer intact. The closed native-value model also includes `Length`, `Ratio`,
`RelativeLength`, `Angle`, `Fraction`, `Color`, `Stroke`, `Dash`, `Insets`, `Mark`, `TextLabel`,
and `MathSymbol`, plus recursively validated arrays and string-keyed dictionaries. Finite
layout, placement, routing, anchor, pattern, mark, and debug choices use their exported enums;
the structured Python reference lists the choices for each option.

`GraphStyleOptions`, `LayoutOptions`, `DrawOptions`, and `PhysicsOptions` cover the corresponding
Linnest surfaces. Layout passes retain their order:

```python
from linnet_py import (
    AUTO,
    Color,
    DebugLevel,
    DrawOptions,
    GraphStyleOptions,
    LayoutAlgorithm,
    LayoutDirection,
    LayoutOptions,
    Length,
    PhysicsOptions,
    RenderConfig,
    Stroke,
    TextLabel,
)

layouts = LayoutOptions(
    algorithm=LayoutAlgorithm.Force,
    direction=LayoutDirection.Right,
    steps=200,
).then(label_steps=80)
graph.render_config = RenderConfig(
    title=TextLabel("Example", size=Length.pt(11)),
    style=GraphStyleOptions(node_label=AUTO),
    layouts=layouts,
    drawing=DrawOptions(
        debug=DebugLevel.Off,
        node_radius=AUTO,
        edge_stroke=Stroke(paint=Color("black"), thickness=Length.pt(0.6)),
    ),
    physics=PhysicsOptions(
        momentum_arrows=True,
        show_momentum=True,
        show_edge_index=True,
    ),
)

output = graph.render("diagram.pdf")
svg = graph.to_svg()
graph  # the final expression in a notebook renders inline
```

The four configuration layers are shipped defaults, the selected diagram-mode preset, the
graph's `render_config`, and a sparse per-call `config` overlay. `render(output, config=None)`
writes PDF, SVG, or PNG according to the output suffix. `to_svg(config=None)` returns SVG text,
and `_repr_svg_()` supports notebooks. These are the only high-level rendering methods; there are
no raw command-line inputs or string-expression escape hatches. Rendering requires Typst 0.15 or
newer and keeps generated inputs and imported modules alive until compilation completes.

`graph.render_config` is a live mutable object: assign its fields directly, or replace the whole
configuration. The nested style, layout, drawing, and physics option objects are immutable values,
and their getters return copies, so replace those fields rather than trying to mutate a nested
object in place. `overlay()` returns a new configuration, and a per-call overlay never mutates the
graph's stored defaults.

Python drawing selectors belong on `GraphStyleOptions`. They may inspect arbitrary `.data`, but
node selectors must return a style dictionary or module value, while edge, source, and sink
selectors return one style layer or an array of layers. Only that validated result crosses the
process boundary. Typst callbacks that require measured geometry use an explicit module function
reference:

```python
from linnet_py import TypstModule

styles = TypstModule.file("styles.typ")
graph.node("in").drawing.label = styles.content("incoming_label")
graph.edge("propagator").drawing.style = styles.function("edge_style").bind(
    emphasis=True
)
graph.edge("propagator").drawing.extensions = {
    "particle-map": styles.value("particle_map"),
}
```

`TypstModule.file` and `TypstModule.package` are the only routes for custom executable Typst.
`value`, `content`, and `function` select a typed export; `call` and `bind` construct typed
applications. The binding never accepts a raw Typst expression and never calls Typst `eval`.

A custom rendering template is a Typst module with one mandatory V1 export:

```typst
#let render(config) = {
  assert(config.schema == "linnet-render-config")
  assert(config.version == 1)
  // Read typed values from config and return the rendered content.
  config.elements.nodes.at(0).label
}
```

Bundled generic and GammaLoop templates use this same contract. Diagram-mode presets retain
particle decoration, momentum arrowheads and index labels, amplitude edge ordering and side
labels, cut-ID-matched cross-section side labels, directional placement, and label-aware sizing.
Explicit element or configuration values override mode-generated values, which override the
template defaults. Final labels and styles are applied before layout so their measured sizes
affect spacing.

=== Explicit DOT codecs

DOT conversion is never implicit because arbitrary Python data and rich drawing values have no
unique DOT representation. Supply a `DotCodec` to `Graph.from_dot`, `Graph.from_dot_set`,
`Graph.from_dot_file`, and `Graph.to_dot`:

```python
from linnet_py import (
    DotCodec,
    DotEdgeData,
    DotHalfEdgeData,
    DotVertexData,
    EdgeValue,
    Graph,
    HalfEdgeValue,
    NodeDrawing,
    NodeValue,
    TextLabel,
)

codec = DotCodec(
    encode_node=lambda value: DotVertexData(
        name=value.data["name"],
        statements={"kind": value.data["kind"]},
    ),
    decode_node=lambda dot: NodeValue(
        data={"name": dot.name, "kind": dot.statements["kind"]},
        drawing=NodeDrawing(label=TextLabel(dot.name or "")),
    ),
    encode_edge=lambda value: DotEdgeData(
        statements={"kind": value.data["kind"]},
    ),
    decode_edge=lambda dot: EdgeValue(
        data={"kind": dot.statements.get("kind")},
    ),
    encode_half_edge=lambda value: DotHalfEdgeData(statement=value.data),
    decode_half_edge=lambda dot: HalfEdgeValue(data=dot.statement),
)

source = r'''digraph demo {
  incoming [kind="external"]
  outgoing [kind="external"]
  incoming -> outgoing [kind="propagator"]
}'''
graph = Graph.from_dot(source, codec)
canonical = graph.to_dot()  # reuses the codec stored by from_dot
```

Codec callbacks receive detached `NodeValue`, `EdgeValue`, or `HalfEdgeValue` snapshots with only
`data` and `drawing`; they never receive graph or index context. `GlobalData` remains separate, as
in the Rust serializer contract, and is not staged during rendering. An explicit codec passed to
`to_dot` overrides the stored codec, and callback exceptions propagate without publishing a
partial graph.

Use `DotCodec.linnest()` to retain canonical, non-executable DOT-representable drawing fields while
dropping arbitrary element `.data`. It rejects Typst module references, calls, and bindings on
both encode and decode; use a custom codec when trusted executable values must be reconstructed
explicitly. Use `DotCodec.topology()` to drop element data and drawing. Both helpers preserve graph
and element names as well as the separate `GlobalData`. Round trips promise canonical semantics
rather than the original whitespace or attribute order. Consult the
#link("reference/python/")[structured Python API] for every drawing field, option, enum, and DOT
data constructor.
]
