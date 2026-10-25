#import "../../shared.typ": boundary, source-link

#let api = [
= Choose a Linnet interface

Use Rust when Linnet is part of an application or library, Python for standalone
graph workflows, and Linnest when a Typst document owns layout and drawing. Each interface
shares the half-edge model while exposing the conventions natural to its language.

== Typst package

Linnest combines Typst functions with a WebAssembly plugin. Follow the
#link("guides/linnest/")[Linnest workflow guide] for a first drawing, then use the focused
#link("reference/typst/")[Typst API reference] for its graph, layout, drawing, subgraph, and
template-composition surfaces. The curve helpers are Kurvst re-exports and stay documented in their
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

== Community graph bindings

The `linnet-py` binding library owns `symbolica.community.graph`. Install the complete
Symbolica Community wheel to use its builders, live graph views, subgraphs, DOT codecs,
layout streams, graph algorithms, and rendering. There is no separate Linnet Python wheel.

`RenderSettings` combines drawing options, layout sequences, and selectors.
`LayoutSettings`, `Stroke`, `DiagramView`, and `DiagramRender` are shared by graphs,
HepKit diagrams, and tensor networks. The notebook browser measures direct MathML labels
before Rust lays out and draws the graph. Layout pass sequences run natively;
unsupported settings raise errors.

```python
from symbolica.community import graph

settings = graph.RenderSettings(
    layouts=graph.LayoutSettings(algorithm=graph.LayoutAlgorithm.Dot),
    drawing=graph.DrawOptions(
        edge_stroke=graph.Stroke(paint=graph.Color("#6f4d85"), thickness=1.2),
    ),
)
```

Call `.render(config=settings)` to prepare a `DiagramView`, then display that same view.
After browser measurement completes, `view.measured()` returns an immutable `DiagramRender`.
Its SVG, HTML, and file exports reuse the snapshot, even after later graph or settings mutations.
HepKit's physics-specific options live in `DiagramStyle`, passed separately as `style`.
For multipage snapshots, use `to_svg_pages()`; `to_svg()` requires exactly one page.

`Graph.generate()` takes `EdgeSignature` values and returns graph/group-size pairs.
`HalfEdge` is a live graph view. Canonicalization and isomorphism compare expression-compatible
payloads; use `node_key`, `edge_key`, and `half_edge_key` for custom data. Canonicalization
preserves payload references and drawing metadata, which does not affect equivalence.

== Python graph workflows

Try the #link("playground/")[live Python playground] for editable DOT and a browser notebook,
or start with the runnable graph example in the #link("quickstart/python/")[Python quickstart].

The Python binding separates arbitrary application data from drawing data. `Node.data`,
`Edge.data`, and `HalfEdge.data` may contain any Python object. Linnet retains those objects by
reference: it does not copy, compare, hash, pickle, stringify, or send them to Typst. The matching
`drawing` properties contain only validated, typed values that the renderer can stage safely.

Build a graph declaratively with node specifications, endpoint specifications, and edge
specifications. The payload classes below belong to the example application; they are not
Linnet types:

```python
from dataclasses import dataclass

from symbolica.community.graph import build, edge, node, sink, source


@dataclass
class UserNodeData:
    state: object


@dataclass
class UserEdgeData:
    interaction: object


@dataclass
class UserHalfEdgeData:
    port: object


incoming_data = UserNodeData(state=object())
outgoing_data = UserNodeData(state=object())
propagator_data = UserEdgeData(interaction=object())

incoming = node("in", data=incoming_data)
outgoing = node("out", data=outgoing_data)
graph = build(
    incoming,
    outgoing,
    edge(
        source(incoming),
        "propagator",
        sink(outgoing),
        data=propagator_data,
    ),
    name="example",
)

assert graph.node("in").data is incoming_data
assert graph.edge("propagator").data is propagator_data
```

`graph.node("in")`, `graph.edge("propagator")`, and `graph.half_edge(0)` return live views.
Names and integer indices are accepted for node and edge lookup. Topology properties such as
`index`, `incidence`, `source`, `sink`, `pair`, and `flow` are read-only, while `data` and every
field on `drawing` are mutable. Endpoint-specific `source(..., data=...)` and
`sink(..., data=...)` payloads are optional; use them only when the application has distinct
half-edge data. Drawing fields include labels and styles as well as placement,
routing, rank, size, bend, decoration, and half-edge anchor information. Put a
template-specific field such as a particle or cut identifier in `drawing.extensions`; built-in
field names are reserved and Linnet assigns no domain meaning to extensions.

`Graph.map(node=..., edge=..., source=..., sink=...)` evaluates the supplied callbacks against
live element views and returns a new graph. An omitted callback preserves the exact Python data
object. Drawing and rendering configuration are copied by value. Reordering nodes or edges and
reversing an edge preserve each data/drawing association, but invalidate every older live view;
using a stale view raises `ReferenceError`.

Incremental edits reuse the declarative constructors instead of introducing a second set of
drawing dictionaries. `add_node` accepts a `NodeSpec`; `add_edge` accepts the same internal or
dangling `EdgeSpec` used by `build`. Declarative endpoints resolve a passed `NodeSpec`, name,
index, or live-node key. Incremental endpoints resolve a current live `Node`, name, index, or a
named `NodeSpec`:

```python
from symbolica.community.graph import Compass, Graph, edge, node, sink, source

graph = Graph()
graph.add_node(node("in", data=UserNodeData(object()), label="incoming"))
graph.add_node(node("out", data=UserNodeData(object())))

propagator = graph.add_edge(edge(
    source("in", data=UserHalfEdgeData(object()), compass=Compass.E),
    "propagator",
    sink(graph.node("out"), data=UserHalfEdgeData(object()), compass=Compass.W),
    data=UserEdgeData(object()),
    label="dependency",
))
```

Each successful insertion is one atomic topology revision and returns a fresh view. Reacquire a
view by name or index after the next edit. Validation failures do not publish the candidate and
therefore leave existing views live. Arbitrary data keeps Python object identity; specs and their
drawing dictionaries are copied by value and remain reusable.

`split_edge(at, replacement, ...)` turns a paired edge into two dangling edges. The selected
half-edge receives the detached replacement `EdgeValue`; the opposite half retains the original
edge value. `connect(source, sink, replacement, ...)` performs the topological inverse, with the
first dangling half-edge becoming the source. The two dangling
edge records are replaced by the supplied value, name, and orientation. Use
`delete(graph.subgraph(...))` for removal; there are no duplicate `remove_node` or `remove_edge`
aliases.

The long-lived topology is Linnet's native `HedgeGraph`; it is not a separate Python adjacency
model. `NodeStore.Vec` is the default dense node store, while `NodeStore.Forest` selects the
identification-preserving forest store. Pass `node_store=` to `Graph`, `build`, or a DOT importer,
inspect `graph.node_store`, and use `graph.to_node_store(...)` for an explicit non-mutating
conversion. Derived graphs preserve their input backend; append and join use the left graph's
backend.

A `Subgraph` is a graph-bound, revision-checked native half-edge selection, augmented only with
explicit zero-crown node IDs so isolated nodes remain first-class members. Construct one from named
or indexed nodes and edges, exact half-edge indices, or live-view predicates, then use the same
object with Linnet's topology algorithms and owning transformations:

```python
from symbolica.community.graph import DirectionBasis

# Explicit selections are unioned. Selecting a node includes its incident crown.
selected = graph.subgraph(nodes=["in"], edges=["propagator"])

# Predicates see the same live views as normal graph access.
propagators = graph.filter(
    edge=lambda value: isinstance(value.data, UserEdgeData),
)

component_copies = [graph.concretize(part) for part in graph.connected_components()]
tree = graph.depth_first_traverse("in", subgraph=selected)
cycles, covered = graph.cycle_basis(selected)

# Boundary operations return composable structural selections.
inside_boundary = graph.internal_boundary(selected)
incident_boundary = graph.boundary(selected)
external = graph.external_half_edges()

# Directed algorithms default to the underlying source/sink flow.
assert graph.is_reachable("in", "out")
underlying_order = graph.topological_order()
superficial_order = graph.topological_order(
    direction=DirectionBasis.Superficial,
)

# Enumeration is explicit about its potentially exponential search space.
circuits = graph.all_cycles(subgraph=selected, max_results=31)
one_bond = graph.find_bond(subgraph=selected, min_size=1, max_size=1)
tadpole_components = graph.tadpoles(["in"])

# The result owns a reduced topology; views into graph remain live.
reduced = graph.transitive_reduction()

for partition in graph.all_cuts(["in"], ["out"]):
    source_region = partition.source_side
    cut_edges = partition.edges
    target_region = partition.target_side
    relative_orientation = partition.boundary.orientation(cut_edges[0])

# Inclusion-minimal connected bipartition boundaries of size one.
bridge_bonds = graph.all_bonds(subgraph=selected, min_size=1, max_size=1)

# concretize() is non-mutating; insertion, split/connect, extract(), delete(),
# contract(), append_mut(), and join_mut() stale old views and selections.
detached = graph.concretize(propagators)
removed = graph.extract(propagators)
```

`Subgraph` supports immutable union, intersection, difference, symmetric difference, complement,
and their Python operators. The binding also exposes connectivity and component queries, DFS and
BFS trees, tree parents/children/ancestors and fundamental cycles, boundary/crown queries, bridges,
cycle bases and bounded circuit enumeration, spanning forests, cut partitions, bond queries,
tadpole components, contraction, disconnected append, and callback-driven joining of dangling
half-edges. `append()` and `join()` return a graph while their
`*_mut()` counterparts update the left graph; the right graph is unchanged. Transformations retain
arbitrary data by Python identity and copy drawing/configuration by value. Callback failures and
name collisions never publish a candidate topology. Join callbacks receive live views, so any
explicit `data` or `drawing` mutations they perform take effect immediately and are not rolled back.
Every cut terminal must be part of a disjoint source or target node group and have an incident
half-edge, while each combined terminal group must have a boundary. An isolated node mixed into
an otherwise connected terminal group is rejected instead of being silently dropped. `all_cuts`
eagerly materializes the native admissible separating
partitions and can be combinatorial. Its result is sorted by structural half-edge indices, while
`CutPartition.boundary.left` and `.right` retain the orientation of each boundary half-edge.
`all_bonds` enumerates inclusion-minimal connected bipartition boundaries in a size range; it is
not a minimum-cardinality or max-flow operation. `find_bond` returns one deterministic native
witness without materializing the complete family. `all_cycles(max_results=...)` checks the full
non-empty cycle-space candidate count before enumeration; `max_results` therefore bounds work even
when only a few combinations are circuits. Tadpole terminals must be distinct, valid, and incident
to at least one half-edge; identification-history aliases of the same structural node are duplicates.

`DirectionBasis.Underlying` is the default for reachability, topological order, and transitive
reduction. It follows the intrinsic half-edge source/sink roles used by GammaLoop momentum routing.
`DirectionBasis.Superficial` applies each edge's drawing-independent `Orientation`: reversed edges
reverse the arc and undirected edges contribute no arc. A direction choice is local to the
algorithm and never changes stored flow or orientation. `transitive_reduction()` returns a new
graph with retained arbitrary payload objects shared by identity; Python transitive closure is
intentionally absent until synthesized edges have an explicit data/drawing factory.
When `is_reachable(..., subgraph=...)` is used, both endpoint nodes must belong to that selection;
an endpoint outside it is a `ValueError`, rather than an implicit unreachable result.

Traversal-tree relationship methods validate that a node belongs to that exact tree revision.
`parent` returns `None` at the root, `children` returns immediate children, and `ancestors` runs from
the immediate parent through the root. `fundamental_cycle` returns `None` for a tree edge and rejects
dangling edges or edges whose endpoints are outside the traversal tree. An `OrientedCut` exposes its
edges, the orientation of each edge relative to the cut, and its winding number against a cycle from
the same graph revision.

=== Typed drawing configuration

Plain Python strings are data, never executable source expressions. `TextLabel` carries text
and `MathSymbol` carries mathematical labels. `AUTO` selects a field's automatic value, and
`INHERIT` leaves the value from the preceding configuration layer intact.
The closed native-value model also includes `Length`, `Ratio`,
`RelativeLength`, `Angle`, `Fraction`, `Color`, `Stroke`, `Dash`, `Insets`, `Mark`, `TextLabel`,
and `MathSymbol`, plus recursively validated arrays and string-keyed dictionaries. Finite
layout, placement, routing, anchor, pattern, mark, and debug choices use their exported enums;
the structured Python reference lists the choices for each option.

`GraphStyleOptions`, `LayoutSettings`, and `DrawOptions` cover the corresponding generic Linnest
surfaces. `DrawingSelectors` maps arbitrary Python data to detached typed drawing values at render
time. Layout passes retain their order:

```python
from symbolica.community.graph import (
    Color,
    DrawingSelectors,
    EdgeDrawing,
    LayoutAlgorithm,
    LayoutDirection,
    LayoutSettings,
    Length,
    RenderSettings,
    Stroke,
)

layouts = LayoutSettings(
    algorithm=LayoutAlgorithm.Force,
    direction=LayoutDirection.Right,
    steps=200,
).then(label_steps=80)
graph.render_config = RenderSettings(
    title="Example",
    layouts=layouts,
    selectors=DrawingSelectors(
        edge=lambda value: EdgeDrawing(
            label="dependency",
            style={
                "stroke": Stroke(
                    paint=Color("navy"),
                    thickness=Length.pt(0.8),
                ),
            },
        ) if isinstance(value.data, UserEdgeData) else None,
    ),
)

view = graph.render()
view  # the notebook browser measures labels and completes this view
```

SVG output preserves the native drawing and adds vertex, edge, and half-edge inspection.
Hover to see an element's identity; click to pin its details without changing the selection. Shift-, Ctrl-, or
Meta-click toggles the element in that displayed selection. Enter and Space perform the same
actions on a focused element; Escape dismisses the details while retaining the selection.
The first and last quarters of arc length pick the source and sink half-edges; the middle
half picks the whole edge. Both end quarters of a dangling edge pick its sole half-edge,
which selects the same subgraph as the whole edge. Half-edge highlighting covers its half
of a paired edge, or the entire dangling edge. Copy the panel's construction into
`graph.subgraph(nodes=[...], edges=[...], half_edges=[...])` to create a Python selection.
Figure selections are browser state and do not mutate the Python graph.
Ctrl/Meta-scroll gently zooms the graph around the pointer; drag to pan. With the graph
focused, `+` and `-` zoom in five-percent steps, arrow keys pan, and `0` fits the complete drawing.
Details appear beside the graph in wide outputs and below it in narrow outputs. Zooming
preserves the selection and pinned details, whose text stays at its original size.
Embedding applications can read the SVG element's `linnetSelection` property or listen for
`linnet-selection-change`; its event detail contains sorted `nodes`, `edges`, and `half_edges` arrays.
For custom overlays, `svg.linnetShade({nodes, edges, half_edges})` returns an SVG group
in the drawing's coordinates. Append it inside `.linnet-viewport` and style its color
and group opacity. It reuses native half-edge paths, compositing overlapping
half-edges and vertices once; it does not change the interactive selection.
The optional second argument controls `width` (edge-band width, default 8),
`padding` (beyond each node, default 2), and `outline` (boundary thickness,
default 0), in drawing coordinates. For example,
`svg.linnetShade({nodes: [0, 1], edges: [2]}, {width: 20, padding: 5, outline: 1.6, fillOpacity: 0.09})`
draws a filled circling with an opaque boundary around the union. `fillOpacity`
controls its interior independently of the boundary. Multiple groups can carry
different colors and padding for nested regions, as in the CFF explorer.
Native momentum arrows retain their edge inspection identity in
`[data-linnet-momentum]` groups, so specialized displays can reuse their curves
and arrowheads without duplicating the renderer.
Use the Python selection API explicitly when applying graph algorithms to those IDs.

A `Subgraph` also renders directly as the final expression of a notebook cell:

```python
region = graph.subgraph(edges=["propagator"])
region
```

The notebook view shows the full owning graph, highlights the selected
half-edges, and draws the complement dotted and muted through Linnest's subgraph drawing mode.
Selecting one half of a paired edge highlights only that half; selected isolated nodes are
highlighted too. The display uses the owner's layout, drawing values, and render configuration
without changing the graph or its stored settings. `region.render(config=None)` prepares a
`DiagramView`; retain it for browser display and later snapshot exports.
Like other selection operations, rendering checks the topology
revision and raises `ReferenceError` for a stale selection.

The graph's `render_config` is overlaid by a sparse per-call `config`.
Preparing a view evaluates selectors once and retains the resulting drawing inputs.
Arbitrary Python `.data` stays in Python; selectors may explicitly convert it to drawing values.
Display the retained view, then capture it in a later notebook cell after the drawing appears:

```python
snapshot = view.measured()
svg = snapshot.to_svg()
snapshot.save("diagram.svg")
snapshot.save("diagram.html")
```

Default `view.snapshot()`, `view.to_svg()` and `view.to_html()` export complete
static output without launching a browser. They reuse completed measurements,
otherwise use optional typst-py outlining when available, otherwise estimated
STIX MathML label bounds. Estimates leave measurement pending. Installing
`anywidget` adds live notebook measurement; earlier snapshots remain immutable.

An unfinished view raises a measurement-pending error only for strict `measured()`.
No browser is started automatically. Scripts can explicitly use
`view.measured(browser="chromium")`; install Playwright and Chromium first.
The measured SVG embeds its font and XHTML `foreignObject` MathML labels.
`snapshot.save_pages(directory, format="svg")` exports every captured SVG page.

Vector-label exports are a separate explicit operation using optional Python
`typst>=0.15`. They do not require a browser:

```python
vector = view.outlined()
vector.save("diagram.pdf")
vector.save("diagram.png")
vector.save("diagram.typ")
```

Outlining measures labels with Typst and runs native layout for those boxes; it does not
replace the browser-measured snapshot. PDF and PNG conversion uses the outlined SVG.
`vector.to_linnest()` emits portable Typst source embedding that SVG.
Default notebook rendering has no Typst compiler dependency. Standalone Typst packages
and the documentation builder's Typst compiler remain independent tools.

The notebook display is generic and model-neutral, including for graphs with dangling edges.
Domain renderers supply their own native scenes and labels. Use typed settings and Python
selectors for customization. `INHERIT` leaves the preceding configuration layer unchanged.

```python
config = RenderSettings(
    layouts=LayoutSettings(impred_steps=50).then(impred_steps=20),
    drawing=DrawOptions(node_fill=Color("teal")),
)
view = graph.render(config=config)
view
```

Layout passes execute natively in order. An ImPrEd pass following another ImPrEd pass retains
the previous geometry instead of constructing a fresh seed. Unsupported executable styles
raise an error instead of falling back to full-document Typst compilation.

`graph.render_config` is a live mutable object: assign its fields directly, or replace the whole
configuration. The nested style, layout, drawing, and selector option objects are immutable
values, and their getters return copies, so replace those fields rather than trying to mutate a
nested object in place. `overlay()` returns a new configuration, and a per-call overlay never
mutates the graph's stored defaults.

Python drawing selectors belong on `RenderSettings.selectors`. They may inspect arbitrary `.data`,
but must return the corresponding detached `NodeDrawing`, `EdgeDrawing`, or `HalfEdgeDrawing`, or
`None` when no defaults apply. Selector fields fill only drawing fields that are absent on the
element, so explicit element drawing remains authoritative. Only the validated drawing snapshot
crosses the drawing boundary. The example above happens to color an application-defined edge
payload; custom particle styling works the same way and is not a special Linnet mode:

```python
from dataclasses import dataclass

from symbolica.community.graph import (
    Color,
    DrawingSelectors,
    EdgeDrawing,
    Length,
    MathSymbol,
    RenderSettings,
    Stroke,
)


@dataclass
class ParticleStyle:
    color: Color
    symbol: str


def particle_drawing(value):
    if not isinstance(value.data, ParticleStyle):
        return None
    return EdgeDrawing(
        label=MathSymbol(value.data.symbol),
        style={"stroke": Stroke(paint=value.data.color, thickness=Length.pt(0.8))},
    )


particle_style_example = RenderSettings(
    selectors=DrawingSelectors(edge=particle_drawing),
)
```

This is ordinary userland settings composition: Linnet neither defines `ParticleStyle` nor
inspects it. GammaLoop's richer particle decorations, momentum annotations, and diagram modes
belong to domain renderers, exposed through HepKit's `DiagramStyle` rather than a
`PhysicsOptions` type in `linnet`.
The #source-link(
  "crates/linnet-py/examples/physics_render_settings.py",
  label: "editable DOT physics notebook",
) loads a model and parses compact or annotated DOT with `FeynmanDiagram.from_dot(model, dot)`.
FeynKit renders the parsed diagram with a shared `RenderSettings`, including momentum labels
from its stored loop momentum basis and interactive SVG hover and selection.
Run it from a checkout using a Python environment with Marimo, Linnet, and the Symbolica host
containing `symbolica.community.hepkit`:
`python -m marimo edit crates/linnet-py/examples/physics_render_settings.py`.
The #source-link("crates/linnet-py/examples/layout_stream.py", label: "streaming layout notebook")
previews the force solver as it runs. Its DOT editor and collapsible sliders restart the
simulation; pause/resume preserves the current state. The viewer creates SVG topology
once, then updates coordinates from `LayoutStream.from_dot(...)` frames. It also runs
in an exported Marimo WebAssembly notebook with the bundled Emscripten wheel.

```python
from symbolica.community.graph import LayoutStream

stream = LayoutStream.from_dot("digraph { a -> b; b -> c; c -> a; }", every=4)
node_names, endpoints = stream.node_names, stream.endpoints
for frame in stream:
    positions = frame.nodes
    edge_positions = frame.edges
```

The first frame has iteration zero; subsequent frames advance up to `every` iterations
until a final `done` frame. Batch size affects observation frequency, not the random
initialization, cooling, or depth-collapse schedule. No history is retained by the
iterator. Frames own their position data; changing an earlier frame's returned arrays
cannot affect later frames. Node names and endpoint indices stay in the same order.
This proof of concept supports a single DOT graph and the force algorithm. Its numeric
geometry includes explicit DOT placement constraints and per-edge spring lengths, but
not browser label measurement, styling, or final label relaxation. It does not alter an
existing Python `Graph` or its payloads.

Custom graph behavior is expressed through Python drawing selectors and native scenes.
Executable Typst module references cannot define graph geometry or styles. Use `TextLabel`
and `MathSymbol` for labels; their sizes are measured before native layout. Standalone
Typst documents can still use Linnest's separate Typst package.

=== Explicit DOT codecs

DOT conversion is never implicit because arbitrary Python data and rich drawing values have no
unique DOT representation. Supply a `DotCodec` to `Graph.from_dot`, `Graph.from_dot_set`,
`Graph.from_dot_file`, and `Graph.to_dot`:

```python
from dataclasses import dataclass

from symbolica.community.graph import (
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


@dataclass
class VertexPayload:
    name: str | None
    kind: str | None


@dataclass
class EdgePayload:
    kind: str | None


@dataclass
class PortPayload:
    statement: str | None

codec = DotCodec(
    encode_node=lambda value: DotVertexData(
        name=value.data.name,
        statements={"kind": value.data.kind or ""},
    ),
    decode_node=lambda dot: NodeValue(
        data=VertexPayload(dot.name, dot.statements.get("kind")),
        drawing=NodeDrawing(label=TextLabel(dot.name or "")),
    ),
    encode_edge=lambda value: DotEdgeData(
        statements={"kind": value.data.kind or ""},
    ),
    decode_edge=lambda dot: EdgeValue(
        data=EdgePayload(dot.statements.get("kind")),
    ),
    encode_half_edge=lambda value: DotHalfEdgeData(statement=value.data.statement),
    decode_half_edge=lambda dot: HalfEdgeValue(data=PortPayload(dot.statement)),
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
