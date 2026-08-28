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

#boundary("Rendering runs through Typst", [
  Linnet can compute layout coordinates, but the supported renderer is not a native Rust drawing
  backend. Clinnet invokes an external Typst 0.15 executable for command-line figure batches;
  `linnet-py` uses the `typst` Python package to compile the same prepared project in-process.
  Use #link("guides/clinnet/")[Clinnet] for batch rendering, or use
  #link("guides/linnest/")[Linnest] when a Typst document owns the final drawing.
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
specifications. The payload classes below belong to the example application; they are not
Linnet types:

```python
from dataclasses import dataclass

from linnet_py import build, edge, node, sink, source


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
from linnet_py import Compass, Graph, edge, node, sink, source

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
from linnet_py import DirectionBasis

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

=== Typed Typst configuration

Plain Python strings are always escaped Typst strings, never source expressions. `None` maps to
Typst `none`, `AUTO` maps to `auto`, and `INHERIT` leaves the value from the preceding
configuration layer intact. The closed native-value model also includes `Length`, `Ratio`,
`RelativeLength`, `Angle`, `Fraction`, `Color`, `Stroke`, `Dash`, `Insets`, `Mark`, `TextLabel`,
and `MathSymbol`, plus recursively validated arrays and string-keyed dictionaries. Finite
layout, placement, routing, anchor, pattern, mark, and debug choices use their exported enums;
the structured Python reference lists the choices for each option.

`GraphStyleOptions`, `LayoutOptions`, and `DrawOptions` cover the corresponding generic Linnest
surfaces. `DrawingSelectors` maps arbitrary Python data to detached typed drawing values at render
time. Layout passes retain their order:

```python
from linnet_py import (
    Color,
    DrawingSelectors,
    EdgeDrawing,
    LayoutAlgorithm,
    LayoutDirection,
    LayoutOptions,
    Length,
    RenderConfig,
    Stroke,
)

layouts = LayoutOptions(
    algorithm=LayoutAlgorithm.Force,
    direction=LayoutDirection.Right,
    steps=200,
).then(label_steps=80)
graph.render_config = RenderConfig(
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

output = graph.render("diagram.pdf")
svg = graph.to_svg()
graph  # the final expression in a notebook renders inline
```

The selected template's defaults are overlaid by the graph's `render_config` and then by a sparse
per-call `config`. `render(output, config=None)` writes PDF, SVG, or PNG according to the output
suffix. `to_svg(config=None)` returns SVG text, and `_repr_svg_()` supports notebooks. These are
the only high-level rendering methods; there are no raw command-line inputs or string-expression
escape hatches. The Python distribution depends on `typst` 0.15.0 and compiles in-process without
looking up or launching a Typst executable. Generated inputs and imported modules remain alive
until compilation completes.

The shipped notebook display is generic and model-neutral, including for graphs with dangling
edges. It never infers amplitude, cross-section, particle, or momentum semantics from topology.
An application that needs those concepts selects its own template and passes that template's
closed native settings through `RenderConfig.template_options`. `INHERIT` leaves the preceding
configuration layer unchanged, while `None` becomes Typst `none`.

```python
domain_template = RenderConfig(
    template="templates/domain-figure.typ",
    template_options={
        "theme": "review",
        "accent": Color("teal"),
    },
)
```

`template_options` is intentionally a string-keyed dictionary because the selected template,
not Linnet, defines that open schema. Its values still use the closed native-value model; it is
not a dictionary of Typst source strings.

`graph.render_config` is a live mutable object: assign its fields directly, or replace the whole
configuration. The nested style, layout, drawing, and selector option objects are immutable
values, and their getters return copies, so replace those fields rather than trying to mutate a
nested object in place. `overlay()` returns a new configuration, and a per-call overlay never
mutates the graph's stored defaults.

Python drawing selectors belong on `RenderConfig.selectors`. They may inspect arbitrary `.data`,
but must return the corresponding detached `NodeDrawing`, `EdgeDrawing`, or `HalfEdgeDrawing`, or
`None` when no defaults apply. Selector fields fill only drawing fields that are absent on the
element, so explicit element drawing remains authoritative. Only the validated drawing snapshot
crosses the process boundary. The example above happens to color an application-defined edge
payload; custom particle styling works the same way and is not a special Linnet mode:

```python
from dataclasses import dataclass

from linnet_py import (
    Color,
    DrawingSelectors,
    EdgeDrawing,
    Length,
    MathSymbol,
    RenderConfig,
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


particle_style_example = RenderConfig(
    selectors=DrawingSelectors(edge=particle_drawing),
)
```

This is ordinary userland settings composition: Linnet neither defines `ParticleStyle` nor
inspects it. GammaLoop's richer particle decorations, momentum annotations, and diagram modes
belong to GammaLoop's own template rather than to a `PhysicsOptions` type in `linnet_py`.
The #source-link(
  "crates/linnet-py/examples/physics_render_settings.py",
  label: "complete Python-only physics theme",
) composes edge and half-edge layers using the same generic API.
Typst callbacks that require measured geometry use an explicit module function reference:

```python
from linnet_py import TypstModule

styles = TypstModule.file("styles.typ")
graph.node("in").drawing.label = styles.content("incoming_label")
graph.edge("propagator").drawing.style = styles.function("edge_style").bind(
    emphasis=True
)
graph.edge("propagator").drawing.extensions = {
    "lane-map": styles.value("lane_map"),
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
  let options = config.at("options", default: (:))
  // Read typed element values and template-owned options, then render content.
  config.elements.nodes.at(0).label
}
```

The bundled generic and GammaLoop templates use this same contract. The generic template ignores
unknown domain concepts. GammaLoop's template owns its mode presets, particle decoration,
momentum arrowheads and index labels, amplitude edge ordering and side labels, cut-ID-matched
cross-section side labels, directional placement, and label-aware sizing. Python selects that
template explicitly and passes its template-specific settings through `template_options`.
Within GammaLoop's template, explicit element or configuration values override generated values,
which override template defaults. Final labels and styles are applied before layout so their
measured sizes affect spacing.

=== Explicit DOT codecs

DOT conversion is never implicit because arbitrary Python data and rich drawing values have no
unique DOT representation. Supply a `DotCodec` to `Graph.from_dot`, `Graph.from_dot_set`,
`Graph.from_dot_file`, and `Graph.to_dot`:

```python
from dataclasses import dataclass

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
