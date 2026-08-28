= GammaLoop Drawing Architecture

#quote(block: true)[
#strong[Status:] Current contract; audited against the implementation on 2026-08-28

The pipeline, ownership, DOT parser fields, typed rendering data, path styles, and subgraph
sections describe the implemented drawing stack. The superseded `*-eval` and
placeholder-interpolation design is preserved separately in
#link("gammaloop-drawing-evaluated-fields-history.typ")[the evaluated-field proposal record].
]

This document describes the drawing path used by GammaLoop and the DOT syntax
that matters at each layer. The important boundary is:

- GammaLoop owns physics graph data, model-to-style policy, and its domain render template.
- Linnet/Linnest owns DOT parsing, half-edge graph shape, layout, and generic drawing primitives.
- Kurvst owns Bezier splitting, trimming, and path patterns.
- User-authored render templates are an optional drawing feature, not part of
  GammaLoop's generated physics interaction.

== Pipeline

GammaLoop drawing uses the same DOT files that describe the Feynman graphs. A
rendering run has these steps:

+ GammaLoop writes DOT graph files.
+ GammaLoop writes a model-specific `edge-style.typ`. This file maps particle
   names to Typst style dictionaries and exports callbacks named
   `source-style`, `sink-style`, and `edge-label`.
+ The drawing templates are extracted to `drawings/templates/`. GammaLoop owns
   its app templates there (`figure.typ`, `grid.typ`, `layout.typ`, `layout-core.typ`), while the
   shared Linnest/Kurvst package files keep their canonical workspace layout
   under `drawings/templates/crates/{linnest,kurvst}/typst/`.
+ Clinnet's shared preparation layer stages the topology DOT and one generated Typst entrypoint
   for each graph. The `linnet` CLI compiles it with an external Typst 0.15 executable.
   `linnet-py` passes the same prepared entrypoint and project root to `typst-py` 0.15.0 for
   in-process compilation; it does not discover or launch an executable.
+ The entrypoint statically imports the selected template, imports each referenced user module
   once, constructs the native configuration dictionary, adds the staged `data-path`, and calls
   the template's mandatory `render(config)` export. Configuration is not flattened into
   `sys.inputs`. Its size does not change Clinnet's process argument list, and the Python path has
   no compiler subprocess arguments at all.
+ The generic rendering path materializes Clinnet's embedded figure/layout templates and the
   embedded Linnest and Kurvst package assets. It has no particle, momentum, amplitude, or
   cross-section mode. The GammaLoop path supplies its generated template bundle instead, binding
   the model-specific `edge-style.typ` callbacks and interpreting GammaLoop-owned options for its
   generic, amplitude, and cross-section presets.
+ The selected `layout.typ` parses the staged DOT topology, attaches the typed node, edge, and
   half-edge records from `config.elements`, and patches half-edge statements, port labels, and
   compass points into the native graph by topology index. It applies final labels and styles with
   `graph.style`, and only then runs the ordered layout passes and `draw`. Label and node
   measurements therefore participate in spacing.
+ `draw` uses ordinary Linnest node and edge styles for a model-neutral graph. GammaLoop's selected
   template installs the particle and momentum callbacks from `edge-style.typ`. It draws edges and
   labels first, then nodes last so nodes sit on top of edges.

The rendering contract has no evaluated-string mode. Generated particle styles are ordinary
Typst dictionaries/functions in `edge-style.typ`, backed by
#link("../../assets/embedded/drawing/templates/physics-edge-style.typ")[GammaLoop's
`physics-edge-style.typ`].
Custom executable Typst enters only through an explicitly imported module; DOT strings and plain
Python strings remain data.

== Components

/ `crates/gammalooprs`: Generates physics DOT and the model-specific `edge-style.typ`. It decides how
  a particle should look: photon wave, gluon coil, scalar dashed line, mass
  dependent thickness, labels, and related policy.

/ `crates/linnet`: Provides the half-edge graph data structure and DOT parser. Its parser turns
  invisible DOT nodes into dangling half-edges and preserves unconsumed
  attributes as statement dictionaries.

/ `crates/linnest`: Provides the Typst-facing wasm plugin and Typst wrappers. It parses DOT into
  graph bytes, lays out node and edge positions, exposes graph queries, and
  provides the domain-neutral CeTZ draw API.

/ `crates/kurvst`: Provides the Typst-facing curve wasm plugin. It splits and trims Bezier
  curves, builds Hobby curves through edge layout points, and generates wave,
  zigzag, and coil path patterns.

/ `crates/clinnet`: Provides the shared external-Typst V1 renderer, embedded generic Clinnet,
  Linnest, and Kurvst assets, and the `linnet` CLI used to batch-render DOT files and assemble
  grid PDFs. Every custom figure template must export `render(config)`.

/ `crates/linnet-py`: Owns Linnet `HedgeGraph` topology through a selectable vector or forest node
  store. It exposes declarative construction, graph-bound subgraphs, topology algorithms and
  transformations,
  typed drawing metadata and render configuration, explicit DOT codecs, and module references.
  Subgraphs retain Linnet's native half-edge set and add only the zero-crown node IDs that the set
  cannot represent, keeping isolated nodes visible to Python graph operations.
  Its rendering methods and notebook SVG display reuse Clinnet's prepared project, then compile it
  in-process through the `typst` Python package. This is a second compiler backend for the same
  rendering contract, not a separate template or staging implementation.

GammaLoop's particle policy is not part of Clinnet's generic embedded template. Python callers
that want GammaLoop particle colors, labels, path decorations, or mode presets select the
GammaLoop `figure.typ` entry point, pass its settings through `RenderConfig.template_options`,
and keep its imported template bundle together. Callers implementing another domain can instead
use generic `DrawingSelectors` to derive typed drawing records from arbitrary Python data.

== Data Ownership

GammaLoop DOT has two kinds of data.

Physics data is read by GammaLoop when a graph is parsed as a physics graph.
Examples are `particle`, `pdg`, `lmb_id`, `num`, `int_id`, and
`overall_factor`.

Generic drawing data is read by Linnest's drawing template. Examples are display labels,
half-edge compass points, and subgraph selections by compass. Particle and cut identifiers are
GammaLoop template extensions rather than Linnest drawing fields.

The Python API makes this boundary explicit. `Node.data`, `Edge.data`, and
`HalfEdge.data` retain arbitrary Python objects by identity; the renderer never inspects or
serializes them. Their corresponding mutable `drawing` views accept only the closed native value
model and template extensions. Python drawing selectors may inspect live element views, but only
their validated results enter `config.elements`.

DOT conversion is a separate, explicit boundary. Python import and export require a `DotCodec`.
The canonical Linnest codec retains its supported drawing subset while dropping arbitrary element
`.data`; the topology codec drops element data and drawing. Both preserve graph and element names
and the separate `GlobalData`. Rendering does not use either codec to guess a representation for
arbitrary `.data`: it stages topology separately and sends typed drawing records through
`config.elements`.

The same DOT file may contain both kinds of data for drawing. If a manually
decorated DOT file is later fed back into GammaLoop's physics parser, keep in
mind that drawing-only fields may produce unknown-attribute warnings unless the
physics parser has explicitly whitelisted them. The drawing templates treat
GammaLoop-only physics fields as ordinary metadata unless a callback chooses to
use them.

== DOT Shape

Use DOT `digraph` syntax:

```dot
digraph demo {
  graph [overall_factor="1"];
  edge [particle="a"];

  ext0 [style=invis];
  ext0 -> v0:0 [id=0, is_cut=0];
  v0:1 -> v1:2 [id=1, particle="d", lmb_id=0];
  v1:3 -> ext1 [id=2, particle="a"];
  ext1 [style=invis];
}
```

Important parser rules:

- `style=invis` on a node marks that node as a dangling external half-edge. It
  is not a normal drawn node in the half-edge graph.
- An edge from an invisible node to a real node becomes an incoming dangling
  half-edge at the real node.
- An edge from a real node to an invisible node becomes an outgoing dangling
  half-edge at the real node.
- An edge between two invisible nodes is invalid.
- `node:port` and `node:port:compass` are preserved as half-edge data. The port
  becomes the hedge id/port label; the compass is used by drawing subgraphs such
  as `subgraph.compass(g, "e")`.
- `edge [key=value]` and `node [key=value]` defaults are merged into individual
  edges and nodes before GammaLoop or Linnest sees the final statements.
- Attribute values are handled as strings after DOT parsing. Quote complex
  expressions and values containing spaces or punctuation.

== GammaLoop DOT Syntax

GammaLoop expects canonical attribute names. Some older aliases may be mentioned
by warnings, but new DOT should use the names below.

=== Graph Attributes

/ `num`: Global numerator factor. Default is `1`.

/ `overall_factor`: Symbolica expression multiplying the graph. Default is `1`.

/ `projector`: Optional projector expression. If omitted, GammaLoop builds the polarization
  projector from external particles.

/ `params`: Semicolon-separated Symbolica expressions used as additional parameters.

/ `group_id`: Optional graph-group id. Graphs with the same group id are evaluated as one
  group.

/ `is_group_master`: Boolean marking the master graph inside a group. If no master is provided,
  GammaLoop chooses one.

Export-only graph attributes include `overall_factor_evaluated`; they are useful
for inspection but are not input knobs.

=== Node Attributes

/ `int_id`: UFO vertex-rule id. If omitted, GammaLoop infers the vertex rule from the
  oriented incident particles when possible.

/ `num`: Explicit vertex numerator. If present, it is used instead of a UFO vertex
  rule.

/ `dod`: Degree of divergence override for the vertex.

/ `name`: Optional semantic name stored on the parsed vertex. The DOT node id itself is
  still the graph topology handle.

GraphViz-only presentation fields such as `label`, `shape`, `style`, `pos`,
`color`, and `fillcolor` are ignored by GammaLoop's physics parser.

=== Edge Attributes

/ `id`: Numeric edge id. Linnet consumes this as the internal edge index rather than
  keeping it as a normal edge statement. Drawing callbacks expose the drawn edge
  index as `eid`.

/ `particle`: Model particle name, for example `"a"`, `"d"`, `"d~"`, `"g"`, `"ghG"`, or
  `"W+"`.

/ `pdg`: Alternative to `particle`; looked up through the model PDG code.

/ `mass`: Symbolica mass expression. With `particle`, this overrides the model mass for
  that edge. Without `particle`, it creates a mass-only scalar-like edge.

/ `dir`: DOT direction/orientation. `forward` means default orientation, `back` means
  reversed, and `none` means undirected. If omitted, GammaLoop derives the
  orientation from the particle.

/ `source`: Half-edge payload for the source half-edge. GammaLoop expects JSON5 when the
  payload carries structured data, for example `source="{ufo_order:2}"`.

/ `sink`: Half-edge payload for the sink half-edge, with the same JSON5 convention as
  `source`.

/ `lmb_id`: Loop-momentum-basis id for a chosen loop edge.

/ `is_cut`: Hedge id used to mark an initial-state cut/external cut.

/ `num`: Explicit edge numerator. The parser localizes `edgeid(...)`, `sourceid(...)`,
  and `sinkid(...)` placeholders to the concrete edge and hedge indices.

/ `dod`: Degree of divergence override for the edge.

/ `name`: Optional semantic edge name.

/ `is_dummy`: Boolean marking a dummy edge. Dummy edges are filtered out of some physics
  operations and vertex matching.

/ `momtrop_edge_power`: Optional Symbolica expression controlling the momentum power used by the
  momtrop sampler. This does not change the graph topology.

/ `vakint_edge_power`: Optional integer controlling the momentum power used by vakint evaluation.
  This does not change the graph topology.

Physics DOT exporters may also write fields such as `lmb_rep`, `dod_autogen`,
`num_autogen`, and `name_autogen`. These are inspection/export metadata, not
normal user input.

== Drawing DOT Syntax

The drawing templates receive the parsed graph after Linnest layout. All
unconsumed DOT statements are available to callbacks as edge or node data.

=== Node Data In Drawing

The `draw` callback data for nodes includes:

- `vid`: zero-based node index in the drawn graph.
- `node`: the full node object.
- `name`: the node name.
- every node statement preserved from DOT.

By default, `draw` uses the node name as the label and computes a circle radius
that fits the label. Users can override this in Typst with `node-label` and
`node-style` callbacks.

=== Edge Data In Drawing

The `draw` callback data for edges includes:

- `eid`: zero-based edge index in the drawn graph.
- `edge`: the full edge object.
- `source-statement`: source half-edge statement, if present.
- `sink-statement`: sink half-edge statement, if present.
- `source-half-edge`: source half-edge object with node, hedge, port, and compass data.
- `sink-half-edge`: sink half-edge object with node, hedge, port, and compass data.
- `orientation`: `default`, `reversed`, or `undirected`.
- `ext`: boolean, true for dangling half-edges.
- every edge statement preserved from DOT.

GammaLoop's generated `edge-style.typ` uses `particle` to look up the default
edge style. A user can add drawing-only fields without affecting the generated
GammaLoop styles.

== Typed Render Configuration And Precedence

`Graph.render_config` stores typed `GraphStyleOptions`, ordered `LayoutOptions`, `DrawOptions`,
generic `DrawingSelectors`, and the selected template's closed native `template_options`. A sparse
configuration supplied to `render` or `to_svg` overlays the rendering layers in this order:

```text
selected-template defaults -> graph configuration -> per-call configuration
```

The `Graph.render_config` object itself is live and mutable, while its nested option objects are
immutable values whose getters reconstruct copies. Replace `style`, `layouts`, `drawing`,
`selectors`, or `template_options` after changing an option. `RenderConfig.overlay` returns a new
configuration, and per-call overlays do not mutate the graph configuration. The template-options
dictionary is open only in its string keys; every value still belongs to the validated native
Typst model.

The shipped Python notebook display is always generic, even when topology has dangling edges. Core
Linnet has no `RenderMode`, `PhysicsOptions`, or topology-based diagram inference. GammaLoop's CLI
and generated bundle select GammaLoop's `figure.typ` and provide its mode and particle settings as
template-owned options. Another application can select another template without teaching Linnet
its domain vocabulary.

The renderer correlates every typed node, edge, and half-edge record with its topology index
before GammaLoop mode generation. Explicit element or configuration values override generated
amplitude/cross-section values, which override template defaults. This preserves particle styles,
momentum arrowheads and index labels, amplitude edge ordering and side labels, cut-ID-matched
cross-section side labels, and directional placement without converting drawing values to DOT
strings.

Generic Python selectors run before staging and may inspect arbitrary element `.data`. They return
detached `NodeDrawing`, `EdgeDrawing`, or `HalfEdgeDrawing` values and fill only fields that the
element did not set explicitly. Selector results never become Typst callbacks and arbitrary
Python objects never enter `config.elements`.

Plain strings are escaped content. `None`, `AUTO`, and `INHERIT` have distinct native meanings;
lengths, colors, strokes, marks, labels, layouts, routing, and the remaining finite choices use
their value classes or enums. A `TypstModule` reference is the only executable extension point:
its `value`, `content`, or `function` export is statically imported, and `call`/`bind` produce a
validated reference. The renderer does not accept raw Typst expressions and does not call
`eval`.

These rules complement the generated particle callback implementation in
#link("../../assets/embedded/drawing/templates/impl/physics-edge-style.typ")[the GammaLoop-owned
implementation]; its public callback API lives in
#link("../../assets/embedded/drawing/templates/physics-edge-style.typ")[the same GammaLoop
template bundle].

== Pattern Style Dictionaries

Kurvst patterns are selected through Typst style dictionaries, not through a
special physics DOT field. A style dictionary may contain:

- `pattern`: `"wave"`, `"zigzag"`, `"coil"`, `"normal"`, or `"curve"`.
- `pattern-amplitude`
- `pattern-wavelength`
- `pattern-phase`
- `pattern-samples-per-period`
- `pattern-coil-longitudinal-scale`
- `pattern-accuracy`

Style dictionaries may also contain path-geometry keys. Linnest consumes these
keys while resolving graph edge layers and delegates the actual path operations
to Kurvst's path-in/path-out helpers:

- `offset`: normal offset for the half-edge path.
- `length`: maximum visible arc length for a centered parallel path.
- `ratio`: maximum visible fraction of the base edge length for a centered
  parallel path.
- `resolve-length`: how to combine `length` and `ratio`.
- `accuracy`: Kurbo fitting tolerance for the parallel path.
- `optimize`: whether Kurbo should optimize the fitted path.
- `offset-side: "label"`: choose the sign of `offset` so the path is on the
  same side as the edge label.

`draw` also accepts `edge-offset`, `edge-length`, `edge-ratio`,
`edge-resolve-length`, `edge-accuracy`, and `edge-optimize` as defaults for both
half edges. Derived paths are computed on the base edge geometry before patterns
and other decorations; node outsets then trim the shifted path, so it remains
shortened at node boundaries. When both length and ratio limits are set,
`resolve-length` decides which visible span to use. For paired edges, Linnest
computes that centered visible interval once on the full edge and then projects
it onto the source and sink halves, so a short layer can cross the source/sink
split instead of being shortened independently on each half. The reusable
operations live in Kurvst (`layer`, `center-outset`, `segments`, and `length`);
Linnest's `draw` layer only maps graph styles and label
positions onto those primitives.

Kurvst follows CeTZ's Rust boundary style: Typst sees explicit CBOR wire
dictionaries and Rust converts those dictionaries to `kurbo::BezPath` internally.
That keeps the public Typst API path-in/path-out without tying path objects to
Kurbo's serde representation.

GammaLoop-generated particle styles produce these dictionaries directly. Python callers express
the same values through typed edge and half-edge drawing styles, or reference a dictionary
exported by an explicit `TypstModule`. A DOT attribute is never interpreted as Typst source.

== Subgraph Drawing

Subgraph shading is a drawing feature. In Typst, callers can construct a
subgraph and pass it to `draw`:

```typ
#let east = subgraph.compass(layed-out, "e")
#draw(layed-out, subgraph: east)
```

Compass subgraphs use half-edge compass data from DOT ports such as `v:0:e` or
from Typst-built half-edge dictionaries such as `(node: a, compass: "e")`.

`draw` shades included half-edges with `subgraph-edge-style`. By default this
is an underlay, so the normal edge style remains visible on top.

== Practical Guidance

For GammaLoop-generated diagrams:

- Put physics data in canonical GammaLoop fields.
- Let GammaLoop generate `edge-style.typ`.
- Select the GammaLoop V1 template when rendering through Python or Clinnet directly.
- Put executable custom styles in an imported Typst module, not in DOT strings.

For manually edited drawing DOT:

- Use `display-label` or `label` as data labels; templates do not expand `{field}` placeholders
  or recognize `label-template`.
- Keep custom styles in the typed render configuration or an imported Typst module instead of a
  DOT string.
- Treat drawing-only DOT as a render artifact if those fields are not yet
  accepted by the GammaLoop physics parser.
