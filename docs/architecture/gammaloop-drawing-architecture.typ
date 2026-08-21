= GammaLoop Drawing Architecture

#quote(block: true)[
#strong[Status:] Current contract; audited against the implementation on 2026-08-21

The pipeline, ownership, DOT parser fields, callback data, path styles, and subgraph sections
describe the implemented drawing stack. The superseded `*-eval` and placeholder-interpolation
design is preserved separately in
#link("gammaloop-drawing-evaluated-fields-history.typ")[the evaluated-field proposal record].
]

This document describes the drawing path used by GammaLoop and the DOT syntax
that matters at each layer. The important boundary is:

- GammaLoop owns physics graph data and model-to-style policy.
- Linnet/Linnest owns DOT parsing, half-edge graph shape, layout, and drawing.
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
+ Clinnet's shared renderer compiles `figure.typ` for each DOT graph by invoking an external
   Typst 0.15 or newer executable. The `linnet` CLI and `linnet-py`'s `DotGraph.render` and
   `DotGraph.to_svg` use this same renderer; it is orchestration around Typst, not a native Rust
   drawing backend.
+ The generic rendering path materializes Clinnet's embedded figure/layout templates and the
   embedded Linnest and Kurvst package assets. The GammaLoop path instead supplies its generated
   figure-template bundle so model-specific `edge-style.typ` callbacks remain available. The
   figure template reads the DOT data through `sys.inputs.data-path` and forwards it to
   `layout.typ`.
+ `layout.typ` binds the extracted Linnest package and generated edge styles;
   `layout-core.typ` then calls `graph.parse`, `layout`, and `draw`.
+ `draw` calls the generated callbacks from `edge-style.typ`, draws edges and
   labels, then draws nodes last so nodes sit on top of edges.

The normal GammaLoop path does not need evaluated strings. Generated particle
styles are ordinary Typst dictionaries/functions in `edge-style.typ`, backed by
#link("../../crates/linnest/typst/src/physics-edge-style.typ")[`physics-edge-style.typ`].
A user can opt into evaluating the recognized `label`, `display-label`,
`source-style`, and `sink-style` string fields at drawing time.

== Components

/ `crates/gammalooprs`: Generates physics DOT and the model-specific `edge-style.typ`. It decides how
  a particle should look: photon wave, gluon coil, scalar dashed line, mass
  dependent thickness, labels, and related policy.

/ `crates/linnet`: Provides the half-edge graph data structure and DOT parser. Its parser turns
  invisible DOT nodes into dangling half-edges and preserves unconsumed
  attributes as statement dictionaries.

/ `crates/linnest`: Provides the Typst-facing wasm plugin and Typst wrappers. It parses DOT into
  graph bytes, lays out node and edge positions, exposes graph queries, and
  provides the CeTZ draw API.

/ `crates/kurvst`: Provides the Typst-facing curve wasm plugin. It splits and trims Bezier
  curves, builds Hobby curves through edge layout points, and generates wave,
  zigzag, and coil path patterns.

/ `crates/clinnet`: Provides the shared external-Typst renderer, embedded generic Clinnet,
  Linnest, and Kurvst assets, and the `linnet` CLI used to batch-render DOT files and assemble
  grid PDFs.

/ `crates/linnet-py`: Provides the Python graph API. Its `DotGraph` rendering methods delegate
  to Clinnet's shared renderer, and notebook SVG display delegates to the same path rather than
  implementing a separate renderer.

GammaLoop's particle policy is not part of Clinnet's generic embedded template. Python callers
that want GammaLoop particle colors, labels, and path decorations must pass the generated
GammaLoop `figure.typ` entry point with the rest of its template bundle kept in place.

== Data Ownership

GammaLoop DOT has two kinds of data.

Physics data is read by GammaLoop when a graph is parsed as a physics graph.
Examples are `particle`, `pdg`, `lmb_id`, `num`, `int_id`, and
`overall_factor`.

Drawing data is read by the Linnest drawing templates. Examples are
`display-label`, `label`, `source-style`, `sink-style`, half-edge
compass points, and subgraph selections by compass.

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

== Current Callback Precedence And Eval Mode

The generated model entry supplies the base source and sink styles. A normal `source-style` or
`sink-style` value is then read from the edge and its corresponding half-edge data. Dictionaries
are accepted directly. String values are ignored in `typst-fields: "plain"` mode and evaluated in
`typst-fields: "eval"` mode. There is no `source-style-eval` or `sink-style-eval` fallback.

The default edge label is selected from edge-data `display-label`/`label`, then top-level
`display-label`/`label`, then the generated particle-map label. The optional `show-momentum`,
`show-edge-index`, `show-half-edge-index`, and `show-particle` controls build an explicit metadata
label instead. There is no `label-template` fallback or placeholder expansion. These rules live in
#link("../../crates/linnest/typst/src/impl/physics-edge-style.typ")[the callback implementation],
while the public options live in
#link("../../crates/linnest/typst/src/physics-edge-style.typ")[the callback API].

The embedded figure still defaults `sys.inputs.typst-fields` to `"plain"`. Opt into executable
Typst strings only for a deliberately hand-authored rendering:

```bash
linnet draw graphs --input typst-fields=eval
```

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

GammaLoop-generated particle styles produce these dictionaries directly. A user
can also produce them through `source-style` or `sink-style` in eval mode:

```dot
a -> b [
  particle="a",
  source-style="(stroke: red + 1pt, pattern: \"wave\", pattern-amplitude: 0.14)",
  sink-style="(stroke: blue + 1pt, pattern: \"coil\", pattern-amplitude: 0.14)"
];
```

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
- Keep `typst-fields` at the default `plain`.
- Do not use the superseded `*-eval` field names.

For manually edited drawing DOT:

- Use `display-label` or `label` for labels; the current callback does not expand
  `{field}` placeholders or recognize `label-template`.
- Use `--input typst-fields=eval` only when a recognized render field contains
  deliberately executable Typst code. Do not use the superseded `*-eval` field names.
- Prefer drawing-only fields such as `source-style` and `sink-style` over
  changing physics fields such as `particle`, `pdg`, or `mass`.
- Treat drawing-only DOT as a render artifact if those fields are not yet
  accepted by the GammaLoop physics parser.
