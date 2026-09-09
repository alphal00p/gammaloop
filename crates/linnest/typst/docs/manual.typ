#import "@preview/tidy:0.4.3"
#import "../src/lib.typ": draw, graph, layout, layouts, subgraph
#import graph: (
  build, dot, edge, edge-data, edges, node, node-data, nodes, parse, sink,
  source, update-edge-data, update-node-data,
)


#set document(title: "Linnest Typst API")
#set page(margin: 22mm)
#set text(size: 10pt, font: "Roboto")

#align(center)[
  #set text(size: 3em, font: "Roboto")
  linnest
  #set text(size: 0.3em, font: "Roboto")

  _a typst interface to the `linnet` crate_

  `linnet` is a Rust crate that provides a graph data structure, graph algorithms,
  several layout algorithms, and a DOT parser built on `dot_parser`.



]

= Linnest Typst API

#outline(depth: 3)
#let crown = $star$

#let linnest-guide = [
  == Linnet

  The `linnet` crate wrapped by this package is built around a half-edge graph data structure.
  This means that instead of a graph being represented as a set of nodes and edges, it is represented as a set of half-edges $H$, and a set of vertices $V$.
  The graph structure is then encoded through two maps.
  The first map, $partial : H --> V$ maps each half-edge to its corresponding vertex.
  The preimage of any vertex $v$ is the set of half-edges that map to it, called the crown of $v$.
  The second map, $iota : H --> H$, is an involution that _glues_ half-edges together to form edges.
  If a half-edge is glued to itself, we call that an external half-edge. This means that linnet graphs are strictly more capable than normal edge and vertex graphs.

  #let g = build({
    node(<a>, label: [$v$])
    node(<b>)
    edge(source(<a>), <a-b>, sink(<b>), label: [e])
    edge(source(<a>), <in-a1>, label: [e])
    edge(source(<a>), <in-a2>, label: [e])
  })
  // #edges(g)
  #context if target() == "paged" {
    figure(draw(layout(g, g-center: 0.005, length-scale: .3)))
  }


  Native half-edges also make subgraphs more granular because they can be encoded
  as sets of half-edges. This directly supports vertex-induced subgraphs: the
  union of the crowns of a set of vertices.


  == Linnest

  Linnest is the Typst and WebAssembly interface to Linnet. It provides layout
  algorithms, DOT parsing, and selected graph algorithms without a separate
  runtime process.

  Graphs can be constructed in two ways: parse a DOT string with `parse`:
  ```typ
  #let g = parse("digraph { a -> b }")
  ```
  or build from edges and nodes with `build`, using a Fletcher-inspired syntax:
  ```typ
  #let g = build({
    node(<a>, label: [$v$])
    node(<b>)
    edge(source(<a>), <a-b>, sink(<b>), label: [e])
    edge(source(<a>), <in-a1>, label: [e])
    edge(source(<a>), <in-a2>, label: [e])
  })```

  In either case, `type(g)` is `dictionary`: graph values wrap an archived Linnet
  graph together with native Typst data. Rust owns topology, layout state,
  statement metadata, and internal opaque payload bytes. User data captured from
  Typst stays in Typst and is merged back into query records.

  The main use case is to place nodes and edges on a canvas with the `layout`
  function and render the result with `draw`.



  - `graph` for construction, parsing, inspection, directed cuts, joins, and graph
    algorithms.
  - `subgraph` for half-edge selection, annotations, and inspection.
  - `layout` for the separate layout pass.
  - `draw` for rendering a laid-out graph object with CeTZ.

  Domain-specific styles are ordinary Typst data and callbacks composed with
  `graph.style` and `draw`; Linnest does not reserve a domain or export a physics
  module.

  === Choose an import path

  Linnest and Kurvst are currently bundled source packages, not Typst Universe
  packages. A Clinnet run writes both package trees below `build/templates/`.
  From a custom template in that directory, import Linnest with:

  ```typ
  #import "crates/linnest/typst/src/lib.typ": draw, graph, layout, subgraph
  ```

  From this repository's `crates/linnest/typst/examples/` directory, the equivalent
  checkout-relative import is `../src/lib.typ`. Keep the package directory and its
  `linnest.wasm` file together when copying it elsewhere. The examples below use
  the checkout-relative form because they are also compiled as repository tests.

  === Minimal Build Example

  ```typ
  #import "../src/lib.typ": draw, graph, layout, subgraph
  #import graph: build, dot, edge, edges, node, nodes, parse, sink, source

  #let g = build({
    node(<a>)
    node(<c>)
    edge(
      source(<a>, compass: "e"),
      <a-c>,
      sink(<c>, compass: "w"),
      label: [a-c],
      statements: (
        color: "0055ff",
        source-color: "d72638",
        sink-color: "1b7f4c",
      ),
    )
  },
    name: "demo",
  )
  #let g = layout(g)
  #let east = subgraph.compass(g, "e")
  #let edge-records = edges(g, subgraph: east)
  #let dot-text = dot(g)
  #let edge-label(edge) = text(fill: rgb("#" + edge.color))[#edge.label]
  #let source-style(edge) = (stroke: rgb("#" + edge.source-color) + 0.5pt)
  #let sink-style(edge) = (stroke: rgb("#" + edge.sink-color) + 0.5pt)
  #context if target() == "paged" {
    draw(
      g,
      subgraph: east,
      edge-label: edge-label,
      edge-label-style: (anchor: "south"),
      source-style: source-style,
      sink-style: sink-style,
    )
  } else {
    [The downloadable PDF renders this CeTZ result. The HTML manual keeps the
    copyable source because Typst's experimental HTML target does not yet emit
    the drawing content.]
  }
  ```

  #import "../src/lib.typ": draw, graph, layout, subgraph
  #import graph: build, dot, edge, edges, node, nodes, parse, sink, source

  #let g = build(
    {
      node(<a>)
      node(<c>)
      edge(
        source(<a>, compass: "e"),
        <a-c>,
        sink(<c>, compass: "w"),
        label: [a-c],
        statements: (
          color: "0055ff",
          source-color: "d72638",
          sink-color: "1b7f4c",
        ),
      )
    },
    name: "demo",
  )
  #let g = layout(g)
  #let east = subgraph.compass(g, "e")
  #let edge-records = edges(g, subgraph: east)
  #let dot-text = dot(g)
  #let edge-label(edge) = text(fill: rgb("#" + edge.color))[#edge.label]
  #let source-style(edge) = (stroke: rgb("#" + edge.source-color) + 0.5pt)
  #let sink-style(edge) = (stroke: rgb("#" + edge.sink-color) + 0.5pt)
  #draw(
    g,
    subgraph: east,
    edge-label: edge-label,
    edge-label-style: (anchor: "south"),
    source-style: source-style,
    sink-style: sink-style,
  )

  === Graph Objects

  Graph values combine archived Linnet topology with native Typst data. The
  #link("reference/typst/graph/#graph-object-api")[focused graph reference]
  documents their constructors, transforms, queries, and update operations.
]

#let graph-objects = [
  === Graph object API

  Graph objects are Typst dictionaries wrapping archived Rust graph bytes plus
  native Typst data arrays for graph, node, edge, source, and sink data. The Rust
  graph may also carry an internal opaque payload, but that payload is used only
  by the Typst wrapper and is not exposed in public records. Build or parse graph
  objects with `graph`, transform graph objects with `layout`, and pass objects
  back to `graph` or `subgraph` for inspection. Subgraphs are also dictionaries:
  they combine an archived half-edge selection with a topology signature and
  Typst-only annotations, rather than exposing raw bytes as public inputs.

  - `graph.parse(input)` parses one or more DOT digraphs and returns an array of
    graph objects. Its `eval-graph-fields`, `eval-node-fields`,
    `eval-edge-fields`, `eval-source-fields`, and `eval-sink-fields` arguments
    are convenience arguments for `graph.eval-fields`.
  - `graph.build(..)` constructs one graph object from a stream of node and edge
    items.
  - `graph.map(graph, ..)` maps graph, node, edge, source, and sink records to
    new native data or structural patches without changing topology. `node` and
    `edge` accept a callback or a dictionary keyed by current record names.
    Each dictionary value is a patch, a callback receiving the full record, or
    `none`. Edge patches may set `spring-length` directly; it is validated and
    stored as the layout statement. Unlisted records remain unchanged; unknown
    names are errors.
  - `graph.cut(graph, left: left, right: right, boundary: patch)` opens a
    weighted directed cut into a new graph, preserving Typst data and recording
    origins. `graph.boundaries(view)` queries its current boundary endpoints.
  - `graph.node-data(graph, <name>)` and `graph.edge-data(graph, <name>)` return
    one named node or edge data.
  - `graph.update-node-data(graph, <name>, data)` and
    `graph.update-edge-data(graph, <name>, data)` update one named node or edge
    data. Direct replacements and callbacks run in Typst with `(data, record)`.
  - `graph.eval-fields(graph, ..)` evaluates selected record fields into data
    entries. It works on parsed and built graph objects.
  - `node(..)` returns a node item.
  - `source(..)` and `sink(..)` return half-edge endpoints.
  - `edge(..)` returns an edge item built from source/sink endpoints and an
    optional edge name.
  - `layout(graph, ..)` runs layout as an explicit
    second step. Its settings are named parameters so calls stay descriptive and
    Tidy can document each field.
  - `draw(graph, ..)` draws a laid-out graph object with CeTZ.
  - `dot(graph)` returns a DOT string for inspection or export.
]

#let graph-concepts = [
  === Graph Specs

  The Typst construction API is half-edge first: create node items, create source
  and sink half-edge endpoints, then pass edge items to `graph.build`.
  `graph.build` accepts both comma-separated items and ordinary Typst code
  blocks. Typst labels such as `<a>`, `<h1>`, and `<e1>` are API names; node and
  edge names are emitted back from `nodes(g)` and `edges(g)` as Typst labels.
  They are resolved before the wire format is sent to the Rust plugin. Numeric `id`
  arguments choose graph indexes or ordering: on nodes, `id` fixes the resulting
  node index and must be unique and in bounds. On nodes, edges, sources, and
  sinks, extra named arguments are captured as opaque Typst data fields:
  `edge(source(<a>, style: (stroke: red + 0.7pt)), sink(<b>), kind: "link")`
  stores `(style: ..)` in the source data and `(kind: "link")` in the edge
  data. These user data values are not sent through the Rust plugin boundary.
  Typst sends an internal opaque payload with correlation data, asks Rust to
  resolve the graph indexes, then stores user data in arrays at those resolved
  graph, node, edge, and half-edge ids. `graph.info`, `graph.nodes`, and
  `graph.edges` merge those native values into the returned records. A captured
  `label` data field on nodes
  and edges is display content used by the default drawing style; use
  `statements: (label: "...")` when a flat metadata label string is needed.
  Statements are flat metadata used by DOT; they cannot nest. Values are scalar
  strings/numbers/booleans. Use data fields for structured Typst data or
  content.

  `graph.build` and `graph.parse` also accept `default-node-data`,
  `default-edge-data`, `default-source-data`, and `default-sink-data`.
  These defaults are merged into the corresponding data; captured data
  fields on nodes, edges, sources, and sinks override the defaults. Drawing
  callbacks decide how those fields affect the result. For example, endpoint
  styles can read `data.style` from the source and sink half-edge records:

  ```typ
  #let g = build({
    node(<a>, label: [a])
    node(<b>, label: [b])
    edge(
      source(<a>, style: (stroke: red + 0.7pt)),
      <e>,
      sink(<b>, style: (stroke: blue + 0.7pt)),
      label: [connection],
      kind: "dependency",
    )
  },
    default-edge-data: (kind: "propagator"),
  )
  #let endpoint-style(edge, side) = {
    let endpoint = edge.at(side + "-half-edge", default: none)
    if endpoint == none { return (:) }
    endpoint.data.at("style", default: (:))
  }
  #draw(
    layout(g),
    source-style: edge => endpoint-style(edge, "source"),
    sink-style: edge => endpoint-style(edge, "sink"),
    edge-label: edge => edge.at("label", default: none),
  )
  ```

  When parsing DOT, the same native data arrays can be filled from selected
  string fields. Rust parses DOT into topology and statement metadata; Typst
  applies defaults and evaluates selected fields afterward. The selected fields
  are evaluated with the record's merged `fields` dictionary in scope. Since node
  names are labels, use `#str(name)` when a name should become visible text:

  ```typ
  #let g = parse(
    "digraph g { a [label=\"A\"]; a -> b [label=\"$p$\", source=\"out\", sink=\"in\"] }",
    eval-node-fields: ("label",),
    eval-edge-fields: ("label",),
    eval-source-fields: ("statement",),
    eval-sink-fields: ("statement",),
  ).first()
  #nodes(g).first().data.label
  #edges(g).first().data.label
  ```

  The same transform can run after construction. This is useful for global edge
  statements that should apply to every edge while still seeing local edge
  fields:

  ```typ
  #let g = build({
    node(<a>)
    node(<b>)
    edge(source(<a>), sink(<b>), statements: (mom: "p"))
  }, default-edge-statements: (display-label: "$#mom$"))
  #let g = graph.eval-fields(g, eval-edge-fields: ("display-label",))
  #edges(g).first().data.at("display-label")
  ```

  ```typ
  #let g = build({
    node(<a>)
    node(<c>)
    edge(source(<a>), <e1>, sink(<c>))
  })
  ```

  Named nodes and edges can be updated after construction without scanning in the
  caller. The update replaces the native data, or it can be a callback receiving
  `(data, record)`:

  ```typ
  #let g = build({
    node(<a>)
    node(<c>)
    edge(source(<a>), <e1>, sink(<c>))
  })
  #let g = update-node-data(g, <a>, (label: [A]))
  #let g = update-edge-data(g, <e1>, (data, edge) => (
    label: [$p$],
    source: edge.source.node,
  ))
  #node-data(g, <a>).label
  #edge-data(g, <e1>).label
  ```

  `source(..)` and `sink(..)` accept a node reference plus optional `name`, `id`,
  `statement`, and `compass`. `name` is a Typst label name for the half-edge;
  `id` is numeric:

  ```typ
  edge(source(<a>, name: <h1>, id: 0), <e1>, sink(<c>, id: 2), label: [a-c])
  ```

  One half-edge creates an external edge. The side is determined by the
  constructor, so there is no public `flow` argument:

  ```typ
  edge(<incoming>, sink(<a>))
  edge(source(<c>), <outgoing>)
  ```

  `graph.edge(..., spring-length: 1.5)` multiplies that edge's existing preferred
  spring length by `1.5` in force/anneal layout. The factor must be positive,
  finite, and dimensionless; it applies to each node-to-control-point spring,
  including the existing doubled rest length for dangling edges. The argument
  defaults to `none`, preserving local or default `spring-length` statements;
  a missing statement means `1`. An explicit argument overrides those statements.
  This is not a guaranteed rendered length, a rendering option, or a change to
  global layout `spring.length`; tree/dot/stable-layered layout is unaffected.

  `graph.build` does not interpolate statement strings on the Rust side. Use
  `graph.eval-fields` or `graph.map` when a default statement should turn into
  Typst content or structured data. The evaluation scope includes the
  record's merged fields, so default edge statements can still refer to local edge
  fields:

  ```typ
  #let g = build({
    node(<a>)
    node(<c>)
    edge(
      source(<a>),
      <a-c>,
      sink(<c>),
      label: [a-c],
      statements: (color: "0055ff", label: "a-c"),
    )
  },
    default-edge-statements: (
      color: "000000",
      display-label: "$#label$",
    ),
  )
  #let g = graph.eval-fields(g, eval-edge-fields: ("display-label",))
  ```
]

#let placement-concepts = [
  === Placements

  `graph.pos` creates a first-class placement. The default `mode: "pin"` turns a
  coordinate into a fixed layout constraint and, for x/y, a drawable position.
  Use `mode: "start"` when the coordinate should only seed the layout, or use
  `graph.pin(value)` / `graph.start(value)` to choose a mode for one axis:

  ```typ
  #let g = build({
    node(<a>, pos: graph.pos(x: -2, y: 0))
    node(<c>, pos: graph.pos(ref: <a>, dx: 4, dy: 0))
    edge(source(<a>), <a-c>, sink(<c>), pos: graph.pos(x: 0, y: 1.2))
  })
  ```

  Auxiliary `z` is force-layout depth in the same numeric layout units as x/y,
  not a Typst length or a rendered coordinate. Both `graph.build` and `graph.map`
  accept `graph.pos(z: graph.pin(2))` for a hard raw-depth pin and
  `graph.pos(z: graph.start(2))` for a movable initial depth. A bare
  `graph.pos(z: 2)` follows `mode`, defaulting to `"pin"`:

  ```typ
  #let g = graph.map(g,
    node: node => (pos: graph.pos(z: graph.pin(2))),
    edge: edge => (pos: graph.pos(z: graph.start(-1))),
  )
  ```

  A z-only placement leaves XY coordinates, constraints, and partial-position
  flags alone, including automatic edge midpoint seeds. Omitting z preserves
  existing depth metadata during unrelated XY patches. `ref`, `dx`, and `dy`
  affect XY only; z groups, relative depth references, `dz`, and non-finite z
  values are not supported.

  Node and edge `statements` expose `"pos-z"` as finite numeric text and
  `"pos-z-mode"` as `"pin"` or `"start"`. These describe the supplied raw pin
  or seed; output `pos` remains a two-coordinate x/y record. The force solver
  multiplies raw depth by a global scale that flattens to zero: a hard raw pin
  remains fixed even when its effective depth is zero. Other layout backends do
  not use z. Auxiliary depth does not change draw order or provide an over/under
  rendering layer.

  `graph.group` links one XY coordinate across several nodes or edge control
  points. A `side` of "+" keeps the coordinate positive, and "-" keeps it negative.
  Use `start` to seed the shared coordinate without fixing it. Matching node and
  edge groups are one degree of freedom during force and annealing layouts, so
  repulsion and springs act on their combined force rather than being reconciled
  afterward. If several members supply starts for one group, their mean initializes
  the shared coordinate.
  GammaLoop external-edge columns use this to keep incoming and outgoing external
  legs on opposite sides while pairing rows by a shared `y` group:

  ```typ
  #let edge-items = (
    edge(
      source(<right-ext>),
      sink(<center>),
      pos: graph.pos(x: graph.group("right", side: "+", start: 4), y: graph.group("edgee0")),
    ),
    edge(
      source(<center>),
      sink(<left-ext>),
      pos: graph.pos(x: graph.group("left", side: "-"), y: graph.group("edgee0")),
    ),
  )
  ```

  Raw DOT input uses the same placement model through the `pos` attribute. The
  standard Graphviz subset is preserved: `pos="x,y"` is a starting coordinate and
  `pos="x,y!"` is pinned. Linnest extends the value with explicit id references
  and axis entries. In axis entries, `!` belongs to that axis: numeric entries
  without `!` are starting coordinates, numeric entries with `!` are fixed
  constraints, and grouped entries must use `!` because groups are constraints.

  ```dot
  digraph {
    a [id=0 pos="0,0!"]
    b [id=1 pos="ref(node:0)+4,0!"]
    a -> b [id=0 pos="ref(node:1)+0,1!"]
    b -> c [id=1 pos="ref(edge:0)+1,0!"]
    c [id=2 pos="x:2!"]
    d [id=3 pos="x:2!,y:1"]
    ext -> a [id=2 pos="x:@-left!,y:@edge0!"]
  }
  ```

  The DOT grammar is intentionally id-based: `ref(node:0)` uses a node `id` and
  `ref(edge:0)` uses an edge `id`. Bare names and implicit edge order are not
  placement references. The `@` syntax denotes grouped coordinate constraints;
  `@+name` and `@-name` keep the grouped coordinate on the positive or negative
  side respectively.

  ```text
  pos="x,y"                   // start x and y
  pos="x,y!"                  // pin x and y
  pos="x:<coord>"             // start numeric x
  pos="x:<coord>!"            // pin x
  pos="y:<coord>!"            // pin y
  pos="x:<coord>!,y:<coord>"  // pin x, start numeric y
  pos="x:<coord>,y:<coord>!"  // start numeric x, pin y
  ```
]

#let drawing-concepts = [
  === Drawing

  Draw styling is Typst-native. Pass dictionaries or callbacks to `draw`; edge
  callbacks receive the merged `scope`, edge statements, source/sink half-edge
  records, and the edge index:

  ```typ
  #let edge-label(edge) = text(fill: rgb("#" + edge.color))[#edge.display-label]
  #let source-style(edge) = (stroke: red + 0.5pt)
  #let sink-style(edge) = (stroke: blue + 0.5pt)
  #draw(layout(g), edge-label: edge-label, source-style: source-style, sink-style: sink-style)
  ```

  Edge labels are drawn with CeTZ content at the layout label position. Use
  `edge-label-style: (anchor: "south")`, or an edge-data callback returning a
  style dictionary, to choose which point of the label is anchored there.

  Marks can follow graph orientation as a first-class draw option. Put the same
  mark layer on both halves and set `mark-orientation: "edge"`; default-oriented
  edges mark the source half, reversed edges mark the sink half with the marker
  flipped, and undirected edges suppress the mark.
  Use `mark-position: "center"` or a numeric ratio from `0` to `1`, then
  `mark-shift` for a signed arc-length adjustment along the derived path. This
  positioning is independent of whether the layer is straight or patterned.
  `"center-if-dangling"` centers a dangling mark but keeps an orientation-selected
  paired mark at the source/sink split point. Heads span a chord between two
  points on the full carrier, not a tangent. Triangle and straight heads place
  their geometric tip and the center of their back on the curve; other marks
  use their declared tip/base anchors. Centered heads straddle the requested
  arc position, so the chord midpoint can lie off the curve. Near endpoints the
  sampling interval moves inward, and short carriers compress the head to fit.
  Longitudinal fitting preserves the head's width and stroke thickness.
  Interior marks overlay the unshortened curve. With shortening enabled, filled
  end heads meet the shaft at the inward contact; straight open heads keep the
  shaft running to their tip.

  ```typ
  #let oriented-arrow = (
    stroke: black + 0.7pt,
    mark: (end: (symbol: ">", fill: black, anchor: "center", shorten-to: auto), scale: 0.75),
    mark-position: "center-if-dangling",
    mark-orientation: "edge",
  )
  #draw(layout(g), source-style: oriented-arrow, sink-style: oriented-arrow)
  ```
]

#let domain-style-concepts = [
  === Domain Styles And Templates

  Linnest's drawing layer is domain-neutral. A caller can map arbitrary node,
  edge, and half-edge data to reusable `draw` callbacks without changing topology
  or adding a domain-specific core option. The following map uses an
  application-defined `kind` field; the names and styles have no meaning to
  Linnest itself:

  ```typ
  #import "../src/lib.typ": draw, graph, layout
  #import graph: build, edge, node, sink, source

  #let kinds = (
    dependency: (stroke: rgb("#315f9f") + 0.8pt),
    event: (
      stroke: rgb("#a24b36") + 0.8pt,
      mark: (end: ">"),
      mark-position: "center",
      mark-orientation: "edge",
    ),
  )
  #let edge-style(edge) = kinds.at(edge.fields.at("kind", default: "dependency"))
  #let g = build({
    node(<queued>, label: [queued])
    node(<running>, label: [running])
    node(<done>, label: [done])
    edge(source(<queued>), <starts>, sink(<running>), kind: "event", label: [starts])
    edge(source(<running>), <needs>, sink(<done>), kind: "dependency", label: [needs])
  })
  #let g = graph.style(g, edge-label: edge => edge.label)
  #draw(layout(g, layout-algo: "stable-layered"), edge-style: edge-style)
  ```

  An edge can sparsely patch the `edge-style` passed to `draw` with
  `style: (...)`, compute the patch with a callback, use `auto` to delegate, or
  use `none` to hide its paint while
  retaining topology. Local style callbacks receive the complete post-layout
  drawing record. Custom templates can package the same policy behind their
  mandatory `render(config)` export; the selected application owns the meaning
  and schema of any extra options.

  Set `edge-offset` on `draw`, or `offset` in a `source-style`/`sink-style`
  dictionary, to draw a fitted parallel path. Style callbacks may also return an
  array of dictionaries; the layers are drawn in order on the same graph edge, so
  parallel strokes do not require duplicate graph edges.

  ```typ
  #let base-style(edge) = (
    stroke: (paint: gray, thickness: 0.7pt, cap: "round"),
  )
  #let offset-style(edge) = (
    offset: 0.18,
    length: 1.4,
    ratio: 0.5,
    resolve-length: "min",
    shift: 0.2,
    stroke: (paint: rgb("#2f6f4e"), thickness: 1pt, cap: "round"),
    label: [active],
    label-side: "left",
    label-gap: 0.12,
  )
  #let source-style(edge) = (base-style(edge), offset-style(edge))
  #let sink-style(edge) = (base-style(edge), offset-style(edge) + (mark: (end: ">")))
  ```

  The parallel path is applied to the base edge geometry before patterns and other
  decorations; node outsets then trim the shifted path, so shifted paths still
  start and end outside fitted node circles. Add `edge-length` or `length` to
  center-trim the shifted path to a fixed arc length, and add `edge-ratio` or
  `ratio` to cap it by a fraction of the full offset path length. `edge-resolve-length`
  / `resolve-length` decides how to combine both limits: `"min"`/`"shorter"` (default), `"max"`/`"longer"`,
  `"length"`/`"fixed"`, `"ratio"`/`"relative"`, `"none"`/`"full"`, or a function
  receiving `(offset-path-length, length, ratio)`.
  Lengths and shortening distances are measured after offsetting, since a curved
  parallel path need not have the same length as its centerline.
  For a finite layer, `shift` is an arc-length displacement along the complete
  offset path: positive values move toward the path end and values that would
  cross an endpoint are clamped. A layer can attach `label` content to its own
  path. `label-shift` (default `0`) moves the label's reference point by signed
  arc length from that derived path's midpoint, clamped to its endpoints:
  `clamp(path-length / 2 + label-shift, 0, path-length)`. Positive values move
  toward the path end, without moving or trimming the layer itself.
  `label-side` chooses `"left"` or `"right"` relative to the local path direction
  at the shifted point; `auto` follows the side selected by ordinary edge-label
  layout.
  With `label-style.anchor` omitted or set to `auto` (also `"auto"`), the label
  is centered and moved until its entire CeTZ box clears the local tangent line
  by `label-gap`. This measures the box, including text bounds, wrapping,
  padding and rotation, not the nearest point of the finite curved shaft.
  An explicit anchor, including `"center"`, instead sits at the shifted reference
  point plus `label-gap` along the chosen normal, without any box-clearance
  correction. Negative gaps are clamped to zero in both modes.
  Other `label-style` fields are forwarded to CeTZ content drawing.
  Automatic clearance is a local tangent-line heuristic, not collision
  avoidance; explicit anchors intentionally allow the box to cross that line.
  In `examples/map-style.typ`, `momentum-label-anchor` selects this behavior:
  omit it or use `auto` for box clearance, or choose e.g. `"east"` for direct
  anchor placement. `momentum-label-gap` supplies the gap in graph units.
  The momentum label uses a separate full, unpainted offset-path layer, so
  `momentum-label-shift` is independent of `momentum-arrow-length`, even when
  arrow and label shifts are equal. It defaults to the requested
  `momentum-arrow-shift`, not the arrow's clamped center. Near endpoints, the
  label point can therefore travel farther than the finite arrow's center;
  changing arrow length never moves the label.

  Import `momentum as mom` from `examples/map-style.typ` for compact, sparse
  options. `mom(side: "right", length: .7, shift: -.4,
  label: (gap: .2, shift: -.75, anchor: "south-west"))` expands to the existing
  `momentum-arrow-*` and `momentum-label-*` fields. Arrow options are `side`,
  `offset`, `length`, and `shift`; label options are `gap`, `shift`, and `anchor`.
  Only supplied options are emitted, so `mom(label: (gap: .2))` preserves an
  inherited arrow offset or label shift. Combine the result with other edge
  patches using dictionary addition, or spread it into `edge(..)`.

  The attached label replaces the ordinary painted edge
  label, while that ordinary label may still supply the pre-layout size used by
  label layout and side selection; attached labels do not add a second collision
  constraint.
  When compatible source and sink layers request a centered mark, Linnest emits
  that mark once on the complete derived layer. The mark therefore follows
  `shift` instead of being centered separately on either half-edge.
  Set `offset-side: "label"` on an offset layer to choose the sign of `offset`
  so the layer is drawn on the same side of the curve as the edge label.

  Paired edges are Kurvst paths split at their edge layout point. Set
  `edge-split-gap` on `draw` to open a centered arc-length gap there, or set
  `split-gap` on an individual source/sink style layer to override the global
  value for only that layer. Wave and coil phases continue across the hidden
  span. `edge-halves` and `to-cetz-edge-halves` accept the same `split-gap`
  control and report the effective gap if either half is too short. This cuts at
  `edge.pos`; it does not detect arbitrary crossings.

  For a known centerline crossing, put `crossing-under: 5` (using the target's
  integer edge id) or `crossing-under: <bridge>` (using its Typst edge name), and
  optionally `crossing-gap` (default `0.55`) on the style layer that should be
  interrupted. Targets are resolved independently of edge order.
  Kurvst locates proper intersections and trims the current path by arc length;
  wave and coil phases continue across every hidden span. This works for
  dangling layers and for paired layers whose source and sink have one continuous
  style. Both paired half styles must name the same target and gap. A cut layer's
  mark is placed once on the original uncut carrier, while only its painted path
  is split; this keeps the mark present without duplicating it on every fragment.
  A cut layer cannot participate in a subgraph underlay. Self, unknown, and invisible references are reported as
  errors. A valid edge pair with no proper interior intersection, including one
  that only shares an endpoint, is left unchanged.

  ```typ
  #let crossed-style(edge) = if edge.eid == 7 {
    (stroke: black, pattern: "coil", crossing-under: 5, crossing-gap: 0.55)
  } else {
    (stroke: black)
  }
  ```

  For dangling edges, `edge-dangling-tangent: "horizontal"` or `"vertical"`
  constrains the tangent at the free edge position while retaining the edge's
  `bend`; the per-layer override is `dangling-tangent`. `auto` keeps the ordinary
  bent route, and `anchor-control-distance` controls the constrained handle
  length. `bend` affects dangling edges only. A paired edge instead follows its
  source node, `edge.pos`, and sink node, with `edge-omega` controlling its Hobby
  curve.

  ```typ
  #draw(
    g,
    edge-split-gap: 0.28,
    edge-dangling-tangent: "horizontal",
    source-style: edge => (stroke: black, split-gap: 0.4),
    sink-style: edge => (stroke: black, split-gap: 0.4),
  )
  ```

  Data defaults are also the global styling hook for all sources, sinks,
  nodes, and edges. More specific data can be added with captured named arguments
  on node, edge, source, and sink items, or by running
  `graph.map`/`graph.eval-fields` after construction. This keeps evaluated Typst
  values in the native data channel instead of adding renderer-specific eval
  fields to the Rust topology spec:

  ```typ
  #let g = build({
    node(<a>)
    node(<c>)
    edge(
      source(<a>),
      <a-c>,
      sink(<c>),
      label: [a-c],
      kind: "highlight",
    )
  },
    name: "demo",
    default-source-data: (style: (stroke: red + 0.5pt)),
    default-sink-data: (style: (stroke: blue + 0.5pt)),
  )
  #let endpoint-style(edge, side) = {
    let endpoint = edge.at(side + "-half-edge", default: none)
    if endpoint == none { return (:) }
    endpoint.data.at("style", default: (:))
  }
  #draw(
    layout(g),
    source-style: edge => endpoint-style(edge, "source"),
    sink-style: edge => endpoint-style(edge, "sink"),
  )
  ```

  `graph.build` also accepts comma-separated items:

  ```typ
  #let g = build(
    node(<a>),
    node(<b>),
    edge(
      source(<a>, compass: "e"),
      <ab>,
      sink(<b>, compass: "w"),
      label: [ab],
      statements: (color: "0055ff", label: "ab"),
    ),
    name: "demo",
    statements: (full_num: "x + y"),
    default-node-statements: (shape: "circle"),
    default-edge-statements: (
      color: "000000",
      display-label: "$#label$",
    ),
  )
  ```
]

#let graph-query-concepts = [
  === Graph Queries

  `graph.info(g)` returns graph metadata. `nodes(g)` returns node records,
  and `edges(g)` returns edge records. Node and edge record `name` values are
  Typst labels when present. Pass `subgraph: sg` to filter nodes
  or edges by a compatible subgraph object. Nodes must be incident to a selected
  half-edge; an edge is included if either half is selected, but its record keeps
  both available endpoints. These queries do not extract or open topology.

  After `graph.cut`, node, edge, and half-edge records carry `origin` relative
  to the immediate input graph. Boundary nodes and dangling cut edges also carry
  `boundary`. Use `graph.boundaries(view)` for a uniform endpoint query, including
  current positions after placement or layout. Auxiliary nodes remain visible in
  `graph.nodes(view)`; filter with `node.boundary == none` to exclude them.

  `graph.join(left, right, key: "statement")` joins matching dangling half edges.
  The key is read from half-edge statements or numeric ids and can be
  `"statement"`, `"compass"`, or `"id"`. Its existing Typst-sidecar data loss is
  not fixed by the cut API; do not use it as a data-preserving inverse of `cut`.

  `graph.cycles(g)` returns subgraph objects for a cycle basis.
  `graph.forests(g)` returns subgraph objects for spanning forests. Both return
  arrays of the same topology-bound dictionaries accepted by `subgraph`, graph
  queries, `layout`, and `draw`, not arrays of raw selection bytes.
]

#let layout-concepts = [
  === Layout Model

  `layout` starts from a traversal-tree placement, then optimizes the positions of
  graph nodes and edge control points. The initial tree spacing is

  $ L = lambda sqrt((W H) / max(n, 1)) $,

  with horizontal spacing $tau_x L$ and vertical spacing $tau_y L$. Here
  $lambda$ is `length-scale`, $W$ is `viewport-w`, $H$ is `viewport-h`, $tau_x$
  is `tree-dx`, and $tau_y$ is `tree-dy`. These fields set the geometry scale for
  both layout modes.

  For deterministic, non-iterative placement, use `layout-algo: "tree"` or
  `layout-algo: "dot"`. `"tree"` places a traversal forest by levels. `"dot"`
  uses a directed layered placement for acyclic inputs, falling back to tree
  placement when the selected edges contain a cycle. Both modes also work on a
  subgraph:

  ```typ
  #let g = parse("digraph partial { a -> b; b -> c; c -> d; d -> a }").at(0)
  #let tree = graph.forests(g).at(0)
  #let g = layout(g, layout-algo: "tree", subgraph: tree)
  #draw(g)
  ```

  Only nodes touched by the selected subgraph are placed. Edge control points are
  then resolved from the current node positions, so edges outside the selected
  subgraph are drawn as straight lines unless their own position constraints say
  otherwise.

  `layout-algo: "force"` and `layout-algo: "anneal"` also accept `subgraph`.
  For these iterative modes, nodes and edge control points outside the selected
  subgraph stay fixed and act as boundary points while the selected subgraph is
  optimized.

  Set `layout-nodes: "fixed"` to keep every node at its current `pos` for this
  layout pass and move only edge control points. With `subgraph`, only edges in
  the selected subgraph are moved; all other edge control points keep their
  current positions. The fixed-node policy is temporary; the returned graph stores
  the resulting coordinates as `pos`, but it does not turn them into persistent
  `pin` constraints. This is useful after one node-placement pass when a later pass
  should route or relax selected edges without disturbing the node layout:

  ```typ
  #let g = parse("digraph partial { a [pos=\"0,0\"]; b [pos=\"4,0\"]; c [pos=\"8,0\"]; a -> b; b -> c }").at(0)
  #let first = subgraph.bits(g, (true, true, false, false))
  #let g = layout(g, layout-algo: "tree", layout-nodes: "fixed", subgraph: first)
  #draw(g)
  ```

  Related controls can be supplied as shallow semantic records. These records
  only name graph-layout concepts; they do not select a drawing domain or style:

  ```typ
  #let common = layouts.options(
    spring: (strength: 8, length: 0.5),
    repulsion: (edge-node: 0.05, dangling: 2),
    constraints: (side-strength: 5),
    labels: (steps: 0),
    solver: (algorithm: "force", steps: 50),
  )
  #let spacious = layouts.options(
    base: common,
    spring: (length: 0.7),
    repulsion: (dangling-centroid: 3),
  )
  #let g = layout(g, ..spacious)
  ```

  `layouts.options` performs a sparse update: the second spring record above
  changes `length` without losing the inherited `strength`. Grouped fields
  override their corresponding flat arguments. Flat arguments remain useful for
  numerical experiments and controls that do not have a semantic grouping.

  The shared spring/charge model uses the following coefficients:

  $ c_("vv") = beta L^3 $
  $ c_("ev") = beta gamma_("ev") L^3 $
  $ c_("ee") = beta gamma_("ee") L^3 $
  $ c_("center") = beta g_("center") $
  $ c_("dangling") = beta gamma_("dangling") L^3 $
  $ c_("dangling-centroid") = beta gamma_("dangling-centroid") L^3 $

  The Typst parameter names are `beta` for $beta$, `gamma-ev` for
  $gamma_("ev")$, `gamma-ee` for $gamma_("ee")$, `g-center` for $g_("center")$,
  `gamma-dangling` for $gamma_("dangling")$, and
  `gamma-dangling-centroid` for $gamma_("dangling-centroid")$. The spring
  stiffness $k$ is `k-spring`, and the softening constant $epsilon$ is `eps`.

  In `layout-algo: "anneal"`, linnest minimizes an energy:

  $ E =
  sum_(i < j) 1/2 c_("vv") / (d(v_i, v_j) + epsilon)
  + sum_(i, e) c_("ev") / (d(v_i, e) + epsilon)
  + sum_((v, e) " incident") 1/2 k (ell_e - d(v, e))^2
  \ + sum_("local edge pairs") 1/2 c_("ee") / (d(e_i, e_j) + epsilon)
  + sum_("dangling pairs") 1/2 c_("dangling") / (d(e_i, e_j) + epsilon)
  + sum_("dangling" e) c_("dangling-centroid") / (d(e, overline(v)) + epsilon)
  + sum_(i) 1/2 c_("center") d(v_i, 0)^2
  + p_("cross") N_("cross") $.

  Here $overline(v)$ is the mean node position. The centroid term pushes every
  dangling endpoint away from that mean; its equal-and-opposite reaction is
  shared over the nodes so it introduces no net force. $p_("cross")$ is
  `crossing-penalty` and $N_("cross")$ is the number of detected edge crossings.
  The quadratic center term pulls nodes toward the origin at every radius.
  `temp`, `step`, `seed`, `steps`, `epochs`, `cool`,
  `accept-floor`, `step-shrink`, and `incremental-energy` belong to this
  simulated annealing mode. `crossing-penalty` is also anneal-only; force mode
  does not currently add a crossing force.

  In `layout-algo: "force"`, linnest applies the direct forces corresponding to
  the same vertex-vertex, edge-vertex, incidence spring, local edge-edge,
  dangling-edge, dangling-centroid, and center terms. `step` is the integration
  step, `delta` clamps per-step movement, `steps` and `epochs` set the iteration
  budget, and `cool` shrinks the step after each epoch. `early-tol` can stop the
  run when movement is small only after the effective depth scale reaches zero;
  small movement or a cooled-to-zero step cannot skip flattening.

  The force-only `depth-scale` (default `1.0`) and `flattening-end` (default `0.5`)
  control auxiliary depth. `depth-scale` must be finite and non-negative;
  `flattening-end` must be finite and in `0..=1`. They can be supplied as flat
  arguments or in `solver: (depth-scale: 1.0, flattening-end: 0.5)`; grouped
  values override flat ones. The old depth-spring controls are not aliases.

  Effective depth is `scale * raw-z`. The integrator starts with auxiliary raw
  depths to break overlaps and smoothly reduces `scale` from `depth-scale` to
  zero over the initial `flattening-end` fraction of the total `steps` × `epochs`
  iteration budget. For normalized progress `u` from 0 to 1 over that interval,
  the smoothstep scale is `depth-scale * (1 - u)^2 * (1 + 2 * u)`. The remaining
  iterations relax projected overlaps with effective depth exactly zero, so
  repulsion cannot be satisfied by invisible separation. `flattening-end: 0`
  or `depth-scale: 0` starts in 2D; `flattening-end: 1` still evaluates the final
  iteration on the exact plane.

  There is a single original iteration budget and cooling schedule: flattening
  neither adds a second phase nor restarts `step`. Unpinned raw depths stay
  bounded by the initial spread and supplied depth magnitudes; hard raw-depth
  pins remain fixed even after their effective depth has flattened to zero.
  XY pins, shared coordinates, and fixed subgraph boundaries apply throughout;
  label layout runs only after this single graph pass. Output positions remain
  2D, and auxiliary z does not alter draw order. Flattening is exact, while force
  convergence remains limited by the iteration budget and movement tolerance.

  `directional-force` is applied in both modes as an extra bias derived from
  pin/port direction constraints.

  After either graph layout mode, labels are relaxed separately. If $L_l$ is the
  label target distance and $q_l$ is the label repulsion strength, then
  $L_l = alpha_l L$ and $q_l = beta_l L^2$. The Typst names are
  `label-length-scale` for $alpha_l$ and `label-charge` for $beta_l$.
  `label-spring` is the spring constant pulling each label toward its target.
  `label-layout: "normal"` uses a perpendicular offset target. With
  `label-layout: "dangling-tangent"`, paired edges still use that perpendicular
  target, but dangling half-edge labels are offset along the edge direction away
  from the attached node. With `label-layout: "fixed-length"`, the label remains
  at distance $L_l$ from the edge point and only rotates around it under repulsive
  forces. `label-steps`, `label-step`, `label-early-tol`, and
  `label-max-delta-scale` control the label relaxation iteration.
]

#let subgraph-concepts = [
  === Subgraphs

  A subgraph is a dictionary containing half-edge membership, a topology
  signature, subgraph-wide `data`, and per-half-edge `hedge-data`. Annotations
  remain Typst values: content and callbacks are not serialized through Wasm.
  Independent selections can annotate the same hedge without changing its graph
  data or each other.

  - `subgraph.select(g, nodes: (), edges: (), source: (), sink: (), hedges: ())`
    selects the union of the supplied references. Node and edge references are
    Typst labels or numeric IDs, supplied in arrays. `nodes` selects incident
    hedges, `edges` selects both available halves, and `source` / `sink` selects
    the named edges' structural halves, independently of drawing orientation.
    `hedges` selects exact numeric hedge IDs. Unknown references and missing
    requested halves are errors; an isolated node contributes no hedges.
  - `subgraph.label(g, label)` constructs a subgraph from a base62 label.
  - `subgraph.bits(g, bits)` constructs a subgraph from a boolean hedge array.
  - `subgraph.compass(g, compass)` selects half edges with a DOT compass point.
  - `subgraph.hedges(sg)`, `subgraph.contains(sg, hedge)`, and
    `subgraph.contains-edge(sg, edge)` inspect membership. The last checks whether
    either half of an edge record is selected.
  - `subgraph.complement(g, sg)` selects all other hedges. It preserves
    subgraph-wide data, but the newly selected hedges have no annotations.

  `subgraph.with-data(g, sg, data: value, hedge: updater)` returns an annotated
  selection. `data: auto` preserves subgraph-wide data and `hedge: none` preserves
  hedge annotations. A non-function `hedge` value replaces every selected
  hedge's annotation. A callback receives its original endpoint record plus
  `edge`, `edge-name`, and `flow`; `data` holds the existing selection annotation
  and `graph-data` holds the graph's half-edge data. Its return value replaces
  the annotation, even when it is `none`. To store a callback as an annotation,
  return it from the updater rather than passing it as the updater itself.
  `subgraph.hedge-data(sg, h)` reads one selected hedge's annotation.

  Graph-aware consumers check a signature of IDs, names, and source/sink
  incidence, not graph identity or placement. Selections survive layout and
  drawing/data changes, but incompatible topology or names are rejected.
  Re-select on a cut view rather than reusing its master's selections.

  *Breaking change:* raw subgraph bytes are no longer accepted by public APIs.
  Use the constructors above, or the wrapped results of `graph.cycles` and
  `graph.forests`, and pass the whole object, not its internal `bytes` field.
  `subgraph.to-label(sg)` exports membership only: rebuilding with
  `subgraph.label(g, label)` binds it to that graph but does not restore annotations
  or the old topology signature. Keep the object to preserve Typst sidecar data;
  a base62 label is not a complete serialization of an annotated selection.

  Passing `subgraph: sg` to `draw` highlights half-edges, while `layout` restricts
  placement or optimization. Neither opens edges; use `graph.cut` for that.

  === Directed Cut Views

  `graph.cut(g, left: left, right: right, boundary: patch)` opens *one weighted
  directed cut* and leaves `g` unchanged. The two selections are its drawing
  sides, not separate cuts or a partition of the vertices. They must be disjoint
  and contain exactly opposite halves of every selected paired edge in `g`.
  Already dangling edges cannot be selected for opening.

  Think of the drawing boundaries as L and R, with a positive seam passage from
  R to L. Selecting an edge's source on R and sink on L follows that positive
  direction; swapping them reverses it. This is independent of superficial
  drawing orientation, including reversed fermion arrows. Side names do not
  place endpoints automatically: `boundary` supplies the view's geometry.

  An annotation `(winding: n)` on either selected half requests a positive integer
  number of same-direction passages; it defaults to one. If both halves specify
  winding, they must agree. Each selected edge becomes a source stub, a sink
  stub, and `n - 1` middle segments. Thus winding two produces three segments
  without separate named seams or a caller-supplied ordering of cut objects.
  Winding is not a net count of cancelling back-and-forth crossings.

  This small example uses the manual's existing imports. It stores content on
  the original edge and endpoints, annotates the cut sides, and opens one edge
  with winding two. The boundary callback places both ordinary free endpoints
  and the generated middle endpoints using the same record interface:

  ```typ
  #let master = graph.build({
    node(<a>, pos: graph.pos(x: 0, y: 0))
    node(<b>, pos: graph.pos(x: 0, y: -2))
    edge(
      <link>,
      source(<a>, caption: [output]),
      sink(<b>, caption: [input]),
      label: [connection],
    )
  })
  #let left = subgraph.with-data(
    master, subgraph.select(master, sink: (<link>,)),
    data: [left boundary],
    hedge: h => (winding: 2, caption: h.graph-data.caption),
  )
  #let right = subgraph.with-data(
    master, subgraph.select(master, source: (<link>,)),
    data: [right boundary],
    hedge: h => (caption: h.graph-data.caption),
  )
  #let view = graph.cut(master, left: left, right: right, boundary: record => {
    let b = record.boundary
    (pos: graph.pos(
      x: if b.side == "left" { -3 } else { 3 },
      y: -2 * b.crossing,
    ))
  })
  #let view = layout(view, layout-algo: "tree")
  #assert(graph.edges(master).len() == 1)
  #assert(graph.edges(view).len() == 3)
  #assert(graph.edge-data(view, <link.1>).label == [connection])
  #assert(graph.boundaries(view).len() == 4)
  #draw(view, edge-label: edge => edge.label)
  ```

  Fragments are named `<link.0>` through `<link.n>` in underlying source-to-sink
  order, irrespective of L/R or drawing arrows. Unnamed edges use the native
  prefix `__linnest_cut_edge_ID`; generated name collisions are errors. Prefer
  explicit edge names when patching a view. Node, edge, and hedge `origin` refer
  to the immediate input graph, including on repeated cuts. Edge origins contain
  the input `edge` ID and `name`, plus `segment` and `winding`; uncut edges have
  `segment: none` and `winding: 0`. Node and hedge origins are input IDs; newly
  generated nodes have `origin: none`.

  Physical graph, node, edge, and half-edge Typst data are remapped, including
  content and callbacks. Middle source hedges inherit the original sink hedge's
  data, and middle sink hedges inherit the original source's data, following the
  boundary side they continue. The cut resets fragment geometry for re-layout,
  rather than copying an old midpoint pin to new endpoints. Apply view-specific
  bends, placement, labels, and `crossing-under` targets with `graph.map`; a
  reference to a split edge needs the appropriate fragment name.

  The `boundary` callback runs only on newly created endpoints, receiving a node
  or dangling-edge record and returning ordinary `graph.map`-style placement/data
  patches. Its `record.boundary` has
  output `edge` and optional `node` IDs, the original side's `hedge`, `side`
  (`"left"` or `"right"`), zero-based `crossing` along underlying flow, `data`
  from that hedge's annotation, `cut-data` from that side's subgraph-wide data,
  and the fragment's `origin` at creation. For example, the code above preserves
  each endpoint caption in `boundary.data.caption`.

  On subsequent cuts, `boundary.origin` and `boundary.hedge` still refer to the
  input to the cut that created the endpoint; `side`, `crossing`, `data`, and
  `cut-data` are retained. Only `boundary.edge` and `boundary.node` remap to
  current anchors. Inherited boundaries retain their placements and annotations
  without invoking the callback again.

  `graph.boundaries(view)` returns these descriptors plus current `pos`, not
  coordinates cached at cut time. A dangling endpoint uses its edge position
  and has `node: none`; a middle endpoint uses its auxiliary node position.
  Query after layout, or inside `draw-after`, to group solved endpoints into
  background boxes with existing CeTZ primitives. Auxiliary boundary nodes are
  automatically zero-sized and unpainted, including their labels and custom
  node-drawing callbacks; they remain addressable anchors for incident edges.
]

#let _show-example-source(code, ..args) = {
  let displayed-code = code
    .text
    .split("\n")
    .filter(line => not line.starts-with(">>>"))
    .map(line => line.trim("<<<", at: start))
    .join("\n")
  block(raw(displayed-code, lang: "typ", block: true))
}
// Typst 0.15 drops Tidy's variable-name stack from HTML, so restore its heading there.
#let _show-variable(variable-doc, style-args) = {
  if target() == "html" {
    heading(variable-doc.name, level: style-args.first-heading-level + 1)
  }
  tidy.styles.default.show-variable(variable-doc, style-args)
}
#let _api-link(prefix, name, display: none) = {
  raw(if display == none { name } else { display }, lang: none)
}

#let _tidy-scope = (
  api-link: _api-link,
  draw: draw,
  graph: graph,
  layout: layout,
  subgraph: subgraph,
)

#let _show-reference(
  source,
  name,
  preamble: "",
  function-aliases: (:),
) = context {
  let style = dictionary(tidy.styles.default)
  let _ = style.insert("show-example", _show-example-source)
  let _ = style.insert("show-variable", _show-variable)
  let docs = tidy.parse-module(
    source,
    name: name,
    preamble: preamble,
    scope: _tidy-scope,
  )
  for (index, function-doc) in docs.functions.enumerate() {
    if function-doc.name in function-aliases {
      function-doc.name = function-aliases.at(function-doc.name)
      docs.functions.at(index) = function-doc
    }
  }
  tidy.show-module(
    docs,
    style: style,
    sort-functions: none,
    enable-cross-references: false,
  )
}

#let graph-reference = _show-reference(
  read("../src/graph.typ"),
  "graph",
  preamble: "#import graph: *\n",
)
#let layout-reference = _show-reference(
  read("../src/layout.typ"),
  "Reference",
  function-aliases: (sequence: "layout-sequence"),
)
#let drawing-reference = _show-reference(read("../src/draw.typ"), "Reference")
#let subgraph-reference = _show-reference(
  read("../src/subgraph.typ"),
  "subgraph",
)

#let manual = [
  #linnest-guide
  #graph-objects
  #graph-concepts
  #placement-concepts
  #drawing-concepts
  #domain-style-concepts
  #graph-query-concepts
  #layout-concepts
  #subgraph-concepts

  == Generated Reference

  #graph-reference
  #layout-reference
  #drawing-reference
  #subgraph-reference
]

#manual
