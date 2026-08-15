// Public drawing API.
//
// Keep the user-facing drawing functions and their documentation here. The
// helper-heavy implementation lives in `impl/draw.typ`.

#import "impl/draw.typ" as _impl

/// Split a laid-out graph edge into source and sink half-edge paths.
///
/// The returned dictionary has `source`, `sink`, and `curve`. The split point is
/// the edge layout point, so the two half-edges join smoothly there.
/// -> dictionary
#let edge-halves(
  /// Edge record returned by `graph.edges(layout(g))`. -> dictionary
  edge,
  /// Node records returned by `graph.nodes(layout(g))`. -> array
  nodes,
  /// Hobby curl used for the source-to-edge and edge-to-sink curves. -> float
  omega: 1.0,
  /// Arc-length trim applied at the source node side. -> int | float
  source-outset: 0,
  /// Arc-length trim applied at the sink node side. -> int | float
  sink-outset: 0,
  /// Arc-length accuracy used while trimming. -> float
  accuracy: 0.001,
) = {
  _impl.edge-halves(
    edge,
    nodes,
    (
      omega: omega,
      source-outset: source-outset,
      sink-outset: sink-outset,
      accuracy: accuracy,
    ),
  )
}

/// Draw the two halves of a laid-out graph edge through CeTZ.
///
/// `source-style` applies from the source node to the edge layout point, and
/// `sink-style` applies from the edge layout point to the sink node.
/// -> content
#let to-cetz-edge-halves(
  /// Edge record returned by `graph.edges(layout(g))`. -> dictionary
  edge,
  /// Node records returned by `graph.nodes(layout(g))`. -> array
  nodes,
  /// Coordinate length for one graph-layout unit. Numbers are interpreted as em.
  /// -> int | float | length | ratio
  unit: 1,
  /// Hobby curl used for the source-to-edge and edge-to-sink curves. -> float
  omega: 1.0,
  /// Arc-length trim applied at the source node side. -> int | float
  source-outset: 0,
  /// Arc-length trim applied at the sink node side. -> int | float
  sink-outset: 0,
  /// Arc-length accuracy used while trimming. -> float
  accuracy: 0.001,
  /// CeTZ style for the source half edge. -> dictionary
  source-style: (:),
  /// CeTZ style for the sink half edge. -> dictionary
  sink-style: (:),
) = {
  _impl.to-cetz-edge-halves(
    edge,
    nodes,
    (
      unit: unit,
      omega: omega,
      source-outset: source-outset,
      sink-outset: sink-outset,
      accuracy: accuracy,
      source-style: source-style,
      sink-style: sink-style,
    ),
  )
}

/// Draw a graph object with CeTZ.
///
/// The graph must already have positions, either from `layout` or from explicit
/// `pos` values passed to graph node or edge items. Paired
/// edges without an explicit edge position use the midpoint of their endpoint
/// nodes. Node and edge Typst style dictionaries or callbacks are evaluated and
/// forwarded to CeTZ.
///
/// #example(`
/// #let positioned = graph.build({
///   graph.node(<left>, label: [left], pos: graph.pos(x: 0, y: 0))
///   graph.node(<right>, label: [right], pos: graph.pos(ref: <left>, dx: 2.5, dy: 0))
///   graph.edge(graph.source(<left>), graph.sink(<right>))
///   graph.edge(graph.source(<right>), <right-out>, pos: graph.pos(ref: <right>, dx: 0.9, dy: 0.7))
/// },
///   name: "positioned",
/// )
/// #draw(positioned, source-style: (stroke: black + 0.7pt), sink-style: (stroke: black + 0.7pt))
/// `,dir:ttb)
///
/// #example(`
/// #let parallel-base-edge = 0
/// #let source-patterns = (
///   "coil",
///   "zigzag",
///   "coil",
///   "wave",
///   "zigzag",
///   "coil",
///   "wave",
///   "zigzag",
/// )
/// #let sink-patterns = (
///   "coil",
///   "zigzag",
///   "coil",
///   "zigzag",
///   "coil",
///   "wave",
///   "coil",
///   "wave",
/// )
/// #let parallel-layer(edge, mark: none) = if edge.eid == parallel-base-edge {
///   (
///     offset: -0.5,
///     length: 1.5,
///     ratio: 0.5,
///     resolve-length: "min",
///     stroke: (paint: rgb("#2f6f4e"), thickness: 1.1pt, cap: "round"),
///     mark: mark,
///   )
/// } else { none }
/// #let stack(base, layer) = if layer == none { base } else { (base, layer) }
/// #let source-stroke(edge) = if edge.eid == parallel-base-edge {
///   (paint: gray, thickness: 0.75pt, cap: "round")
/// } else {
///   (paint: red, thickness: 1.2pt, cap: "round")
/// }
/// #let sink-stroke(edge) = if edge.eid == parallel-base-edge {
///   (paint: gray, thickness: 0.75pt, cap: "round")
/// } else {
///   (paint: blue, thickness: 1.2pt, cap: "round", dash: "dotted")
/// }
/// #let source-style(edge) = {
///   let base = (
///     stroke: source-stroke(edge),
///     pattern: source-patterns.at(edge.eid),
///     pattern-amplitude: 0.18,
///     pattern-wavelength: 0.55,
///     pattern-coil-longitudinal-scale: 1.6,
///   )
///   stack(base, parallel-layer(edge))
/// }
/// #let sink-style(edge) = {
///   let base = (
///     stroke: sink-stroke(edge),
///     pattern: sink-patterns.at(edge.eid),
///     pattern-amplitude: 0.18,
///     pattern-wavelength: 0.55,
///     pattern-coil-longitudinal-scale: 1.6,
///   )
///   stack(base, parallel-layer(edge, mark: (end: (symbol: "straight"), scale: 0.75)))
/// }
/// #let g = graph.build({
///   graph.node(<a>, label: [a hi])
///   graph.node(<c>)
///   graph.node(<d>)
///   graph.node(<e>)
///   graph.edge(graph.source(<a>), <ac>, graph.sink(<c>), pos: graph.pos(x: 0, y: 0.75, mode: "pin"))
///   graph.edge(graph.source(<c>), <ca>, graph.sink(<a>, compass: "e"))
///   graph.edge(graph.source(<c>), <cd>, graph.sink(<d>, compass: "e"))
///   graph.edge(graph.source(<e>), <ed>, graph.sink(<d>, compass: "e"))
///   graph.edge(graph.source(<e>), <ea>, graph.sink(<a>, compass: "e"))
///   graph.edge(graph.source(<d>), <d-out>)
///   graph.edge(graph.source(<e>, compass: "e"), <e-out>)
///   graph.edge(graph.source(<a>), <a-out>)
/// },
///   name: "demo",
/// )
/// #let layed-out = layout(g)
/// #let east = subgraph.compass(layed-out, "e")
/// #draw(layed-out, subgraph: east, source-style: source-style, sink-style: sink-style)
/// `,dir:ttb)
///
/// #example(`
/// #let p = graph.build({
///   graph.node(<pa>, label: [a], pos: graph.pos(y: 0, mode: "pin"))
///   graph.node(<pc>, label: [c], pos: graph.pos(y: 0, mode: "pin"))
///   graph.node(<pc1>, label: [c], pos: graph.pos(y: 0, mode: "pin"))
///   graph.node(<pc2>, label: [c], pos: graph.pos(y: 0, mode: "pin"))
///   graph.edge(graph.source(<pa>), <e0>, graph.sink(<pc>))
///   graph.edge(graph.source(<pc2>), <e1>, graph.sink(<pc1>))
/// },
///   name: "parallel demo",
/// )
/// #let parallel-edge-style(edge) = (
///     offset: -0.5,
///     length: 1.4,
///     ratio: 0.5,
///     resolve-length: "min",
///     stroke: (paint: rgb("#2f6f4e"), thickness: 1.1pt, cap: "round"),
///   )
///
/// #let focused-base-style(edge) = (
///   stroke: (paint: gray, thickness: 0.7pt, cap: "round"),
///   pattern: "coil",
///   pattern-amplitude: 0.14,
///   pattern-wavelength: 0.45,
///   pattern-coil-longitudinal-scale: 1.5,
/// )
/// #let focused-source-style(edge) = (focused-base-style(edge), parallel-edge-style(edge))
/// #let focused-sink-style(edge) = (focused-base-style(edge), parallel-edge-style(edge) + (
///   mark: (end: ">"),
/// ))
/// #draw(layout(p, g-center: 0.004, label-steps: 0,), unit: 1.25, source-style: focused-source-style, sink-style: focused-sink-style)
/// `,dir:ttb)
///
/// #example(`
/// #let g = graph.build({
///   graph.node(<a>, pos: graph.pos(x: 0, y: 0, mode: "pin"))
///   graph.node(<c>, pos: graph.pos(x: 3, y: 0, mode: "pin"))
///   graph.node(<d>, pos: graph.pos(x: 3, y: -1, mode: "pin"))
///   graph.edge(graph.source(<a>), <ac>, graph.sink(<c>), orientation: "default")
///   graph.edge(graph.source(<a>), <ad>, graph.sink(<d>), orientation: "reversed")
/// },
///   name: "oriented marks",
/// )
/// #let arrow = (
///   end: (symbol: ">", fill: black, anchor: "center", shorten-to: auto),
///   scale: 0.75,
/// )
/// #let oriented-arrow = (
///   stroke: black + 0.7pt,
///   mark: arrow,
///   mark-position: "center-if-dangling",
///   mark-orientation: "edge",
/// )
/// #draw(layout(g), unit: 1.4, source-style: oriented-arrow, sink-style: oriented-arrow)
/// `,dir:ttb)
/// -> content
#let draw(
  /// Graph object with positions from `layout` or explicit graph API `pos` fields. -> bytes
  graph,
  /// Additional values merged into node and edge callback dictionaries.
  /// -> dictionary
  scope: (:),
  /// Coordinate length for one graph-layout unit. `auto` uses any unit stored
  /// by `graph.style`, falling back to `1em`. Numbers are interpreted as em.
  /// -> auto | int | float | length | ratio
  unit: auto,
  /// Optional title displayed above the diagram. Use `auto` for the graph name.
  /// -> none | auto | content | string
  title: none,
  /// Optional subgraph or array of subgraphs whose half-edges are shaded.
  /// Array entries may be raw subgraphs or records like
  /// `(subgraph: sg, edge-style: (stroke: red + 2pt))`.
  /// -> none | bytes | array
  subgraph: none,
  /// Debug level. `1` enables CeTZ canvas debug; `2` also marks edge positions.
  /// -> bool | int
  debug: false,
  /// Default CeTZ node radius. Use `auto` to fit the node label. -> auto | int | float | array
  node-radius: auto,
  /// Minimum radius used when `node-radius` is `auto`. -> int | float
  node-min-radius: 0.16,
  /// Extra canvas-unit padding around labels when `node-radius` is `auto`.
  /// -> int | float
  node-label-padding: 0.08,
  /// Default CeTZ node fill. -> any
  node-fill: white,
  /// Default CeTZ node stroke. -> any
  node-stroke: black,
  /// Edge clearance from node centers. `auto` uses each node circle radius.
  /// Increase this when labels extend beyond the circle. -> auto | int | float
  node-outset: auto,
  /// Default node-label style forwarded to `cetz.draw.content`. -> dictionary
  node-label-style: (:),
  /// Node style dictionary or callback. A callback receives node data. -> dictionary | function | none
  node-style: (:),
  /// Node label content or callback. `auto` uses the node name. -> auto | content | string | function | none
  node-label: auto,
  /// CeTZ-compatible node drawing callback. `auto` draws the default circular
  /// node and label. A callback receives `(node, box)` and should return CeTZ
  /// draw elements; `box` contains `name`, `center`, `width`, `height`, `unit`,
  /// `label`, `label-style`, `style`, `radius`, and `node`.
  /// -> auto | function
  draw-node: auto,
  /// Default CeTZ edge stroke. -> any
  edge-stroke: 0.1em,
  /// Default normal offset for edge paths. Applied to the base edge geometry
  /// before patterns; node outsets then trim the shifted path. -> int | float
  edge-offset: 0,
  /// Maximum visible arc length for centered parallel edge paths. `none` keeps
  /// the full shifted path. -> none | int | float
  edge-length: none,
  /// Maximum visible fraction of the base edge length for centered parallel edge
  /// paths. Combined with `edge-length` according to `edge-resolve-length`.
  /// -> none | int | float
  edge-ratio: none,
  /// Resolve `edge-length` and `edge-ratio`. Accepted string
  /// values are `"min"`/`"shorter"`, `"max"`/`"longer"`, `"length"`/`"fixed"`,
  /// `"ratio"`/`"relative"`, or `"none"`/`"full"`. A function receives
  /// `(base-length, length, ratio)`. -> string | function
  edge-resolve-length: "min",
  /// Arc-length accuracy for fitted parallel edge paths. -> float
  edge-accuracy: 0.001,
  /// Let Kurbo optimize the fitted parallel path. -> bool
  edge-optimize: true,
  /// Source half-edge style dictionary, array of layer dictionaries, or
  /// callback. `mark-position: "center-if-dangling"` keeps an end marker at the
  /// paired-edge split point while centering it on dangling half edges.
  /// `mark-orientation: "edge"` makes a mark follow `edge.orientation` instead
  /// of raw path direction; reversed edges move the mark to the sink half and
  /// flip it, while undirected edges suppress it.
  /// `source-anchor` may be a CeTZ anchor name such as `"north"` or `"south"`
  /// to route this endpoint from a measured node-box anchor. By default,
  /// anchored paired edges use two smooth cubic halves through the edge layout
  /// point while preserving the endpoint anchor tangents. Set
  /// `route: "hobby-through"` to instead build a single Hobby spline through
  /// the source anchor, source guide, edge point, sink guide, and sink anchor,
  /// then split source/sink styling at the edge point. `route: "direct"` keeps
  /// the same anchored cubic routing but suppresses the default edge-position
  /// Hobby route.
  /// `route-points: "through"` also threads any layout-provided half-edge route
  /// points through that same Hobby path.
  /// `route: "straight-through"` draws the two straight force springs from
  /// source to edge position and from edge position to sink.
  /// -> dictionary | array | function | none
  source-style: (:),
  /// Sink half-edge style dictionary, array of layer dictionaries, or callback.
  /// A callback receives edge data. `mark-position: "center-if-dangling"` has
  /// the same dangling-edge behavior as for `source-style`, and
  /// `mark-orientation: "edge"` participates in the same orientation-aware mark
  /// placement.
  /// `sink-anchor` may be a CeTZ anchor name such as `"north"` or `"south"`
  /// to route this endpoint into a measured node-box anchor.
  /// `route: "direct"` has the same meaning as in `source-style`.
  /// -> dictionary | array | function | none
  sink-style: (:),
  /// Edge label content or callback. A callback receives edge data.
  /// -> content | string | function | none
  edge-label: none,
  /// Edge-label style forwarded to `cetz.draw.content`; use `anchor` to choose
  /// which point of the label is placed at the layout label position. A
  /// callback receives edge data.
  /// -> dictionary | function | none
  edge-label-style: (:),
  /// Hobby curl used at the endpoints of paired edge curves. -> float
  edge-omega: 1.0,
  /// Optional style key for anchored source/sink routes. Set
  /// `anchor-control-distance` in `source-style` or `sink-style` to override
  /// the automatic guide distance used by cubic anchored routes.
  /// -> auto | int | float
  /// Arc-length accuracy for trimming edge curves at node outsets. -> float
  edge-trim-accuracy: 0.001,
  /// CeTZ canvas padding. -> none | int | float | array | dictionary
  padding: 0.4,
  /// Radius for edge-position markers shown at `debug >= 2`. -> int | float
  debug-edge-radius: 0.08,
  /// Fill for edge-position markers shown at `debug >= 2`. -> any
  debug-edge-fill: rgb("#ff9f1c"),
  /// Stroke for edge-position markers shown at `debug >= 2`. -> any
  debug-edge-stroke: rgb("#d72638") + 0.35pt,
  /// Label fill for edge-position markers shown at `debug >= 2`. -> any
  debug-edge-label-fill: rgb("#7a1020"),
  /// Default CeTZ style used to shade half-edges included in `subgraph`.
  /// Individual subgraph records may override this with `edge-style`.
  /// -> dictionary
  subgraph-edge-style: (stroke: rgb("#ffd166") + 4.5pt),
  /// Draw subgraph shading below the normal half-edge style. -> bool
  subgraph-edge-underlay: true,
) = {
  _impl.draw(
    graph,
    (
      scope: scope,
      unit: unit,
      title: title,
      subgraph: subgraph,
      debug: debug,
      node-radius: node-radius,
      node-min-radius: node-min-radius,
      node-label-padding: node-label-padding,
      node-fill: node-fill,
      node-stroke: node-stroke,
      node-outset: node-outset,
      node-label-style: node-label-style,
      node-style: node-style,
      node-label: node-label,
      draw-node: draw-node,
      edge-stroke: edge-stroke,
      edge-offset: edge-offset,
      edge-length: edge-length,
      edge-ratio: edge-ratio,
      edge-resolve-length: edge-resolve-length,
      edge-accuracy: edge-accuracy,
      edge-optimize: edge-optimize,
      source-style: source-style,
      sink-style: sink-style,
      edge-label: edge-label,
      edge-label-style: edge-label-style,
      edge-omega: edge-omega,
      edge-trim-accuracy: edge-trim-accuracy,
      padding: padding,
      debug-edge-radius: debug-edge-radius,
      debug-edge-fill: debug-edge-fill,
      debug-edge-stroke: debug-edge-stroke,
      debug-edge-label-fill: debug-edge-label-fill,
      subgraph-edge-style: subgraph-edge-style,
      subgraph-edge-underlay: subgraph-edge-underlay,
    ),
  )
}
