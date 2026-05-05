#import "@preview/cetz:0.3.4" as cetz
#import "curve.typ" as curve-api
#import "graph.typ" as graph-api
#import "subgraph.typ" as subgraph-api

#let _point(p) = (p.x, p.y)

#let _canvas-length(unit) = {
  if type(unit) in (int, float) {
    unit * 1em
  } else {
    unit
  }
}

#let _radius-outset(radius) = {
  if type(radius) == array {
    calc.max(..radius)
  } else {
    radius
  }
}

#let _line-segment(start, end) = {
  let ctrl-a = (
    x: start.x + (end.x - start.x) / 3,
    y: start.y + (end.y - start.y) / 3,
  )
  let ctrl-b = (
    x: start.x + (end.x - start.x) * 2 / 3,
    y: start.y + (end.y - start.y) * 2 / 3,
  )
  (
    start: start,
    ctrl-a: ctrl-a,
    ctrl-b: ctrl-b,
    end: end,
  )
}

#let _pattern-style-keys = (
  "pattern",
  "pattern-amplitude",
  "pattern-wavelength",
  "pattern-phase",
  "pattern-samples-per-period",
  "pattern-coil-longitudinal-scale",
  "pattern-accuracy",
)

#let _parallel-style-keys = (
  "parallel-distance",
  "parallel-length",
  "parallel-ratio",
  "parallel-accuracy",
  "parallel-optimize",
)

#let _without-pattern-style(style) = {
  let clean = style
  for key in _pattern-style-keys {
    if clean.keys().contains(key) {
      let _ = clean.remove(key)
    }
  }
  clean
}

#let _draw-style(style) = {
  let clean = _without-pattern-style(style)
  for key in _parallel-style-keys {
    if clean.keys().contains(key) {
      let _ = clean.remove(key)
    }
  }
  clean
}

#let _style-value(style, key, default) = style.at(key, default: default)

#let _call(value, data) = {
  if type(value) == function {
    value(data)
  } else {
    value
  }
}

#let _style(value, data) = {
  let value = _call(value, data)
  if value == none { (:) } else { value }
}

#let _content(value, data, default: none) = {
  if value == auto {
    default
  } else {
    let value = _call(value, data)
    if type(value) == str { [#value] } else { value }
  }
}

#let _pattern-name(style) = _style-value(style, "pattern", none)

#let _has-pattern(style) = {
  let pattern = _pattern-name(style)
  pattern != none and pattern != "normal" and pattern != "curve"
}

#let _has-parallel(style) = {
  let distance = _style-value(style, "parallel-distance", 0)
  distance != none and distance != 0
}

#let _same-parallel-geometry(source-style, sink-style) = {
  let same = _style-value(source-style, "parallel-distance", 0) == _style-value(sink-style, "parallel-distance", 0)
  same = same and _style-value(source-style, "parallel-length", none) == _style-value(sink-style, "parallel-length", none)
  same = same and _style-value(source-style, "parallel-ratio", none) == _style-value(sink-style, "parallel-ratio", none)
  same = same and _style-value(source-style, "parallel-accuracy", 0.001) == _style-value(sink-style, "parallel-accuracy", 0.001)
  same = same and _style-value(source-style, "parallel-optimize", true) == _style-value(sink-style, "parallel-optimize", true)
  same
}

#let _same-pattern-geometry(source-style, sink-style) = {
  let same = _has-pattern(source-style)
  same = same and _pattern-name(source-style) == _pattern-name(sink-style)
  same = same and _same-parallel-geometry(source-style, sink-style)
  same = (
    same and _style-value(source-style, "pattern-amplitude", 0.1) == _style-value(sink-style, "pattern-amplitude", 0.1)
  )
  same = (
    same
      and _style-value(source-style, "pattern-wavelength", 1.0) == _style-value(sink-style, "pattern-wavelength", 1.0)
  )
  same = same and _style-value(source-style, "pattern-phase", 0) == _style-value(sink-style, "pattern-phase", 0)
  same = (
    same
      and _style-value(source-style, "pattern-samples-per-period", 16)
        == _style-value(sink-style, "pattern-samples-per-period", 16)
  )
  same = (
    same
      and _style-value(source-style, "pattern-coil-longitudinal-scale", 1.25)
        == _style-value(sink-style, "pattern-coil-longitudinal-scale", 1.25)
  )
  same = (
    same
      and _style-value(source-style, "pattern-accuracy", 0.001) == _style-value(sink-style, "pattern-accuracy", 0.001)
  )
  same
}

#let _parallel-path(segment, style, start-outset: 0, end-outset: 0) = {
  curve-api.parallel-segment(
    segment,
    distance: _style-value(style, "parallel-distance", 0),
    start-outset: start-outset,
    end-outset: end-outset,
    accuracy: _style-value(style, "parallel-accuracy", 0.001),
    optimize: _style-value(style, "parallel-optimize", true),
  )
}

#let _segment-length(segment, accuracy: 0.001) = {
  curve-api.parallel-segment(segment, accuracy: accuracy, optimize: false).length
}

#let _segments-length(segments, accuracy: 0.001) = {
  let length = 0
  for segment in segments {
    length = length + _segment-length(segment, accuracy: accuracy)
  }
  length
}

#let _parallel-target-length(base-length, style) = {
  let limits = ()
  let length = _style-value(style, "parallel-length", none)
  if length != none and length > 0 {
    limits.push(length)
  }
  let ratio = _style-value(style, "parallel-ratio", none)
  if ratio != none and ratio > 0 {
    limits.push(base-length * ratio)
  }
  if limits.len() == 0 { none } else { calc.min(..limits) }
}

#let _parallel-center-outset(base-length, style, source-outset: 0, sink-outset: 0) = {
  let target = _parallel-target-length(base-length, style)
  if target == none {
    0
  } else {
    let visible-length = calc.max(0, base-length - source-outset - sink-outset)
    calc.max(0, (visible-length - target) / 2)
  }
}

#let _path-segments(path) = {
  let segments = ()
  if path.keys().contains("curves") and path.curves.len() > 0 {
    for segment in path.curves {
      segments.push(segment)
    }
  } else if path.keys().contains("points") and path.points.len() > 1 {
    for i in range(0, path.points.len() - 1) {
      segments.push(_line-segment(path.points.at(i), path.points.at(i + 1)))
    }
  }
  segments
}

#let _trim-segments(segments, start-outset: 0, end-outset: 0, accuracy: 0.001) = {
  let trimmed = ()
  for (index, segment) in segments.enumerate() {
    let piece = segment
    if index == 0 and start-outset != 0 {
      piece = curve-api.trim-segment(piece, start-outset: start-outset, accuracy: accuracy)
    }
    if index == segments.len() - 1 and end-outset != 0 {
      piece = curve-api.trim-segment(piece, end-outset: end-outset, accuracy: accuracy)
    }
    trimmed.push(piece)
  }
  trimmed
}

#let _geometry-segments(segment, style, start-outset: 0, end-outset: 0, accuracy: 0.001) = {
  if _has-parallel(style) {
    _path-segments(_parallel-path(segment, style, start-outset: start-outset, end-outset: end-outset))
  } else {
    _trim-segments((segment,), start-outset: start-outset, end-outset: end-outset, accuracy: accuracy)
  }
}

#let _edge-geometry-halves(edge, nodes, source-style, sink-style, omega: 1.0, source-outset: 0, sink-outset: 0, accuracy: 0.001) = {
  if _has-parallel(source-style) or _has-parallel(sink-style) {
    let base = curve-api.edge-halves(edge, nodes, omega: omega, accuracy: accuracy)
    let base-length = _segments-length(base.curve.segments, accuracy: accuracy)
    let source-center-outset = _parallel-center-outset(base-length, source-style, source-outset: source-outset, sink-outset: sink-outset)
    let source-geometry = _geometry-segments(
      base.curve.segments.at(0),
      source-style,
      start-outset: source-outset + source-center-outset,
      accuracy: accuracy,
    )
    let sink-geometry = if _same-parallel-geometry(source-style, sink-style) {
      _geometry-segments(
        base.curve.segments.at(1),
        source-style,
        end-outset: sink-outset + source-center-outset,
        accuracy: accuracy,
      )
    } else {
      let sink-center-outset = _parallel-center-outset(base-length, sink-style, source-outset: source-outset, sink-outset: sink-outset)
      _geometry-segments(
        base.curve.segments.at(1),
        sink-style,
        end-outset: sink-outset + sink-center-outset,
        accuracy: accuracy,
      )
    }
    (
      source: source-geometry,
      sink: sink-geometry,
      curve: base.curve,
    )
  } else {
    let trimmed = curve-api.edge-halves(
      edge,
      nodes,
      omega: omega,
      source-outset: source-outset,
      sink-outset: sink-outset,
      accuracy: accuracy,
    )
    (
      source: (trimmed.source,),
      sink: (trimmed.sink,),
      curve: trimmed.curve,
    )
  }
}

#let _pattern-path(segment, style, phase: auto, anchor-start: true, anchor-end: true) = {
  curve-api.pattern-segment(
    segment,
    pattern: _pattern-name(style),
    amplitude: _style-value(style, "pattern-amplitude", 0.1),
    wavelength: _style-value(style, "pattern-wavelength", 1.0),
    phase: if phase == auto { _style-value(style, "pattern-phase", 0) } else { phase },
    samples-per-period: _style-value(style, "pattern-samples-per-period", 16),
    coil-longitudinal-scale: _style-value(style, "pattern-coil-longitudinal-scale", 1.25),
    anchor-start: anchor-start,
    anchor-end: anchor-end,
    accuracy: _style-value(style, "pattern-accuracy", 0.001),
  )
}

#let _segments-elements(segments, style, phase: auto, anchor-start: true, anchor-end: true) = {
  let elements = ()
  let length = 0
  if _has-pattern(style) {
    let current-phase = if phase == auto { _style-value(style, "pattern-phase", 0) } else { phase }
    let wavelength = _style-value(style, "pattern-wavelength", 1.0)
    for (index, segment) in segments.enumerate() {
      let piece = _pattern-path(
        segment,
        style,
        phase: current-phase,
        anchor-start: anchor-start and index == 0,
        anchor-end: anchor-end and index == segments.len() - 1,
      )
      elements.push(curve-api.cetz-pattern(piece, .._draw-style(style)))
      length = length + piece.length
      current-phase = current-phase + 2 * calc.pi * piece.length / wavelength
    }
  } else {
    for segment in segments {
      elements.push(curve-api.cetz-bezier(segment, .._draw-style(style)))
    }
  }
  (elements: elements, length: length)
}

#let _segment-elements(segment, style, phase: auto, anchor-start: true, anchor-end: true) = {
  _segments-elements(_geometry-segments(segment, style), style, phase: phase, anchor-start: anchor-start, anchor-end: anchor-end)
}

#let _pattern-edge-halves(halves, source-style, sink-style) = {
  let elements = ()
  if _same-pattern-geometry(source-style, sink-style) {
    let source = _segments-elements(halves.source, source-style, anchor-end: false)
    let wavelength = _style-value(source-style, "pattern-wavelength", 1.0)
    let phase = _style-value(source-style, "pattern-phase", 0)
    let sink-phase = phase + 2 * calc.pi * source.length / wavelength
    for element in source.elements {
      elements.push(element)
    }
    for element in _segments-elements(halves.sink, sink-style, phase: sink-phase, anchor-start: false).elements {
      elements.push(element)
    }
  } else {
    for element in _segments-elements(halves.source, source-style).elements {
      elements.push(element)
    }
    for element in _segments-elements(halves.sink, sink-style).elements {
      elements.push(element)
    }
  }
  elements
}

#let _pattern-line(start, end, style) = {
  _segment-elements(_line-segment(start, end), style).elements
}

#let _node-outset(style, node-outset) = {
  if node-outset == auto {
    _radius-outset(style.at("radius", default: 0))
  } else {
    node-outset
  }
}

#let _fit-node-radius(ctx, label, minimum, padding) = {
  if label == none {
    minimum
  } else {
    let size = measure(text(top-edge: "cap-height", bottom-edge: "bounds", label))
    let width = calc.abs(size.width / ctx.length) + 2 * padding
    let height = calc.abs(size.height / ctx.length) + 2 * padding
    calc.max(minimum, calc.sqrt(width * width + height * height) / 2)
  }
}

#let _node-radius(ctx, label, style, minimum, padding) = {
  let radius = style.at("radius", default: auto)
  if radius == auto {
    _fit-node-radius(ctx, label, minimum, padding)
  } else {
    radius
  }
}

#let _debug-level(debug) = {
  if type(debug) == bool {
    if debug { 1 } else { 0 }
  } else {
    debug
  }
}

#let _in-subgraph(hedges, endpoint) = {
  hedges != none and endpoint != none and hedges.contains(endpoint.hedge)
}

/// Draw a laid-out graph object with CeTZ.
///
/// The graph must already have positions, so call `layout` first. Node and edge
/// Typst style dictionaries or callbacks are evaluated and forwarded to CeTZ.
///
/// #example(`
/// #let b = graph.builder(name: "demo")
/// #let parallel-base-edge = 0
/// #let parallel-edge = 1
/// #let source-patterns = (
///   none,
///   none,
///   "zigzag",
///   "coil",
///   "wave",
///   "zigzag",
///   "coil",
///   "wave",
///   "zigzag",
/// )
/// #let sink-patterns = (
///   none,
///   none,
///   "zigzag",
///   "coil",
///   "zigzag",
///   "coil",
///   "wave",
///   "coil",
///   "wave",
/// )
/// #let parallel-style(edge) = if edge.eid == parallel-edge {
///   (
///     parallel-distance: 0.16,
///     parallel-length: 1.1,
///     parallel-ratio: 0.5,
///   )
/// } else { (:) }
/// #let source-stroke(edge) = if edge.eid == parallel-base-edge {
///   (paint: gray, thickness: 0.75pt, cap: "round")
/// } else if edge.eid == parallel-edge {
///   (paint: rgb("#2f6f4e"), thickness: 1.1pt, cap: "round")
/// } else {
///   (paint: red, thickness: 1.2pt, cap: "round")
/// }
/// #let sink-stroke(edge) = if edge.eid == parallel-base-edge {
///   (paint: gray, thickness: 0.75pt, cap: "round")
/// } else if edge.eid == parallel-edge {
///   (paint: rgb("#2f6f4e"), thickness: 1.1pt, cap: "round")
/// } else {
///   (paint: blue, thickness: 1.2pt, cap: "round", dash: "dotted")
/// }
/// #let source-style(edge) = (
///   stroke: source-stroke(edge),
///   pattern: source-patterns.at(edge.eid),
///   pattern-amplitude: 0.18,
///   pattern-wavelength: 0.55,
///   pattern-coil-longitudinal-scale: 1.6,
/// ) + parallel-style(edge)
/// #let sink-style(edge) = (
///   stroke: sink-stroke(edge),
///   pattern: sink-patterns.at(edge.eid),
///   pattern-amplitude: 0.18,
///   pattern-wavelength: 0.55,
///   pattern-coil-longitudinal-scale: 1.6,
/// ) + parallel-style(edge)
/// #let (node: a, builder: b) = graph.node(b, name: "a hi")
/// #let (node: c, builder: b) = graph.node(b, name: "c")
/// #let (node: d, builder: b) = graph.node(b, name: "d")
/// #let (node: e, builder: b) = graph.node(b, name: "e")
/// #let b = graph.edge(b, source: (node: a), sink: (node: c), statements: (pin: "x:0,y:0.75"))
/// #let b = graph.edge(b, source: (node: a), sink: (node: c), statements: (pin: "x:0,y:0.75"))
/// #let b = graph.edge(b, source: (node: c), sink: (node: a, compass: "e"))
/// #let b = graph.edge(b, source: (node: c), sink: (node: d, compass: "e"))
/// #let b = graph.edge(b, source: (node: e), sink: (node: d, compass: "e"))
/// #let b = graph.edge(b, source: (node: e), sink: (node: a, compass: "e"))
/// #let b = graph.edge(b, source: (node: d), sink: none)
/// #let b = graph.edge(b, source: (node: e, compass: "e"), sink: none)
/// #let b = graph.edge(b, source: (node: a), sink: none)
/// #let layed-out = layout(graph.finish(b))
/// #let east = subgraph.compass(layed-out, "e")
/// #draw(layed-out, subgraph: east, source-style: source-style, sink-style: sink-style)
/// `,dir:ttb)
///
/// #example(`
/// #let p = graph.builder(name: "parallel demo")
/// #let main-line-edge = 0
/// #let parallel-edge = 1
/// #let (node: pa, builder: p) = graph.node(p, name: "a", statements: (pin: "x:-2,y:0"))
/// #let (node: pc, builder: p) = graph.node(p, name: "c", statements: (pin: "x:2,y:0"))
/// #let p = graph.edge(p, source: (node: pa), sink: (node: pc), statements: (pin: "x:0,y:1.1"))
/// #let p = graph.edge(p, source: (node: pa), sink: (node: pc), statements: (pin: "x:0,y:1.1"))
/// #let parallel-edge-style(edge) = if edge.eid == parallel-edge {
///   (
///     parallel-distance: 0.18,
///     parallel-length: 1.4,
///     parallel-ratio: 0.5,
///   )
/// } else { (:) }
/// #let focused-style(edge) = (
///   stroke: if edge.eid == parallel-edge {
///     (paint: rgb("#2f6f4e"), thickness: 1.1pt, cap: "round")
///   } else {
///     (paint: gray, thickness: 0.7pt, cap: "round")
///   },
///   pattern: if edge.eid == main-line-edge { "coil" } else { none },
///   pattern-amplitude: 0.14,
///   pattern-wavelength: 0.45,
///   pattern-coil-longitudinal-scale: 1.5,
/// ) + parallel-edge-style(edge)
/// #draw(layout(graph.finish(p), steps: 1, epochs: 1, label-steps: 0), unit: 1.25, node-radius: 0.24, source-style: focused-style, sink-style: focused-style)
/// `,dir:ttb)
/// -> content
#let draw(
  /// Laid-out graph object returned by `layout`. -> bytes
  graph,
  /// Additional values merged into node and edge callback dictionaries.
  /// -> dictionary
  scope: (:),
  /// Coordinate length for one graph-layout unit. Numbers are interpreted as em.
  /// -> int | float | length | ratio
  unit: 1,
  /// Optional title displayed above the diagram. Use `auto` for the graph name.
  /// -> none | auto | content | string
  title: none,
  /// Optional subgraph whose half-edges are shaded. -> none | bytes
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
  /// Default CeTZ edge stroke. -> any
  edge-stroke: 0.1em,
  /// Default normal offset for edge paths. Applied to the base edge geometry
  /// before patterns; node outsets then trim the shifted path. -> int | float
  edge-parallel-distance: 0,
  /// Maximum visible arc length for centered parallel edge paths. `none` keeps
  /// the full shifted path. -> none | int | float
  edge-parallel-length: none,
  /// Maximum visible fraction of the base edge length for centered parallel edge
  /// paths. Combined with `edge-parallel-length` by taking the shorter length.
  /// -> none | int | float
  edge-parallel-ratio: none,
  /// Arc-length accuracy for fitted parallel edge paths. -> float
  edge-parallel-accuracy: 0.001,
  /// Let Kurbo optimize the fitted parallel path. -> bool
  edge-parallel-optimize: true,
  /// Source half-edge style dictionary or callback. A callback receives edge data.
  /// -> dictionary | function | none
  source-style: (:),
  /// Sink half-edge style dictionary or callback. A callback receives edge data.
  /// -> dictionary | function | none
  sink-style: (:),
  /// Edge label content or callback. A callback receives edge data.
  /// -> content | string | function | none
  edge-label: none,
  /// Hobby curl used at the endpoints of paired edge curves. -> float
  edge-omega: 1.0,
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
  /// CeTZ style used to shade half-edges included in `subgraph`. -> dictionary
  subgraph-edge-style: (stroke: rgb("#ffd166") + 4.5pt),
  /// Draw subgraph shading below the normal half-edge style. -> bool
  subgraph-edge-underlay: true,
) = {
  let info = graph-api.info(graph)
  let debug-level = _debug-level(debug)
  let diagram-content = cetz.canvas(
    length: _canvas-length(unit),
    debug: debug-level >= 1,
    padding: padding,
    {
      cetz.draw.get-ctx(ctx => {
        let nodes = graph-api.nodes(graph)
        let edges = graph-api.edges(graph)
        let elements = ()
        let node-elements = ()
        let node-outsets = ()
        let debug-level = _debug-level(debug)
        let subgraph-hedges = if subgraph == none { none } else { subgraph-api.hedges(subgraph) }

        for (i, v) in nodes.enumerate() {
          let pos = v.pos
          let node-data = scope + v.statements + (
            vid: i,
            node: v,
            name: v.name,
          )
          let default-label = if v.name == none { none } else { [#v.name] }
          let label = _content(node-label, node-data, default: default-label)
          let node-style-value = (
            (
              radius: node-radius,
              fill: node-fill,
              stroke: node-stroke,
            )
              + _style(node-style, node-data)
          )
          let node-style = (
            node-style-value
              + (
                radius: _node-radius(ctx, label, node-style-value, node-min-radius, node-label-padding),
              )
          )
          node-elements.push(cetz.draw.circle(_point(pos), name: "n" + str(i), ..node-style))
          node-outsets.push(_node-outset(node-style, node-outset))

          if label != none {
            node-elements.push(cetz.draw.content(
              _point(pos),
              label,
              padding: 0,
              ..node-label-style,
            ))
          }
        }

        for (i, e) in edges.enumerate() {
          let start = e.source
          let end = e.sink
          let source = if start == none { none } else { start.statement }
          let sink = if end == none { none } else { end.statement }
          let ext = start == none or end == none
          let data = e.statements
          let orientation = e.orientation
          let edge-data = scope + data + (
            eid: i,
            edge: e,
            source: source,
            sink: sink,
            source-endpoint: start,
            sink-endpoint: end,
            orientation: orientation,
            ext: ext,
          )
          let source-in-subgraph = _in-subgraph(subgraph-hedges, start)
          let sink-in-subgraph = _in-subgraph(subgraph-hedges, end)

          let geometry-style = (
            parallel-distance: edge-parallel-distance,
            parallel-length: edge-parallel-length,
            parallel-ratio: edge-parallel-ratio,
            parallel-accuracy: edge-parallel-accuracy,
            parallel-optimize: edge-parallel-optimize,
          )
          let source-style-value = (stroke: edge-stroke) + geometry-style + _style(source-style, edge-data)
          let sink-style-value = (stroke: edge-stroke) + geometry-style + _style(sink-style, edge-data)
          let source-draw-style = if source-in-subgraph and not subgraph-edge-underlay {
            source-style-value + subgraph-edge-style
          } else {
            source-style-value
          }
          let sink-draw-style = if sink-in-subgraph and not subgraph-edge-underlay {
            sink-style-value + subgraph-edge-style
          } else {
            sink-style-value
          }
          let ev-label = _content(edge-label, edge-data, default: none)

          if start != none and end != none {
            let halves = _edge-geometry-halves(
              e,
              nodes,
              source-style-value,
              sink-style-value,
              omega: edge-omega,
              source-outset: node-outsets.at(start.node),
              sink-outset: node-outsets.at(end.node),
              accuracy: edge-trim-accuracy,
            )
            if source-in-subgraph and subgraph-edge-underlay {
              for element in _segments-elements(halves.source, _without-pattern-style(source-style-value) + subgraph-edge-style).elements {
                elements.push(element)
              }
            }
            if sink-in-subgraph and subgraph-edge-underlay {
              for element in _segments-elements(halves.sink, _without-pattern-style(sink-style-value) + subgraph-edge-style).elements {
                elements.push(element)
              }
            }
            for element in _pattern-edge-halves(halves, source-draw-style, sink-draw-style) {
              elements.push(element)
            }
          } else if start != none {
            let line-start = curve-api.outset-point(
              nodes.at(start.node).pos,
              e.pos,
              distance: node-outsets.at(start.node),
            )
            if source-in-subgraph and subgraph-edge-underlay {
              for element in _pattern-line(line-start, e.pos, _without-pattern-style(source-style-value) + subgraph-edge-style) {
                elements.push(element)
              }
            }
            for element in _pattern-line(line-start, e.pos, source-draw-style) {
              elements.push(element)
            }
          } else if end != none {
            let line-end = curve-api.outset-point(
              nodes.at(end.node).pos,
              e.pos,
              distance: node-outsets.at(end.node),
            )
            if sink-in-subgraph and subgraph-edge-underlay {
              for element in _pattern-line(e.pos, line-end, _without-pattern-style(sink-style-value) + subgraph-edge-style) {
                elements.push(element)
              }
            }
            for element in _pattern-line(e.pos, line-end, sink-draw-style) {
              elements.push(element)
            }
          }

          if ev-label != none {
            let label-pos = if e.label-pos == none { e.pos } else { e.label-pos }
            elements.push(cetz.draw.content(_point(label-pos), ev-label, padding: 0))
          }

          if debug-level >= 2 {
            elements.push(cetz.draw.circle(
              _point(e.pos),
              radius: debug-edge-radius,
              fill: debug-edge-fill,
              stroke: debug-edge-stroke,
            ))
            elements.push(cetz.draw.content(
              _point(e.pos),
              text(size: 0.65em, fill: debug-edge-label-fill)[e#i],
              padding: 0,
            ))
          }
        }

        for element in node-elements {
          elements.push(element)
        }

        for element in elements {
          element
        }
      })
    },
  )

  if title == none {
    diagram-content
  } else {
    let shown-title = if title == auto { info.name } else { title }
    grid(
      align: center + top,
      gutter: 1em,
      [#shown-title],
      diagram-content,
    )
  }
}
