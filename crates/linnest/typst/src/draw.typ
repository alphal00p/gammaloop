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

#let _line(start, end, style) = cetz.draw.line(_point(start), _point(end), ..style)

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

#let _without-pattern-style(style) = {
  let clean = style
  for key in _pattern-style-keys {
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

#let _same-pattern-geometry(source-style, sink-style) = {
  let same = _has-pattern(source-style)
  same = same and _pattern-name(source-style) == _pattern-name(sink-style)
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

#let _pattern-curve(segment, style) = {
  if _has-pattern(style) {
    curve-api.cetz-pattern(_pattern-path(segment, style), .._without-pattern-style(style))
  } else {
    curve-api.cetz-bezier(segment, .._without-pattern-style(style))
  }
}

#let _pattern-edge-halves(halves, source-style, sink-style) = {
  let elements = ()
  if _same-pattern-geometry(source-style, sink-style) {
    let source-path = _pattern-path(halves.source, source-style, anchor-end: false)
    let wavelength = _style-value(source-style, "pattern-wavelength", 1.0)
    let phase = _style-value(source-style, "pattern-phase", 0)
    let sink-phase = phase + 2 * calc.pi * source-path.length / wavelength
    elements.push(curve-api.cetz-pattern(source-path, .._without-pattern-style(source-style)))
    elements.push(curve-api.cetz-pattern(
      _pattern-path(halves.sink, sink-style, phase: sink-phase, anchor-start: false),
      .._without-pattern-style(sink-style),
    ))
  } else {
    elements.push(_pattern-curve(halves.source, source-style))
    elements.push(_pattern-curve(halves.sink, sink-style))
  }
  elements
}

#let _pattern-line(start, end, style) = {
  let pattern = _style-value(style, "pattern", none)
  let draw-style = _without-pattern-style(style)
  if pattern == none or pattern == "normal" or pattern == "curve" {
    _line(start, end, draw-style)
  } else {
    _pattern-curve(_line-segment(start, end), style)
  }
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
/// #let source-patterns = (
///   none,
///   "wave",
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
///   "wave",
///   "zigzag",
///   "coil",
///   "zigzag",
///   "coil",
///   "wave",
///   "coil",
///   "wave",
/// )
/// #let source-style(edge) = (
///   stroke: (paint: red, thickness: 1.2pt, cap: "round"),
///   pattern: source-patterns.at(edge.eid),
///   pattern-amplitude: 0.18,
///   pattern-wavelength: 0.55,
///   pattern-coil-longitudinal-scale: 1.6,
/// )
/// #let sink-style(edge) = (
///   stroke: (paint: blue, thickness: 1.2pt, cap: "round", dash: "dotted"),
///   pattern: sink-patterns.at(edge.eid),
///   pattern-amplitude: 0.18,
///   pattern-wavelength: 0.55,
///   pattern-coil-longitudinal-scale: 1.6,
/// )
/// #let (node: a, builder: b) = graph.node(b, name: "a hi")
/// #let (node: c, builder: b) = graph.node(b, name: "c")
/// #let (node: d, builder: b) = graph.node(b, name: "d")
/// #let (node: e, builder: b) = graph.node(b, name: "e")
/// #let b = graph.edge(b, source: (node: a), sink: (node: c, compass: "e"))
/// #let b = graph.edge(b, source: (node: a, compass: "e"), sink: (node: c,  compass: "e"))
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

          let source-style-value = (stroke: edge-stroke) + _style(source-style, edge-data)
          let sink-style-value = (stroke: edge-stroke) + _style(sink-style, edge-data)
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
            let halves = curve-api.edge-halves(
              e,
              nodes,
              omega: edge-omega,
              source-outset: node-outsets.at(start.node),
              sink-outset: node-outsets.at(end.node),
              accuracy: edge-trim-accuracy,
            )
            if source-in-subgraph and subgraph-edge-underlay {
              elements.push(curve-api.cetz-bezier(halves.source, ..subgraph-edge-style))
            }
            if sink-in-subgraph and subgraph-edge-underlay {
              elements.push(curve-api.cetz-bezier(halves.sink, ..subgraph-edge-style))
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
              elements.push(_line(line-start, e.pos, subgraph-edge-style))
            }
            elements.push(_pattern-line(line-start, e.pos, source-draw-style))
          } else if end != none {
            let line-end = curve-api.outset-point(
              nodes.at(end.node).pos,
              e.pos,
              distance: node-outsets.at(end.node),
            )
            if sink-in-subgraph and subgraph-edge-underlay {
              elements.push(_line(e.pos, line-end, subgraph-edge-style))
            }
            elements.push(_pattern-line(e.pos, line-end, sink-draw-style))
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
