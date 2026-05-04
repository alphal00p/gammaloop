#import "@preview/cetz:0.3.4" as cetz
#import "curve.typ" as curve-api
#import "graph.typ" as graph-api
#import "subgraph.typ" as subgraph-api

#let _eval-dict(data, name, scope) = {
  let dict = data.at(name, default: "(:)")
  if dict == none {
    dict = "(:)"
  }
  eval(dict, scope: scope + data, mode: "code")
}

#let _eval-content(data, name, scope) = {
  let content = data.at(name, default: "[]")
  if content == none {
    content = "[]"
  }
  eval(content, scope: scope + data, mode: "code")
}

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
    ctrl_a: ctrl-a,
    ctrl_b: ctrl-b,
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

#let _pattern-name(style) = _style-value(style, "pattern", none)

#let _has-pattern(style) = {
  let pattern = _pattern-name(style)
  pattern != none and pattern != "normal" and pattern != "curve"
}

#let _same-pattern-geometry(source-style, sink-style) = {
  let same = _has-pattern(source-style)
  same = same and _pattern-name(source-style) == _pattern-name(sink-style)
  same = same and _style-value(source-style, "pattern-amplitude", 0.1) == _style-value(sink-style, "pattern-amplitude", 0.1)
  same = same and _style-value(source-style, "pattern-wavelength", 1.0) == _style-value(sink-style, "pattern-wavelength", 1.0)
  same = same and _style-value(source-style, "pattern-phase", 0) == _style-value(sink-style, "pattern-phase", 0)
  same = same and _style-value(source-style, "pattern-samples-per-period", 16) == _style-value(sink-style, "pattern-samples-per-period", 16)
  same = same and _style-value(source-style, "pattern-coil-longitudinal-scale", 1.25) == _style-value(sink-style, "pattern-coil-longitudinal-scale", 1.25)
  same = same and _style-value(source-style, "pattern-accuracy", 0.001) == _style-value(sink-style, "pattern-accuracy", 0.001)
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
/// DOT `eval_*` statements are evaluated in `scope` and forwarded to CeTZ.
///
/// ```example
/// #let b = graph.builder(
///   name: "demo",
///   edge-statements: (
///     eval_source: "(stroke: red + 0.5pt)",
///     eval_sink: "(stroke: blue + 0.5pt)",
///   ),
/// )
/// #let (node: a, builder: b) = graph.node(b, name: "a")
/// #let (node: c, builder: b) = graph.node(b, name: "c")
/// #let b = graph.edge(b, source: (node: a), sink: (node: c))
/// #draw(layout(graph.finish(b)))
/// ```
/// -> content
#let draw(
  /// Laid-out graph object returned by `layout`. -> bytes
  graph,

  /// Additional symbols available while evaluating DOT `eval_*` statements.
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

  /// Default CeTZ node radius. Use `auto` to fit the node label.
  /// -> auto | int | float | array
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

  /// Default CeTZ edge stroke. -> any
  edge-stroke: 0.1em,

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
  subgraph-edge-style: (stroke: rgb("#ffd166") + 0.35em),

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
        let subgraph-hedges = if subgraph == none { none } else { subgraph-api.hedges(subgraph) }

        for (i, v) in nodes.enumerate() {
          let pos = v.pos
          let ev = v.eval
          if ev == none {
            ev = "(:)"
          }
          let label = if v.name == none { none } else { [#v.name] }
          let node-style = (
            radius: node-radius,
            fill: node-fill,
            stroke: node-stroke,
          ) + eval(ev, scope: scope + (vid: i), mode: "code")
          let node-style = node-style + (
            radius: _node-radius(ctx, label, node-style, node-min-radius, node-label-padding),
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
          let local-scope = scope + (orientation: orientation) + (eid: i) + (ext: ext)
          let source-in-subgraph = _in-subgraph(subgraph-hedges, start)
          let sink-in-subgraph = _in-subgraph(subgraph-hedges, end)

          let source-style = (stroke: edge-stroke) + _eval-dict(data, "eval_source", local-scope)
          let sink-style = (stroke: edge-stroke) + _eval-dict(data, "eval_sink", local-scope)
          let source-draw-style = if source-in-subgraph and not subgraph-edge-underlay {
            source-style + subgraph-edge-style
          } else {
            source-style
          }
          let sink-draw-style = if sink-in-subgraph and not subgraph-edge-underlay {
            sink-style + subgraph-edge-style
          } else {
            sink-style
          }
          let ev-label = _eval-content(
            data,
            "eval_label",
            local-scope + (sink: sink) + (source: source) + data,
          )

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

          let label-pos = if e.label_pos == none { e.pos } else { e.label_pos }
          elements.push(cetz.draw.content(_point(label-pos), ev-label, padding: 0))

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
    grid(align: center + top, gutter: 1em, [#shown-title], diagram-content)
  }
}
