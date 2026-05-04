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

#let _node-outset(style, node-outset) = {
  if node-outset == auto {
    style.at("radius", default: 0)
  } else {
    node-outset
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
///     eval_source: "(stroke: red + 1.5pt)",
///     eval_sink: "(stroke: blue + 1.5pt)",
///   ),
/// )
/// #let (node: a, builder: b) = graph.node(b, name: "a")
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
/// #let layed-out = layout(graph.finish(b), beta: 50,gamma-ee:0.1,gamma-dangling:5,gamma-ev:0.01,seed: 2,epochs:30, steps: 30, g-center: 0.005, length-scale: 0.5)
/// #let east = subgraph.compass(layed-out, "e")
/// #draw(layed-out,node-radius: 1,subgraph:east)
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
  /// Default CeTZ node radius. -> int | float
  node-radius: 0.16,
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
  subgraph-edge-style: (stroke: rgb("#ffd166") + 4.5pt),
  /// Draw subgraph shading below the normal half-edge style. -> bool
  subgraph-edge-underlay: true,
) = {
  let info = graph-api.info(graph)
  let nodes = graph-api.nodes(graph)
  let edges = graph-api.edges(graph)
  let elements = ()
  let node-elements = ()
  let node-outsets = ()
  let debug-level = _debug-level(debug)
  let subgraph-hedges = if subgraph == none { none } else { subgraph-api.hedges(subgraph) }

  for (i, v) in nodes.enumerate() {
    let pos = v.pos
    let ev = v.eval
    if ev == none {
      ev = "(:)"
    }
    let node-style = (
      (
        radius: node-radius,
        fill: node-fill,
        stroke: node-stroke,
      )
        + eval(ev, scope: scope + (vid: i), mode: "code")
    )
    node-elements.push(cetz.draw.circle(_point(pos), name: "n" + str(i), ..node-style))
    node-outsets.push(_node-outset(node-style, node-outset))

    if v.name != none {
      node-elements.push(cetz.draw.content(
        _point(pos),
        [#v.name],
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
      elements.push(curve-api.cetz-bezier(halves.source, ..source-draw-style))
      elements.push(curve-api.cetz-bezier(halves.sink, ..sink-draw-style))
    } else if start != none {
      let line-start = curve-api.outset-point(
        nodes.at(start.node).pos,
        e.pos,
        distance: node-outsets.at(start.node),
      )
      if source-in-subgraph and subgraph-edge-underlay {
        elements.push(_line(line-start, e.pos, subgraph-edge-style))
      }
      elements.push(_line(line-start, e.pos, source-draw-style))
    } else if end != none {
      let line-end = curve-api.outset-point(
        nodes.at(end.node).pos,
        e.pos,
        distance: node-outsets.at(end.node),
      )
      if sink-in-subgraph and subgraph-edge-underlay {
        elements.push(_line(e.pos, line-end, subgraph-edge-style))
      }
      elements.push(_line(e.pos, line-end, sink-draw-style))
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

  let diagram-content = cetz.canvas(
    length: _canvas-length(unit),
    debug: debug-level >= 1,
    padding: padding,
    {
      for element in elements {
        element
      }
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
