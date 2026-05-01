#import "@preview/fletcher:0.5.8" as fletcher: diagram, edge, node
#import "graph.typ" as graph-api

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

/// Draw a laid-out graph object with Fletcher.
///
/// The graph must already have positions, so call `layout` first. Node and edge
/// DOT `eval_*` statements are evaluated in `scope` and forwarded to Fletcher.
///
/// ```example
/// #let b = graph.builder(name: "demo")
/// #let (node: a, builder: b) = graph.node(b, name: "a")
/// #let (node: c, builder: b) = graph.node(b, name: "c")
/// #let (node: d, builder: b) = graph.node(b, name: "d")
/// #let b = graph.edge(b, source: (node: a), sink: (node: c))
///  #let layed-out = layout(graph.finish(b), seed: 2, steps: 2,g-center:500000,length-scale:0.01)
/// //#graph.info(layed-out)
///
/// #graph.nodes(layed-out)
/// //#draw(layed-out)
/// ```
/// -> content
#let draw(
  /// Laid-out graph object returned by `layout`. -> bytes
  graph,
  /// Additional symbols available while evaluating DOT `eval_*` statements.
  /// -> dictionary
  scope: (:),
  /// Coordinate scale applied to node and edge positions. -> int | float
  unit: 1,
  /// Optional title displayed above the diagram. Use `auto` for the graph name.
  /// -> none | auto | content | string
  title: none,
  /// Fletcher diagram debug flag. -> int
  debug: 0,
  /// Default Fletcher node shape. -> any
  node-shape: circle,
  /// Default Fletcher node fill. -> any
  node-fill: black,
  /// Default Fletcher edge stroke. -> any
  edge-stroke: 0.1em,
  /// Fletcher spacing option passed to `diagram`. -> any
  spacing: 2em,
) = {
  let info = graph-api.info(graph)
  let nodes = graph-api.nodes(graph)
  let edges = graph-api.edges(graph)
  let elements = ()
  let node-by-index = (:)

  for (i, v) in nodes.enumerate() {
    node-by-index.insert(str(i), v)
    let pos = v.pos
    let ev = v.eval
    if ev == none {
      ev = "(:)"
    }
    let node-label = if v.name == none { [] } else { [#v.name] }
    elements.push(node(
      pos: (pos.x * unit, pos.y * unit),
      name: label(str(i)),
      label: node-label,
      ..eval(ev, scope: scope + (vid: i), mode: "code"),
      layer: 2,
    ))
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

    let ev-sink = _eval-dict(data, "eval_sink", local-scope)
    let ev-label = _eval-content(
      data,
      "eval_label",
      local-scope + (sink: sink) + (source: source) + data,
    )
    let ev-source = _eval-dict(data, "eval_source", local-scope)

    let bend = e.bend
    if bend == none {
      bend = 0.
    }

    let end-marker = label("em" + str(i))
    let (end-node, end-node-pos) = if end != none {
      let nodelab = label(str(end.node))
      if start != none and nodelab == label(str(start.node)) {
        bend = bend + 2
      }
      elements.push(edge(
        vertices: ((e.pos.x * unit, e.pos.y * unit), nodelab),
        bend: bend * 0.5rad,
        ..ev-sink,
      ))
      (nodelab, node-by-index.at(str(end.node)).pos)
    } else {
      let lab = label("exte" + str(i))
      elements.push(node(
        (e.pos.x * unit, e.pos.y * unit),
        name: lab,
        outset: -5mm,
        radius: 5mm,
        fill: none,
      ))
      (lab, e.pos)
    }

    let start-marker = label("sm" + str(i))
    let (start-node, start-node-pos) = if start != none {
      let nodelab = label(str(start.node))
      elements.push(edge(
        vertices: (nodelab, (e.pos.x * unit, e.pos.y * unit)),
        bend: bend * 0.5rad,
        ..ev-source,
      ))
      (nodelab, node-by-index.at(str(start.node)).pos)
    } else {
      let lab = label("exts" + str(i))
      elements.push(node(
        (e.pos.x * unit, e.pos.y * unit),
        name: lab,
        outset: -5mm,
        radius: 5mm,
        fill: none,
      ))
      (lab, e.pos)
    }

    let bend-scale = 1 + calc.abs(bend / calc.pi)
    let marker-outset = (
      calc.sqrt(
        calc.pow(start-node-pos.x - end-node-pos.x, 2) + calc.pow(start-node-pos.y - end-node-pos.y, 2),
      )
        * 2.5
        * bend-scale
        * unit
    )

    elements.push(node(
      pos: (start-node-pos.x * unit, start-node-pos.y * unit),
      name: start-marker,
      label: [],
      outset: marker-outset * 1em,
    ))
    elements.push(node(
      pos: (end-node-pos.x * unit, end-node-pos.y * unit),
      name: end-marker,
      label: [],
      outset: marker-outset * 1em,
    ))

    let label-pos = if e.label_pos == none { e.pos } else { e.label_pos }
    elements.push(node(
      (label-pos.x * unit, label-pos.y * unit),
      ev-label,
      inset: 0mm,
      snap: false,
      name: label("e" + str(i)),
      fill: none,
    ))
  }

  let diagram-content = diagram(
    debug: debug,
    node-shape: node-shape,
    node-fill: node-fill,
    edge-stroke: edge-stroke,
    spacing: spacing,
    ..elements,
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
