// Internal drawing implementation. Public users should import `draw.typ`.

#import "@preview/cetz:0.5.1" as cetz
#import "../curve.typ" as curve-api
#import "../graph.typ" as graph-api
#import "../subgraph.typ" as subgraph-api

#let _point-x(p) = if type(p) == array { p.at(0) } else { p.x }
#let _point-y(p) = if type(p) == array { p.at(1) } else { p.y }
#let _point(p) = (_point-x(p), _point-y(p))

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

#let _edge-geometry-defaults = (
  offset: 0,
  length: none,
  ratio: none,
  resolve-length: "min",
  accuracy: 0.001,
  optimize: true,
  offset-side: none,
)

#let _mark-defaults = (
  mark-position: "end",
  mark-orientation: "path",
  mark-direction: "forward",
)

#let _pattern-defaults = (
  pattern: none,
  pattern-amplitude: 0.1,
  pattern-wavelength: 1.0,
  pattern-phase: 0,
  pattern-samples-per-period: 16,
  pattern-coil-longitudinal-scale: 1.25,
  pattern-accuracy: 0.001,
)

#let _without-keys(style, keys) = {
  let clean = style
  for key in keys {
    if clean.keys().contains(key) {
      let _ = clean.remove(key)
    }
  }
  clean
}

#let _without-pattern-style(style) = _without-keys(style, _pattern-defaults.keys())

#let _draw-style(style) = {
  let clean = _without-pattern-style(style)
  clean = _without-keys(clean, _edge-geometry-defaults.keys())
  _without-keys(clean, _mark-defaults.keys())
}

#let _without-mark-style(style) = {
  let clean = style
  if clean.keys().contains("mark") {
    let _ = clean.remove("mark")
  }
  clean
}

#let _start-mark-entry(entry) = if type(entry) == dictionary {
  entry
} else {
  (symbol: entry)
}

#let _backward-mark(mark) = if mark == none {
  none
} else if type(mark) == dictionary {
  let clean = mark
  let entry = none
  if clean.keys().contains("end") {
    entry = clean.end
    let _ = clean.remove("end")
  } else if clean.keys().contains("symbol") {
    entry = clean.symbol
    let _ = clean.remove("symbol")
  } else if clean.keys().contains("start") {
    entry = clean.start
    let _ = clean.remove("start")
  }
  if entry == none {
    mark
  } else {
    clean + (start: _start-mark-entry(entry))
  }
} else {
  (start: _start-mark-entry(mark))
}

#let _style-value(style, key) = {
  if _edge-geometry-defaults.keys().contains(key) {
    style.at(key, default: _edge-geometry-defaults.at(key))
  } else if _pattern-defaults.keys().contains(key) {
    style.at(key, default: _pattern-defaults.at(key))
  } else if _mark-defaults.keys().contains(key) {
    style.at(key, default: _mark-defaults.at(key))
  } else {
    style.at(key, default: none)
  }
}

#let _edge-geometry(style) = _edge-geometry-defaults + style

#let _pattern-style(style) = _pattern-defaults + style

#let _style-offset(style) = {
  _style-value(style, "offset")
}

#let _resolve-length-method(style) = {
  _style-value(style, "resolve-length")
}

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

#let _style-layers(value, data) = {
  let value = _call(value, data)
  if value == none {
    ((:),)
  } else if type(value) == array {
    value
  } else {
    (value,)
  }
}

#let _style-layer(layers, index) = {
  if index < layers.len() {
    layers.at(index)
  } else {
    none
  }
}

#let _as-content(value) = if value == none {
  none
} else if type(value) == str or type(value) == label {
  [#str(value)]
} else {
  value
}

#let _content(value, data, default) = {
  if value == auto {
    default
  } else {
    _as-content(_call(value, data))
  }
}

#let _data-label(record) = {
  let data = record.at("data", default: none)
  if type(data) == dictionary {
    data.at("label", default: none)
  } else {
    none
  }
}

#let _data-fields(record) = {
  let data = record.at("data", default: none)
  if type(data) == dictionary { data } else { (:) }
}

#let _pattern-name(style) = _style-value(style, "pattern")

#let _has-pattern(style) = {
  let pattern = _pattern-name(style)
  pattern != none and pattern != "normal" and pattern != "curve"
}

#let _has-mark(style) = style.at("mark", default: none) != none

#let _mark-position(style) = _style-value(style, "mark-position")

#let _mark-orientation(style) = _style-value(style, "mark-orientation")

#let _mark-direction(style) = _style-value(style, "mark-direction")

#let _dangling-mark-style(style) = {
  if _mark-position(style) == "center-if-dangling" {
    style + (mark-position: "center")
  } else {
    style
  }
}

#let _edge-has-source(edge) = edge.at("source-half-edge", default: none) != none

#let _edge-has-sink(edge) = edge.at("sink-half-edge", default: none) != none

#let _edge-oriented-mark-half(edge) = {
  let has-source = _edge-has-source(edge)
  let has-sink = _edge-has-sink(edge)
  if has-source and has-sink {
    if edge.at("orientation", default: "default") == "reversed" { "sink" } else { "source" }
  } else if has-source {
    "source"
  } else if has-sink {
    "sink"
  } else {
    none
  }
}

#let _orient-mark-style(style, edge, half) = {
  if style == none or not _has-mark(style) or _mark-orientation(style) != "edge" {
    style
  } else {
    let orientation = edge.at("orientation", default: "default")
    let oriented-half = _edge-oriented-mark-half(edge)
    if orientation == "undirected" or half != oriented-half {
      _without-mark-style(style)
    } else if orientation == "reversed" {
      style + (
        mark: _backward-mark(style.mark),
        mark-direction: "backward",
      )
    } else {
      style + (mark-direction: "forward")
    }
  }
}

#let _same-layer-geometry(source-style, sink-style) = {
  let same = _style-offset(source-style) == _style-offset(sink-style)
  same = (
    same and _style-value(source-style, "length") == _style-value(sink-style, "length")
  )
  same = same and _style-value(source-style, "ratio") == _style-value(sink-style, "ratio")
  same = same and _resolve-length-method(source-style) == _resolve-length-method(sink-style)
  same = (
    same and _style-value(source-style, "accuracy") == _style-value(sink-style, "accuracy")
  )
  same = (
    same and _style-value(source-style, "optimize") == _style-value(sink-style, "optimize")
  )
  same = same and _style-value(source-style, "offset-side") == _style-value(sink-style, "offset-side")
  same
}

#let _same-pattern-geometry(source-style, sink-style) = {
  let same = _has-pattern(source-style)
  same = same and _pattern-name(source-style) == _pattern-name(sink-style)
  same = same and _same-layer-geometry(source-style, sink-style)
  same = (
    same and _style-value(source-style, "pattern-amplitude") == _style-value(sink-style, "pattern-amplitude")
  )
  same = (
    same
      and _style-value(source-style, "pattern-wavelength") == _style-value(sink-style, "pattern-wavelength")
  )
  same = same and _style-value(source-style, "pattern-phase") == _style-value(sink-style, "pattern-phase")
  same = (
    same
      and _style-value(source-style, "pattern-samples-per-period")
        == _style-value(sink-style, "pattern-samples-per-period")
  )
  same = (
    same
      and _style-value(source-style, "pattern-coil-longitudinal-scale")
        == _style-value(sink-style, "pattern-coil-longitudinal-scale")
  )
  same = (
    same and _style-value(source-style, "pattern-accuracy") == _style-value(sink-style, "pattern-accuracy")
  )
  same
}

#let _center-outset(base-length, style, source-outset, sink-outset) = {
  let geometry = _edge-geometry(style)
  curve-api.center-outset(
    base-length,
    length: geometry.length,
    ratio: geometry.ratio,
    resolve-length: geometry.resolve-length,
    start-outset: source-outset,
    end-outset: sink-outset,
  )
}

#let _clamp(value, low, high) = calc.min(high, calc.max(low, value))

#let _segment-length(segment, accuracy) = {
  curve-api.length(curve-api.from-cubic(segment), accuracy: accuracy)
}

#let _visible-half-outsets(
  base-length,
  half-start,
  half-length,
  style,
  source-outset,
  sink-outset,
) = {
  let center-outset = _center-outset(
    base-length,
    style,
    source-outset,
    sink-outset,
  )
  let visible-start = source-outset + center-outset
  let visible-end = base-length - sink-outset - center-outset
  let local-start = _clamp(visible-start - half-start, 0, half-length)
  let local-end = _clamp(visible-end - half-start, 0, half-length)
  if local-end <= local-start {
    none
  } else {
    (
      start: local-start,
      end: half-length - local-end,
    )
  }
}

#let _layer(segment, style, start-outset, end-outset, label-pos, center-outset) = {
  let geometry = _edge-geometry(style)
  let side-point = if geometry.offset-side == "label" { label-pos } else { none }
  let length = geometry.length
  let ratio = geometry.ratio
  if center-outset != auto {
    start-outset = start-outset + center-outset
    end-outset = end-outset + center-outset
    length = none
    ratio = none
  }
  curve-api.layer(
    curve-api.from-cubic(segment),
    offset: geometry.offset,
    length: length,
    ratio: ratio,
    resolve-length: geometry.resolve-length,
    start-outset: start-outset,
    end-outset: end-outset,
    side-point: if side-point == none { none } else { _point(side-point) },
    accuracy: geometry.accuracy,
    optimize: geometry.optimize,
  )
}

#let _geometry-segments(
  segment,
  style,
  start-outset,
  end-outset,
  label-pos,
  center-outset,
) = {
  curve-api.segments(_layer(
    segment,
    style,
    start-outset,
    end-outset,
    label-pos,
    center-outset,
  ))
}

#let _half-geometry(segment, style, outsets, label-pos) = {
  if outsets == none {
    ()
  } else {
    _geometry-segments(
      segment,
      style,
      outsets.start,
      outsets.end,
      label-pos,
      0,
    )
  }
}
#let edge-halves(edge, nodes, options) = {
  let omega = options.omega
  let source-outset = options.source-outset
  let sink-outset = options.sink-outset
  let accuracy = options.accuracy
  if edge.source == none or edge.sink == none {
    panic("edge-halves currently requires a paired edge with source and sink")
  }

  let points = (
    _point(nodes.at(edge.source.node).pos),
    _point(edge.pos),
    _point(nodes.at(edge.sink.node).pos),
  )
  let split = curve-api.split-through(
    points,
    omega: omega,
    start-outset: source-outset,
    end-outset: sink-outset,
    accuracy: accuracy,
  )

  (
    source: split.parts.at(0),
    sink: split.parts.at(1),
    curve: split.curve,
  )
}
#let to-cetz-edge-halves(
  edge,
  nodes,
  options,
) = {
  let halves = edge-halves(
    edge,
    nodes,
    (
      omega: options.omega,
      source-outset: options.source-outset,
      sink-outset: options.sink-outset,
      accuracy: options.accuracy,
    ),
  )
  curve-api.to-cetz(halves.source, unit: options.unit, ..options.source-style)
  curve-api.to-cetz(halves.sink, unit: options.unit, ..options.sink-style)
}

#let _edge-geometry-halves(
  edge,
  nodes,
  source-style,
  sink-style,
  options,
) = {
  let omega = options.omega
  let source-outset = options.source-outset
  let sink-outset = options.sink-outset
  let accuracy = options.accuracy
  let label-pos = options.label-pos
  let base = edge-halves(edge, nodes, (
    omega: omega,
    source-outset: 0,
    sink-outset: 0,
    accuracy: accuracy,
  ))
  let base-segments = curve-api.segments(base.curve)
  let source-segment = base-segments.at(0)
  let sink-segment = base-segments.at(1)
  let source-length = _segment-length(source-segment, accuracy)
  let sink-length = _segment-length(sink-segment, accuracy)
  let base-length = source-length + sink-length
  let source-outsets = _visible-half-outsets(
    base-length,
    0,
    source-length,
    source-style,
    source-outset,
    sink-outset,
  )
  let source-geometry = _half-geometry(source-segment, source-style, source-outsets, label-pos)
  let sink-geometry = if _same-layer-geometry(source-style, sink-style) {
    let sink-outsets = _visible-half-outsets(
      base-length,
      source-length,
      sink-length,
      source-style,
      source-outset,
      sink-outset,
    )
    _half-geometry(sink-segment, source-style, sink-outsets, label-pos)
  } else {
    let sink-outsets = _visible-half-outsets(
      base-length,
      source-length,
      sink-length,
      sink-style,
      source-outset,
      sink-outset,
    )
    _half-geometry(sink-segment, sink-style, sink-outsets, label-pos)
  }
  (
    source: source-geometry,
    sink: sink-geometry,
    curve: base.curve,
  )
}

#let _pattern(segment, style, phase, anchor-start, anchor-end) = {
  let pattern-style = _pattern-style(style)
  curve-api.pattern(
    curve-api.from-cubic(segment),
    pattern: pattern-style.pattern,
    amplitude: pattern-style.pattern-amplitude,
    wavelength: pattern-style.pattern-wavelength,
    phase: if phase == auto { pattern-style.pattern-phase } else { phase },
    samples-per-period: pattern-style.pattern-samples-per-period,
    coil-longitudinal-scale: pattern-style.pattern-coil-longitudinal-scale,
    anchor-start: anchor-start,
    anchor-end: anchor-end,
    accuracy: pattern-style.pattern-accuracy,
  )
}

#let _bezier-element(segment, style) = {
  cetz.draw.bezier(
    _point(segment.start),
    _point(segment.end),
    _point(segment.control-start),
    _point(segment.control-end),
    .._draw-style(style),
  )
}

#let _segments-elements(segments, style, phase, anchor-start, anchor-end) = {
  let elements = ()
  let length = 0
  if _has-pattern(style) {
    let current-phase = if phase == auto { _style-value(style, "pattern-phase") } else { phase }
    let wavelength = _style-value(style, "pattern-wavelength")
    for (index, segment) in segments.enumerate() {
      let piece = _pattern(
        segment,
        style,
        current-phase,
        anchor-start and index == 0,
        anchor-end and index == segments.len() - 1,
      )
      elements.push(curve-api.to-cetz(piece, .._draw-style(style)))
      let piece-length = _segment-length(segment, _style-value(style, "pattern-accuracy"))
      length = length + piece-length
      current-phase = current-phase + 2 * calc.pi * piece-length / wavelength
    }
  } else {
    let mark-index = if _mark-direction(style) == "backward" { 0 } else { segments.len() - 1 }
    for (index, segment) in segments.enumerate() {
      if _has-mark(style) and index == mark-index {
        elements.push(_bezier-element(segment, style))
      } else if _has-mark(style) {
        elements.push(_bezier-element(segment, _without-mark-style(style)))
      } else {
        elements.push(curve-api.to-cetz(curve-api.from-cubic(segment), .._draw-style(style)))
      }
    }
  }
  (elements: elements, length: length)
}

#let _segment-elements(segment, style, phase, anchor-start, anchor-end, label-pos) = {
  let path = _layer(segment, style, 0, 0, label-pos, auto)
  if _has-mark(style) and _mark-position(style) == "center" and not _has-pattern(style) {
    let length = curve-api.length(path, accuracy: _style-value(style, "accuracy"))
    let first = curve-api.trim(path, end-outset: length / 2, accuracy: _style-value(style, "accuracy"))
    let second = curve-api.trim(path, start-outset: length / 2, accuracy: _style-value(style, "accuracy"))
    let backward = _mark-direction(style) == "backward"
    let first-style = if backward { _without-mark-style(style) } else { style }
    let second-style = if backward { style } else { _without-mark-style(style) }
    let first-elements = _segments-elements(
      curve-api.segments(first),
      first-style,
      phase,
      anchor-start,
      false,
    ).elements
    let second-elements = _segments-elements(
      curve-api.segments(second),
      second-style,
      phase,
      false,
      anchor-end,
    ).elements
    let pieces = ()
    if backward {
      for element in first-elements {
        pieces.push(element)
      }
      for element in second-elements {
        pieces.push(element)
      }
    } else {
      for element in second-elements {
        pieces.push(element)
      }
      for element in first-elements {
        pieces.push(element)
      }
    }
    (elements: pieces, length: length)
  } else {
    _segments-elements(
      curve-api.segments(path),
      style,
      phase,
      anchor-start,
      anchor-end,
    )
  }
}

#let _pattern-edge-halves(halves, source-style, sink-style) = {
  let elements = ()
  if _same-pattern-geometry(source-style, sink-style) {
    let source = _segments-elements(halves.source, source-style, auto, true, false)
    let wavelength = _style-value(source-style, "pattern-wavelength")
    let phase = _style-value(source-style, "pattern-phase")
    let sink-phase = phase + 2 * calc.pi * source.length / wavelength
    let sink-elements = _segments-elements(halves.sink, sink-style, sink-phase, false, true).elements
    if _has-mark(source-style) and not _has-mark(sink-style) {
      for element in sink-elements {
        elements.push(element)
      }
      for element in source.elements {
        elements.push(element)
      }
    } else {
      for element in source.elements {
        elements.push(element)
      }
      for element in sink-elements {
        elements.push(element)
      }
    }
  } else {
    let source-elements = _segments-elements(halves.source, source-style, auto, true, true).elements
    let sink-elements = _segments-elements(halves.sink, sink-style, auto, true, true).elements
    if _has-mark(source-style) and not _has-mark(sink-style) {
      for element in sink-elements {
        elements.push(element)
      }
      for element in source-elements {
        elements.push(element)
      }
    } else {
      for element in source-elements {
        elements.push(element)
      }
      for element in sink-elements {
        elements.push(element)
      }
    }
  }
  elements
}

#let _pattern-line(start, end, style, label-pos) = {
  _segment-elements(curve-api.line-segment(_point(start), _point(end)), style, auto, true, true, label-pos).elements
}

#let _bent-line-segment(start, end, bend) = {
  let start = _point(start)
  let end = _point(end)
  if bend == none or bend == 0 {
    return curve-api.line-segment(start, end)
  }

  let dx = _point-x(end) - _point-x(start)
  let dy = _point-y(end) - _point-y(start)
  let length = calc.sqrt(dx * dx + dy * dy)
  if length == 0 {
    return curve-api.line-segment(start, end)
  }

  let amount = calc.sin(bend) * length * 0.5
  let control = (
    (_point-x(start) + _point-x(end)) / 2 - dy / length * amount,
    (_point-y(start) + _point-y(end)) / 2 + dx / length * amount,
  )
  curve-api.segments(curve-api.quad(start, control, end)).first()
}

#let _pattern-dangling(start, end, bend, style, label-pos) = {
  _segment-elements(_bent-line-segment(start, end, bend), style, auto, true, true, label-pos).elements
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

#let _in-subgraph(hedges, half-edge) = {
  hedges != none and half-edge != none and hedges.contains(half-edge.hedge)
}

#let _node-pos(node) = {
  if node.pos == none {
    panic("draw: graph node has no position; call layout(...) or pass pos: to graph.node/build")
  }
  node.pos
}

#let _midpoint(a, b) = ((_point-x(a) + _point-x(b)) / 2, (_point-y(a) + _point-y(b)) / 2)

#let _edge-pos(edge, nodes) = {
  if edge.pos != none {
    edge.pos
  } else if edge.source != none and edge.sink != none {
    _midpoint(_node-pos(nodes.at(edge.source.node)), _node-pos(nodes.at(edge.sink.node)))
  } else {
    panic("draw: dangling graph edge has no position; call layout(...) or pass pos: to graph.edge/build")
  }
}
#let draw(graph, options) = {
  let scope = options.scope
  let unit = options.unit
  let title = options.title
  let subgraph = options.subgraph
  let debug = options.debug
  let node-radius = options.node-radius
  let node-min-radius = options.node-min-radius
  let node-label-padding = options.node-label-padding
  let node-fill = options.node-fill
  let node-stroke = options.node-stroke
  let node-outset = options.node-outset
  let node-label-style = options.node-label-style
  let node-style = options.node-style
  let node-label = options.node-label
  let edge-stroke = options.edge-stroke
  let edge-offset = options.edge-offset
  let edge-length = options.edge-length
  let edge-ratio = options.edge-ratio
  let edge-resolve-length = options.edge-resolve-length
  let edge-accuracy = options.edge-accuracy
  let edge-optimize = options.edge-optimize
  let source-style = options.source-style
  let sink-style = options.sink-style
  let edge-label = options.edge-label
  let edge-label-style = options.edge-label-style
  let edge-omega = options.edge-omega
  let edge-trim-accuracy = options.edge-trim-accuracy
  let padding = options.padding
  let debug-edge-radius = options.debug-edge-radius
  let debug-edge-fill = options.debug-edge-fill
  let debug-edge-stroke = options.debug-edge-stroke
  let debug-edge-label-fill = options.debug-edge-label-fill
  let subgraph-edge-style = options.subgraph-edge-style
  let subgraph-edge-underlay = options.subgraph-edge-underlay
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
          let pos = _node-pos(v)
          let node = v + (pos: pos)
          let record-label = _data-label(v)
          let statement-label = v.statements.at("label", default: none)
          let data-label = if record-label == none { statement-label } else { record-label }
          let node-data = (
            scope
              + v.statements
              + _data-fields(v)
              + (
                vid: i,
                node: node,
                name: node.name,
                data: v.at("data", default: none),
                label: data-label,
              )
          )
          let default-label-value = node-data.at("label", default: node.name)
          let default-label = _as-content(default-label-value)
          let label = _content(node-label, node-data, default-label)
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
          let edge-pos = _edge-pos(e, nodes)
          let edge = e + (pos: edge-pos)
          let source-half-edge = edge.source
          let sink-half-edge = edge.sink
          let source-statement = if source-half-edge == none { none } else { source-half-edge.statement }
          let sink-statement = if sink-half-edge == none { none } else { sink-half-edge.statement }
          let ext = source-half-edge == none or sink-half-edge == none
          let data = edge.statements
          let record-label = _data-label(edge)
          let statement-label = data.at("label", default: none)
          let data-label = if record-label == none { statement-label } else { record-label }
          let orientation = edge.orientation
          let label-pos = if edge.label-pos == none { edge.pos } else { edge.label-pos }
          let edge-data = (
            scope
              + data
              + _data-fields(edge)
              + (
                eid: edge.edge,
                edge: edge,
                source-statement: source-statement,
                sink-statement: sink-statement,
                source-half-edge: source-half-edge,
                sink-half-edge: sink-half-edge,
                data: edge.at("data", default: none),
                label: data-label,
                label-pos: label-pos,
                orientation: orientation,
                ext: ext,
              )
          )
          let source-in-subgraph = _in-subgraph(subgraph-hedges, source-half-edge)
          let sink-in-subgraph = _in-subgraph(subgraph-hedges, sink-half-edge)

          let geometry-style = _edge-geometry-defaults + (
            offset: edge-offset,
            length: edge-length,
            ratio: edge-ratio,
            resolve-length: edge-resolve-length,
            accuracy: edge-accuracy,
            optimize: edge-optimize,
          )
          let source-style-layers = _style-layers(source-style, edge-data)
          let sink-style-layers = _style-layers(sink-style, edge-data)
          let layer-count = calc.max(source-style-layers.len(), sink-style-layers.len())
          let ev-label = _content(edge-label, edge-data, _as-content(data-label))
          let edge-label-style = _style(edge-label-style, edge-data)

          for layer-index in range(0, layer-count) {
            let source-layer = _style-layer(source-style-layers, layer-index)
            let sink-layer = _style-layer(sink-style-layers, layer-index)
            let source-style-value = if source-layer == none {
              none
            } else {
              (stroke: edge-stroke) + geometry-style + source-layer
            }
            let sink-style-value = if sink-layer == none {
              none
            } else {
              (stroke: edge-stroke) + geometry-style + sink-layer
            }
            let source-draw-style = if source-style-value == none {
              none
            } else if source-in-subgraph and not subgraph-edge-underlay {
              source-style-value + subgraph-edge-style
            } else {
              source-style-value
            }
            source-draw-style = _orient-mark-style(source-draw-style, edge-data, "source")
            let sink-draw-style = if sink-style-value == none {
              none
            } else if sink-in-subgraph and not subgraph-edge-underlay {
              sink-style-value + subgraph-edge-style
            } else {
              sink-style-value
            }
            sink-draw-style = _orient-mark-style(sink-draw-style, edge-data, "sink")

            if source-style-value == none and sink-style-value == none {
              ()
            } else if source-half-edge != none and sink-half-edge != none {
              let source-geometry-style = if source-style-value == none { sink-style-value } else { source-style-value }
              let sink-geometry-style = if sink-style-value == none { source-style-value } else { sink-style-value }
              let halves = _edge-geometry-halves(
                edge,
                nodes,
                source-geometry-style,
                sink-geometry-style,
                (
                  omega: edge-omega,
                  source-outset: node-outsets.at(source-half-edge.node),
                  sink-outset: node-outsets.at(sink-half-edge.node),
                  accuracy: edge-trim-accuracy,
                  label-pos: label-pos,
                ),
              )
              if source-style-value != none and source-in-subgraph and subgraph-edge-underlay {
                for element in (
                  _segments-elements(
                    halves.source,
                    _without-mark-style(_without-pattern-style(source-style-value)) + subgraph-edge-style,
                    auto,
                    true,
                    true,
                  ).elements
                ) {
                  elements.push(element)
                }
              }
              if sink-style-value != none and sink-in-subgraph and subgraph-edge-underlay {
                for element in (
                  _segments-elements(
                    halves.sink,
                    _without-mark-style(_without-pattern-style(sink-style-value)) + subgraph-edge-style,
                    auto,
                    true,
                    true,
                  ).elements
                ) {
                  elements.push(element)
                }
              }
              if source-style-value != none and sink-style-value != none {
                for element in _pattern-edge-halves(halves, source-draw-style, sink-draw-style) {
                  elements.push(element)
                }
              } else if source-style-value != none {
                for element in _segments-elements(halves.source, source-draw-style, auto, true, true).elements {
                  elements.push(element)
                }
              } else if sink-style-value != none {
                for element in _segments-elements(halves.sink, sink-draw-style, auto, true, true).elements {
                  elements.push(element)
                }
              }
	            } else if source-half-edge != none {
	              let line-start = curve-api.outset-point(
	                _point(_node-pos(nodes.at(source-half-edge.node))),
	                _point(edge.pos),
	                distance: node-outsets.at(source-half-edge.node),
	              )
              let bend = edge.at("bend", default: none)
              if source-style-value != none and source-in-subgraph and subgraph-edge-underlay {
                for element in _pattern-dangling(
                  line-start,
                  edge.pos,
                  bend,
                  _without-mark-style(_without-pattern-style(_dangling-mark-style(source-style-value))) + subgraph-edge-style,
                  label-pos,
                ) {
                  elements.push(element)
                }
              }
              if source-style-value != none {
                for element in _pattern-dangling(line-start, edge.pos, bend, _dangling-mark-style(source-draw-style), label-pos) {
                  elements.push(element)
                }
              }
	            } else if sink-half-edge != none {
	              let line-end = curve-api.outset-point(
	                _point(_node-pos(nodes.at(sink-half-edge.node))),
	                _point(edge.pos),
	                distance: node-outsets.at(sink-half-edge.node),
	              )
              let bend = edge.at("bend", default: none)
              if sink-style-value != none and sink-in-subgraph and subgraph-edge-underlay {
                for element in _pattern-dangling(
                  edge.pos,
                  line-end,
                  bend,
                  _without-mark-style(_without-pattern-style(_dangling-mark-style(sink-style-value))) + subgraph-edge-style,
                  label-pos,
                ) {
                  elements.push(element)
                }
              }
              if sink-style-value != none {
                for element in _pattern-dangling(edge.pos, line-end, bend, _dangling-mark-style(sink-draw-style), label-pos) {
                  elements.push(element)
                }
              }
            }
          }

          if ev-label != none {
            elements.push(cetz.draw.content(_point(label-pos), ev-label, padding: 0, ..edge-label-style))
          }

          if debug-level >= 2 {
            elements.push(cetz.draw.circle(
              _point(edge.pos),
              radius: debug-edge-radius,
              fill: debug-edge-fill,
              stroke: debug-edge-stroke,
            ))
            elements.push(cetz.draw.content(
              _point(edge.pos),
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
