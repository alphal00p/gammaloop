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

#let _content(value, data, default: none) = {
  if value == auto {
    default
  } else {
    let value = _call(value, data)
    if type(value) == str { [#value] } else { value }
  }
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

#let _center-outset(base-length, style, source-outset: 0, sink-outset: 0) = {
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

#let _segment-length(segment, accuracy: 0.001) = {
  curve-api.path-length(curve-api.cubic-path(..segment, accuracy: accuracy), accuracy: accuracy)
}

#let _visible-half-outsets(
  base-length,
  half-start,
  half-length,
  style,
  source-outset: 0,
  sink-outset: 0,
) = {
  let center-outset = _center-outset(
    base-length,
    style,
    source-outset: source-outset,
    sink-outset: sink-outset,
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

#let _layer-path(segment, style, start-outset: 0, end-outset: 0, label-pos: none, center-outset: auto) = {
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
  curve-api.layer-path(
    curve-api.cubic-path(..segment, accuracy: geometry.accuracy),
    offset: geometry.offset,
    length: length,
    ratio: ratio,
    resolve-length: geometry.resolve-length,
    start-outset: start-outset,
    end-outset: end-outset,
    side-point: side-point,
    accuracy: geometry.accuracy,
    optimize: geometry.optimize,
  )
}

#let _geometry-segments(
  segment,
  style,
  start-outset: 0,
  end-outset: 0,
  accuracy: 0.001,
  label-pos: none,
  center-outset: auto,
) = {
  curve-api.path-segments(_layer-path(
    segment,
    style,
    start-outset: start-outset,
    end-outset: end-outset,
    label-pos: label-pos,
    center-outset: center-outset,
  ))
}

#let _half-geometry(segment, style, outsets, accuracy: 0.001, label-pos: none) = {
  if outsets == none {
    ()
  } else {
    _geometry-segments(
      segment,
      style,
      start-outset: outsets.start,
      end-outset: outsets.end,
      accuracy: accuracy,
      label-pos: label-pos,
      center-outset: 0,
    )
  }
}

#let _edge-geometry-halves(
  edge,
  nodes,
  source-style,
  sink-style,
  omega: 1.0,
  source-outset: 0,
  sink-outset: 0,
  accuracy: 0.001,
  label-pos: none,
) = {
  let base = curve-api.edge-halves(edge, nodes, omega: omega, accuracy: accuracy)
  let source-segment = base.curve.curves.at(0)
  let sink-segment = base.curve.curves.at(1)
  let source-length = _segment-length(source-segment, accuracy: accuracy)
  let sink-length = _segment-length(sink-segment, accuracy: accuracy)
  let base-length = source-length + sink-length
  let source-outsets = _visible-half-outsets(
    base-length,
    0,
    source-length,
    source-style,
    source-outset: source-outset,
    sink-outset: sink-outset,
  )
  let source-geometry = _half-geometry(source-segment, source-style, source-outsets, accuracy: accuracy, label-pos: label-pos)
  let sink-geometry = if _same-layer-geometry(source-style, sink-style) {
    let sink-outsets = _visible-half-outsets(
      base-length,
      source-length,
      sink-length,
      source-style,
      source-outset: source-outset,
      sink-outset: sink-outset,
    )
    _half-geometry(sink-segment, source-style, sink-outsets, accuracy: accuracy, label-pos: label-pos)
  } else {
    let sink-outsets = _visible-half-outsets(
      base-length,
      source-length,
      sink-length,
      sink-style,
      source-outset: source-outset,
      sink-outset: sink-outset,
    )
    _half-geometry(sink-segment, sink-style, sink-outsets, accuracy: accuracy, label-pos: label-pos)
  }
  (
    source: source-geometry,
    sink: sink-geometry,
    curve: base.curve,
  )
}

#let _pattern-path(segment, style, phase: auto, anchor-start: true, anchor-end: true) = {
  let pattern-style = _pattern-style(style)
  curve-api.pattern-path(
    curve-api.cubic-path(..segment, accuracy: pattern-style.pattern-accuracy),
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
    _point(segment.ctrl-a),
    _point(segment.ctrl-b),
    .._draw-style(style),
  )
}

#let _segments-elements(segments, style, phase: auto, anchor-start: true, anchor-end: true) = {
  let elements = ()
  let length = 0
  if _has-pattern(style) {
    let current-phase = if phase == auto { _style-value(style, "pattern-phase") } else { phase }
    let wavelength = _style-value(style, "pattern-wavelength")
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
    let mark-index = if _mark-direction(style) == "backward" { 0 } else { segments.len() - 1 }
    for (index, segment) in segments.enumerate() {
      if _has-mark(style) and index == mark-index {
        elements.push(_bezier-element(segment, style))
      } else if _has-mark(style) {
        elements.push(_bezier-element(segment, _without-mark-style(style)))
      } else {
        elements.push(curve-api.cetz-path(curve-api.cubic-path(..segment), .._draw-style(style)))
      }
    }
  }
  (elements: elements, length: length)
}

#let _segment-elements(segment, style, phase: auto, anchor-start: true, anchor-end: true, label-pos: none) = {
  let path = _layer-path(segment, style, label-pos: label-pos)
  if _has-mark(style) and _mark-position(style) == "center" and not _has-pattern(style) {
    let length = curve-api.path-length(path, accuracy: _style-value(style, "accuracy"))
    let first = curve-api.trim-path(path, end-outset: length / 2, accuracy: _style-value(style, "accuracy"))
    let second = curve-api.trim-path(path, start-outset: length / 2, accuracy: _style-value(style, "accuracy"))
    let backward = _mark-direction(style) == "backward"
    let first-style = if backward { _without-mark-style(style) } else { style }
    let second-style = if backward { style } else { _without-mark-style(style) }
    let first-elements = _segments-elements(
      curve-api.path-segments(first),
      first-style,
      phase: phase,
      anchor-start: anchor-start,
      anchor-end: false,
    ).elements
    let second-elements = _segments-elements(
      curve-api.path-segments(second),
      second-style,
      phase: phase,
      anchor-start: false,
      anchor-end: anchor-end,
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
      curve-api.path-segments(path),
      style,
      phase: phase,
      anchor-start: anchor-start,
      anchor-end: anchor-end,
    )
  }
}

#let _pattern-edge-halves(halves, source-style, sink-style) = {
  let elements = ()
  if _same-pattern-geometry(source-style, sink-style) {
    let source = _segments-elements(halves.source, source-style, anchor-end: false)
    let wavelength = _style-value(source-style, "pattern-wavelength")
    let phase = _style-value(source-style, "pattern-phase")
    let sink-phase = phase + 2 * calc.pi * source.length / wavelength
    let sink-elements = _segments-elements(halves.sink, sink-style, phase: sink-phase, anchor-start: false).elements
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
    let source-elements = _segments-elements(halves.source, source-style).elements
    let sink-elements = _segments-elements(halves.sink, sink-style).elements
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

#let _pattern-line(start, end, style, label-pos: none) = {
  _segment-elements(curve-api.line-segment(start, end), style, label-pos: label-pos).elements
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

/// Draw a laid-out graph object with CeTZ.
///
/// The graph must already have positions, so call `layout` first. Node and edge
/// Typst style dictionaries or callbacks are evaluated and forwarded to CeTZ.
///
/// #example(`
/// #let b = graph.builder(name: "demo")
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
/// #let (node: a, builder: b) = graph.node(b, name: "a hi")
/// #let (node: c, builder: b) = graph.node(b, name: "c")
/// #let (node: d, builder: b) = graph.node(b, name: "d")
/// #let (node: e, builder: b) = graph.node(b, name: "e")
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
/// #let (node: pa, builder: p) = graph.node(p, name: "a", statements: (pin: "y:0"))
/// #let (node: pc, builder: p) = graph.node(p, name: "c", statements: (pin: "y:0"))
/// #let (node: pc1, builder: p) = graph.node(p, name: "c", statements: (pin: "y:0"))
/// #let (node: pc2, builder: p) = graph.node(p, name: "c", statements: (pin: "y:0"))
/// #let p = graph.edge(p, source: (node: pa), sink: (node: pc))
/// #let p = graph.edge(p, source: (node: pc2), sink: (node: pc1))
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
/// #draw(layout(graph.finish(p), g-center: 0.004, label-steps: 0,), unit: 1.25, source-style: focused-source-style, sink-style: focused-sink-style)
/// `,dir:ttb)
///
/// #example(`
/// #let b = graph.builder(name: "oriented marks")
/// #let (node: a, builder: b) = graph.node(b, name: "a", statements: (pin: "x:0,y:0"))
/// #let (node: c, builder: b) = graph.node(b, name: "c", statements: (pin: "x:3,y:0"))
/// #let (node: d, builder: b) = graph.node(b, name: "d", statements: (pin: "x:3,y:-1"))
/// #let b = graph.edge(b, source: (node: a), sink: (node: c), orientation: "default")
/// #let b = graph.edge(b, source: (node: a), sink: (node: d), orientation: "reversed")
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
/// #draw(layout(graph.finish(b)), unit: 1.4, source-style: oriented-arrow, sink-style: oriented-arrow)
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
  edge-offset: _edge-geometry-defaults.offset,
  /// Maximum visible arc length for centered parallel edge paths. `none` keeps
  /// the full shifted path. -> none | int | float
  edge-length: _edge-geometry-defaults.length,
  /// Maximum visible fraction of the base edge length for centered parallel edge
  /// paths. Combined with `edge-length` according to `edge-resolve-length`.
  /// -> none | int | float
  edge-ratio: _edge-geometry-defaults.ratio,
  /// Resolve `edge-length` and `edge-ratio`. Accepted string
  /// values are `"min"`/`"shorter"`, `"max"`/`"longer"`, `"length"`/`"fixed"`,
  /// `"ratio"`/`"relative"`, or `"none"`/`"full"`. A function receives
  /// `(base-length, length, ratio)`. -> string | function
  edge-resolve-length: _edge-geometry-defaults.resolve-length,
  /// Arc-length accuracy for fitted parallel edge paths. -> float
  edge-accuracy: _edge-geometry-defaults.accuracy,
  /// Let Kurbo optimize the fitted parallel path. -> bool
  edge-optimize: _edge-geometry-defaults.optimize,
  /// Source half-edge style dictionary, array of layer dictionaries, or
  /// callback. `mark-position: "center-if-dangling"` keeps an end marker at the
  /// paired-edge split point while centering it on dangling half edges.
  /// `mark-orientation: "edge"` makes a mark follow `edge.orientation` instead
  /// of raw path direction; reversed edges move the mark to the sink half and
  /// flip it, while undirected edges suppress it.
  /// -> dictionary | array | function | none
  source-style: (:),
  /// Sink half-edge style dictionary, array of layer dictionaries, or callback.
  /// A callback receives edge data. `mark-position: "center-if-dangling"` has
  /// the same dangling-edge behavior as for `source-style`, and
  /// `mark-orientation: "edge"` participates in the same orientation-aware mark
  /// placement.
  /// -> dictionary | array | function | none
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
          let node-data = (
            scope
              + v.statements
              + (
                vid: i,
                node: v,
                name: v.name,
              )
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
          let source-half-edge = e.source
          let sink-half-edge = e.sink
          let source-statement = if source-half-edge == none { none } else { source-half-edge.statement }
          let sink-statement = if sink-half-edge == none { none } else { sink-half-edge.statement }
          let ext = source-half-edge == none or sink-half-edge == none
          let data = e.statements
          let orientation = e.orientation
          let edge-data = (
            scope
              + data
              + (
                eid: i,
                edge: e,
                source-statement: source-statement,
                sink-statement: sink-statement,
                source-half-edge: source-half-edge,
                sink-half-edge: sink-half-edge,
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
          let ev-label = _content(edge-label, edge-data, default: none)
          let label-pos = if e.label-pos == none { e.pos } else { e.label-pos }

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
                e,
                nodes,
                source-geometry-style,
                sink-geometry-style,
                omega: edge-omega,
                source-outset: node-outsets.at(source-half-edge.node),
                sink-outset: node-outsets.at(sink-half-edge.node),
                accuracy: edge-trim-accuracy,
                label-pos: label-pos,
              )
              if source-style-value != none and source-in-subgraph and subgraph-edge-underlay {
                for element in (
                  _segments-elements(
                    halves.source,
                    _without-mark-style(_without-pattern-style(source-style-value)) + subgraph-edge-style,
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
                for element in _segments-elements(halves.source, source-draw-style).elements {
                  elements.push(element)
                }
              } else if sink-style-value != none {
                for element in _segments-elements(halves.sink, sink-draw-style).elements {
                  elements.push(element)
                }
              }
            } else if source-half-edge != none {
              let line-start = curve-api.outset-point(
                nodes.at(source-half-edge.node).pos,
                e.pos,
                distance: node-outsets.at(source-half-edge.node),
              )
              if source-style-value != none and source-in-subgraph and subgraph-edge-underlay {
                for element in _pattern-line(
                  line-start,
                  e.pos,
                  _without-mark-style(_without-pattern-style(_dangling-mark-style(source-style-value))) + subgraph-edge-style,
                  label-pos: label-pos,
                ) {
                  elements.push(element)
                }
              }
              if source-style-value != none {
                for element in _pattern-line(line-start, e.pos, _dangling-mark-style(source-draw-style), label-pos: label-pos) {
                  elements.push(element)
                }
              }
            } else if sink-half-edge != none {
              let line-end = curve-api.outset-point(
                nodes.at(sink-half-edge.node).pos,
                e.pos,
                distance: node-outsets.at(sink-half-edge.node),
              )
              if sink-style-value != none and sink-in-subgraph and subgraph-edge-underlay {
                for element in _pattern-line(
                  e.pos,
                  line-end,
                  _without-mark-style(_without-pattern-style(_dangling-mark-style(sink-style-value))) + subgraph-edge-style,
                  label-pos: label-pos,
                ) {
                  elements.push(element)
                }
              }
              if sink-style-value != none {
                for element in _pattern-line(e.pos, line-end, _dangling-mark-style(sink-draw-style), label-pos: label-pos) {
                  elements.push(element)
                }
              }
            }
          }

          if ev-label != none {
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
