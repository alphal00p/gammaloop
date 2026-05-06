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

#let _lerp(a, b, t) = a + (b - a) * t

#let _bezier-point(segment, t) = {
  let ab = (
    x: _lerp(segment.start.x, segment.ctrl-a.x, t),
    y: _lerp(segment.start.y, segment.ctrl-a.y, t),
  )
  let bc = (
    x: _lerp(segment.ctrl-a.x, segment.ctrl-b.x, t),
    y: _lerp(segment.ctrl-a.y, segment.ctrl-b.y, t),
  )
  let cd = (
    x: _lerp(segment.ctrl-b.x, segment.end.x, t),
    y: _lerp(segment.ctrl-b.y, segment.end.y, t),
  )
  let abc = (
    x: _lerp(ab.x, bc.x, t),
    y: _lerp(ab.y, bc.y, t),
  )
  let bcd = (
    x: _lerp(bc.x, cd.x, t),
    y: _lerp(bc.y, cd.y, t),
  )
  (
    x: _lerp(abc.x, bcd.x, t),
    y: _lerp(abc.y, bcd.y, t),
  )
}

#let _bezier-tangent(segment, t) = {
  let ab = (
    x: _lerp(segment.start.x, segment.ctrl-a.x, t),
    y: _lerp(segment.start.y, segment.ctrl-a.y, t),
  )
  let bc = (
    x: _lerp(segment.ctrl-a.x, segment.ctrl-b.x, t),
    y: _lerp(segment.ctrl-a.y, segment.ctrl-b.y, t),
  )
  let cd = (
    x: _lerp(segment.ctrl-b.x, segment.end.x, t),
    y: _lerp(segment.ctrl-b.y, segment.end.y, t),
  )
  let abc = (
    x: _lerp(ab.x, bc.x, t),
    y: _lerp(ab.y, bc.y, t),
  )
  let bcd = (
    x: _lerp(bc.x, cd.x, t),
    y: _lerp(bc.y, cd.y, t),
  )
  (
    x: bcd.x - abc.x,
    y: bcd.y - abc.y,
  )
}

#let _offset-on-label-side(style, segment, label-pos) = {
  let offset = style.at("offset", default: 0)
  if label-pos == none or offset == none or offset == 0 {
    offset
  } else {
    let mid = _bezier-point(segment, 0.5)
    let tangent = _bezier-tangent(segment, 0.5)
    let cross = tangent.x * (label-pos.y - mid.y) - tangent.y * (label-pos.x - mid.x)
    let sign = if cross < 0 { -1 } else { 1 }
    sign * calc.abs(offset)
  }
}

#let _mark-style-on-label-side(style, segment, label-pos) = {
  style + (offset: _offset-on-label-side(style, segment, label-pos))
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
  "offset",
  "length",
  "ratio",
  "resolve-length",
  "accuracy",
  "optimize",
)

#let _decoration-style-keys = ("parallel-mark",)

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
  for key in _decoration-style-keys {
    if clean.keys().contains(key) {
      let _ = clean.remove(key)
    }
  }
  clean
}

#let _without-mark-style(style) = {
  let clean = style
  if clean.keys().contains("mark") {
    let _ = clean.remove("mark")
  }
  clean
}

#let _style-value(style, key, default) = style.at(key, default: default)

#let _style-offset(style) = {
  style.at("offset", default: 0)
}

#let _resolve-length-method(style) = {
  style.at("resolve-length", default: "min")
}

#let _resolve-parallel-target(base-length, fixed, relative, method) = {
  if type(method) == function {
    method((base-length: base-length, length: fixed, ratio: relative))
  } else if fixed == none {
    relative
  } else if relative == none {
    fixed
  } else if method in ("max", "longer") {
    calc.max(fixed, relative)
  } else if method in ("length", "fixed") {
    fixed
  } else if method in ("ratio", "relative") {
    relative
  } else if method in ("none", "full") {
    none
  } else {
    calc.min(fixed, relative)
  }
}

#let _parallel-mark-style(style) = {
  let mark = if style == none { none } else { _style-value(style, "parallel-mark", none) }
  if mark == none {
    none
  } else {
    (
      stroke: mark.at("stroke", default: style.at("stroke", default: (paint: black, thickness: 0.45pt, cap: "round"))),
      mark: mark.at("mark", default: (end: ">")),
      offset: mark.at("offset", default: 0),
      length: mark.at("length", default: none),
      ratio: mark.at("ratio", default: none),
      resolve-length: mark.at("resolve-length", default: _resolve-length-method(style)),
      accuracy: mark.at("accuracy", default: _style-value(style, "accuracy", 0.001)),
      optimize: mark.at("optimize", default: _style-value(style, "optimize", true)),
    )
  }
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

#let _pattern-name(style) = _style-value(style, "pattern", none)

#let _has-pattern(style) = {
  let pattern = _pattern-name(style)
  pattern != none and pattern != "normal" and pattern != "curve"
}

#let _has-parallel(style) = {
  let offset = _style-offset(style)
  offset != none and offset != 0
}

#let _has-mark(style) = _style-value(style, "mark", none) != none

#let _same-parallel-geometry(source-style, sink-style) = {
  let same = _style-offset(source-style) == _style-offset(sink-style)
  same = (
    same and _style-value(source-style, "length", none) == _style-value(sink-style, "length", none)
  )
  same = same and _style-value(source-style, "ratio", none) == _style-value(sink-style, "ratio", none)
  same = same and _resolve-length-method(source-style) == _resolve-length-method(sink-style)
  same = (
    same
      and _style-value(source-style, "accuracy", 0.001) == _style-value(sink-style, "accuracy", 0.001)
  )
  same = (
    same
      and _style-value(source-style, "optimize", true) == _style-value(sink-style, "optimize", true)
  )
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
  curve-api.parallel-path(
    curve-api.cubic-path(..segment, accuracy: _style-value(style, "accuracy", 0.001)),
    distance: _style-offset(style),
    start-outset: start-outset,
    end-outset: end-outset,
    accuracy: _style-value(style, "accuracy", 0.001),
    optimize: _style-value(style, "optimize", true),
  )
}

#let _segment-length(segment, accuracy: 0.001) = {
  curve-api.cubic-path(..segment, accuracy: accuracy).length
}

#let _segments-length(segments, accuracy: 0.001) = {
  let length = 0
  for segment in segments {
    length = length + _segment-length(segment, accuracy: accuracy)
  }
  length
}

#let _parallel-target-length(base-length, style) = {
  let length = _style-value(style, "length", none)
  let ratio = _style-value(style, "ratio", none)
  let fixed = if length != none and length > 0 { length } else { none }
  let relative = if ratio != none and ratio > 0 { base-length * ratio } else { none }
  _resolve-parallel-target(base-length, fixed, relative, _resolve-length-method(style))
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
  } else if path.keys().contains("segments") and path.segments.len() > 0 {
    for segment in path.segments {
      segments.push(_line-segment(segment.start, segment.end))
    }
  } else if path.keys().contains("points") and path.points.len() > 1 {
    for i in range(0, path.points.len() - 1) {
      segments.push(_line-segment(path.points.at(i), path.points.at(i + 1)))
    }
  }
  segments
}

#let _trim-path-curves(segments, start-outset: 0, end-outset: 0, accuracy: 0.001) = {
  let trimmed = ()
  for (index, segment) in segments.enumerate() {
    let piece = segment
    if index == 0 and start-outset != 0 {
      piece = curve-api
        .trim-path(curve-api.cubic-path(..piece, accuracy: accuracy), start-outset: start-outset, accuracy: accuracy)
        .curves
        .at(0)
    }
    if index == segments.len() - 1 and end-outset != 0 {
      piece = curve-api
        .trim-path(curve-api.cubic-path(..piece, accuracy: accuracy), end-outset: end-outset, accuracy: accuracy)
        .curves
        .at(0)
    }
    trimmed.push(piece)
  }
  trimmed
}

#let _geometry-segments(segment, style, start-outset: 0, end-outset: 0, accuracy: 0.001) = {
  if _has-parallel(style) {
    _path-segments(_parallel-path(segment, style, start-outset: start-outset, end-outset: end-outset))
  } else {
    _trim-path-curves((segment,), start-outset: start-outset, end-outset: end-outset, accuracy: accuracy)
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
) = {
  if _has-parallel(source-style) or _has-parallel(sink-style) {
    let base = curve-api.edge-halves(edge, nodes, omega: omega, accuracy: accuracy)
    let base-length = _segments-length(base.curve.curves, accuracy: accuracy)
    let source-center-outset = _parallel-center-outset(
      base-length,
      source-style,
      source-outset: source-outset,
      sink-outset: sink-outset,
    )
    let source-geometry = _geometry-segments(
      base.curve.curves.at(0),
      source-style,
      start-outset: source-outset + source-center-outset,
      accuracy: accuracy,
    )
    let sink-geometry = if _same-parallel-geometry(source-style, sink-style) {
      _geometry-segments(
        base.curve.curves.at(1),
        source-style,
        end-outset: sink-outset + source-center-outset,
        accuracy: accuracy,
      )
    } else {
      let sink-center-outset = _parallel-center-outset(
        base-length,
        sink-style,
        source-outset: source-outset,
        sink-outset: sink-outset,
      )
      _geometry-segments(
        base.curve.curves.at(1),
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
      source: _path-segments(trimmed.source),
      sink: _path-segments(trimmed.sink),
      curve: trimmed.curve,
    )
  }
}

#let _pattern-path(segment, style, phase: auto, anchor-start: true, anchor-end: true) = {
  curve-api.pattern-path(
    curve-api.cubic-path(..segment, accuracy: _style-value(style, "pattern-accuracy", 0.001)),
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
    for (index, segment) in segments.enumerate() {
      if _has-mark(style) and index == segments.len() - 1 {
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

#let _segment-elements(segment, style, phase: auto, anchor-start: true, anchor-end: true) = {
  _segments-elements(
    _geometry-segments(segment, style),
    style,
    phase: phase,
    anchor-start: anchor-start,
    anchor-end: anchor-end,
  )
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

#let _parallel-mark-half-elements(
  segments,
  style,
  center-outset: 0,
  trim-start: false,
  trim-end: false,
  label-pos: none,
) = {
  let elements = ()
  if style != none and segments.len() > 0 {
    let geometry = ()
    for (index, segment) in segments.enumerate() {
      let start-outset = if trim-start and index == 0 { center-outset } else { 0 }
      let end-outset = if trim-end and index == segments.len() - 1 { center-outset } else { 0 }
      let piece-style = _mark-style-on-label-side(style, segment, label-pos)
      for piece in _geometry-segments(
        segment,
        piece-style,
        start-outset: start-outset,
        end-outset: end-outset,
        accuracy: _style-value(piece-style, "accuracy", 0.001),
      ) {
        geometry.push(piece)
      }
    }
    for element in _segments-elements(geometry, style).elements {
      elements.push(element)
    }
  }
  elements
}

#let _parallel-mark-edge-elements(halves, source-style, sink-style, label-pos: none) = {
  let source-mark = _parallel-mark-style(source-style)
  let sink-mark = _parallel-mark-style(sink-style)
  let mark-style = if sink-mark != none { sink-mark } else { source-mark }
  if mark-style == none {
    ()
  } else {
    let base-segments = halves.source + halves.sink
    let center-outset = _parallel-center-outset(
      _segments-length(base-segments, accuracy: _style-value(mark-style, "accuracy", 0.001)),
      mark-style,
    )
    _parallel-mark-half-elements(
      base-segments,
      mark-style,
      center-outset: center-outset,
      trim-start: true,
      trim-end: true,
      label-pos: label-pos,
    )
  }
}

#let _pattern-line(start, end, style) = {
  _segment-elements(_line-segment(start, end), style).elements
}

#let _parallel-mark-line-elements(start, end, style, label-pos: none) = {
  let mark-style = _parallel-mark-style(style)
  if mark-style == none {
    ()
  } else {
    let segment = _line-segment(start, end)
    let mark-style = _mark-style-on-label-side(mark-style, segment, label-pos)
    let center-outset = _parallel-center-outset(
      _segment-length(segment, accuracy: _style-value(mark-style, "accuracy", 0.001)),
      mark-style,
    )
    let geometry = _geometry-segments(
      segment,
      mark-style,
      start-outset: center-outset,
      end-outset: center-outset,
      accuracy: _style-value(mark-style, "accuracy", 0.001),
    )
    _segments-elements(geometry, mark-style).elements
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
/// #let p = graph.edge(p, source: (node: pa), sink: (node: pc))
/// #let parallel-edge-style(edge) = (
///     offset: 0.18,
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
/// #draw(layout(graph.finish(p),  label-steps: 0,), unit: 1.25, source-style: focused-source-style, sink-style: focused-sink-style)
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
  /// callback. A callback receives edge data. -> dictionary | array | function | none
  source-style: (:),
  /// Sink half-edge style dictionary, array of layer dictionaries, or callback.
  /// A callback receives edge data. -> dictionary | array | function | none
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
          let start = e.source
          let end = e.sink
          let source = if start == none { none } else { start.statement }
          let sink = if end == none { none } else { end.statement }
          let ext = start == none or end == none
          let data = e.statements
          let orientation = e.orientation
          let edge-data = (
            scope
              + data
              + (
                eid: i,
                edge: e,
                source: source,
                sink: sink,
                source-endpoint: start,
                sink-endpoint: end,
                orientation: orientation,
                ext: ext,
              )
          )
          let source-in-subgraph = _in-subgraph(subgraph-hedges, start)
          let sink-in-subgraph = _in-subgraph(subgraph-hedges, end)

          let geometry-style = (
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
            let sink-draw-style = if sink-style-value == none {
              none
            } else if sink-in-subgraph and not subgraph-edge-underlay {
              sink-style-value + subgraph-edge-style
            } else {
              sink-style-value
            }

            if source-style-value == none and sink-style-value == none {
              ()
            } else if start != none and end != none {
              let source-geometry-style = if source-style-value == none { sink-style-value } else { source-style-value }
              let sink-geometry-style = if sink-style-value == none { source-style-value } else { sink-style-value }
              let halves = _edge-geometry-halves(
                e,
                nodes,
                source-geometry-style,
                sink-geometry-style,
                omega: edge-omega,
                source-outset: node-outsets.at(start.node),
                sink-outset: node-outsets.at(end.node),
                accuracy: edge-trim-accuracy,
              )
              if source-style-value != none and source-in-subgraph and subgraph-edge-underlay {
                for element in (
                  _segments-elements(
                    halves.source,
                    _without-pattern-style(source-style-value) + subgraph-edge-style,
                  ).elements
                ) {
                  elements.push(element)
                }
              }
              if sink-style-value != none and sink-in-subgraph and subgraph-edge-underlay {
                for element in (
                  _segments-elements(
                    halves.sink,
                    _without-pattern-style(sink-style-value) + subgraph-edge-style,
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
              for element in _parallel-mark-edge-elements(
                halves,
                source-style-value,
                sink-style-value,
                label-pos: label-pos,
              ) {
                elements.push(element)
              }
            } else if start != none {
              let line-start = curve-api.outset-point(
                nodes.at(start.node).pos,
                e.pos,
                distance: node-outsets.at(start.node),
              )
              if source-style-value != none and source-in-subgraph and subgraph-edge-underlay {
                for element in _pattern-line(
                  line-start,
                  e.pos,
                  _without-pattern-style(source-style-value) + subgraph-edge-style,
                ) {
                  elements.push(element)
                }
              }
              if source-style-value != none {
                for element in _pattern-line(line-start, e.pos, source-draw-style) {
                  elements.push(element)
                }
                for element in _parallel-mark-line-elements(
                  line-start,
                  e.pos,
                  source-style-value,
                  label-pos: label-pos,
                ) {
                  elements.push(element)
                }
              }
            } else if end != none {
              let line-end = curve-api.outset-point(
                nodes.at(end.node).pos,
                e.pos,
                distance: node-outsets.at(end.node),
              )
              if sink-style-value != none and sink-in-subgraph and subgraph-edge-underlay {
                for element in _pattern-line(
                  e.pos,
                  line-end,
                  _without-pattern-style(sink-style-value) + subgraph-edge-style,
                ) {
                  elements.push(element)
                }
              }
              if sink-style-value != none {
                for element in _pattern-line(e.pos, line-end, sink-draw-style) {
                  elements.push(element)
                }
                for element in _parallel-mark-line-elements(e.pos, line-end, sink-style-value, label-pos: label-pos) {
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
