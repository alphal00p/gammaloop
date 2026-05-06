#import "@preview/cetz:0.3.4" as cetz

#let _plugin = plugin("kurvst.wasm")

#let _point(p, unit: 1) = (p.x * unit, p.y * unit)

/// Default geometry options for derived path layers.
///
/// These defaults are shared by @layer-path and drawing packages that build on
/// Kurvst. The dictionary shape intentionally mirrors the public function
/// arguments so callers can merge user styles into it and unpack the result.
/// -> dictionary
#let layer-defaults = (
  offset: 0,
  length: none,
  ratio: none,
  resolve-length: "min",
  start-outset: 0,
  end-outset: 0,
  side-point: none,
  accuracy: 0.001,
  optimize: true,
)

#let _sampled-pattern(samples-per-period, offset) = {
  let samples = calc.max(1, samples-per-period)
  range(0, samples + 1).map(index => {
    let at = index / samples
    let theta = 2 * calc.pi * at
    let point = offset(theta)
    (at: at, x: point.x, y: point.y)
  })
}

/// A smooth sinusoidal path pattern.
///
/// ```example
/// #let pattern = kurvst.wave(samples-per-period: 24)
/// #let points = pattern.points.map(point => (
///   x: point.at * 3.1 + point.x * 0.18,
///   y: 1 + point.y * 0.32,
/// ))
/// #native-scene({
///   native-polyline(((x: 0.0, y: 1.0), (x: 3.15, y: 1.0)), stroke: rgb("#bbbbbb") + 0.35pt)
///   native-polyline(points, stroke: rgb("#1b7f4c") + 0.75pt)
/// })
/// ```
///
/// -> dictionary
#let wave(samples-per-period: 16) = (
  kind: "points",
  name: "wave",
  interpolation: "smooth",
  points: _sampled-pattern(samples-per-period, theta => (x: 0, y: calc.sin(theta))),
)

/// A straight-segment triangular path pattern.
///
/// ```example
/// #let pattern = kurvst.zigzag()
/// #let points = pattern.points.map(point => (
///   x: point.at * 3.1 + point.x * 0.18,
///   y: 1 + point.y * 0.32,
/// ))
/// #native-scene({
///   native-polyline(((x: 0.0, y: 1.0), (x: 3.15, y: 1.0)), stroke: rgb("#bbbbbb") + 0.35pt)
///   native-polyline(points, stroke: rgb("#1b7f4c") + 0.75pt)
/// })
/// ```
///
/// -> dictionary
#let zigzag() = (
  kind: "points",
  name: "zigzag",
  interpolation: "linear",
  points: (
    (at: 0, x: 0, y: 0),
    (at: 0.25, x: 0, y: 1),
    (at: 0.75, x: 0, y: -1),
    (at: 1, x: 0, y: 0),
  ),
)

/// A smooth coil path pattern.
///
/// ```example
/// #let pattern = kurvst.coil(samples-per-period: 24, longitudinal-scale: 1.6)
/// #let points = pattern.points.map(point => (
///   x: point.at * 3.1 + point.x * 0.18,
///   y: 1 + point.y * 0.32,
/// ))
/// #native-scene({
///   native-polyline(((x: 0.0, y: 1.0), (x: 3.15, y: 1.0)), stroke: rgb("#bbbbbb") + 0.35pt)
///   native-polyline(points, stroke: rgb("#1b7f4c") + 0.75pt)
/// })
/// ```
///
/// -> dictionary
#let coil(samples-per-period: 16, longitudinal-scale: 1.25) = (
  kind: "points",
  name: "coil",
  interpolation: "smooth",
  endpoint-ramp: true,
  points: _sampled-pattern(
    samples-per-period,
    theta => (x: longitudinal-scale * calc.cos(theta), y: calc.sin(theta)),
  ),
)

#let _resolve-pattern(pattern, samples-per-period: 16, coil-longitudinal-scale: 1.25) = {
  if type(pattern) != str {
    pattern
  } else {
    let name = pattern.trim()
    if name == "wave" or name == "sine" or name == "sin" {
      wave(samples-per-period: samples-per-period)
    } else if name == "zigzag" or name == "zig-zag" or name == "triangle" {
      zigzag()
    } else if name == "coil" or name == "helix" or name == "spring" {
      coil(samples-per-period: samples-per-period, longitudinal-scale: coil-longitudinal-scale)
    } else {
      panic("Unsupported path pattern: " + pattern)
    }
  }
}

#let outset-point(from, toward, distance: 0) = {
  let dx = toward.x - from.x
  let dy = toward.y - from.y
  let length = calc.sqrt(dx * dx + dy * dy)
  if distance == 0 or length == 0 {
    from
  } else {
    let applied = calc.min(distance, length * 0.45)
    (
      x: from.x + dx / length * applied,
      y: from.y + dy / length * applied,
    )
  }
}

#let _path-value(path) = {
  if type(path) == dictionary and path.keys().contains("path") {
    path.path
  } else {
    path
  }
}

/// Build a cubic segment dictionary for a straight line.
///
/// -> dictionary
#let line-segment(start, end) = {
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

/// Evaluate a cubic segment at parameter `t`.
///
/// -> dictionary
#let cubic-point(segment, t) = {
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

/// Evaluate the tangent of a cubic segment at parameter `t`.
///
/// -> dictionary
#let cubic-tangent(segment, t) = {
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

/// Return drawable cubic segments for any Kurvst path dictionary.
///
/// Line segments and sampled point paths are converted to equivalent cubic
/// line segments so downstream renderers can use one code path.
/// -> array
#let path-segments(path) = {
  let segments = ()
  if path.keys().contains("curves") and path.curves.len() > 0 {
    for segment in path.curves {
      segments.push(segment)
    }
  } else if path.keys().contains("segments") and path.segments.len() > 0 {
    for segment in path.segments {
      segments.push(line-segment(segment.start, segment.end))
    }
  } else if path.keys().contains("points") and path.points.len() > 1 {
    for i in range(0, path.points.len() - 1) {
      segments.push(line-segment(path.points.at(i), path.points.at(i + 1)))
    }
  }
  segments
}

/// Compute the arc length of a path dictionary.
///
/// If the path already carries a `length` field, that value is returned.
/// Otherwise the length is computed from the drawable segments.
/// -> int | float
#let path-length(path, accuracy: 0.001) = {
  if path.keys().contains("length") {
    path.length
  } else {
    path-segments(path).map(segment => cubic-path(..segment, accuracy: accuracy).length).sum()
  }
}

/// Resolve a fixed and relative visible path length.
///
/// `length` is interpreted as a fixed arc length. `ratio` is interpreted as a
/// fraction of `base-length`. `method` may be `"min"`/`"shorter"`,
/// `"max"`/`"longer"`, `"length"`/`"fixed"`, `"ratio"`/`"relative"`,
/// `"none"`/`"full"`, or a function receiving
/// `(base-length, length, ratio)`.
/// -> none | int | float
#let _resolve-length-target(base-length, length: none, ratio: none, method: "min") = {
  let fixed = if length != none and length > 0 { length } else { none }
  let relative = if ratio != none and ratio > 0 { base-length * ratio } else { none }
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

#let resolve-length(base-length, length: none, ratio: none, method: "min") = {
  _resolve-length-target(base-length, length: length, ratio: ratio, method: method)
}

/// Compute the symmetric trim needed to center a shorter path layer.
///
/// -> int | float
#let center-outset(
  base-length,
  length: none,
  ratio: none,
  resolve-length: "min",
  start-outset: 0,
  end-outset: 0,
) = {
  let target = _resolve-length-target(base-length, length: length, ratio: ratio, method: resolve-length)
  if target == none {
    0
  } else {
    let visible-length = calc.max(0, base-length - start-outset - end-outset)
    calc.max(0, (visible-length - target) / 2)
  }
}

#let _side-segment(path) = {
  let segments = path-segments(path)
  if segments.len() == 0 { none } else { segments.at(calc.quo(segments.len(), 2)) }
}

#let _offset-toward-side-point(path, offset, side-point) = {
  if side-point == none or offset == none or offset == 0 {
    offset
  } else {
    let segment = _side-segment(path)
    if segment == none {
      offset
    } else {
      let mid = cubic-point(segment, 0.5)
      let tangent = cubic-tangent(segment, 0.5)
      let cross = tangent.x * (side-point.y - mid.y) - tangent.y * (side-point.x - mid.x)
      let sign = if cross < 0 { -1 } else { 1 }
      sign * calc.abs(offset)
    }
  }
}

/// Build a path dictionary from one cubic Bezier.
///
/// ```example
/// #let path = kurvst.cubic-path(..demo-segment)
/// #native-scene({
///   native-cubics(path.curves, stroke: rgb("#1b7f4c") + 0.8pt)
/// })
/// ```
///
/// -> dictionary
#let cubic-path(start: none, end: none, ctrl-a: none, ctrl-b: none, accuracy: 0.001) = {
  cbor(_plugin.curve_cubic_path(cbor.encode((
    start: start,
    end: end,
    ctrl-a: ctrl-a,
    ctrl-b: ctrl-b,
    accuracy: accuracy,
  ))))
}

/// Trim a path by curve distance from either end.
///
/// The returned dictionary is another path dictionary with `path`, `points`,
/// `curves`, `segments`, and `length`.
///
/// ```example
/// #let base = kurvst.cubic-path(..demo-segment)
/// #let trimmed = kurvst.trim-path(base, start-outset: 0.35, end-outset: 0.3)
/// #native-scene({
///   native-cubics(base.curves, stroke: rgb("#c8c8c8") + 0.45pt)
///   native-cubics(trimmed.curves, stroke: rgb("#d72638") + 1pt)
/// })
/// ```
///
/// -> dictionary
#let trim-path(path, start-outset: 0, end-outset: 0, accuracy: 0.001) = {
  cbor(_plugin.curve_trim_path(cbor.encode((
    path: _path-value(path),
    start-outset: start-outset,
    end-outset: end-outset,
    accuracy: accuracy,
  ))))
}

/// Build an open Hobby curve through `start`, `through`, and `end`.
///
/// The returned dictionary is a path dictionary. For this three-point variant,
/// `curves` contains the two cubic halves meeting at `through`.
///
/// ```example
/// #let path = kurvst.hobby-through(demo-start, demo-through, demo-end, omega: 1.2)
/// #native-scene({
///   native-cubics(path.curves, stroke: rgb("#1b7f4c") + 0.8pt)
///   native-dot(demo-through, fill: black)
/// })
/// ```
///
/// -> dictionary
#let hobby-through(start, through, end, omega: 1.0, accuracy: 0.001) = {
  cbor(_plugin.curve_hobby_through(cbor.encode((
    start: start,
    through: through,
    end: end,
    omega: omega,
    accuracy: accuracy,
  ))))
}

/// Build an open Hobby spline through two or more points.
///
/// The returned dictionary is a path dictionary backed by Kurvst's explicit
/// wire path. Rust converts that wire path to Kurbo's `BezPath` internally, so
/// it can be passed directly to @parallel-path, @pattern-path, @trim-path, or
/// @cetz-path.
///
/// ```example
/// #let spline = kurvst.hobby-spline((
///   (x: 0.0, y: 0.0),
///   (x: 0.9, y: 0.8),
///   (x: 1.8, y: -0.3),
///   (x: 2.8, y: 0.4),
/// ))
/// #let shifted = kurvst.parallel-path(spline, distance: 0.14)
/// #native-scene({
///   native-cubics(spline.curves, stroke: rgb("#c8c8c8") + 0.45pt)
///   native-cubics(shifted.curves, stroke: rgb("#1b7f4c") + 0.8pt)
/// })
/// ```
///
/// -> dictionary
#let hobby-spline(points, omega: 1.0, accuracy: 0.001) = {
  cbor(_plugin.curve_hobby_spline(cbor.encode((
    points: points,
    omega: omega,
    accuracy: accuracy,
  ))))
}

/// Generate a sampled 1D pattern along a path.
///
/// `pattern` may be `wave()`, `zigzag()`, `coil()`, a compatible string name,
/// or a point pattern:
/// `(kind: "points", interpolation: "linear" or "smooth", points: ((at: 0, x: 0, y: 0), ...))`.
/// The whole input path is sampled continuously, so pattern phase does not
/// restart at cubic segment boundaries.
///
/// ```example
/// #let spline = kurvst.hobby-spline((
///   (x: 0.0, y: 0.0),
///   (x: 0.9, y: 0.8),
///   (x: 1.8, y: -0.3),
///   (x: 2.8, y: 0.4),
/// ))
/// #let path = kurvst.pattern-path(
///   spline,
///   pattern: kurvst.wave(samples-per-period: 24),
///   amplitude: 0.12,
///   wavelength: 0.55,
/// )
/// #native-scene({
///   native-cubics(spline.curves, stroke: rgb("#c8c8c8") + 0.45pt)
///   native-cubics(path.curves, stroke: rgb("#1b7f4c") + 0.8pt)
/// })
/// ```
///
/// -> dictionary
#let pattern-path(
  path,
  pattern: "wave",
  amplitude: 0.1,
  wavelength: 1.0,
  phase: 0,
  samples-per-period: 16,
  coil-longitudinal-scale: 1.25,
  anchor-start: true,
  anchor-end: true,
  accuracy: 0.001,
) = {
  let pattern = _resolve-pattern(
    pattern,
    samples-per-period: samples-per-period,
    coil-longitudinal-scale: coil-longitudinal-scale,
  )
  cbor(_plugin.curve_pattern_path(cbor.encode((
    path: _path-value(path),
    pattern: pattern,
    amplitude: amplitude,
    wavelength: wavelength,
    phase: phase,
    samples-per-period: samples-per-period,
    coil-longitudinal-scale: coil-longitudinal-scale,
    anchor-start: anchor-start,
    anchor-end: anchor-end,
    accuracy: accuracy,
  ))))
}

/// Generate a parallel path for a path.
///
/// Pass any path returned by this module, such as @hobby-spline. The result is a
/// drawable path dictionary with the fitted parallel path and its explicit wire
/// path.
///
/// ```example
/// #let spline = kurvst.hobby-spline((
///   (x: 0.0, y: 0.0),
///   (x: 0.9, y: 0.8),
///   (x: 1.8, y: -0.3),
///   (x: 2.8, y: 0.4),
/// ))
/// #let parallel = kurvst.parallel-path(spline, distance: 0.14)
/// #native-scene({
///   native-cubics(spline.curves, stroke: rgb("#c8c8c8") + 0.45pt)
///   native-cubics(parallel.curves, stroke: rgb("#1b7f4c") + 0.8pt)
/// })
/// ```
///
/// -> dictionary
#let parallel-path(
  path,
  distance: 0,
  start-outset: 0,
  end-outset: 0,
  accuracy: 0.001,
  optimize: true,
) = {
  cbor(_plugin.curve_parallel_path(cbor.encode((
    path: _path-value(path),
    distance: distance,
    start-outset: start-outset,
    end-outset: end-outset,
    accuracy: accuracy,
    optimize: optimize,
  ))))
}

/// Build a derived visible path layer.
///
/// `layer-path` combines the common operations needed by drawing packages:
/// optional side-aware offsetting, endpoint trimming, and centered shortening by
/// a fixed `length`, a relative `ratio`, or both. The return value is a normal
/// Kurvst path dictionary that can be passed to @pattern-path, @parallel-path,
/// @trim-path, or @cetz-path.
///
/// ```example
/// #let base = kurvst.hobby-spline((
///   (x: 0.0, y: 0.0),
///   (x: 0.9, y: 0.8),
///   (x: 1.8, y: -0.3),
///   (x: 2.8, y: 0.4),
/// ))
/// #let layer = kurvst.layer-path(
///   base,
///   offset: 0.16,
///   length: 1.6,
///   ratio: 0.5,
///   resolve-length: "min",
/// )
/// #native-scene({
///   native-cubics(base.curves, stroke: rgb("#c8c8c8") + 0.45pt)
///   native-cubics(layer.curves, stroke: rgb("#1b7f4c") + 0.8pt)
/// })
/// ```
///
/// -> dictionary
#let layer-path(
  path,
  offset: 0,
  length: none,
  ratio: none,
  resolve-length: "min",
  start-outset: 0,
  end-outset: 0,
  side-point: none,
  accuracy: 0.001,
  optimize: true,
) = {
  let distance = _offset-toward-side-point(path, offset, side-point)
  let trim = center-outset(
    path-length(path, accuracy: accuracy),
    length: length,
    ratio: ratio,
    resolve-length: resolve-length,
    start-outset: start-outset,
    end-outset: end-outset,
  )
  let start-outset = start-outset + trim
  let end-outset = end-outset + trim
  if distance == none or distance == 0 {
    trim-path(path, start-outset: start-outset, end-outset: end-outset, accuracy: accuracy)
  } else {
    parallel-path(
      path,
      distance: distance,
      start-outset: start-outset,
      end-outset: end-outset,
      accuracy: accuracy,
      optimize: optimize,
    )
  }
}

#let _cetz-args(segment, unit: 1) = (
  _point(segment.start, unit: unit),
  _point(segment.end, unit: unit),
  _point(segment.ctrl-a, unit: unit),
  _point(segment.ctrl-b, unit: unit),
)

#let _cetz-bezier(segment, unit: 1, ..style) = {
  cetz.draw.bezier(.._cetz-args(segment, unit: unit), ..style)
}

#let _draw-merged-curves(segments, unit: 1, ..style) = {
  if segments.len() > 0 {
    cetz.draw.merge-path({
      for segment in segments {
        _cetz-bezier(segment, unit: unit)
      }
    }, ..style)
  } else {
    ()
  }
}

/// Draw a path dictionary through CeTZ.
///
/// ```example
/// #let path = kurvst.parallel-path(kurvst.cubic-path(..demo-segment), distance: 0.18)
/// #cetz.canvas({
///   kurvst.cetz-path(path, stroke: rgb("#355c9a") + 0.55pt)
/// })
/// ```
///
/// -> content
#let cetz-path(path, unit: 1, ..style) = {
  let style = style.named()
  if path.keys().contains("curves") and path.curves.len() > 0 {
    _draw-merged-curves(path.curves, unit: unit, ..style)
  } else if path.keys().contains("points") and path.points.len() > 1 {
    cetz.draw.line(..path.points.map(point => _point(point, unit: unit)), ..style)
  } else {
    ()
  }
}

/// Draw a sampled path pattern.
///
/// ```example
/// #let path = kurvst.pattern-path(
///   kurvst.cubic-path(..demo-segment),
///   pattern: kurvst.coil(longitudinal-scale: 1.6),
///   amplitude: 0.12,
///   wavelength: 0.55,
/// )
/// #cetz.canvas({
///   kurvst.cetz-pattern(path, stroke: rgb("#355c9a") + 0.55pt)
/// })
/// ```
///
/// -> content
#let cetz-pattern(path, unit: 1, ..style) = cetz-path(path, unit: unit, ..style)

/// Split a laid-out graph edge into source and sink half-edge paths.
///
/// The returned dictionary has `source`, `sink`, and `curve`. The split point is
/// the edge layout point, so the two half-edges join smoothly there.
///
/// ```example
/// #let halves = kurvst.edge-halves(
///   demo-edge,
///   demo-nodes,
///   source-outset: 0.2,
///   sink-outset: 0.2,
/// )
/// #native-scene({
///   native-cubics(halves.source.curves, stroke: rgb("#d72638") + 0.9pt)
///   native-cubics(halves.sink.curves, stroke: rgb("#355c9a") + 0.9pt)
///   native-dot(demo-through, fill: black)
/// })
/// ```
///
/// -> dictionary
#let edge-halves(edge, nodes, omega: 1.0, source-outset: 0, sink-outset: 0, accuracy: 0.001) = {
  if edge.source == none or edge.sink == none {
    panic("edge-halves currently requires a paired edge with source and sink")
  }

  let start = nodes.at(edge.source.node).pos
  let end = nodes.at(edge.sink.node).pos
  let split = hobby-through(start, edge.pos, end, omega: omega, accuracy: accuracy)
  let source = trim-path(cubic-path(..split.curves.at(0), accuracy: accuracy), start-outset: source-outset, accuracy: accuracy)
  let sink = trim-path(cubic-path(..split.curves.at(1), accuracy: accuracy), end-outset: sink-outset, accuracy: accuracy)

  (
    source: source,
    sink: sink,
    curve: split,
  )
}

/// Draw the two halves of a laid-out graph edge through CeTZ.
///
/// `source-style` applies from the source node to the edge layout point, and
/// `sink-style` applies from the edge layout point to the sink node.
///
/// ```example
/// #cetz.canvas({
///   kurvst.cetz-edge-halves(
///     demo-edge,
///     demo-nodes,
///     source-outset: 0.2,
///     sink-outset: 0.2,
///     source-style: (stroke: rgb("#d72638") + 0.6pt),
///     sink-style: (stroke: rgb("#355c9a") + 0.6pt),
///   )
/// })
/// ```
///
/// -> content
#let cetz-edge-halves(
  edge,
  nodes,
  unit: 1,
  omega: 1.0,
  source-outset: 0,
  sink-outset: 0,
  accuracy: 0.001,
  source-style: (:),
  sink-style: (:),
) = {
  let halves = edge-halves(
    edge,
    nodes,
    omega: omega,
    source-outset: source-outset,
    sink-outset: sink-outset,
    accuracy: accuracy,
  )
  cetz-path(halves.source, unit: unit, ..source-style)
  cetz-path(halves.sink, unit: unit, ..sink-style)
}
