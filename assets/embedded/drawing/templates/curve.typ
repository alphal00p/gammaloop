#import "@preview/cetz:0.3.4" as cetz

#let _plugin = plugin("kurvst.wasm")

#let _point(p, unit: 1) = (p.x * unit, p.y * unit)

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
/// The returned dictionary is a path dictionary backed by a serialized Kurbo
/// `BezPath`, so it can be passed directly to @parallel-path, @pattern-path,
/// @trim-path, or @cetz-path.
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

/// Generate a sampled 1D pattern along a serialized Kurbo `BezPath`.
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

/// Generate a parallel path for a serialized Kurbo `BezPath`.
///
/// Pass any path returned by this module, such as @hobby-spline. The result is a
/// drawable path dictionary with the fitted parallel path and its serialized
/// `BezPath`.
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
