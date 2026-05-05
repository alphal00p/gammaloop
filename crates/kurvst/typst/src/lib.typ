#import "@preview/cetz:0.3.4" as cetz

#let _plugin = plugin("../kurvst.wasm")

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
/// #native-pattern-object(kurvst.wave(samples-per-period: 24))
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
/// #native-pattern-object(kurvst.zigzag())
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
/// #native-pattern-object(kurvst.coil(samples-per-period: 24, longitudinal-scale: 1.6))
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

/// Move `from` toward `toward` by `distance`, capped before `toward`.
///
/// ```example
/// #let moved = kurvst.outset-point(demo-start, demo-end, distance: 0.65)
/// #native-outset-demo(demo-start, demo-end, moved)
/// ```
///
/// -> dictionary
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

/// Split a cubic Bezier curve into two cubic Bezier segments.
///
/// Returns a dictionary with `curve` and `segments`. Each segment has `start`,
/// `ctrl-a`, `ctrl-b`, and `end` points.
///
/// ```example
/// #let split = kurvst.split-cubic(
///   demo-segment.start,
///   demo-segment.end,
///   demo-segment.ctrl-a,
///   demo-segment.ctrl-b,
///   t: 0.45,
/// )
/// #native-split-demo(split)
/// ```
///
/// -> dictionary
#let split-cubic(start, end, ctrl-a, ctrl-b, t: 0.5) = {
  cbor(_plugin.curve_split_cubic(cbor.encode((
    start: start,
    end: end,
    ctrl-a: ctrl-a,
    ctrl-b: ctrl-b,
    t: t,
  ))))
}

/// Trim a cubic Bezier by curve distance from either end.
///
/// Returns the remaining cubic segment after splitting by arc length.
///
/// ```example
/// #let trimmed = kurvst.trim-cubic(
///   demo-segment.start,
///   demo-segment.end,
///   demo-segment.ctrl-a,
///   demo-segment.ctrl-b,
///   start-outset: 0.35,
///   end-outset: 0.3,
/// )
/// #native-trim-demo(demo-segment, trimmed)
/// ```
///
/// -> dictionary
#let trim-cubic(start, end, ctrl-a, ctrl-b, start-outset: 0, end-outset: 0, accuracy: 0.001) = {
  cbor(_plugin.curve_trim_cubic(cbor.encode((
    start: start,
    end: end,
    ctrl-a: ctrl-a,
    ctrl-b: ctrl-b,
    start-outset: start-outset,
    end-outset: end-outset,
    accuracy: accuracy,
  ))))
}

/// Trim a cubic segment returned by @split-cubic, @split-through, or
/// @hobby-through by curve distance from either end.
///
/// ```example
/// #let half = kurvst.split-cubic(
///   demo-segment.start,
///   demo-segment.end,
///   demo-segment.ctrl-a,
///   demo-segment.ctrl-b,
/// ).segments.at(0)
/// #let trimmed = kurvst.trim-segment(half, start-outset: 0.25, end-outset: 0.15)
/// #native-trim-demo(half, trimmed)
/// ```
///
/// -> dictionary
#let trim-segment(segment, start-outset: 0, end-outset: 0, accuracy: 0.001) = {
  trim-cubic(
    segment.start,
    segment.end,
    segment.ctrl-a,
    segment.ctrl-b,
    start-outset: start-outset,
    end-outset: end-outset,
    accuracy: accuracy,
  )
}

/// Build a quadratic Bezier through `start`, `through`, and `end`, promote it
/// to a cubic, then split it at `t`.
///
/// This is useful for graph edges whose layout point should lie on the smooth
/// curve between source and sink.
///
/// ```example
/// #let split = kurvst.split-through(demo-start, demo-through, demo-end)
/// #native-split-demo(split)
/// ```
///
/// -> dictionary
#let split-through(start, through, end, t: 0.5) = {
  cbor(_plugin.curve_split_quadratic_through(cbor.encode((
    start: start,
    through: through,
    end: end,
    t: t,
  ))))
}

/// Build an open Hobby curve through `start`, `through`, and `end`.
///
/// Returns a dictionary with two cubic Bezier `segments`, split at `through`.
///
/// ```example
/// #let split = kurvst.hobby-through(demo-start, demo-through, demo-end, omega: 1.2)
/// #native-split-demo(split)
/// ```
///
/// -> dictionary
#let hobby-through(start, through, end, omega: 1.0) = {
  cbor(_plugin.curve_hobby_through(cbor.encode((
    start: start,
    through: through,
    end: end,
    omega: omega,
  ))))
}

/// Generate a sampled 1D pattern along a cubic Bezier path.
///
/// `pattern` may be `wave()`, `zigzag()`, `coil()`, a compatible string name,
/// or a point pattern:
/// `(kind: "points", interpolation: "linear" or "smooth", points: ((at: 0, x: 0, y: 0), ...))`.
/// The returned dictionary contains `points`, smooth `curves` for smooth
/// patterns, and straight `segments` in graph coordinates.
///
/// ```example
/// #let path = kurvst.pattern-cubic(
///   demo-segment.start,
///   demo-segment.end,
///   demo-segment.ctrl-a,
///   demo-segment.ctrl-b,
///   pattern: kurvst.wave(samples-per-period: 24),
///   amplitude: 0.14,
///   wavelength: 0.55,
/// )
/// #native-pattern-path(path)
/// ```
///
/// -> dictionary
#let pattern-cubic(
  start,
  end,
  ctrl-a,
  ctrl-b,
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
  cbor(_plugin.curve_pattern_cubic(cbor.encode((
    start: start,
    end: end,
    ctrl-a: ctrl-a,
    ctrl-b: ctrl-b,
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

/// Generate a sampled 1D pattern along a cubic segment returned by this module.
///
/// ```example
/// #let path = kurvst.pattern-segment(
///   demo-segment,
///   pattern: kurvst.zigzag(),
///   amplitude: 0.13,
///   wavelength: 0.55,
/// )
/// #native-pattern-path(path)
/// ```
///
/// -> dictionary
#let pattern-segment(
  segment,
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
  pattern-cubic(
    segment.start,
    segment.end,
    segment.ctrl-a,
    segment.ctrl-b,
    pattern: pattern,
    amplitude: amplitude,
    wavelength: wavelength,
    phase: phase,
    samples-per-period: samples-per-period,
    coil-longitudinal-scale: coil-longitudinal-scale,
    anchor-start: anchor-start,
    anchor-end: anchor-end,
    accuracy: accuracy,
  )
}

/// Generate a parallel path for a cubic Bezier.
///
/// Positive `distance` offsets to the left normal of the cubic direction;
/// negative values offset to the right. The result contains `points`, fitted
/// cubic `curves`, straight `segments`, and `length`.
///
/// ```example
/// #let parallel = kurvst.parallel-cubic(
///   demo-segment.start,
///   demo-segment.end,
///   demo-segment.ctrl-a,
///   demo-segment.ctrl-b,
///   distance: 0.18,
///   start-outset: 0.35,
///   end-outset: 0.35,
/// )
/// #native-pattern-path(parallel, base: demo-segment)
/// ```
///
/// -> dictionary
#let parallel-cubic(
  start,
  end,
  ctrl-a,
  ctrl-b,
  distance: 0,
  start-outset: 0,
  end-outset: 0,
  accuracy: 0.001,
  optimize: true,
) = {
  cbor(_plugin.curve_parallel_cubic(cbor.encode((
    start: start,
    end: end,
    ctrl-a: ctrl-a,
    ctrl-b: ctrl-b,
    distance: distance,
    start-outset: start-outset,
    end-outset: end-outset,
    accuracy: accuracy,
    optimize: optimize,
  ))))
}

/// Generate a parallel path for a cubic segment returned by this module.
///
/// ```example
/// #let parallel = kurvst.parallel-segment(demo-segment, distance: -0.18)
/// #native-pattern-path(parallel, base: demo-segment)
/// ```
///
/// -> dictionary
#let parallel-segment(segment, distance: 0, start-outset: 0, end-outset: 0, accuracy: 0.001, optimize: true) = {
  parallel-cubic(
    segment.start,
    segment.end,
    segment.ctrl-a,
    segment.ctrl-b,
    distance: distance,
    start-outset: start-outset,
    end-outset: end-outset,
    accuracy: accuracy,
    optimize: optimize,
  )
}

/// Convert a cubic segment returned by @split-cubic or @split-through into
/// positional arguments accepted by `cetz.draw.bezier`.
///
/// ```example
/// #cetz.canvas({
///   cetz.draw.bezier(..kurvst.cetz-args(demo-segment), stroke: black + 0.6pt)
/// })
/// ```
///
/// -> array
#let cetz-args(segment, unit: 1) = (
  _point(segment.start, unit: unit),
  _point(segment.end, unit: unit),
  _point(segment.ctrl-a, unit: unit),
  _point(segment.ctrl-b, unit: unit),
)

/// Draw one cubic segment through CeTZ.
///
/// ```example
/// #cetz.canvas({
///   kurvst.cetz-bezier(demo-segment, stroke: rgb("#d72638") + 0.6pt)
/// })
/// ```
///
/// -> content
#let cetz-bezier(segment, unit: 1, ..style) = {
  cetz.draw.bezier(..cetz-args(segment, unit: unit), ..style)
}

#let _draw-merged-curves(segments, unit: 1, ..style) = {
  if segments.len() > 0 {
    cetz.draw.merge-path({
      for segment in segments {
        cetz-bezier(segment, unit: unit)
      }
    }, ..style)
  } else {
    ()
  }
}

/// Draw a path returned by @pattern-cubic, @pattern-segment, @parallel-cubic,
/// or @parallel-segment.
///
/// ```example
/// #let path = kurvst.parallel-segment(demo-segment, distance: 0.18)
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
  } else if path.points.len() > 1 {
    cetz.draw.line(..path.points.map(point => _point(point, unit: unit)), ..style)
  } else {
    ()
  }
}

/// Draw a sampled pattern returned by @pattern-cubic or @pattern-segment.
///
/// ```example
/// #let path = kurvst.pattern-segment(
///   demo-segment,
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

/// Split a laid-out graph edge into source and sink half-edge cubic segments.
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
/// #native-edge-halves-demo(halves)
/// ```
///
/// -> dictionary
#let edge-halves(edge, nodes, omega: 1.0, source-outset: 0, sink-outset: 0, accuracy: 0.001) = {
  if edge.source == none or edge.sink == none {
    panic("edge-halves currently requires a paired edge with source and sink")
  }

  let start = nodes.at(edge.source.node).pos
  let end = nodes.at(edge.sink.node).pos
  let split = hobby-through(start, edge.pos, end, omega: omega)
  let source = trim-segment(split.segments.at(0), start-outset: source-outset, accuracy: accuracy)
  let sink = trim-segment(split.segments.at(1), end-outset: sink-outset, accuracy: accuracy)

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
  cetz-bezier(halves.source, unit: unit, ..source-style)
  cetz-bezier(halves.sink, unit: unit, ..sink-style)
}
