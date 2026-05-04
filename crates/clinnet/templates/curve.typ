#import "@preview/cetz:0.3.4" as cetz

#let _plugin = plugin("./linnest-curve.wasm")

#let _point(p, unit: 1) = (p.x * unit, p.y * unit)

/// Move `from` toward `toward` by `distance`, capped before `toward`.
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
/// `ctrl_a`, `ctrl_b`, and `end` points.
/// -> dictionary
#let split-cubic(start, end, ctrl-a, ctrl-b, t: 0.5) = {
  cbor(_plugin.curve_split_cubic(cbor.encode((
    start: start,
    end: end,
    ctrl_a: ctrl-a,
    ctrl_b: ctrl-b,
    t: t,
  ))))
}

/// Trim a cubic Bezier by curve distance from either end.
///
/// Returns the remaining cubic segment after splitting by arc length.
/// -> dictionary
#let trim-cubic(start, end, ctrl-a, ctrl-b, start-outset: 0, end-outset: 0, accuracy: 0.001) = {
  cbor(_plugin.curve_trim_cubic(cbor.encode((
    start: start,
    end: end,
    ctrl_a: ctrl-a,
    ctrl_b: ctrl-b,
    start_outset: start-outset,
    end_outset: end-outset,
    accuracy: accuracy,
  ))))
}

/// Trim a cubic segment returned by @split-cubic, @split-through, or
/// @hobby-through by curve distance from either end.
/// -> dictionary
#let trim-segment(segment, start-outset: 0, end-outset: 0, accuracy: 0.001) = {
  trim-cubic(
    segment.start,
    segment.end,
    segment.ctrl_a,
    segment.ctrl_b,
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
/// -> dictionary
#let hobby-through(start, through, end, omega: 1.0) = {
  cbor(_plugin.curve_hobby_through(cbor.encode((
    start: start,
    through: through,
    end: end,
    omega: omega,
  ))))
}

/// Convert a cubic segment returned by @split-cubic or @split-through into
/// positional arguments accepted by `cetz.draw.bezier`.
/// -> array
#let cetz-args(segment, unit: 1) = (
  _point(segment.start, unit: unit),
  _point(segment.end, unit: unit),
  _point(segment.ctrl_a, unit: unit),
  _point(segment.ctrl_b, unit: unit),
)

/// Draw one cubic segment through CeTZ.
/// -> content
#let cetz-bezier(segment, unit: 1, ..style) = {
  cetz.draw.bezier(..cetz-args(segment, unit: unit), ..style)
}

/// Split a laid-out graph edge into source and sink half-edge cubic segments.
///
/// The returned dictionary has `source`, `sink`, and `curve`. The split point is
/// the edge layout point, so the two half-edges join smoothly there.
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
