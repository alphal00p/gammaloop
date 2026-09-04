#import "../src/lib.typ": graph

#let massless = 0.5mm
#let massive = 1pt
#let dashed = (0.1em, 0.45em)
#let edge-stroke = (paint: black, thickness: massless, cap: "round")
#let straight = (stroke: edge-stroke)
#let wave = (
  stroke: edge-stroke,
  pattern: "wave",
  pattern-amplitude: 0.20,
  pattern-wavelength: 0.50,
)
#let coil = (
  stroke: edge-stroke,
  pattern: "coil",
  pattern-amplitude: 0.25,
  pattern-wavelength: 0.60,
  pattern-coil-longitudinal-scale: 1.60,
)
#let scalar = (
  source: (stroke: edge-stroke + (thickness: massive, dash: dashed)),
  sink: (stroke: edge-stroke + (thickness: massive, dash: dashed)),
)

#let fermion-mark = (
  end: (
    symbol: ">",
    fill: black,
    stroke: black + 0.2pt,
    anchor: "center",
    shorten-to: auto,
  ),
  scale: 1.20,
)

#let fermion = (source: straight, sink: straight, fermion: true)
#let photon = (source: wave, sink: wave)
#let gluon = (source: coil, sink: coil)

#let particle-map = (
  "a": photon,
  "mu+": fermion + (label: [$mu^+$]),
  "mu-": fermion + (label: [$mu^-$]),
  "e+": fermion + (label: [$e^+$]),
  "e-": fermion + (label: [$e^-$]),
  "q": fermion + (label: [$q$]),
  "d": fermion,
  "t": fermion + (label: [$t$]),
  "photon": photon,
  "g": gluon,
  "gluon": gluon,
  "fermion": fermion + (label: [$f$]),
  "scalar": scalar + (label: [$s$]),
  "ghG": scalar + (label: [$s$]),
)

#let _field(record, key, default) = {
  let data = record.at("data", default: none)
  let statements = record.at("statements", default: (:))
  if type(data) == dictionary and data.keys().contains(key) {
    data.at(key)
  } else if type(statements) == dictionary and statements.keys().contains(key) {
    statements.at(key)
  } else {
    record.at(key, default: default)
  }
}

#let _number(record, key, default) = {
  let value = _field(record, key, default)
  if type(value) in (int, float) { value } else { float(str(value).trim("\"")) }
}

#let _particle(edge) = str(_field(edge, "particle", "")).trim("\"")
#let _entry(edge) = particle-map.at(
  _particle(edge),
  default: (source: straight, sink: straight),
)
#let _has-half(edge, half) = {
  edge.at(half + "-half-edge", default: edge.at(half, default: none)) != none
}
#let _momentum-mark-half(edge) = if _has-half(edge, "sink") { "sink" } else {
  "source"
}

#let _particle-style(edge, half) = {
  let entry = _entry(edge)
  let style = entry.at(half, default: entry.source)
  let route = _field(edge, "route", none)
  if route != none { style += (route: route) }
  if entry.at("fermion", default: false) {
    (
      style
        + (
          mark: fermion-mark,
          mark-position: "center-if-dangling",
          mark-orientation: "edge",
        )
    )
  } else {
    style
  }
}

#let momentum-mark = (
  end: (
    symbol: ")>",
    fill: black,
    stroke: black + 0.2pt,
    anchor: "center",
    shorten-to: auto,
  ),
  scale: 0.90,
)

#let _momentum-style(edge, half) = {
  let length = _number(edge, "momentum-arrow-length", 1.70)
  let shift = _number(edge, "momentum-arrow-shift", 0)
  let paired = _has-half(edge, "source") and _has-half(edge, "sink")
  // Opposite per-half lengths translate the visible interval along a paired
  // edge: positive shifts toward the sink and negative shifts toward source.
  let shifted-length = if paired {
    length + (if half == "source" { -2 * shift } else { 2 * shift })
  } else {
    length
  }
  let style = (
    offset: _number(edge, "momentum-arrow-offset", 0.62),
    offset-side: "label",
    length: calc.max(0.02, shifted-length),
    ratio: none,
    resolve-length: "length",
    stroke: (paint: black, thickness: 1.2pt, cap: "round"),
  )
  let route = _field(edge, "route", none)
  if route != none { style += (route: route) }
  if half == _momentum-mark-half(edge) {
    style + (mark: momentum-mark)
  } else {
    style
  }
}

#let _half-style(edge, half) = (
  _particle-style(edge, half),
  _momentum-style(edge, half),
)

#let source-style(edge) = _half-style(edge, "source")
#let sink-style(edge) = _half-style(edge, "sink")
#let edge-label(edge) = [$p_#(edge.eid)$]
#let edge-label-style(edge) = (anchor: "center", padding: 0.05)

#let _hidden(node) = _field(node, "hidden", false) == true
#let node-label(node) = if _hidden(node) { none } else { [$n_#(node.vid)$] }
#let node-style(node) = if _hidden(node) {
  (radius: 0, fill: none, stroke: none)
} else {
  (fill: white, stroke: black + 0.6pt)
}

// Measure labels and nodes before layout. Draw reuses the stored callbacks.
#let style(graph_, unit: 1.35) = graph.style(
  graph_,
  unit: unit,
  node-label: node-label,
  node-style: node-style,
  edge-label: edge-label,
  edge-label-style: edge-label-style,
)

#let _x(point) = if type(point) == array { point.at(0) } else { point.x }
#let _y(point) = if type(point) == array { point.at(1) } else { point.y }
#let _add(left, right) = (_x(left) + _x(right), _y(left) + _y(right))
#let _sub(left, right) = (_x(left) - _x(right), _y(left) - _y(right))
#let _scale(point, value) = (_x(point) * value, _y(point) * value)
#let _dot(left, right) = _x(left) * _x(right) + _y(left) * _y(right)
#let _unit(point) = {
  let length = calc.sqrt(_dot(point, point))
  if length <= 1e-9 { (1, 0) } else { _scale(point, 1 / length) }
}

#let _edge-tangent(edge, nodes) = {
  let source = edge.at("source", default: none)
  let sink = edge.at("sink", default: none)
  if source != none and sink != none {
    _sub(nodes.at(sink.node).pos, nodes.at(source.node).pos)
  } else if source != none {
    _sub(edge.pos, nodes.at(source.node).pos)
  } else {
    _sub(nodes.at(sink.node).pos, edge.pos)
  }
}

#let _edge-center(edge, nodes) = {
  let source = edge.at("source", default: none)
  let sink = edge.at("sink", default: none)
  if source != none and sink != none {
    edge.pos
  } else {
    let node = if source != none { source.node } else { sink.node }
    _scale(_add(nodes.at(node).pos, edge.pos), 0.5)
  }
}

// Put the complete label box beyond the offset momentum shaft, with a fixed
// gap. Its center follows any along-edge momentum-arrow-shift as well. Keeping
// the collision-selected side avoids moving several nearby labels together.
#let position-labels(graph_, gap: 0.32) = {
  let nodes = graph.nodes(graph_)
  graph.map(graph_, edge: edge => {
    let center = _edge-center(edge, nodes)
    let label = edge.at("label-pos", default: none)
    if label == none { label = center }
    let tangent = _unit(_edge-tangent(edge, nodes))
    let normal = (-_y(tangent), _x(tangent))
    let displacement = _sub(label, center)
    let radial = _dot(displacement, normal)
    let side = if radial < 0 { -1 } else { 1 }
    let paired = _has-half(edge, "source") and _has-half(edge, "sink")
    let shift = if paired { _number(edge, "momentum-arrow-shift", 0) } else {
      0
    }
    let label-extent = (
      (
        calc.abs(_x(normal)) * _number(edge, "label-width", 0)
          + calc.abs(_y(normal)) * _number(edge, "label-height", 0)
      )
        / 2
    )
    let clearance = (
      calc.abs(_number(edge, "momentum-arrow-offset", 0.62))
        + label-extent
        + gap
    )
    let shifted = _add(
      center,
      _add(
        _scale(tangent, shift),
        _scale(normal, side * clearance),
      ),
    )
    (label-pos: shifted)
  })
}
