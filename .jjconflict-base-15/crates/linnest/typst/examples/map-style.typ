#let massless = 0.5pt
#let massive = 1pt
#let edge-stroke = (paint: black, thickness: massless, cap: "round")

#let fermion-mark = (
  end: (
    symbol: ">",
    fill: black,
    stroke: black + 0.2pt,
    anchor: "center",
    shorten-to: auto,
  ),
  scale: .5,
)
#let momentum-mark = (
  end: (
    symbol: "straight",
    fill: black,
    stroke: black + 0.4pt,
    anchor: "center",
    shorten-to: auto,
  ),
  scale: 0.50,
)

#let fermion = (
  stroke: edge-stroke,
  mark: fermion-mark,
  mark-position: "center-if-dangling",
  mark-orientation: "edge",
)
#let photon = (
  stroke: edge-stroke,
  pattern: "wave",
  pattern-amplitude: 0.10,
  pattern-wavelength: 0.50,
)
#let gluon = (
  stroke: edge-stroke,
  pattern: "coil",
  pattern-amplitude: 0.15,
  pattern-wavelength: 0.60,
  pattern-coil-longitudinal-scale: 1.60,
)
#let scalar = (
  stroke: edge-stroke + (thickness: massive, dash: (0.1em, 0.45em)),
)
#let particles = (
  "a": photon,
  "photon": photon,
  "g": gluon,
  "gluon": gluon,
  "scalar": scalar,
  "ghG": scalar,
)

// Sparse edge options: omitted settings preserve inherited defaults. Arrow
// shift and label.shift remain independent; only their field names are expanded.
#let momentum(label: (:), ..arrow) = {
  assert(arrow.pos().len() == 0, message: "momentum: expected named options")
  assert(type(label) == dictionary, message: "momentum: label must be a dictionary")
  let fields = (:)
  for (kind, options, keys) in (
    ("arrow", arrow.named(), ("side", "offset", "length", "shift")),
    ("label", label, ("gap", "shift", "anchor")),
  ) {
    for (key, value) in options {
      assert(key in keys, message: "momentum: unknown " + kind + " option " + key)
      fields.insert("momentum-" + kind + "-" + key, value)
    }
  }
  fields
}

#let _value(element, key, default) = element.fields.at(key, default: default)
#let _text(element, key, default) = str(_value(element, key, default)).trim(
  "\"",
)
#let _number(element, key, default) = {
  let value = _value(element, key, default)
  if type(value) in (int, float) { value } else { float(str(value).trim("\"")) }
}
#let _enabled(element, key) = (
  _value(element, key, false) in (true, "true", "\"true\"")
)
#let _route(edge) = {
  let route = _value(edge, "route", none)
  if route == none { (:) } else { (route: str(route).trim("\"")) }
}

#let _particle-layer(edge) = {
  let particle = _text(edge, "particle", "d")
  let style = particles.at(particle, default: fermion) + _route(edge)
  if _enabled(edge, "cut") {
    style.split-gap = _number(edge, "cut-gap", 0.55)
  }
  let under = _value(edge, "crossing-under", none)
  if under != none {
    style.crossing-under = under
    style.crossing-gap = _number(edge, "crossing-gap", 0.55)
  }
  // Crossing cuts split only the paint, so one carrier mark can move freely.
  style.mark-shift = _number(edge, "fermion-arrow-shift", 0)
  style
}

#let _momentum-layers(edge) = {
  let shift = _number(edge, "momentum-arrow-shift", 0)
  let label-shift = _number(edge, "momentum-label-shift", shift)
  let anchor = _value(edge, "momentum-label-anchor", auto)
  let offset = _number(edge, "momentum-arrow-offset", 0.62)
  let side = _value(edge, "momentum-arrow-side", auto)
  let side = if side == auto { "auto" } else { str(side).trim("\"") }
  assert(
    side in ("auto", "left", "right"),
    message: "momentum-arrow-side must be auto, left, or right",
  )
  // Explicit sides are relative to source -> sink and shared with the label carrier.
  // The gap is in graph units: to the box for auto, or to an explicit anchor.
  let geometry = (
    _route(edge)
      + (
        offset: if side == "auto" { offset } else {
          calc.abs(offset) * if side == "left" { 1 } else { -1 }
        },
        offset-side: if side == "auto" { "label" } else { none },
        label-side: if side == "auto" { auto } else { side },
        label-gap: _number(edge, "momentum-label-gap", 0.45),
        label-style: (anchor: if type(anchor) == str { anchor.trim("\"") } else { anchor }),
      )
  )
  let arrow = (
    geometry
      + (
        length: _number(edge, "momentum-arrow-length", 1.0),
        shift: shift,
        ratio: none,
        resolve-length: "length",
        stroke: (paint: black, thickness: 0.4pt, cap: "round"),
        mark: momentum-mark,
        mark-position: "end",
        mark-orientation: "path",
      )
  )
  // Auto clears the complete label box; explicit anchors use only the gap.
  // The label defaults to the requested arrow shift, but follows a point on the
  // full invisible path, so its endpoint clamps never depend on arrow length.
  (
    arrow,
    geometry + (
      length: none,
      ratio: none,
      resolve-length: "none",
      shift: 0,
      stroke: none,
      mark: none,
      label: edge.momentum,
      label-shift: label-shift,
    ),
  )
}

#let edge-style(edge) = (
  _particle-layer(edge),
  .._momentum-layers(edge),
)

#let node-style(node) = if _enabled(node, "hidden") {
  (radius: 0, fill: none, stroke: none)
} else {
  (radius: 0.18, fill: white, stroke: edge-stroke)
}

// The hidden ordinary label participates in layout and selects the automatic side;
// an explicit momentum-arrow-side moves the visible label and momentum shaft together.
#let graph-style = (
  unit: 1.35,
  node-label: none,
  node-style: node-style,
  edge-label: edge => hide([$p_(#edge.eid)$]),
  edge-label-style: (anchor: "center", padding: 0.05),
)
#let draw-style = (
  edge-style: edge-style,
  padding: 1.5,
  edge-dangling-tangent: "horizontal",
)
