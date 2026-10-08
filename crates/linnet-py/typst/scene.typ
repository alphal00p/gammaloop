// Domain scene primitives shared by advanced graph rendering.
#import "@preview/cetz:0.5.1" as cetz

#let hatch = tiling(size: (4pt, 4pt), line(start: (0pt, 4pt), end: (4pt, 0pt), stroke: 0.4pt + gray))

#let draw-node(node, box) = {
  let center = box.center
  if node.at("scene-rectangular", default: false) {
    cetz.draw.rect(
      (center.at(0) - box.width / 2, center.at(1) - box.height / 2),
      (center.at(0) + box.width / 2, center.at(1) + box.height / 2),
      name: box.name, fill: box.style.fill, stroke: box.style.stroke,
    )
  } else {
    cetz.draw.circle(center, name: box.name, ..box.style)
  }
  if box.label != none {
    cetz.draw.content(center, box.label, padding: 0, ..box.label-style)
  }
}

#let annotation-layer(edge, half) = {
  if not edge.at("scene-momentum", default: false) { return (:) }
  let marked-half = if edge.at("sink-half-edge", default: none) == none { "source" } else { "sink" }
  let arrow = (
    offset: edge.scene-momentum-offset, length: 1.4, ratio: 0.5, resolve-length: "min",
    stroke: (paint: rgb("#333333"), thickness: 1pt, cap: "round"),
    pattern: none,
    mark: if half == marked-half { (end: (symbol: "straight", scale: 0.8)) } else { none },
    mark-position: 1, mark-orientation: "path",
  )
  (
    (:), arrow,
    (label-only: true, label: auto, label-path: arrow, label-gap: 0.1,
     label-slide: true, label-side: auto),
  )
}
#let source-style(edge) = annotation-layer(edge, "source")
#let sink-style(edge) = annotation-layer(edge, "sink")

// Scene styles are defaults. Explicit graph options override them, while
// element selectors remain the final layer in the common renderer.
#let scene-style(kind, defaults, style, draw, element) = {
  let style = if style == none { (:) } else { style }
  let draw = if draw == none { (:) } else { draw }
  let result = defaults
  for field in if kind == "node" { ("radius", "fill", "stroke") } else { ("stroke",) } {
    let key = kind + "-" + field
    if key in draw { result.insert(field, draw.at(key)) }
  }
  if kind == "node" {
    let configured = draw.at("node-style", default: style.at("node-style", default: (:)))
    let configured = if type(configured) == function { configured(element) } else { configured }
    if configured != none { result += configured }
  }
  result
}
