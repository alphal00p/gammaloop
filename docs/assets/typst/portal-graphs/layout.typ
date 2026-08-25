#import "../../../../crates/linnest/typst/src/lib.typ": (
  draw, graph, layout as apply-layout, physics,
)
#import "@preview/cetz:0.5.1" as cetz
#import "../../../../assets/embedded/drawing/templates/layout-core.typ": (
  bind-layout,
)
#import "edge-style.typ" as edge-style
#import "../theme.typ": palette

#let portal-node-style(node) = if node.at("hidden", default: false) in (true, "true") {
  (radius: 0, fill: none, stroke: none)
} else {
  (:)
}

// About-page diagrams may tag an invisible, fixed-position node with one
// editable cubic cut contour. Linnest still owns and draws the physics graph;
// this CeTZ node renderer adds the smooth contour afterward so it sits above
// the propagators without baking generated paths into the SVG source.
#let _plain(value) = str(value).trim("\"")

#let _point(value) = {
  let coordinates = _plain(value).split(",").map(value => float(value.trim()))
  (coordinates.at(0), coordinates.at(1))
}

#let portal-cut-curve(node) = {
  let kind = _plain(node.at("cut_curve"))
  let paint = if kind == "blue" { palette.cut-blue } else { palette.cut-red }
  cetz.draw.bezier(
    _point(node.at("cut_start")),
    _point(node.at("cut_end")),
    _point(node.at("cut_control_start")),
    _point(node.at("cut_control_end")),
    stroke: (paint: paint, thickness: 3.1pt, cap: "round"),
  )
}

#let portal-draw-node(node, box) = {
  if node.at("cut_curve", default: none) != none {
    portal-cut-curve(node)
  } else if not (node.at("hidden", default: false) in (true, "true")) {
    cetz.draw.circle(box.center, name: box.name, ..box.style)
    if box.label != none {
      cetz.draw.content(box.center, box.label, padding: 0, ..box.label-style)
    }
  }
}

// Website adapter: the graph algorithm is exactly GammaLoop's save-dot core;
// only the transparent, theme-aware presentation differs.
#let layout = bind-layout(
  draw: draw,
  graph: graph,
  apply-layout: apply-layout,
  physics: physics,
  edge-style: edge-style,
  diagram-options: (
    title: none,
    node-fill: none,
    node-stroke: palette.ink + 1.45pt,
    node-style: portal-node-style,
  ),
)

#let layout-with-cut-curves = bind-layout(
  draw: draw,
  graph: graph,
  apply-layout: apply-layout,
  physics: physics,
  edge-style: edge-style,
  diagram-options: (
    title: none,
    node-fill: none,
    node-stroke: palette.ink + 1.7pt,
    node-style: portal-node-style,
    draw-node: portal-draw-node,
  ),
)
