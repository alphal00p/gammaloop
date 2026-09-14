#import "../src/lib.typ": draw, graph, layout, layouts, subgraph
#import graph: *
#import "map-style.typ" as feynman
#import feynman: momentum as mom
#import "@preview/cetz:0.5.1" as cetz

// Native font and stroke sizes are independent of the graph coordinate unit.
#let graph-style = feynman.graph-style(unit: 2.6mm, line-width: 0.5pt)
#let diagram-style = (
  font-size: 6pt,
  cut-line-width: 0.8pt,
  endpoint-box: (
    padding: (x: 0.5, y: 0.4),
    radius: 0.25,
    fills: (p1: blue.lighten(85%), p2: orange.lighten(85%)),
  ),
  padding: 0.3,
)
#let draw-initials = false

// One stroke per opening: original, p1 replaced, p2 replaced, both replaced.
#let initial-cut-styles = (
  (paint: blue.darken(45%), dash: "dotted"),
  (paint: blue.darken(15%), dash: "dashed"),
  (paint: blue.lighten(10%), dash: (4pt, 1.5pt, 1pt, 1.5pt)),
  (paint: blue.lighten(30%), dash: (4pt, 1.5pt, 1pt, 1.5pt, 1pt, 1.5pt)),
).map(style => (thickness: 1pt) + style)

#let base-layout = layouts.options(
  spring: (strength: 18, length: 0.15),
  repulsion: (
    strength: 12,
    centering: 0.005,
    edge-node: .29,
    edge-edge: .1,
    dangling: 14.75,
    dangling-centroid: 50,
  ),
  constraints: (side-strength: 10),
  labels: (steps: 100, distance: 0.9, spring: 12, repulsion: 8),
  solver: (
    steps: 50,
    epochs: 50,
    step: 0.81,
    max-movement: 0.4,
    cooling: 0.85,
    depth-scale: 1,
    flattening-end: 0.5,
  ),
)
#let routed-layout = layouts.options(
  spring: (strength: 40, length: 0.2),
  repulsion: (
    strength: 10,
    centering: 0.005,
    edge-node: .39,
    edge-edge: .75,
    dangling: 5.5,
    dangling-centroid: 30,
  ),
)
#let in-x = group("in", side: "-", start: -7)
#let out-x = group("out", side: "+", start: 7)
#let top = group("top", start: 5)
#let mid = group("mid", start: 0)
#let bot = group("bot", start: -5)
#let mid2 = group("mid2", start: 2)
#let edge-data = (particle: "d") + mom(offset: .4)
#let vertices = { for name in (<a>, <b>, <c>, <d>) { node(name) } }
#let compact = {
  node(<compact-top>, hidden: true, pos: pos(x: start(0), y: top))
  node(<compact-bottom>, hidden: true, pos: pos(x: start(0), y: bot))
  edge(
    <compact>,
    source(<compact-top>),
    sink(<compact-bottom>),
    momentum: [],
    style: none,
    pos: pos(x: pin(0)),
  )
}

// `cut-x` is in graph units; auto puts the cut just right of center.
#let diagram(
  g,
  options: base-layout,
  cut-x: auto,
  cut-y: auto,
  initial-cut: none,
  draw-after: none,
) = context {
  draw(
    layout(graph.style(g, ..graph-style), ..options),
    ..feynman.draw-style,
    padding: diagram-style.padding,
    draw-after: (g, bounds) => {
      // Group signed external momenta explicitly, including the through-gluon's
      // hidden endpoints. Shading follows the solved positions, not the layout seeds.
      let edges = graph.edges(g)
      let endpoints = graph
        .boundaries(g)
        .map(b => (
          b
            + (
              group: b.data.crossings.at(b.crossing).group,
            )
        ))
      // Directly built dangling edges can provide grouping without cut provenance.
      endpoints += edges
        .filter(e => (
          e.boundary == none
            and (e.source == none or e.sink == none)
            and type(e.data) == dictionary
            and "group" in e.data
        ))
        .map(e => (group: e.data.group, side: e.data.side, pos: e.pos))
      let box-style = diagram-style.endpoint-box
      cetz.draw.on-layer(-1, {
        for (momentum, fill) in box-style.fills {
          for side in ("left", "right") {
            let points = endpoints
              .filter(b => b.side == side and b.group == momentum)
              .map(b => b.pos)
            if points.len() > 0 {
              let xs = points.map(p => p.x)
              let ys = points.map(p => p.y)
              cetz.draw.rect(
                (
                  calc.min(..xs) - box-style.padding.x,
                  calc.min(..ys) - box-style.padding.y,
                ),
                (
                  calc.max(..xs) + box-style.padding.x,
                  calc.max(..ys) + box-style.padding.y,
                ),
                radius: box-style.radius,
                fill: fill,
                stroke: none,
              )
            }
          }
        }
      })
      // Exclude the invisible compactification spring and its anchor nodes
      // when locating the cut; the renderer supplies the visible drawing bounds.
      let nodes = graph
        .nodes(g)
        .filter(n => (
          n.boundary == none
            and (
              n.data == none or not n.data.at("hidden", default: false)
            )
        ))
      let xs = nodes.map(n => n.pos.x)
      let x = if cut-x == auto {
        calc.min(..xs) + 0.6 * (calc.max(..xs) - calc.min(..xs))
      } else { cut-x }
      let y = if cut-y == auto { 1.5 } else { cut-y }
      bounds = (
        left: bounds.left - box-style.padding.x,
        right: bounds.right + box-style.padding.x,
        top: bounds.top + y,
        bottom: bounds.bottom - y,
        width: bounds.width + 2 * box-style.padding.x,
        height: bounds.height + 2 * y,
      )
      cetz.draw.line((x, bounds.top), (x, bounds.bottom), stroke: (
        paint: red,
        thickness: diagram-style.cut-line-width,
        dash: "dashed",
      ))
      if initial-cut != none and draw-initials {
        for side in ("left", "right") {
          let side-xs = endpoints.filter(b => b.side == side).map(b => b.pos.x)
          if side-xs.len() > 0 {
            let x = side-xs.first()
            cetz.draw.line(
              (x, bounds.top),
              (x, bounds.bottom),
              stroke: initial-cut-styles.at(initial-cut),
            )
          }
        }
      }
      if type(draw-after) == function { draw-after(g, bounds) } else {
        draw-after
      }
    },
  )
}
