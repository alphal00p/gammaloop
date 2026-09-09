#import "../src/lib.typ": draw, graph, layouts, subgraph

#set page(width: auto, height: auto, margin: 10pt)
#set text(size: 8pt)

#let basketball-src = ```dot
digraph basketball {
  node [num = "1"]
  edge [particle=scalar_1]
  e [style=invis]

  e -> A:1 [id=5]
  B:0 -> e [id=4]
  A -> B [id=0 lmb_id=0]
  A -> B [id=1 lmb_id=1]
  A -> B [id=2 lmb_id=2]
  A -> B [id=3]
}
```

#let base = graph.parse(basketball-src.text).first()

#let node-label(node) = node.at("name", default: none)

#let node-style(node) = (
  radius: 0.18,
  fill: white,
  stroke: rgb("#555b66") + 0.45pt,
)

#let curved-edge-style(edge) = (
  stroke: rgb("#8090a8") + 0.5pt,
  route: "edge-pos",
)

#let straight-edge-style(edge) = (
  stroke: rgb("#8090a8") + 0.5pt,
  route: "straight-through",
)

#let shared-layout = (
  layout-algo: "force",
  seed: 7,
  epochs: 52,
  steps: 60,
  viewport-w: 4.0,
  viewport-h: 2.8,
  label-steps: 0,
  g-center: 0.0,
)

#let variants = (
  (
    name: "default",
    note: "prompt z collapse",
    settings: (:),
  ),
  (
    name: "old-z-schedule",
    note: "delayed z collapse",
    settings: (
      z-spring: 0.05,
      z-spring-growth: 1.3,
    ),
  ),
  (
    name: "selected",
    note: "corrected zoom tuning",
    settings: (
      beta: 150.0,
      k-spring: 36.0,
      gamma-ee: 0.003,
      gamma-ev: 0.1,
      gamma-dangling: 1.0,
      length-scale: 0.78,
    ),
  ),
  (
    name: "balanced-3.0",
    note: "sweep target: wider nodes",
    settings: (
      beta: 150.0,
      k-spring: 48.0,
      gamma-ee: 0.01,
      gamma-ev: 0.0,
      gamma-dangling: 1.0,
      length-scale: 0.68,
    ),
  ),
  (
    name: "balanced-3.1",
    note: "slightly longer springs",
    settings: (
      beta: 150.0,
      k-spring: 72.0,
      gamma-ee: 0.01,
      gamma-ev: 0.0,
      gamma-dangling: 1.0,
      length-scale: 0.72,
    ),
  ),
  (
    name: "wide-3.9",
    note: "longer rest length",
    settings: (
      beta: 220.0,
      k-spring: 24.0,
      gamma-ee: 0.005,
      gamma-ev: 0.0,
      gamma-dangling: 1.0,
      length-scale: 0.78,
    ),
  ),
  (
    name: "wide-4.1",
    note: "very long rest length",
    settings: (
      beta: 60.0,
      k-spring: 72.0,
      gamma-ee: 0.005,
      gamma-ev: 0.0,
      gamma-dangling: 1.0,
      length-scale: 0.95,
    ),
  ),
)

#let layout-variant(variant) = layouts.layout(base, ..shared-layout, ..variant.settings)

#let distance(a, b) = calc.sqrt(calc.pow(a.x - b.x, 2) + calc.pow(a.y - b.y, 2))

#let graph-summary(name, g) = {
  let nodes = graph.nodes(g)
  let edges = graph.edges(g)
  let a = nodes.find(node => node.name == <A>)
  let b = nodes.find(node => node.name == <B>)
  let ab = (x: b.pos.x - a.pos.x, y: b.pos.y - a.pos.y)
  let ab-len = distance(a.pos, b.pos)
  let project(edge) = {
    let ap = (x: edge.pos.x - a.pos.x, y: edge.pos.y - a.pos.y)
    (
      id: edge.edge,
      pos: edge.pos,
      t: (ap.x * ab.x + ap.y * ab.y) / calc.pow(ab-len, 2),
      normal: (-ap.x * ab.y + ap.y * ab.x) / ab-len,
    )
  }
  (
    name: name,
    node-distance: ab-len,
    nodes: nodes.map(node => (name: str(node.name), pos: node.pos)),
    edges: edges.map(project),
  )
}

#let render-variant(variant, straight: false) = {
  let g = layout-variant(variant)
  let selected = subgraph.label(g, "GS")
  let nodes = graph.nodes(g)
  let a = nodes.find(node => node.name == <A>)
  let b = nodes.find(node => node.name == <B>)
  let node-distance = distance(a.pos, b.pos)
  let route-name = if straight { "straight springs" } else { "smooth edge-pos" }
  let edge-style = if straight { straight-edge-style } else { curved-edge-style }
  [
    #block(width: 9.2cm)[
      #text(weight: "bold")[#variant.name]
      #h(0.8em)
      #text(fill: gray)[#variant.note; #route-name]
      #linebreak()
      #text(size: 0.85em, fill: gray)[A-B: #(calc.round(node-distance * 100) / 100)]
      #draw(
        g,
        unit: 1.0,
        subgraph: selected,
        node-label: node-label,
        node-style: node-style,
        source-style: edge-style,
        sink-style: edge-style,
        edge-label: none,
        node-min-radius: 0.16,
        node-label-padding: 0.05,
        padding: 0.32,
        debug: 2,
        debug-edge-radius: 0.055,
        subgraph-edge-style: (stroke: rgb("#e4504f") + 1.7pt),
      )
    ]
  ]
}

#let layouts = variants.map(variant => {
  let g = layout-variant(variant)
  graph-summary(variant.name, g)
})

#metadata(layouts) <linnest-basketball-force>

= Basketball Force Layout

The force model does produce midpoint edge controls for this graph once the
temporary z offsets are collapsed promptly. The old z schedule kept those hidden
3D offsets active for too long: repulsion used 3D distances, while the visible
layout only drew x/y coordinates. `beta` does repel the two nodes, but it also
scales edge-node, edge-edge, dangling, and center charges, so the useful tuning
ratios are `k-spring` versus `beta * gamma-*`. With the dimensional scaling
fixed, `length-scale` behaves as a zoom-like parameter: it changes the natural
size while preserving the useful force ratios. The selected sweep point keeps
the internal edge controls near their incidence midpoint and spreads them across
the middle of the A-B segment.

= Smoothed Curves

#grid(
  columns: 2,
  gutter: 12pt,
  row-gutter: 10pt,
  ..variants.map(render-variant),
)

= Straight Force Springs

This uses the same layout positions, but draws each paired edge as the two
straight springs that the force system actually optimizes: source to edge point,
then edge point to sink. That makes the old z schedule's behavior more direct:
parallel edge control points repel into a crossed diamond when the hidden z
coordinate is still carrying much of the separation.

#grid(
  columns: 2,
  gutter: 12pt,
  row-gutter: 10pt,
  ..variants.map(variant => render-variant(variant, straight: true)),
)
