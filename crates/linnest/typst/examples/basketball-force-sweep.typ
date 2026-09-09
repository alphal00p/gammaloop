#import "../src/lib.typ": draw, graph, layouts


#set page(width: auto, height: auto, margin: 10pt)

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

#let shared-layout = (
  layout-algo: "force",
  seed: 7,
  epochs: 52,
  steps: 60,
  viewport-w: 4.0,
  viewport-h: 2.8,
  label-steps: 0,
  g-center: 0.0,
  z-spring: 2.0,
  z-spring-growth: 1.0,
)

#let distance(a, b) = calc.sqrt(calc.pow(a.x - b.x, 2) + calc.pow(a.y - b.y, 2))

#let summarize(settings) = {
  let g = layouts.layout(base, ..shared-layout, ..settings)
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
      t: (ap.x * ab.x + ap.y * ab.y) / calc.pow(ab-len, 2),
      normal: (-ap.x * ab.y + ap.y * ab.x) / ab-len,
    )
  }
  (
    settings: settings,
    node-distance: ab-len,
    edges: edges.map(project),
    drawn: draw(g),
  )
}

#let combos = ()
#for beta in (60.0, 80.0, 100.0, 120.0, 150.0, 180.0, 220.0){
  for k-spring in (12.0, 24.0, 36.0, 48.0, 72.0, 96.0){
    for gamma-ee in (0.003,0.03,0.3){// 0.005, 0.0075, 0.01, 0.015, 0.02)
      for length-scale in (0.68,){// 0.72, 0.78, 0.95)
        for gamma-ev in (0.1,){// 0.1)
          for gamma-dangling in (0.1,){// 1.0)
            combos.push((
              beta: beta,
              k-spring: k-spring,
              gamma-ee: gamma-ee,
              gamma-ev: gamma-ev,
              gamma-dangling: gamma-dangling,
              length-scale: length-scale,
            ))
          }
        }
      }
    }
  }
}
#let summ = combos.map(summarize);
#metadata(summ.map(x => (settings: x.settings, node-distance: x.node-distance, edges: x.edges))) <linnest-basketball-force-sweep>

#summ.map(x=>[
#x.settings

  #x.drawn

]).join()
= Basketball Force Sweep

This document exists for `typst query` experiments. It emits force-layout
summaries for a grid over `beta`, `k-spring`, `gamma-ee`, `gamma-ev`,
`gamma-dangling`, and `length-scale`.
