#import "xbox-common.typ": *

#set page(height: auto, width: auto, margin: 2mm)
#set text(size: diagram-style.font-size)
#show math.equation: set text(size: math-font-size)

#let soft-layout = layouts.options(base: routed-layout, spring: (length: .2))
#let edge-data = (particle: "d", spring-length: .3, show-momentum: false)

#let in-top = pos(x: in-x, y: top, z: pin(0))
#let in-mid = pos(x: in-x, y: mid, z: pin(0))
#let in-bot = pos(x: in-x, y: bot, z: pin(0))
#let out-top = pos(x: out-x, y: top, z: pin(0))
#let out-mid = pos(x: out-x, y: mid, z: pin(0))
#let out-bot = pos(x: out-x, y: bot, z: pin(0))

#let opened = graph.build(
  default-edge-data: edge-data + (show-momentum: true) + mom(offset: .4, label: (gap: .15)),
  vertices,
  {
    edge(<D1.1>, sink(<a>), momentum: [$p_1 - k$], pos: in-bot, group: "p1", side: "left")
    edge(<D2>, source(<a>), sink(<b>), show-momentum: false)
    edge(
      <D3.1>,
      sink(<b>),
      group: "p2",
      side: "left",
      orientation: "reversed",
      momentum: [$p_2-k$],
      pos: in-top,
      ..mom(side: "left", label: (gap: .2, shift: 1)),
    )
    edge(
      <D4>,
      source(<b>),
      sink(<c>),
      show-momentum: false,
      particle: "a",
      spring-length: .1,
      ..mom(side: "right", label: (gap: .2, shift: .5)),
    )
    edge(
      <D3.2>,
      source(<c>),
      group: "p2",
      side: "right",
      orientation: "reversed",
      momentum: [$p_1-k$],
      ..mom(side: "right", shift:1,label: (gap: .2,shift:0.3)),
      pos: out-bot,
    )
    edge(<D5>, sink(<c>), source(<d>), show-momentum: false)
    edge(<D1.2>, source(<d>), orientation: "reversed", momentum: [$p_2-k$], pos: out-top, group: "p1", side: "right",..mom(side: "left", label: (gap: .2, shift: -1.3)))
    edge(<D6.1>, sink(<a>), momentum: [$k$], particle: "g", pos: in-mid,
      group: ("p1", "p2"), side: "left",
      ..mom(side: "left", label: (gap: .2)))
    edge(<D6.2>, source(<d>), momentum: [$k$], particle: "g", pos: out-mid,
      group: ("p1", "p2"), side: "right",..mom(side: "right", label: (gap: .2, shift: 0.3)))
  },
  // ,compact
)


#let soft = graph.build(
  default-edge-data: edge-data,
  vertices,
  {
    edge(<D2>, source(<a>), sink(<b>), momentum: [$p_1+k$])
    edge(
      <D3.1>,
      sink(<b>),
      source(<d>),
      orientation: "reversed",
      momentum: [$p_2-k$],
    )
    edge(
      <D4>,
      source(<b>),
      sink(<c>),
      momentum: [$p_(12)$],
      particle: "a",
      spring-length: .1,
    )
    edge(<top>,sink(<a>), source(<c>), momentum: [$p_2-k$])
    edge(<bottom>, source(<d>), sink(<c>), momentum: [$p_1+k$], pos: out-top)
    // edge(<D1.2>, sink(<d>), momentum: [$p_1$], pos: out-top)
    edge(<D6>, sink(<a>), source(<d>), momentum: [$k$], particle: "g")
  },
  // ,compact
)

#let soft-diagram = context {
  set text(size: diagram-style.font-size)
  show math.equation: set text(size: diagram-style.font-size)
  let selected = subgraph.select(soft, edges: (<D2>, <D3.1>, <D6>, <top>, <bottom>))
  let options = layouts.options(base: soft-layout, solver: (subgraph-mode: "isolated"))
  let g = layout(graph.style(soft, ..graph-style), ..options, subgraph: selected)
  // Position the closing edges relative to the solved internal topology.
  let points = (graph.nodes(g, subgraph: selected) + graph.edges(g, subgraph: selected)).map(item => item.pos)
  let center = (
    x: (calc.min(..points.map(p => p.x)) + calc.max(..points.map(p => p.x))) / 2,
    y: (calc.min(..points.map(p => p.y)) + calc.max(..points.map(p => p.y))) / 2,
  )
  let g = graph.map(g, edge: (
    "D4": (pos: pos(x: center.x+2 , y: center.y - 2)),
  ))
  // Hidden momentum labels should not enlarge the final canvas.
  let g = graph.style(g, ..(graph-style + (edge-label: none)))
  draw(g, ..feynman.draw-style, padding: diagram-style.padding,
    draw-after: (g, bounds) => {


      let x = center.x+1.5
      let y = center.y - 1
      let radius = 0.25
      let diagonal = radius / calc.sqrt(2)
      let stroke = (
        paint: red.transparentize(30%),
        thickness: diagram-style.cut-line-width,
        cap: "round",
      )
      cetz.draw.rotate(-30deg, origin: (x, y))
      cetz.draw.line((x, y + radius), (x, y + 3), stroke: stroke)
      let stroke = stroke + (paint: red.lighten(30%))
      cetz.draw.circle((x, y), radius: radius, fill: none, stroke: stroke)
      cetz.draw.line((x - diagonal, y - diagonal), (x + diagonal, y + diagonal), stroke: stroke)
      cetz.draw.line((x - diagonal, y + diagonal), (x + diagonal, y - diagonal), stroke: stroke)
      

      
    },
  )
}

#context {
  let delimiter-size = measure(soft-diagram).height
  $
  #diagram(opened, options: layouts.options(base: soft-layout, spring: (length: .2)), cut-y: .5, cut-x: .3) = op("disc")_(p_1^2) op("disc")_(p_2^2)lr(
    (#h(-1mm)
      lr(#box(inset: 0mm, baseline: 45%, soft-diagram) |, size: #delimiter-size)_"soft"),
    size: #delimiter-size,
  ) + cal(O)(Lambda^2/q_"h"^2)
  $
}
