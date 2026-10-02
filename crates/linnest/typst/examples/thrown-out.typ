#import "xbox-common.typ": *

#set page(height: auto, width: auto, margin: 2mm)
#set text(size: diagram-style.font-size)
#show math.equation: set text(size: math-font-size)

#let p1-left = (pos: pos(x: in-x, y: top), group: "p1", side: "left")
#let p1-right = (pos: pos(x: out-x, y: top), group: "p1", side: "right")
#let p2-left = (pos: pos(x: in-x, y: bot), group: "p2", side: "left")
#let p2-right = (pos: pos(x: out-x, y: bot), group: "p2", side: "right")
#let k-left = (pos: pos(x: in-x, y: mid), group: "p2", side: "left")
#let k-right = (pos: pos(x: out-x, y: mid), group: "p2", side: "right")

#let direct = {
  let g = graph.build(
    default-edge-data: edge-data + (spring-length: .1),
    vertices,
    {
      edge(<D1.1>, sink(<a>), momentum: [$p_1$], ..p1-left)
      edge(<D2>, source(<a>), sink(<b>), momentum: [$p_1+k$])
      edge(
        <D3.1>,
        sink(<b>),
        orientation: "reversed",
        momentum: [$p_2-k$],
        ..p2-left,
      )
      edge(
        <D4>,
        source(<b>),
        sink(<c>),
        momentum: [$p_(12)$],
        particle: "a",
        spring-length: .1,
      )
      edge(
        <D3.2>,
        source(<c>),
        orientation: "reversed",
        momentum: [$p_2-k$],
        ..p2-right,
      )
      edge(<D5>, source(<c>), sink(<d>), momentum: [$p_1+k$])
      edge(<D1.2>, source(<d>), momentum: [$p_1$], ..p1-right)
      edge(<D6.1>, sink(<a>), momentum: [$k$], particle: "g", ..k-left)
      edge(<D6.2>, source(<d>), momentum: [$k$], particle: "g", ..k-right)
    },
    // ,compact
  )
  graph.map(g, edge: (
    "D1.1": mom(side: "left", label: (gap: .1)),
    D2: (show-momentum: false),
    "D3.1": mom(side: "left", label: (shift: -1.7, gap: .01, anchor: "east")),
    D4: (show-momentum: false),
    "D3.2": mom(side: "left", label: (gap: .2, shift: .6, anchor: "west")),
    D5: (show-momentum: false),
    "D1.2": mom(side: "left", label: (gap: .1)),
    "D6.1": mom(side: "left", label: (gap: .1)),
    "D6.2": mom(side: "left", label: (gap: .2)),
  ))
}

#let crossed = {
  let g = graph.build(default-edge-data: edge-data, vertices, {
    edge(<D1.1>, sink(<a>), momentum: [$p_1$], ..p1-left)
    edge(<D2>, source(<a>), sink(<b>), momentum: [$p_1-k$])
    edge(
      <D3.1>,
      sink(<b>),
      orientation: "reversed",
      momentum: [$p_2-k$],
      ..p2-left,
    )
    edge(<D4>, source(<b>), sink(<c>), momentum: [$p_(12) - 2k$], particle: "a")
    edge(<D3.2>, sink(<c>), momentum: [$p_1$], ..p1-right)
    edge(<D5>, source(<c>), sink(<d>), momentum: [$p_2-2k$])
    edge(<D1.2>, source(<d>), momentum: [$p_2-k$], ..p2-right)
    edge(<D6.1>, source(<a>), momentum: [$k$], particle: "g", ..k-right)
    edge(<D6.2>, sink(<d>), momentum: [$k$], particle: "g", ..k-left)
  })
  graph.map(
    g,
    node: (b: (pos: pos(y: mid)), c: (pos: pos(y: mid))),
    edge: (
      "D1.1": mom(side: "left", length: .8, label: (
        shift: 1,
        gap: .1,
       
      )),
      D2: (show-momentum: false),
      "D3.1": (spring-length: 1.5, crossing-under: <D6.2>, crossing-gap: .7)
        + mom(side: "right", length: .8, label: (gap: .1)),
      D4: (spring-length: 1.5, show-momentum: false),
      "D3.2": (spring-length: 1.5, crossing-under: <D6.1>, crossing-gap: .7) + mom(side: "right", length: .8, label: (gap: .1)),
      D5: (show-momentum: false),
      "D1.2": mom(side: "left", shift: 1, length: .8, label: (
        shift: 2,
        gap: .2,
        anchor: "south-west",
      )),
      "D6.1": (bend: .44, spring-length: 2.5)
        + mom(side: "left", shift: 2, length: .8, label: (gap: .2)),
      "D6.2": (bend: -.54, spring-length: 2.5)
        + mom(side: "left", shift: -2, length: .8, label: (gap: .2)),
    ),
  )
}

#align(center, grid(
  columns: 3,
  align: center + horizon,
  diagram(direct, options: routed-layout, cut-y: .5, cut-x: .3),
  {h(3mm); text("or",size:10pt); h(3mm)},
  diagram(crossed, options: routed-layout, cut-y: .5, cut-x: .3),
))
