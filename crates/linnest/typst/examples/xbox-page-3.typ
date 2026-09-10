#import "xbox-common.typ": *

#set page(height: auto, width: auto, margin: 2mm)
#set text(size: diagram-style.font-size)

#let soft-layout = layouts.options(base: routed-layout, spring: (length: .12))

#let in-top = pos(x: in-x, y: top, z: pin(0))
#let in-mid = pos(x: in-x, y: mid, z: pin(0))
#let in-bot = pos(x: in-x, y: bot, z: pin(0))
#let out-top = pos(x: out-x, y: top, z: pin(0))
#let out-mid = pos(x: out-x, y: mid, z: pin(0))
#let out-bot = pos(x: out-x, y: bot, z: pin(0))

#let opened = graph.build(
  default-edge-data: edge-data + (spring-length: .3),
  vertices,
  {
    edge(<D1.1>, sink(<a>), momentum: [], pos: in-bot)
    edge(<D2>, source(<a>), sink(<b>), momentum: [$p_1+k$])
    edge(
      <D3.1>,
      sink(<b>),
      orientation: "reversed",
      momentum: [$p_2-k$],
      pos: in-top,
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
      pos: out-bot,
    )
    edge(<D5>, sink(<c>), source(<d>), momentum: [$p_1+k$])
    edge(<D1.2>, sink(<d>), momentum: [$p_1$], pos: out-top)
    edge(<D6.1>, sink(<a>), momentum: [$k$], particle: "g", pos: in-mid)
    edge(<D6.2>, source(<d>), momentum: [$k$], particle: "g", pos: out-mid)
  },
  // ,compact
)

#let soft = graph.build(
  default-edge-data: edge-data + (spring-length: .3),
  vertices,
  {
    edge(<D1.1>, sink(<a>), momentum: [], pos: in-bot)
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
    edge(<D3.2>, source(<c>), momentum: [$p_2-k$], pos: out-bot)
    edge(<D5>, sink(<c>), momentum: [$p_1+k$], pos: out-top)
    edge(<D5.1>, source(<d>), momentum: [$p_1+k$], pos: in-top)
    // edge(<D1.2>, sink(<d>), momentum: [$p_1$], pos: out-top)
    edge(<D6>, sink(<a>), source(<d>), momentum: [$k$], particle: "g")
  },
  // ,compact
)

$
  #diagram(opened, options: soft-layout, cut-y: .5, cut-x: .3, show-momentum: false) = op("disc")_(p_1^2) op("disc")_(p_2^2)lr(
    (#h(-3mm)
      #box(inset: -2mm, baseline: 45%, diagram(soft, options: soft-layout, cut-y: .5, cut-x: 1., show-momentum: false)))|
  )_"soft"
$
