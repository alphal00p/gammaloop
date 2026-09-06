#set page(width: auto, height: auto, margin: 0pt)

#import "crates/linnest/typst/src/lib.typ": draw, graph, layout, layouts
#import graph: *

#let close(a, b, epsilon: 1e-6) = calc.abs(a - b) < epsilon
#let same-pos(a, b) = close(a.x, b.x) and close(a.y, b.y)

// Semantic layout groups should be equivalent to their flat spelling.
#let semantic-base = layouts.options(spring: (strength: 8, length: 0.5))
#let semantic-patch = layouts.options(base: semantic-base, spring: (
  length: 0.7,
))
#assert(semantic-base.spring.length == 0.5)
#assert(semantic-patch.spring == (strength: 8, length: 0.7))

#let layout-base = graph.build({
  node(<l0>)
  node(<l1>)
  node(<l2>)
  edge(source(<l0>), sink(<l1>))
  edge(source(<l1>), sink(<l2>))
})
#let flat-layout = layout(
  layout-base,
  steps: 2,
  seed: 17,
  k-spring: 8,
  length-scale: 0.5,
  beta: 4,
  label-steps: 0,
)
#let semantic-layout = layout(
  layout-base,
  spring: (strength: 8, length: 0.5),
  repulsion: (strength: 4),
  labels: (steps: 0),
  solver: (steps: 2, seed: 17),
)
#let flat-pos = graph.nodes(flat-layout).map(node => node.pos)
#let semantic-pos = graph.nodes(semantic-layout).map(node => node.pos)
#assert(
  range(flat-pos.len()).all(index => same-pos(
    flat-pos.at(index),
    semantic-pos.at(index),
  )),
)

// Semantic rank constraints override the corresponding flat option and are
// visible in the public node positions.
#let ranked = layout(
  layout-base,
  rank-same: ((0, 2),),
  constraints: (same-rank: ((0, 1),)),
  labels: (steps: 0),
  solver: (algorithm: "stable-layered"),
)
#let ranked-pos = graph.nodes(ranked).map(node => node.pos)
#assert(close(ranked-pos.at(0).y, ranked-pos.at(1).y))
#assert(not close(ranked-pos.at(0).y, ranked-pos.at(2).y))

// Grouped coordinates remain shared after force layout, including a group
// shared between nodes and an edge position.
#let grouped-base = graph.build({
  node(
    <g0>,
    pos: pos(
      x: group("left", side: "-", start: -2),
      y: group("lane", start: -1),
    ),
  )
  node(
    <g1>,
    pos: pos(
      x: group("left", side: "-", start: -2),
      y: start(2),
    ),
  )
  node(
    <g2>,
    pos: pos(
      x: group("right", side: "+", start: 2),
      y: group("lane", start: 1),
    ),
  )
  edge(
    source(<g0>),
    sink(<g2>),
    pos: pos(x: start(0), y: group("lane")),
  )
  edge(source(<g1>), sink(<g2>))
})
#let grouped = layout(
  grouped-base,
  labels: (steps: 0),
  solver: (algorithm: "force", steps: 12, seed: 17),
)
#let grouped-pos = graph.nodes(grouped).map(node => node.pos)
#let grouped-edge-pos = graph.edges(grouped).first().pos
#assert(close(grouped-pos.at(0).x, grouped-pos.at(1).x))
#assert(
  close(grouped-pos.at(0).y, grouped-pos.at(2).y)
    and close(grouped-pos.at(0).y, grouped-edge-pos.y),
)
#assert(grouped-pos.at(0).x <= 1e-6 and grouped-pos.at(2).x >= -1e-6)

#let invisible-node = (radius: 0, fill: none, stroke: none)
#let arrow-color = rgb("#d119e6")
#let shaft-color = rgb("#16a34a")
#let shaft-reference-color = rgb("#0ea5e9")
#let crossing-color = rgb("#dc2626")
#let crossing-mark-color = rgb("#7e22ce")
#let crossing-reference-color = rgb("#94a3b8")
#let crossing-reference-mark-color = rgb("#9333ea")
#let target-color = rgb("#2563eb")
#let label-left-color = rgb("#ea580c")
#let label-right-color = rgb("#0891b2")
#let label-reference-color = rgb("#475569")
#let label-left-stroke-color = rgb("#c2410c")
#let label-right-stroke-color = rgb("#0e7490")
#let default-color = rgb("#65a30d")
#let callback-color = rgb("#f97316")
#let callback-fallback-color = rgb("#84cc16")
#let local-color = rgb("#ca8a04")
#let endpoint-color = rgb("#0f766e")
#let endpoint-fallback-color = rgb("#4d7c0f")
#let hidden-color = rgb("#be123c")

// A finite shifted layer keeps an end mark on its painted endpoint.
#let arrow-graph = graph.build({
  node(<arrow-a>, pos: pos(x: 0, y: 0, mode: "pin"))
  node(<arrow-b>, pos: pos(x: 6, y: 0, mode: "pin"))
  edge(
    source(<arrow-a>),
    sink(<arrow-b>),
    pos: pos(x: 3, y: 0, mode: "pin"),
  )
})
#let arrow-style = (
  (
    stroke: shaft-reference-color + 0.4pt,
    length: 2,
    resolve-length: "length",
  ),
  (
    stroke: shaft-color + 0.8pt,
    length: 2,
    resolve-length: "length",
    shift: 1,
    mark: (
      end: (
        symbol: ">",
        fill: arrow-color,
        stroke: arrow-color,
        anchor: "center",
        shorten-to: auto,
      ),
    ),
    mark-position: "end",
  ),
)

// A cut edge has a measurable gap but retains exactly one end mark. The
// reference edge selects the same target without intersecting it.
#let crossing-graph = graph.build({
  node(<crossing-left>, pos: pos(x: 0, y: 0, mode: "pin"))
  node(<crossing-right>, pos: pos(x: 6, y: 0, mode: "pin"))
  node(<crossing-down>, pos: pos(x: 2, y: -2, mode: "pin"))
  node(<crossing-up>, pos: pos(x: 2, y: 2, mode: "pin"))
  node(<reference-left>, pos: pos(x: 0, y: 3, mode: "pin"))
  node(<reference-right>, pos: pos(x: 6, y: 3, mode: "pin"))
  edge(
    source(<crossing-left>),
    sink(<crossing-right>),
    style: (
      stroke: crossing-color + 0.8pt,
      crossing-under: <crossing-target>,
      crossing-gap: 1,
      mark: (
        end: (
          symbol: ">",
          fill: crossing-mark-color,
          stroke: crossing-mark-color,
          anchor: "center",
          shorten-to: auto,
        ),
      ),
      mark-position: "end",
    ),
  )
  edge(
    source(<crossing-down>),
    <crossing-target>,
    sink(<crossing-up>),
    style: (stroke: target-color + 0.8pt),
  )
  edge(
    source(<reference-left>),
    sink(<reference-right>),
    style: (
      stroke: crossing-reference-color + 0.8pt,
      crossing-under: <crossing-target>,
      crossing-gap: 1,
      mark: (
        end: (
          symbol: ">",
          fill: crossing-reference-mark-color,
          stroke: crossing-reference-mark-color,
          anchor: "center",
          shorten-to: auto,
        ),
      ),
      mark-position: "end",
    ),
  )
})

// Attached labels are emitted once per common source/sink layer and follow
// that layer's arc-length shift.
#let labels-graph = graph.build({
  node(<labels-left>, pos: pos(x: 0, y: 0, mode: "pin"))
  node(<labels-right>, pos: pos(x: 6, y: 0, mode: "pin"))
  edge(
    source(<labels-left>),
    sink(<labels-right>),
    pos: pos(x: 3, y: 0, mode: "pin"),
  )
})
#let label-square(color) = box(width: 3pt, height: 3pt, fill: color)
#let label-style = (
  (stroke: label-reference-color + 0.3pt),
  (
    stroke: label-left-stroke-color + 0.3pt,
    length: 1,
    resolve-length: "length",
    shift: -1,
    label: label-square(label-left-color),
    label-gap: 0.1,
    label-side: "left",
  ),
  (
    stroke: label-right-stroke-color + 0.3pt,
    length: 1,
    resolve-length: "length",
    shift: 1,
    label: label-square(label-right-color),
    label-gap: 0.1,
    label-side: "left",
  ),
)

// Element callbacks, `auto`, `none`, and endpoint overrides are tested by the
// colors that reach the rendered output. The callback assertions also exercise
// its documented graph record.
#let precedence-graph = graph.build({
  node(<callback-a>, pos: pos(x: 0, y: 0, mode: "pin"))
  node(<callback-b>, pos: pos(x: 4, y: 0, mode: "pin"))
  node(<auto-a>, pos: pos(x: 0, y: 2, mode: "pin"))
  node(<auto-b>, pos: pos(x: 4, y: 2, mode: "pin"))
  node(<hidden-a>, pos: pos(x: 0, y: 4, mode: "pin"))
  node(<hidden-b>, pos: pos(x: 4, y: 4, mode: "pin"))
  node(<endpoint-a>, pos: pos(x: 0, y: 6, mode: "pin"))
  node(<endpoint-b>, pos: pos(x: 4, y: 6, mode: "pin"))
  edge(
    source(<callback-a>, role: "producer"),
    sink(<callback-b>, role: "consumer"),
    kind: "callback",
    style: edge => {
      assert(edge.eid == 0 and edge.fields.kind == "callback")
      assert(edge.source-half-edge.data.role == "producer")
      assert(edge.sink-half-edge.data.role == "consumer")
      (stroke: callback-color + 0.8pt)
    },
  )
  edge(source(<auto-a>), sink(<auto-b>), style: auto)
  edge(source(<hidden-a>), sink(<hidden-b>), style: none)
  edge(
    source(<endpoint-a>),
    sink(<endpoint-b>),
    style: (stroke: local-color + 0.8pt),
  )
})
#let precedence-base(edge) = (
  stroke: (
    callback-fallback-color,
    default-color,
    hidden-color,
    endpoint-fallback-color,
  ).at(edge.eid)
    + 0.8pt,
)
#let precedence-endpoint(edge) = if edge.eid == 3 {
  (stroke: endpoint-color + 0.8pt)
} else {
  (:)
}

#let bare-draw = (
  unit: 10pt,
  padding: 1,
  node-label: none,
  node-style: invisible-node,
  node-outset: 0,
)
#stack(
  spacing: 12pt,
  draw(arrow-graph, edge-style: arrow-style, ..bare-draw),
  draw(crossing-graph, ..bare-draw),
  draw(labels-graph, edge-style: label-style, ..bare-draw),
  draw(
    precedence-graph,
    edge-style: precedence-base,
    source-style: precedence-endpoint,
    sink-style: precedence-endpoint,
    ..bare-draw,
  ),
)
