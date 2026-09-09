#set page(width: auto, height: auto, margin: 0pt)

#import "crates/linnest/typst/src/lib.typ": draw, graph, layout, layouts
#import graph: *
#import "map-style.typ" as feynman
#import "@preview/cetz:0.5.1" as cetz
#import "curved-arrow-behavior.typ": curved-arrow-behavior
#curved-arrow-behavior
#import "weighted-cut-behavior.typ": weighted-cut-behavior
#weighted-cut-behavior
#import "named-map-behavior.typ"

#let close(a, b, epsilon: 1e-6) = calc.abs(a - b) < epsilon
#let same-pos(a, b) = close(a.x, b.x) and close(a.y, b.y)

// Auxiliary depth stays separate from XY, including through graph.map patches.
#let depth-base = graph.build({
  node(<depth-origin>, pos: pos(x: pin(0), y: pin(0), z: pin(2)))
  edge(source(<depth-origin>), pos: pos(x: pin(0), y: start(1), z: pin(0)))
})
#let depth-node = graph.nodes(depth-base).first()
#let depth-edge = graph.edges(depth-base).first()
#assert(depth-node.statements.at("pos-z") == "2")
#assert(depth-node.statements.at("pos-z-mode") == "pin")
#assert(depth-edge.statements.at("pos-z") == "0")
#let depth-patched = graph.map(depth-base, node: node => (
  ..node,
  pos: pos(z: start(-3)),
))
#let patched-node = graph.nodes(depth-patched).first()
#assert(patched-node.pos == depth-node.pos)
#assert(graph.dot(depth-patched).contains("pin=\"x:0,y:0\""))
#assert(patched-node.statements.at("pos-z") == "-3")
#assert(patched-node.statements.at("pos-z-mode") == "start")
#let depth-options = (
  viewport-w: 1,
  viewport-h: 1,
  spring: (strength: 1, length: 1),
  repulsion: (strength: 0),
  labels: (steps: 0),
  solver: (
    steps: 2, epochs: 1, step: 0.1, cooling: 1,
    depth-scale: 1, flattening-end: 1,
  ),
)
#let depth-solved = layout(depth-base, ..depth-options)
#let flat-solved = layout(depth-base, ..depth-options, solver: (
  ..depth-options.solver, depth-scale: 0,
))
#assert(
  graph.edges(depth-solved).first().pos.y < graph.edges(flat-solved).first().pos.y,
  message: "nonzero pinned auxiliary depth must affect the early XY forces",
)
#assert(graph.nodes(depth-solved).first().pos == (x: 0, y: 0))
#assert(graph.edges(depth-solved).first().pos.x == 0)
#assert(graph.nodes(depth-solved).first().statements.at("pos-z") == "2")
#assert(graph.nodes(depth-solved).first().statements.at("pos-z-mode") == "pin")

// Edge spring lengths are statements, not user data; none preserves inheritance.
#let spring-items = {
  node(<spring-a>)
  node(<spring-b>)
  edge(source(<spring-a>), sink(<spring-b>))
  edge(source(<spring-a>), sink(<spring-b>), spring-length: none)
  edge(
    source(<spring-a>), sink(<spring-b>),
    statements: (spring-length: 0.75),
  )
  edge(
    source(<spring-a>), sink(<spring-b>),
    spring-length: 1.5,
    statements: (spring-length: 3),
    kind: "spring-data",
  )
  edge(source(<spring-a>), spring-length: 1)
}
#let spring-records = graph.edges(graph.build(spring-items))
#assert(spring-records.map(edge => edge.statements.at(
  "spring-length", default: none,
)) == (none, none, "0.75", "1.5", "1"))
#assert(spring-records.at(3).data == (kind: "spring-data"))
#let inherited-spring-records = graph.edges(graph.build(
  spring-items,
  default-edge-statements: (spring-length: 2),
))
#assert(inherited-spring-records.map(edge => edge.statements.at(
  "spring-length",
)) == ("2", "2", "0.75", "1.5", "1"))

// With other forces disabled, dangling endpoints approach their scaled rest lengths.
#let spring-probe = graph.build({
  node(<spring-origin>, pos: pos(x: pin(0), y: pin(0)))
  for factor in (0.5, 2) {
    edge(
      source(<spring-origin>),
      spring-length: factor,
      pos: pos(x: start(2), y: pin(0)),
    )
  }
})
#for algorithm in ("force", "anneal") {
  let positions = graph.edges(layout(
    spring-probe,
    viewport-w: 2,
    viewport-h: 2,
    spring: (strength: 1, length: 0.5),
    repulsion: (strength: 0),
    labels: (steps: 0),
    solver: (
      algorithm: algorithm,
      steps: 100,
      epochs: 30,
      step: 0.1,
      cooling: 0.95,
      temperature: 0.001,
    ),
  )).map(edge => edge.pos)
  for (index, expected) in (1, 4).enumerate() {
    assert(close(calc.abs(positions.at(index).x), expected, epsilon: 0.1))
    assert(positions.at(index).y == 0)
  }
}

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

// The opened xbox has three freely ordered y groups, not pinned rows. A 3D
// layout alone left cut and d1 only 0.20764 units apart after projection.
#let opened-planar-base = graph.build({
  let in-x = group("in", side: "-", start: -7)
  let out-x = group("out", side: "+", start: 7)
  let d1 = group("d1", start: -5)
  let d2 = group("d2", start: 5)
  let cut = group("cut", start: 0)
  node(<opened-a>, pos: pos(x: start(0)))
  node(<opened-b>, pos: pos(x: start(-2)))
  node(<opened-c>, pos: pos(x: start(-2)))
  node(<opened-d>, pos: pos(x: start(2)))
  edge(source(<opened-a>), sink(<opened-c>))
  edge(
    source(<opened-b>), sink(<opened-a>),
    pos: pos(x: group("lower-inner", side: "+", start: 0.5)),
  )
  edge(sink(<opened-b>), pos: pos(x: in-x, y: d1))
  edge(sink(<opened-c>), pos: pos(x: in-x, y: d2))
  edge(sink(<opened-d>), pos: pos(x: out-x, y: d2))
  edge(source(<opened-d>), pos: pos(x: out-x, y: d1))
  edge(source(<opened-a>), sink(<opened-d>), spring-length: 0.1)
  edge(source(<opened-b>), pos: pos(x: out-x, y: cut))
  edge(sink(<opened-c>), pos: pos(x: in-x, y: cut))
})
#let opened-planar-options = layouts.options(
  spring: (strength: 18, length: 0.25),
  repulsion: (
    strength: 12,
    centering: 0.0005,
    edge-node: 0.29,
    edge-edge: 0.1,
    dangling: 14.75,
    dangling-centroid: 25,
  ),
  constraints: (side-strength: 10),
  labels: (steps: 0),
  solver: (
    algorithm: "force",
    seed: 2,
    steps: 50,
    epochs: 50,
    step: 0.81,
    max-movement: 0.4,
    cooling: 0.85,
    depth-scale: 1,
    flattening-end: 0.5,
  ),
)
#let opened-planar = layout(opened-planar-base, ..opened-planar-options)
#let opened-planar-edges = graph.edges(opened-planar).map(edge => edge.pos)
#for (left, right) in ((2, 5), (3, 4), (7, 8)) {
  assert(close(opened-planar-edges.at(left).y, opened-planar-edges.at(right).y))
}
#for (first, second) in ((2, 3), (2, 8), (4, 5), (4, 7)) {
  assert(close(opened-planar-edges.at(first).x, opened-planar-edges.at(second).x))
}
#for row in (2, 3) {
  assert(
    calc.abs(opened-planar-edges.at(7).y - opened-planar-edges.at(row).y) > 1,
    message: "opened xbox: projected cut/row separation must exceed one unit",
  )
}
#assert(opened-planar-edges.at(1).x >= 0)
#assert(opened-planar-edges.at(2).x <= 0)
#assert(opened-planar-edges.at(4).x >= 0)
#let opened-planar-positions = (
  graph.nodes(opened-planar) + graph.edges(opened-planar)
).map(record => record.pos)
#for point in opened-planar-positions {
  assert(calc.abs(point.x) < calc.inf and calc.abs(point.y) < calc.inf)
}
#let opened-planar-repeat = layout(opened-planar-base, ..opened-planar-options)
#assert(opened-planar-positions == (
  graph.nodes(opened-planar-repeat) + graph.edges(opened-planar-repeat)
).map(record => record.pos))

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

// The shorthand emits only supplied fields and preserves their native values.
#assert.eq(feynman.momentum(), (:))
#assert.eq(feynman.momentum(label: (:)), (:))
#for (actual, expected) in (
  (feynman.momentum(side: auto), (momentum-arrow-side: auto)),
  (feynman.momentum(side: "left"), (momentum-arrow-side: "left")),
  (feynman.momentum(offset: -.4), (momentum-arrow-offset: -.4)),
  (feynman.momentum(length: "1.4"), (momentum-arrow-length: "1.4")),
  (feynman.momentum(shift: 1.5), (momentum-arrow-shift: 1.5)),
  (feynman.momentum(label: (gap: .2)), (momentum-label-gap: .2)),
  (feynman.momentum(label: (shift: -.75)), (momentum-label-shift: -.75)),
  (feynman.momentum(label: (anchor: "south-west")), (momentum-label-anchor: "south-west")),
) {
  assert.eq(actual, expected)
}
#let compact-momentum = feynman.momentum(
  side: "right", offset: -.4, length: .7, shift: -.4,
  label: (gap: .2, shift: -.75, anchor: "south-west"),
)
#let expanded-momentum = (
  momentum-arrow-side: "right", momentum-arrow-offset: -.4,
  momentum-arrow-length: .7, momentum-arrow-shift: -.4,
  momentum-label-gap: .2, momentum-label-shift: -.75,
  momentum-label-anchor: "south-west",
)
#assert.eq(compact-momentum, expanded-momentum)
#assert.eq(
  feynman.edge-style((momentum: [$k-p_2$], fields: compact-momentum)),
  feynman.edge-style((momentum: [$k-p_2$], fields: expanded-momentum)),
)
#let shorthand-graph = graph.build(
  default-edge-data: feynman.momentum(offset: .4, label: (shift: 1)),
  { node(<shorthand-node>); edge(<shorthand-edge>, source(<shorthand-node>), ..feynman.momentum(label: (gap: .2))) },
)
#let shorthand-graph = graph.map(shorthand-graph, edge: (
  "shorthand-edge": feynman.momentum(side: auto, label: (gap: 0)),
))
#assert.eq(graph.edges(shorthand-graph).first().data,
  feynman.momentum(side: auto, offset: .4, label: (gap: 0, shift: 1)),
)

// Explicit momentum sides override the old layout-label side for both the arrow
// and the full label carrier, including equal and independent requested shifts.
#let momentum-graph = graph.build({
  node(<momentum-a>, pos: pos(x: 0, y: 0))
  node(<momentum-b>, pos: pos(x: 6, y: 0))
  for side in ("left", "right") {
    edge(
      source(<momentum-a>), sink(<momentum-b>),
      pos: pos(x: 3, y: 0),
      label-pos: (3, 2),
      orientation: "reversed",
      momentum: label-square(if side == "left" { rgb("#a21caf") } else { rgb("#0369a1") }),
      momentum-arrow-side: side,
      momentum-arrow-offset: -0.8,
      momentum-arrow-shift: 0.5,
      momentum-label-shift: if side == "left" { 0.5 } else { -1 },
      route: "straight-through",
    )
  }
})
#for side in (auto, "auto", "left", "right") {
  for (arrow-shift, label-shift) in (
    (0, 0), (0, 1), (0.5, 0.5), (0.5, none), (-100, none), (100, none),
  ) {
    let fields = (
      momentum-arrow-side: side,
      momentum-arrow-shift: arrow-shift,
    ) + if label-shift == none { (:) } else {
      (momentum-label-shift: label-shift)
    }
    let layers = feynman.edge-style((momentum: [], fields: fields))
    assert(layers.len() == 3)
    let arrow = layers.at(1)
    let label = layers.last()
    assert(arrow.shift == arrow-shift and arrow.at("label", default: none) == none)
    assert.eq(
      (label.length, label.ratio, label.resolve-length), (none, none, "none"),
    )
    assert(label.shift == 0 and label.stroke == none and label.mark == none)
    assert(label.label == [])
    assert.eq(
      label.label-shift,
      if label-shift == none { arrow-shift } else { label-shift },
    )
    // Arrow length never changes the full label carrier, even for equal shifts.
    for length in (0.2, 2, 100) {
      let customized = feynman.edge-style((
        momentum: [],
        fields: fields + (momentum-arrow-length: length),
      ))
      assert(customized == layers.enumerate().map(((index, layer)) => {
        if index == 1 { layer + (length: length) } else { layer }
      }))
    }
    let automatic = side in (auto, "auto")
    assert(arrow.offset-side == if automatic { "label" } else { none })
    assert(arrow.label-side == if automatic { auto } else { side })
    assert(arrow.offset == if side == "right" { -0.62 } else { 0.62 })
    assert(label.offset == arrow.offset and label.label-side == arrow.label-side)
    assert(label.offset-side == arrow.offset-side)
    assert(arrow.label-gap == 0.45 and label.label-gap == 0.45)
    assert(arrow.label-style.anchor == auto and label.label-style.anchor == auto)
    // Anchor overrides reach both arrow and label layers without changing
    // shaft geometry.
    for anchor in (
      auto, "auto", "\"auto\"", "center", "\"center\"", "east",
      "south-west", "north-east", "\"south-west\"",
    ) {
      let customized = feynman.edge-style((
        momentum: [],
        fields: fields + (momentum-label-anchor: anchor),
      ))
      assert(customized == layers.enumerate().map(((index, layer)) => {
        if index == 0 { layer } else {
          layer + (label-style: (
            anchor: if anchor == auto { auto } else { anchor.trim("\"") },
          ))
        }
      }))
    }
    // Gap overrides must reach both arrow and label layers without moving the arrow.
    for gap in (0, 0.15, 0.8, "0.2") {
      let customized = feynman.edge-style((
        momentum: [],
        fields: fields + (momentum-label-gap: gap),
      ))
      assert(customized == layers.enumerate().map(((index, layer)) => {
        if index == 0 { layer } else { layer + (label-gap: float(gap)) }
      }))
    }
  }
}

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
// Overlays cover both nodes and the edge, and extend beyond the diagram in
// both axes. Horizontal cubics reuse the edge SVG span probes.
#let after-graph = layout(graph.build({
  node(<after-left>, pos: pos(x: -2, y: 1, mode: "pin"))
  node(<after-right>, pos: pos(x: 2, y: 1, mode: "pin"))
  edge(
    source(<after-left>), sink(<after-right>),
    pos: pos(x: 0, y: 1, mode: "pin"),
  )
}), labels: (steps: 0), solver: (steps: 2, seed: 17))
#let after-options = bare-draw + (
  padding: 0.25,
  node-style: (radius: 0.5, fill: rgb("#1e3a8a"), stroke: none),
  edge-style: (stroke: rgb("#b45309") + 0.5pt),
)
#let after-static = draw(
  after-graph,
  ..after-options,
  unit: 10pt,
  draw-after: {
    cetz.draw.rect((-3, 0), (3, 2), fill: rgb("#f472b6"), stroke: none)
    cetz.draw.bezier(
      (-2, 1), (2, 1), (0, 1), (0, 1),
      stroke: rgb("#9f1239") + 0.8pt,
    )
  },
)
#let after-callback = draw(
  after-graph,
  ..after-options,
  unit: 20pt,
  node-style: after-options.node-style + (fill: rgb("#172554")),
  edge-style: (stroke: rgb("#92400e") + 0.5pt),
  draw-after: g => {
    assert(g == after-graph, message: "draw-after must receive the input graph")
    let nodes = graph.nodes(g)
    let edges = graph.edges(g)
    assert(nodes.len() == 2 and edges.len() == 1)
    let left = nodes.first().pos
    let right = nodes.last().pos
    let middle = edges.first().pos
    assert(same-pos(left, (x: -2, y: 1)))
    assert(same-pos(right, (x: 2, y: 1)))
    assert(same-pos(middle, (x: 0, y: 1)))
    cetz.draw.rect(
      (left.x - 1, left.y - 1), (right.x + 1, right.y + 1),
      fill: rgb("#22d3ee"), stroke: none,
    )
    cetz.draw.bezier(
      (left.x, left.y), (right.x, right.y),
      (middle.x, middle.y), (middle.x, middle.y),
      stroke: rgb("#6d28d9") + 0.8pt,
    )
  },
)
#context {
  for (unit, panel) in ((10pt, after-static), (20pt, after-callback)) {
    let baseline = measure(draw(after-graph, ..after-options, unit: unit))
    for empty in (none, ()) {
      assert(baseline == measure(draw(
        after-graph, ..after-options, unit: unit, draw-after: empty,
      )), message: "empty draw-after must preserve diagram bounds")
    }
    assert(close(baseline.width / unit, 5.5))
    assert(close(baseline.height / unit, 1.5))
    let overlaid = measure(panel)
    assert(close(overlaid.width / unit, 6.5), message: "draw-after must extend canvas width")
    assert(close(overlaid.height / unit, 2.5), message: "draw-after must extend canvas height")
  }
}

#stack(
  spacing: 12pt,
  draw(arrow-graph, edge-style: arrow-style, ..bare-draw),
  draw(crossing-graph, ..bare-draw),
  draw(labels-graph, edge-style: label-style, ..bare-draw),
  draw(momentum-graph, edge-style: edge => {
    let layers = feynman.edge-style(edge)
    layers.at(0).stroke = rgb("#78716c") + 0.3pt
    layers.at(1).stroke = if edge.fields.momentum-arrow-side == "left" {
      rgb("#86198f") + 0.8pt
    } else { rgb("#075985") + 0.8pt }
    layers
  }, ..bare-draw),
  draw(
    precedence-graph,
    edge-style: precedence-base,
    source-style: precedence-endpoint,
    sink-style: precedence-endpoint,
    ..bare-draw,
  ),
  after-static,
  after-callback,
)
