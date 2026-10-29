#import "../../../linnest/typst/src/graph.typ" as graph
#import "../../../linnest/typst/src/render/layout.typ" as renderer

#let original = graph.build({
  graph.node(<left>, pos: graph.pos(x: 0, y: 0, mode: "pin"), label: [left])
  graph.node(<right>, pos: graph.pos(x: 4, y: 0, mode: "pin"))
  graph.edge(graph.source(<left>), <middle>, graph.sink(<right>),
    pos: graph.pos(x: 2, y: 0, mode: "pin"), label: [before], keep: 7)
})

// Empty indexed collections and omitted callbacks preserve the complete graph.
#assert(graph.map(original) == original)
#for elements in (none, (:), (graph: none, nodes: (), edges: (:), hedges: none),
  (nodes: (none,), edges: ("0": none))) {
  assert(renderer.attach-elements(original, elements) == original)
}
#let cleared = renderer.attach-elements(original, (
  nodes: ((label: none),), edges: ("0": (label: none)),
))
#assert(graph.nodes(cleared).first().data.label == none)
#assert(graph.edges(cleared).first().data.label == none)
#assert(graph.edges(cleared).first().data.keep == 7)
#assert(graph.nodes(original).first().data.label == [left])

// Half-edge callbacks see the original edge record, including inherited fields,
// even when the same map changes that edge's native data and coordinates.
#let mapped = graph.map(original,
  graph: record => (checked: true),
  node: node => (label: none),
  edge: edge => (label: [after], pos: graph.pos(y: 1)),
  source: half => {
    assert(half.edge.data.label == [before])
    assert(half.fields.keep == 7)
    assert(half.edge.pos.y == 0)
    (statement: "source", observed: half.fields.keep)
  },
  sink: half => {
    assert(half.edge.data.label == [before])
    assert(half.edge.pos.y == 0)
    (statement: "sink", observed: half.fields.keep)
  },
)
#let edge = graph.edges(mapped).first()
#assert(graph.info(mapped).data.checked)
#assert(graph.nodes(mapped).first().data.label == none)
#assert(edge.data.label == [after])
#assert(edge.pos.x == 2 and edge.pos.y == 1)
#assert(edge.source.statement == "source" and edge.sink.statement == "sink")
#assert(edge.source.data.observed == 7 and edge.sink.data.observed == 7)

// Explicit zero constraints and offsets still run; explicit none is a no-op.
#let inactive = graph.map(original,
  node: _ => (rank: none, minimum-size: none, maximum-size: none),
  edge: _ => (minimum-length: none, label-offset: none),
)
#assert(renderer._apply-element-layout-constraints(inactive) == inactive)
#assert(renderer._apply-label-offsets(inactive) == inactive)
#let constrained = graph.map(original,
  node: node => if node.node == 0 {
    (rank: 0, minimum-size: 2)
  } else { (maximum-size: 0) },
  edge: _ => (minimum-length: 0, label-offset: 0.5,
    label-pos: graph.pos(x: 2, y: 0), statements: (layout-label-gap: 0.2)),
)
#let constrained = renderer._apply-element-layout-constraints(constrained)
#let nodes = graph.nodes(constrained)
#assert(nodes.first().statements.at("layout-rank") == "0")
#assert(nodes.first().statements.at("layout-width") == "2")
#assert(nodes.last().statements.at("layout-width") == "0")
#assert(graph.edges(constrained).first().statements.minlen == "0")
#let offset = renderer._apply-label-offsets(constrained)
#let edge = graph.edges(offset).first()
#assert(edge.label-pos.x == 2 and edge.label-pos.y == 0.5)
#assert(calc.abs(float(edge.statements.at("layout-label-gap")) - 0.7) < 1e-12)
#let zero = graph.map(constrained, edge: _ => (label-offset: 0))
#assert(renderer._apply-label-offsets(zero) == zero)

[Preparation preserves graph data, constraints, and label offsets.]
