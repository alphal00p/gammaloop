#import "../../../linnest/typst/src/lib.typ": graph
#import "../../../linnest/typst/src/render/layout.typ" as renderer
#import "../../../../assets/embedded/drawing/templates/layout-core.typ": autogen-external-edge-fields
#import "../../../../assets/embedded/drawing/templates/physics-edge-style.typ" as physics
#set page(width: auto, height: auto, margin: 5pt)

// Explicit coordinates can put incoming legs on the right and outgoing legs
// on the left. Label boxes must face outward according to their actual geometry.
#let cases = (
  (incoming: true, dx: 4, dy: 0, anchor: "west"),
  (incoming: false, dx: -4, dy: 0, anchor: "east"),
  (incoming: true, dx: -4, dy: 0, anchor: "east"),
  (incoming: false, dx: 4, dy: 0, anchor: "west"),
  (incoming: true, dx: 0, dy: 4, anchor: "south"),
  (incoming: false, dx: 0, dy: -4, anchor: "north"),
  (incoming: true, dx: 4, dy: 0, anchor: "north", explicit: true),
)
#let original = graph.build({
  for (index, case) in cases.enumerate() {
    let node = label("node-" + str(index))
    let y = -10 * index
    graph.node(node, pos: graph.pos(x: graph.pin(0), y: graph.pin(y), z: graph.pin(2)))
    graph.edge(
      if case.incoming { graph.sink(node) } else { graph.source(node) },
      label("external-" + str(index)),
      pos: graph.pos(x: graph.pin(case.dx), y: graph.pin(y + case.dy), z: graph.pin(7)),
      label-pos: (case.dx * 1.15, y + case.dy * 1.15),
      particle: "photon", expected-anchor: case.anchor,
      ..if case.at("explicit", default: false) { (label-anchor: case.anchor) } else { (:) },
    )
  }
  graph.edge(graph.source(<node-0>), graph.sink(<node-1>), <internal>,
    pos: graph.pos(x: graph.pin(1), y: graph.pin(-5), z: graph.pin(9)),
    particle: "fermion")
})
#let prepared = autogen-external-edge-fields(original, graph: graph)
#let before = graph.edges(original)
#for edge in graph.edges(prepared) {
  assert(edge.pos == before.at(edge.edge).pos)
  for key in ("pos-x-set", "pos-y-set") {
    assert(edge.at(key) == before.at(edge.edge).at(key))
  }
  assert(edge.statements.at("pos-z") == if edge.name == <internal> { "9" } else { "0" })
  assert(edge.statements.at("pos-z-mode") == "pin")
}
#assert(graph.nodes(prepared).all(node => node.statements.at("pos-z") == "2"))

#let styles = physics.style(momentum-arrows: true, show-particle: false)
#context renderer.layout-graph((
  style: styles + (node-label: none, edge-label-style: record => {
    let result = (styles.edge-label-style)(record)
    // Measurement has no solved node aliases; anchoring must leave text sizing
    // usable there and resolve only when the renderer supplies the geometry.
    if record.at("source-node", default: none) != none or record.at("sink-node", default: none) != none {
      if record.ext {
        assert(result.anchor == record.at("expected-anchor"), message: repr(record.eid))
        let half-style = if record.at("source-node", default: none) == none {
          styles.sink-style
        } else { styles.source-style }
        let layers = half-style(record)
        assert(layers.len() == 2)
        assert(not layers.any(layer => layer.at("label-only", default: false)))
        // An explicit momentum shift still requests path-relative placement.
        assert(half-style(record + (momentum-label-shift: 0.1)).len() == 3)
      }
    } else if record.ext and not cases.at(record.eid).at("explicit", default: false) {
      assert(result == (:))
    }
    result
  }),
  layouts: none,
  draw: (title: none, draw-after: (g, bounds) => {
    assert(bounds.width > 0 and bounds.height > 0)
    assert(graph.edges(g).len() == 8)
    ()
  }),
), prepared)
