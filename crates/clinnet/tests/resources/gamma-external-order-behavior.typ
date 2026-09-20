#import "crates/linnest/typst/src/lib.typ": draw, graph, layout
#import "gamma-physics-edge-style.typ" as physics
#import "gamma-layout-core.typ": autogen-external-edge-fields
#set page(width: auto, height: auto, margin: 5pt)

// Edge IDs deliberately disagree with the physical half-edge order on both
// sides. Internal half-edges fill the gaps. Relative ranks, not raw IDs,
// determine Y coordinates.
#let amplitude = (
  graph
    .parse(
      ```dot
      digraph {
        ext [style=invis];
        a [pos="-2,0!"]; b [pos="2,0!"];
        ext -> a:8 [id=0 particle=fermion];
        ext -> a:1 [id=3 particle=fermion];
        ext -> a:4 [id=5 particle=fermion];
        b:9 -> ext [id=1 particle=fermion];
        b:0 -> ext [id=2 particle=fermion];
        b:5 -> ext [id=4 particle=fermion];
        a -> b [id=6 particle=fermion];
        a -> b [id=7 particle=fermion];
      }
      ```.text,
    )
    .first()
)
#let prepared = autogen-external-edge-fields(amplitude, graph: graph)
#let edges = graph.edges(prepared)
#for (eid, hid, y) in ((3, 1, 10), (5, 4, 0), (0, 8, -10)) {
  assert(edges.at(eid).sink.hedge == hid)
  assert(edges.at(eid).pos.y == y)
}
#for (eid, hid, y) in ((2, 0, 10), (4, 5, 0), (1, 9, -10)) {
  assert(edges.at(eid).source.hedge == hid)
  assert(edges.at(eid).pos.y == y)
}

// Existing coordinates remain authoritative per axis, while every external
// depth is replaced by a zero pin. Internal depths and node depths are retained.
#let positioned = graph.map(
  amplitude,
  node: node => (pos: graph.pos(z: graph.pin(2))),
  edge: edge => (
    pos: graph.pos(
      z: graph.pin(if edge.edge < 6 { 7 } else { 9 }),
      ..if edge.edge == 0 { (x: graph.pin(-30)) } else if edge.edge == 1 {
        (y: graph.pin(17))
      } else if edge.edge == 2 { (x: graph.pin(30), y: graph.pin(-17)) } else {
        (:)
      },
    ),
  ),
)
#let positioned = autogen-external-edge-fields(positioned, graph: graph)
#assert(graph.edges(positioned).at(0).pos.y == -10)
#let positioned = layout(
  positioned,
  solver: (algorithm: "force", seed: 42, steps: 8, epochs: 1, depth-scale: 2),
  labels: (steps: 0),
)
#let edges = graph.edges(positioned)
#assert(edges.at(0).pos.x == -30)
#assert(edges.at(1).pos.y == 17)
#assert(edges.at(2).pos.x == 30 and edges.at(2).pos.y == -17)
#for edge in edges {
  assert(edge.statements.at("pos-z") == if edge.edge < 6 { "0" } else { "9" })
  assert(edge.statements.at("pos-z-mode") == "pin")
}
#assert(graph.nodes(positioned).all(node => node.statements.at("pos-z") == "2"))

// Numeric cut tags, not their spelling, incoming edge IDs, or incoming half-edge
// IDs, determine cross-section rows. Both halves of each cut occupy the same row.
#let cross-section = (
  graph
    .parse(
      ```dot
      digraph {
        ext [style=invis];
        a [pos="-2,0!"]; b [pos="2,0!"];
        ext -> a:0 [id=0 is_cut=10 particle=fermion];
        ext -> a:6 [id=1 is_cut=2 particle=fermion];
        b:1 -> ext [id=2 is_cut=2 particle=fermion];
        b:4 -> ext [id=3 is_cut=10 particle=fermion];
        b:5 -> ext [id=4 is_cut=99 particle=fermion];
        a -> b [id=5 particle=fermion];
      }
      ```.text,
    )
    .first()
)
#let cross-section = autogen-external-edge-fields(
  cross-section,
  graph: graph,
  match-field: "is_cut",
)
#let edges = graph.edges(cross-section)
#for (eid, y) in ((0, -5), (1, 5), (2, 5), (3, -5)) {
  assert(edges.at(eid).pos.y == y)
}
// A cut tag without an incoming counterpart does not acquire automatic XY
// placement. The unconditional external depth pin still applies.
#assert(not edges.at(4).at("pos-x-set") and not edges.at(4).at("pos-y-set"))
#for eid in range(5) {
  assert(edges.at(eid).statements.at("pos-z") == "0")
  assert(edges.at(eid).statements.at("pos-z-mode") == "pin")
}
#let cross-section = layout(
  cross-section,
  solver: (algorithm: "force", seed: 42, steps: 8, epochs: 1, depth-scale: 2),
  labels: (steps: 0),
)
#for (eid, y) in ((0, -5), (1, 5), (2, 5), (3, -5)) {
  assert(graph.edges(cross-section).at(eid).pos.y == y)
}

// Placement ranks never renumber q_(eid), including paired and unmatched cuts.
#let styles = physics.style(momentum-arrows: true, show-particle: false)
#for g in (prepared, positioned, cross-section) {
  for edge in graph.edges(g) {
    let label = repr((styles.edge-label)((eid: edge.edge, data: edge.data)))
    assert(label.contains("base: [q]"))
    assert(label.contains("b: [" + str(edge.edge) + "]"))
  }
  context draw(graph.style(g, ..styles), title: none)
}


// FeynKit supplies native graphs to the same pipeline as DOT. Amplitudes use
// half-edge order, whereas cross sections match native numeric sewing IDs.
#import "gamma-layout-core.typ" as physics-layout
#import "crates/linnest/typst/src/render/layout.typ" as renderer
#context for is-cross-section in (false, true) {
  let native = graph.build({
    graph.node(<a>)
    graph.node(<b>)
    graph.edge(<in0>, graph.sink(<a>), particle: "fermion", is_cut: 10)
    graph.edge(<in1>, graph.sink(<a>), particle: "fermion", is_cut: 2)
    graph.edge(graph.source(<b>), <out0>, particle: "fermion", is_cut: 2)
    graph.edge(graph.source(<b>), <out1>, particle: "fermion", is_cut: 10)
    graph.edge(graph.source(<a>), <internal>, graph.sink(<b>), particle: "photon")
  })
  physics-layout.layout(
    native,
    graph: graph,
    renderer: (
      attach-elements: renderer.attach-elements,
      layout-graph: (config, g) => {
        let edges = graph.edges(g)
        let expected = if is-cross-section { (-5, 5, 5, -5) } else { (5, -5, 5, -5) }
        for (i, y) in expected.enumerate() {
          assert(edges.at(i).pos.y == y)
          assert(edges.at(i).statements.at("pos-z-mode") == "pin")
          assert(edges.at(i).statements.at("pos-z") == "0")
        }
        renderer.layout-graph(config, g)
      },
    ),
    physics: physics,
    edge-style: (map: (:), default-edge: physics.default-edge),
    amplitude-mode: not is-cross-section,
    cross-section-mode: is-cross-section,
  )
}
