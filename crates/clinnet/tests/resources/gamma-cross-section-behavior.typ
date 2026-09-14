#import "crates/linnest/typst/src/lib.typ": graph, layout, draw
#import "crates/linnest/typst/src/render/layout.typ" as renderer
#import "gamma-physics-edge-style.typ" as physics
#import "gamma-layout-core.typ": autogen-external-edge-fields, bind-layout
#set page(width: auto, height: auto, margin: 0pt)

#let input = read("epemttbar.dot")
#let original = graph.parse(input).first()
#let gamma-layout = bind-layout(
  graph: graph, renderer: renderer, physics: physics,
  edge-style: (default-edge: physics.default-edge, map: physics.default-map),
)

// Exercise the actual Gamma entrypoint: automatic cross-section detection must
// enable placement, pair by inherited invisible-node cut tags, and pin each row.
#context gamma-layout(input, auto-mode: true,
  edge-style-options: (momentum-arrows: true, show-particle: false),
  layout-passes: ((solver: (seed: 42, steps: 1200, depth-scale: 2), labels: (steps: 2)),),
  diagram-options: (title: none, draw-after: (g, _) => {
    let nodes = graph.nodes(g)
    let edges = graph.edges(g)
    let left = calc.min(..nodes.map(node => node.pos.x))
    let right = calc.max(..nodes.map(node => node.pos.x))
    let labels = graph.info(g).data.at("linnest-style").at("edge-label")
    for (eid, y) in ((0, -5), (1, 5), (2, 5), (3, -5)) {
      let edge = edges.at(eid)
      assert(calc.abs(edge.pos.y - y) < 1e-9)
      assert(edge.statements.at("pos-z") == "0")
      assert(edge.statements.at("pos-z-mode") == "pin")
      if eid < 2 { assert(edge.pos.x > right) }
      else { assert(edge.pos.x < left) }
      // Cut matching never renumbers EMR momentum labels: equal tags still
      // display distinct q_0/q_3 and q_1/q_2 labels.
      let label = repr(labels((eid: eid, data: edge.data)))
      assert(label.contains("base: [q]"))
      assert(label.contains("b: [" + str(eid) + "]"))
    }
    assert(edges.at(0).pos.x == edges.at(1).pos.x)
    assert(edges.at(2).pos.x == edges.at(3).pos.x)
    assert(edges.at(0).statements.at("is_cut") == edges.at(3).statements.at("is_cut"))
    assert(edges.at(1).statements.at("is_cut") == edges.at(2).statements.at("is_cut"))
    ()
  }),
)

// A provided coordinate is authoritative without suppressing the other axis's
// automatic placement. Include a complete XY override and an internal depth pin.
#let partial = graph.map(original, edge: edge => {
  if edge.edge == 0 { (pos: graph.pos(y: -17, z: graph.start(3))) }
  else if edge.edge == 1 { (pos: graph.pos(x: 25, y: 11, z: graph.pin(4))) }
  else if edge.edge == 2 { (pos: graph.pos(x: -30, z: graph.pin(5))) }
  else if edge.edge == 4 { (pos: graph.pos(z: graph.pin(9))) }
})
#let before = graph.edges(partial)
#assert(not before.at(0).at("pos-x-set") and before.at(0).at("pos-y-set"))
#assert(before.at(2).at("pos-x-set") and not before.at(2).at("pos-y-set"))
#let prepared = autogen-external-edge-fields(partial, graph: graph, match-field: "is_cut")
#let prepared-edges = graph.edges(prepared)
#for eid in (0, 1, 2, 3) {
  assert(prepared-edges.at(eid).at("pos-x-set"))
  assert(prepared-edges.at(eid).at("pos-y-set"))
}
#let positioned = layout(prepared,
  solver: (seed: 42, steps: 80, depth-scale: 2), labels: (steps: 0))
#let edges = graph.edges(positioned)
#assert(edges.at(0).pos.y == -17)
#assert(edges.at(0).pos.x > 0)
#assert(edges.at(1).pos == before.at(1).pos)
#assert(edges.at(2).pos.x == -30)
#assert(edges.at(2).pos.y == 5)
#assert(edges.at(3).pos.y == -5)
#for eid in (0, 1, 2, 3) {
  assert(edges.at(eid).statements.at("pos-z") == "0")
  assert(edges.at(eid).statements.at("pos-z-mode") == "pin")
}
#assert(edges.at(4).statements.at("pos-z") == "9")
#assert(edges.at(4).statements.at("pos-z-mode") == "pin")
#draw(positioned, title: none)

// Public Gamma layout keeps outside labels close to the free endpoint, with
// extra room only when momentum is shown. Explicit layout distances still win.
#let spacing-input = ```dot
digraph {
  a [pos="0,0!"];
  exte0 [style=invis];
  exte1 [style=invis];
  exte0 -> a [id=0 particle="e-" is_cut="42" pos="-8,0!"];
  a -> exte1 [id=1 particle="e-" is_cut="42" pos="8,0!"];
}
```.text
#for cross-section in (false, true) {
  for (arrows, momentum, distance, expected) in (
    (true, auto, none, 0.25),
    (true, false, none, 0.10),
    (false, auto, none, 0.10),
    (false, true, none, 0.25),
    (true, false, 0.6, 0.6),
  ) {
    context gamma-layout(spacing-input,
      amplitude-mode: not cross-section, cross-section-mode: cross-section,
      edge-style-options: (momentum-arrows: arrows, show-momentum: momentum),
      layout-passes: ((viewport-w: 1, viewport-h: 1, spring: (length: 1),
        solver: (steps: 0), labels: (steps: 1, repulsion: 0,
          ..if distance == none { (:) } else { (distance: distance) })),),
      diagram-options: (title: none, draw-after: (g, _) => {
        for edge in graph.edges(g) {
          let sign = if edge.source == none { -1 } else { 1 }
          assert(calc.abs(edge.label-pos.x - edge.pos.x - sign * expected) < 1e-9)
          assert(calc.abs(edge.label-pos.y - edge.pos.y) < 1e-9)
        }
        ()
      }),
    )
  }
}
