#import "../../../linnest/typst/src/lib.typ": draw, graph, layout, subgraph
#import "../../../linnest/typst/src/render/layout.typ" as renderer
#import "../../../../assets/embedded/drawing/templates/layout-core.typ": (
  autogen-external-edge-fields, bind-layout,
)
#import "../../../../assets/embedded/drawing/templates/physics-edge-style.typ" as physics
#import "../../../../docs/assets/typst/portal-graphs/edge-style.typ" as model-style

#set page(width: auto, height: auto)

#let original = graph.build(
  {
    graph.node(
      <a>,
      pos: graph.pos(x: 0, y: 0, z: graph.pin(2)),
      keep: [node data],
    )
    graph.node(<b>, pos: graph.pos(x: 4, y: 0))
    graph.edge(
      graph.sink(<a>),
      <incoming>,
      pos: graph.pos(x: -3, y: 1, z: graph.pin(4)),
      particle: "photon",
      is_cut: "paired",
    )
    graph.edge(
      graph.source(<a>),
      graph.sink(<b>),
      <internal>,
      pos: graph.pos(x: 2, y: 0, z: graph.start(3)),
      particle: "fermion",
    )
    graph.edge(
      graph.source(<b>),
      <outgoing>,
      pos: graph.pos(x: 7, y: -1, z: graph.start(5)),
      particle: "scalar",
      is_cut: "unmatched",
    )
  },
  data: (keep: [graph data]),
)
#let g = renderer.attach-elements(original, (
  graph: (added: true),
  nodes: ((label: [native label], node-style: _ => (radius: 0.2)),),
  edges: ((:), (spring-length: 7, token: [edge data]), (:)),
  hedges: ((anchor: "west", token: [half data]),),
))
#assert(graph.info(g).data.keep == [graph data])
#assert(graph.info(g).data.added)
#assert(graph.nodes(g).find(node => node.name == <a>).data.keep == [node data])
#assert(
  graph.nodes(g).find(node => node.name == <a>).data.label == [native label],
)
#assert(
  type(graph.nodes(g).find(node => node.name == <a>).data.at("node-style"))
    == function,
)
#assert(
  graph
    .edges(g)
    .find(edge => edge.name == <internal>)
    .statements
    .at("spring-length")
    == "7",
)
#assert(
  graph.edges(g).find(edge => edge.name == <internal>).data.token
    == [edge data],
)
#assert(
  graph.edges(g).find(edge => edge.name == <incoming>).sink.data.token
    == [half data],
)
#let selected = subgraph.select(g, edges: (<internal>,))

// Incoming and outgoing legs keep their explicit x/y coordinates even when
// their cut tags are unmatched. External preparation changes only their depth.
#let prepared = autogen-external-edge-fields(
  g,
  graph: graph,
  match-field: "is_cut",
  place: false,
)
#for name in (<incoming>, <outgoing>) {
  let before = graph.edges(g).find(edge => edge.name == name)
  let after = graph.edges(prepared).find(edge => edge.name == name)
  assert(after.pos == before.pos)
  assert(after.statements.at("pos-z") == "0")
  assert(after.statements.at("pos-z-mode") == "pin")
}
#assert(
  graph
    .edges(prepared)
    .find(edge => edge.name == <internal>)
    .statements
    .at("pos-z")
    == "3",
)
#assert(
  graph
    .edges(prepared)
    .find(edge => edge.name == <internal>)
    .statements
    .at("pos-z-mode")
    == "start",
)
#assert(
  graph.nodes(prepared).find(node => node.name == <a>).statements.at("pos-z")
    == "2",
)
#assert(
  subgraph.hedges(selected)
    == subgraph.hedges(subgraph.select(prepared, edges: (<internal>,))),
)

#let automatic = autogen-external-edge-fields(
  graph.build({
    graph.node(<v>)
    graph.edge(graph.sink(<v>), <in-first>)
    graph.edge(graph.sink(<v>), <in-second>)
    graph.edge(graph.source(<v>), <out>)
  }),
  graph: graph,
)
#let automatic-edges = graph.edges(automatic)
#assert(automatic-edges.at(0).pos.y > automatic-edges.at(1).pos.y)
#for edge in automatic-edges {
  assert(edge.statements.at("pos-z") == "0")
  assert(edge.statements.at("pos-z-mode") == "pin")
}

#let callbacks = physics.style(momentum-arrows: true, show-particle: false)
#for eid in (0, 1, 2, 7) {
  let record = (eid: eid, data: (:), ext: eid != 1)
  let label = (callbacks.edge-label)(record)
  assert(repr(label).contains("q"))
  assert(repr(label).contains(str(eid)))
}
#let labelled = (eid: 9, data: (label: [user label]))
#assert((callbacks.edge-label)(labelled) == [user label])

// Match the generated GammaLoop wrapper contract: the explicit map is spread
// once, then merged above the generated and model-neutral maps.
#let generated = model-style.style(
  map: ("a": (source: (stroke: red + 0.5pt), sink: (stroke: red + 0.5pt))),
  momentum-arrows: true,
)
#let sample = (
  eid: 1,
  particle: "a",
  data: (:),
  orientation: "default",
  source-half-edge: (hedge: 0),
  sink-half-edge: (hedge: 1),
)
#assert((generated.source-style)(sample).first().stroke == red + 0.5pt)

#context {
  let styled = graph.style(
    prepared,
    node-label: none,
    ..physics.style(momentum-arrows: true),
  )
  renderer.layout-graph(
    (
      layouts: none,
      draw: (
        title: none,
        draw-after: (result, bounds) => {
          let styles = graph.info(result).data.at("linnest-style")
          assert(type(styles.at("source-style")) == function)
          assert(
            (styles.at("edge-label"))((eid: 7, data: (:)))
              == (callbacks.edge-label)((eid: 7, data: (:))),
          )
          assert(bounds.width > 0)
        },
      ),
    ),
    styled,
  )
  let positioned = layout(
    styled,
    solver: (steps: 8, depth-scale: 2, flattening-end: 0.8),
    labels: (steps: 3),
  )
  for name in (<incoming>, <outgoing>) {
    let edge = graph.edges(positioned).find(edge => edge.name == name)
    assert(edge.statements.at("pos-z") == "0")
    assert(edge.statements.at("pos-z-mode") == "pin")
    assert(
      edge.pos == graph.edges(prepared).find(edge => edge.name == name).pos,
    )
  }
  draw(positioned, title: none)
  renderer.layout-graph(
    (
      style: physics.style(momentum-arrows: true),
      layouts: ((solver: (steps: 2)), (solver: (steps: 2))),
      draw: (title: none, subgraph: selected),
    ),
    prepared,
  )
}

#let gamma-layout = bind-layout(
  graph: graph,
  renderer: renderer,
  physics: physics,
  edge-style: model-style,
)
#context gamma-layout(
  "digraph demo { v; i [style=invis]; o [style=invis]; i -> v [id=0, particle=a]; v -> o [id=1, particle=d]; }",
  edge-style-options: (momentum-arrows: true),
  layout-passes: ((steps: 3, label-steps: 2),),
)
