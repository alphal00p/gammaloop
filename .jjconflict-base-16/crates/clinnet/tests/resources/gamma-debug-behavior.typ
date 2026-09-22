#import "crates/linnest/typst/src/lib.typ": graph
#import "crates/linnest/typst/src/render/layout.typ" as renderer
#import "gamma-physics-edge-style.typ" as physics
#import "gamma-layout-core.typ": bind-render
#set page(width: auto, height: auto)

#let render = bind-render(
  graph: graph,
  renderer: renderer,
  physics: physics,
  edge-style: (default-edge: physics.default-edge, map: physics.default-map),
)

// Observe the ordinary Linnest ID content after actual rendering. Measurement
// alone must not create IDs, and explicit user styles can suppress the defaults.
#for mode in ("plain", "debug", "override") {
  show math.attach: it => {
    [#metadata((mode, repr(it.base), repr(it.b))) <gamma-debug-id>]
    it
  }
  context render((
    data-path: "gamma-debug.dot",
    options: (debug: mode != "plain"),
    layouts: ((steps: 0, label-steps: 0),),
    style: if mode == "override" { (node-label: [custom]) } else { (:) },
    draw: if mode == "override" { (show-half-edge-ids: false) } else { (:) },
  ))
}
#context {
  let ids = query(<gamma-debug-id>).map(it => it.value)
  assert(ids.filter(it => it.at(0) == "plain").len() == 0)
  assert(ids.filter(it => it.at(0) == "override").len() == 0)
  let nodes = ids.filter(it => it.at(0) == "debug" and it.at(1) == "[n]")
  let hedges = ids.filter(it => it.at(0) == "debug" and it.at(1) == "[h]")
  assert(nodes.map(it => it.at(2)).sorted() == ("[0]", "[1]"))
  assert(hedges.map(it => it.at(2)).sorted() == ("[0]", "[1]", "[2]"))
}
