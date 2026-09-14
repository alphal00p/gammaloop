#import "crates/linnest/typst/src/lib.typ": graph, draw
#import "gamma-physics-edge-style.typ" as physics
#set page(width: auto, height: auto, margin: 0pt)

#context {
  let mode = sys.inputs.at("momentum-side-case", default: "full-above")
  let blue = rgb("#2563eb")
  let red = rgb("#dc2626")
  let green = rgb("#16a34a")
  let cases = (
    (bend: 2, kind: "paired", side: auto, expected: "left"),
    (bend: -2, kind: "paired", side: auto, expected: "right"),
    (bend: 0, kind: "paired", side: auto, expected: "left"),
    (bend: 0, kind: "incoming", side: auto, expected: "left"),
    (bend: 0, kind: "outgoing", side: auto, expected: "left"),
    (bend: 2, kind: "reversed", side: auto, expected: "left"),
    (bend: 2, kind: "paired", side: "right", expected: "right"),
    (bend: -2, kind: "paired", side: "left", expected: "left"),
    (bend: 0, kind: "negative", side: auto, expected: "right"),
  )
  let g = graph.build({
    for (index, case) in cases.enumerate() {
      let row = -8 * index
      let a = label("a" + str(index))
      let b = label("b" + str(index))
      graph.node(a, pos: graph.pos(x: 0, y: row))
      graph.node(b, pos: graph.pos(x: 8, y: row))
      graph.edge(
        ..if case.kind == "incoming" { () } else { (graph.source(a),) },
        ..if case.kind == "outgoing" { () } else { (graph.sink(b),) },
        pos: graph.pos(
          x: if case.kind == "incoming" { 0 }
            else if case.kind == "outgoing" { 8 } else { 4 },
          y: row + case.bend,
        ),
        // Cross the entire arc with the native label position. Neither it nor
        // generated label content may select the momentum arrow's side.
        label-pos: graph.pos(x: 4, y: row + if mode.ends-with("above") { 6 } else { -6 }),
        orientation: if case.kind == "reversed" { "reversed" } else { "default" },
        particle: "probe",
        momentum-arrow-side: case.side,
        ..if case.kind == "negative" { (momentum-arrow-offset: -0.35) } else { (:) },
      )
    }
  })
  let line = (stroke: blue + 0.5pt)
  let callbacks = physics.style(
    map: (probe: (
      source: line, sink: line, label: [particle],
      fermion-arrow: true,
      fermion-arrow-mark: (end: (symbol: "triangle", fill: green, stroke: none), scale: 0.5),
    )),
    momentum-arrows: true,
    show-momentum: mode.starts-with("full"),
    show-particle: not mode.starts-with("none"),
    momentum-arrow-stroke: (paint: red, thickness: 0.4pt, cap: "round"),
    momentum-arrow-mark: (end: (symbol: "triangle", fill: red, stroke: none)),
  )
  for half in ("source", "sink") {
    let callback = callbacks.at(half + "-style")
    callbacks.insert(half + "-style", edge => {
      let layers = callback(edge)
      let expected = cases.at(edge.eid).expected
      for layer in layers.slice(1) {
        assert(layer.offset * if expected == "left" { 1 } else { -1 } > 0)
        assert(layer.at("offset-side", default: none) == none)
        assert(layer.label-side == expected)
      }
      layers
    })
  }
  let edge-label = callbacks.edge-label
  callbacks.edge-label = edge => {
    let value = edge-label(edge)
    assert((value == none) == mode.starts-with("none"))
    if mode.starts-with("full") { assert(repr(value).contains("base: [q]")) }
    value
  }
  g = graph.style(g, unit: 10pt, node-label: none,
    node-style: (radius: 0, stroke: none, fill: none), ..callbacks)
  draw(g, title: none, node-outset: 0)
}
