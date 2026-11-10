#import "crates/linnest/typst/src/lib.typ": graph, draw, curve
#import "gamma-physics-edge-style.typ" as physics
#import "gamma-layout-core.typ": autogen-external-edge-fields
#set page(width: auto, height: auto, margin: 0pt)

#assert.eq(physics.momentum-arrow-defaults.offset, 0.35)
#assert.eq(physics.momentum-arrow-defaults.length, 1.4)
#assert.eq(physics.momentum-arrow-defaults.ratio, 0.5)
#assert.eq(physics.fermion-arrow-mark, curve.mark.triangle(length: 0.21cm, width: 0.1575cm))
#assert.eq(physics.momentum-arrow-defaults.mark, curve.mark.straight(length: 0.16cm, width: 0.12cm))
#assert.eq(physics.coil.pattern-natural-endpoints, true)
#assert.eq(physics.coil.pattern-amplitude, 0.15)
#assert.eq(physics.coil.pattern-wavelength, 0.45)
#assert.eq(physics.coil.pattern-coil-longitudinal-scale, 1.4)

#context {
  let mode = sys.inputs.at("momentum-case", default: "short")
  let outside = mode.starts-with("outside-")
  let incoming = mode == "incoming" or mode.starts-with("outside-incoming")
  let outgoing = (mode == "outgoing" or mode.starts-with("outside-outgoing")
    or mode.starts-with("outside-unmatched"))
  let blue = rgb("#2563eb")
  let green = rgb("#16a34a")
  let red = rgb("#dc2626")
  let orange = rgb("#ea580c")
  let explicit-label = box(width: 4pt, height: 3pt, fill: orange, stroke: none)
  let g = graph.build({
    graph.node(<a>, pos: graph.pos(x: 0, y: 0))
    graph.node(<b>, pos: graph.pos(x: 8, y: 0))
    graph.edge(
      ..if incoming { () } else { (graph.source(<a>),) },
      ..if outgoing { () } else { (graph.sink(<b>),) },
      pos: graph.pos(
        x: if incoming { 0 } else if outgoing { 8 } else { 4 },
        y: if mode.starts-with("curved") { 2 } else { 0 },
      ),
      label-pos: graph.pos(
        x: if outside { if incoming { -1 } else { 9 } } else { 4 },
        y: if outside { 0 } else { 2 },
      ),
      orientation: if mode == "reverse" { "reversed" }
        else if mode == "undirected" { "undirected" } else { "default" },
      particle: "probe", is_cut: "unmatched",
      ..if outside { (:) } else { (label: explicit-label) },
    )
  })
  if outside {
    g = autogen-external-edge-fields(g, graph: graph, place: false,
      match-field: if mode.starts-with("outside-unmatched") { "is_cut" } else { none })
  }
  let line = (stroke: blue + 0.5pt)
  let callbacks = physics.style(
    map: (probe: (
      source: line, sink: line, label: [generated particle label],
      fermion-arrow: true,
      fermion-arrow-mark: curve.mark.triangle(
        length: 0.1cm, width: 0.075cm, fill: green, stroke: none,
      ),
    )),
    momentum-arrows: mode != "ordinary" and not mode.ends-with("-ordinary"),
    show-momentum: true,
    // Keep this geometry probe independent of the public default offset and gap.
    momentum-arrow-offset: 0.62,
    momentum-label-gap: 0.45,
    momentum-arrow-length: if mode.ends-with("long") { 2.4 } else { 0.6 },
    momentum-arrow-shift: if mode == "arrow-shift" { 1 }
      else if outside and mode != "outside-incoming-shift" { 0 } else { 0.5 },
    momentum-label-shift: if mode == "label-shift" { -1 }
      else if mode == "outside-outgoing-label-shift" { 0 }
      else if mode == "label-default" or outside { auto } else { 0.5 },
    momentum-label-anchor: if mode == "anchor-center" { "center" }
      else if mode == "outside-outgoing-anchor" { "west" } else { auto },
    momentum-arrow-side: if mode == "right" { "right" } else { "left" },
    momentum-arrow-stroke: (paint: red, thickness: 0.4pt, cap: "round"),
    momentum-arrow-mark: if mode == "no-mark" { none }
      else if mode == "default-mark" { auto }
      else { curve.mark.triangle(
        length: 0.16cm, width: 0.12cm, fill: red, stroke: none,
      ) },
  )
  if outside {
    let generated-label = callbacks.edge-label
    callbacks.edge-label = edge => {
      let label = generated-label(edge)
      assert(repr(label).contains("generated particle label"))
      assert(repr(label).contains("q") and repr(label).contains(str(edge.eid)))
      box(fill: orange, stroke: none, label)
    }
  }
  g = graph.style(g, unit: 10pt, node-label: none,
    node-style: (radius: 0, stroke: none, fill: none), ..callbacks)
  draw(g, title: none, node-outset: 0,
    edge-style: if mode == "fallback" { ((:),) } else { (:) },
  )
}
