#import "crates/linnest/typst/src/lib.typ": graph, draw, subgraph
#import graph: node, edge, source, sink, pos
#set page(width: auto, height: auto, margin: 0pt)

#context {
  let mode = sys.inputs.at("style-case", default: "defaults")
  let unplaced = mode in ("carrier-auto", "carrier-left", "carrier-auto-right", "carrier-right")
  let stage = (
    defaults: 0, "draw-edge": 1, "record-edge": 2, "draw-endpoint": 3,
    "record-endpoint": 4, "half-endpoint": 5, delegate: 0, hidden: 0,
  ).at(mode, default: -1)
  let blue = rgb("#2563eb")
  let green = rgb("#16a34a")
  let orange = rgb("#ea580c")
  let cyan = rgb("#0891b2")
  let purple = rgb("#9333ea")
  let red = rgb("#dc2626")
  let g = graph.build({
    node(<a>, pos: pos(x: 0, y: 0))
    node(<b>, pos: pos(x: 6, y: 0))
    edge(
      source(<a>, style: if stage >= 5 { (stroke: red) } else { auto }),
      sink(<b>),
      pos: pos(x: 3, y: if unplaced { 1 } else { 0 }),
      label-pos: if unplaced { none } else { pos(x: 3, y: 2) },
      style: if mode == "hidden" { none }
        else if stage >= 2 { (stroke: orange) } else { auto },
      source-style: if stage >= 4 { (stroke: purple) } else { auto },
    )
    if mode.starts-with("carrier") or mode == "ordinary" {
      edge(source(<a>), pos: pos(x: -2, y: 1),
        label-pos: if unplaced { none } else { pos(x: -1, y: 2) })
      edge(sink(<b>), pos: pos(x: 8, y: 1),
        label-pos: if unplaced { none } else { pos(x: 7, y: 2) })
    }
  })
  let base = (stroke: blue + 0.5pt)
  let endpoints = if mode.starts-with("carrier") {
    (base, (
      label-only: true, label: auto,
      offset: if mode.ends-with("right") { -0.6 } else { 0.6 },
      offset-side: if mode in ("carrier-auto", "carrier-auto-right") { "label" } else { none },
      label-side: if mode == "carrier-left" { "left" }
        else if mode == "carrier-right" { "right" } else { auto },
      label-gap: 0.4, label-shift: 0.25,
    ))
  } else if mode == "decoration" {
    (base, (stroke: purple + 0.5pt, offset: 0.6))
  } else { base }
  if mode == "decoration" {
    g = graph.map(g, edge: _ => (decoration: (
      (pattern: none, stroke: green + 0.5pt),
      (pattern: "wave", stroke: cyan + 0.5pt, offset: -0.6),
    )))
  }
  let label = if mode.starts-with("carrier") or mode == "ordinary" {
    box(width: 3pt, height: 3pt, fill: orange, stroke: none)
  } else { none }
  let shared = _ => (stroke: black)
  g = graph.style(g, node-label: none, node-style: (radius: 0),
    edge-style: shared,
    source-style: if mode == "carrier-sink" { base }
      else if mode == "carrier-source" { endpoints.last() } else { endpoints },
    sink-style: if mode == "carrier-source" { base }
      else if mode == "carrier-sink" { endpoints.last() } else { endpoints },
    edge-label: label,
  )
  let stored = graph.info(g).data.at("linnest-style")
  assert(stored.at("edge-style") == shared)
  if mode not in ("carrier-source", "carrier-sink") {
    assert(stored.at("source-style") == endpoints)
    assert(stored.at("sink-style") == endpoints)
  }
  let configured = if mode == "carrier-hostile" {
    ((:), (stroke: red + 4pt, fill: cyan, mark: "triangle", pattern: "wave"))
  } else if mode == "carrier-removed" {
    ((:),)
  } else if mode == "delegate" {
    auto
  } else if stage >= 1 {
    (stroke: green)
  } else { (:) }
  draw(g, unit: 10pt, node-outset: 0, edge-style: configured,
    source-style: if stage >= 3 { (stroke: cyan) } else { none },
    subgraph: if mode in ("carrier-clean", "carrier-hostile") {
      subgraph.select(g, edges: (0, 1, 2))
    } else { none },
    subgraph-edge-style: (stroke: green + 2pt),
  )
}
