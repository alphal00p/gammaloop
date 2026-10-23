#import "../lib.typ": draw, graph, layouts, subgraph
#import "@preview/cetz:0.5.1" as cetz
#import "layout.typ" as shared


#let diagram-unit = 1.0
#let operator-node-padding = .86
#let leaf-node-padding = 0.46
#let leaf-node-corner-radius = 0.08
#let operator-node-radius = 0.82
#let operator-node-stroke = 0.25pt + rgb("#666666")
#let leaf-node-stroke = 0.22pt + rgb("#aeb4bd")
#let tree-edge-stroke = rgb("#555555")
#let non-tree-source-edge-stroke = rgb("#7f95b8")
#let non-tree-sink-edge-stroke = rgb("#d8dce3")
#let non-tree-anchor-control-distance = 0.01
#let edge-label-fill = white.transparentize(50%)
#let edge-label-radius = 0.08em
#let edge-label-inset = (x: 0.18em, y: 0.06em)
#let tree-edge-label-style = (
  padding: 0.8,
  frame: none,
  fill: none,
  stroke: none,
)

#let edge-depth(edge, depths) = {
  let depth = none
  let source = edge.at("source-half-edge", default: none)
  let sink = edge.at("sink-half-edge", default: none)
  if source == none and edge.keys().contains("source") {
    source = edge.source
  }
  if sink == none and edge.keys().contains("sink") {
    sink = edge.sink
  }
  for half-edge in (source, sink) {
    if half-edge != none {
      let value = depths.at(str(half-edge.node), default: 0)
      if depth == none or value < depth {
        depth = value
      }
    }
  }
  if depth == none { 0 } else { depth }
}

#let record-label-value(record, fallback: none) = {
  let data = record.at("data", default: (:))
  let value = if type(data) == dictionary {
    data.at("value", default: data.at("slot", default: none))
  } else { none }
  if value != none { return value }
  if value == none { value = record.at("label", default: none) }
  let statements = record.at("statements", default: none)
  if value == none and statements != none {
    value = statements.at("label", default: none)
  }
  if value == none {
    value = fallback
  }
  if value == none {
    none
  } else {
    let value = str(value).trim("\"")
    if value == "" and fallback != none { str(fallback) } else { value }
  }
}

#let node-label-value(node) = record-label-value(node, fallback: node.at("name", default: none))

#let display-leaf-label(value) = {
  if value.starts-with("L:") or value.starts-with("T:") {
    value.slice(2)
  } else if value.starts-with("S:") {
    let rest = value.slice(2)
    let parts = rest.split(":")
    if parts.len() > 1 {
      parts.slice(1).join(":")
    } else {
      rest
    }
  } else {
    value
  }
}

#let classified-leaf-kind(value) = {
  if value.starts-with("S:") {
    "scalar"
  } else if value.starts-with("L:") {
    "library"
  } else if ("T:", "TS:", "TT:", "TTS:").any(prefix => value.starts-with(prefix)) {
    "tensor"
  } else {
    "leaf"
  }
}

#let node-kind(record) = {
  let value = record-label-value(record)
  if value == "∏" {
    "product"
  } else if value == "∑" {
    "sum"
  } else {
    let data = record.at("data", default: (:))
    if type(data) == dictionary and data.keys().contains("kind") {
      data.kind
    } else { classified-leaf-kind(value) }
  }
}

#let typed-leaf-kinds = ("scalar", "library", "tensor")

#let is-leaf-node(record) = node-kind(record) in ("leaf", ..typed-leaf-kinds)

#let tree-hedge-set(tree) = {
  let hedges = (:)
  for hedge in subgraph.hedges(tree) {
    hedges.insert(str(hedge), true)
  }
  hedges
}

#let edge-endpoint(edge, key) = {
  let half = edge.at(key + "-half-edge", default: none)
  if half == none { edge.at(key, default: none) } else { half }
}

#let tree-edge-selected(hedges, edge) = {
  let source = edge-endpoint(edge, "source")
  let sink = edge-endpoint(edge, "sink")
  let source-selected = source != none and hedges.at(str(source.hedge), default: false)
  let sink-selected = sink != none and hedges.at(str(sink.hedge), default: false)
  source-selected or sink-selected
}

#let node-label(node) = {
  let data = node.at("data", default: (:))
  let math-source = data.at("label-typst", default: none)
  if math-source != none { return eval("$" + math-source + "$", mode: "markup") }
  let value = node-label-value(node)
  if value == none {
    none
  } else if value == "∏" {
    $product$
  } else if value == "∑" {
    $sum$
  } else {
    let data = node.at("data", default: (:))
    let label = if type(data) == dictionary and data.keys().contains("value") {
      value
    } else { display-leaf-label(value) }
    text(weight: "bold")[#label]
  }
}

#let tree-node-label-style = node => {
  if node-kind(node) in ("product", "sum") {
    (padding: operator-node-padding, frame: none)
  } else if node-kind(node) in typed-leaf-kinds {
    (padding: leaf-node-padding, frame: none)
  } else {
    (padding: leaf-node-padding, frame: "rect", fill: white, stroke: none)
  }
}

#let tree-node-shape-style = node => {
  if node-kind(node) in ("product", "sum") {
    (radius: operator-node-radius)
  } else {
    (:)
  }
}

#let tree-edge-label(edge) = {
  let math-source = edge.at("data", default: (:)).at("label-typst", default: none)
  let value = if math-source != none {
    eval("$" + math-source + "$", mode: "markup")
  } else { record-label-value(edge) }
  if value == none or value == "" {
    none
  } else {
    box(
      inset: edge-label-inset,
      radius: edge-label-radius,
      fill: edge-label-fill,
      stroke: none,
      text(size: 0.72em, fill: rgb("#737985"))[#value],
    )
  }
}

#let typed-leaf-fill(kind) = {
  if kind == "scalar" {
    rgb("#fde8e8")
  } else if kind == "library" {
    rgb("#e7f6e9")
  } else if kind == "tensor" {
    rgb("#e7f0ff")
  } else {
    none
  }
}

#let draw-tree-node(node, box) = {
  let kind = node-kind(node)
  let center = box.center
  let half-width = box.width / 2
  let half-height = box.height / 2
  if kind == "product" {
    cetz.draw.circle(
      center,
      radius: calc.max(box.width, box.height) / 2,
      name: box.name,
      fill: white,
      stroke: operator-node-stroke,
    )
  } else if kind == "sum" {
    cetz.draw.polygon(
      center,
      4,
      angle: 0deg,
      radius: calc.max(box.width, box.height) / 2,
      name: box.name,
      fill: white,
      stroke: operator-node-stroke,
    )
  } else if kind in typed-leaf-kinds {
    cetz.draw.rect(
      (center.at(0) - half-width, center.at(1) - half-height),
      (center.at(0) + half-width, center.at(1) + half-height),
      radius: leaf-node-corner-radius,
      name: box.name,
      fill: typed-leaf-fill(kind),
      stroke: leaf-node-stroke,
    )
  } else {
    cetz.draw.rect(
      (center.at(0) - half-width, center.at(1) - half-height),
      (center.at(0) + half-width, center.at(1) + half-height),
      name: box.name,
      fill: white,
      stroke: none,
    )
  }
  if box.label != none {
    cetz.draw.content(center, box.label, padding: 0, ..box.label-style)
  }
}

#let edge-stroke(edge, tree-hedges, depths) = {
  if tree-edge-selected(tree-hedges, edge.edge) {
    let depth = edge-depth(edge, depths)
    tree-edge-stroke + 1pt / (1 + depth / 4)
  } else {
    non-tree-sink-edge-stroke + 0.38pt
  }
}

#let source-half-edge-stroke(edge, tree-hedges, depths) = {
  if tree-edge-selected(tree-hedges, edge.edge) {
    edge-stroke(edge, tree-hedges, depths)
  } else {
    non-tree-source-edge-stroke + 0.72pt
  }
}

#let sink-half-edge-stroke(edge, tree-hedges, depths) = {
  edge-stroke(edge, tree-hedges, depths)
}

#let tree-endpoint-anchor(edge, role, depths) = {
  if role == "source" { "north" } else { "south" }
}

#let non-tree-endpoint-anchor(node, role) = {
  if node != none and is-leaf-node(node) {
    "south"
  } else {
    auto
  }
}

#let non-tree-route-exit-statements(nodes, edge) = {
  let statements = (:)
  let source = edge.at("source", default: none)
  let source-node = if source == none { none } else { nodes.at(source.node) }
  let source-anchor = non-tree-endpoint-anchor(source-node, "source")
  if source-anchor != auto {
    statements.insert("source-route-exit", source-anchor)
  }

  let sink = edge.at("sink", default: none)
  let sink-node = if sink == none { none } else { nodes.at(sink.node) }
  let sink-anchor = non-tree-endpoint-anchor(sink-node, "sink")
  if sink-anchor != auto {
    statements.insert("sink-route-exit", sink-anchor)
  }

  if statements.len() == 0 { none } else { (statements: statements) }
}

#let source-edge-style(tree-hedges, depths) = edge => {
  let is-tree = tree-edge-selected(tree-hedges, edge.edge)
  (
    stroke: source-half-edge-stroke(edge, tree-hedges, depths),
    route: if is-tree { "direct" } else { "hobby-through" },
    route-points: if is-tree { "ignore" } else { "through" },
    source-anchor: if is-tree {
      tree-endpoint-anchor(edge, "source", depths)
    } else {
      non-tree-endpoint-anchor(edge.source-node, "source")
    },
    anchor-control-distance: if is-tree { auto } else { non-tree-anchor-control-distance },
  )
}

#let sink-edge-style(tree-hedges, depths) = edge => {
  let is-tree = tree-edge-selected(tree-hedges, edge.edge)
  (
    stroke: sink-half-edge-stroke(edge, tree-hedges, depths),
    route: if is-tree { "direct" } else { "hobby-through" },
    route-points: if is-tree { "ignore" } else { "through" },
    sink-anchor: if is-tree {
      tree-endpoint-anchor(edge, "sink", depths)
    } else {
      non-tree-endpoint-anchor(edge.sink-node, "sink")
    },
    anchor-control-distance: if is-tree { auto } else { non-tree-anchor-control-distance },
  )
}

#let non-tree-edge-label(tree-hedges) = edge => {
  if tree-edge-selected(tree-hedges, edge.edge) {
    none
  } else {
    tree-edge-label(edge)
  }
}

// Draw a network graph, distinguishing the expression tree from contraction
// edges. Both native callers and the DOT example supply an explicit subgraph.
#let render(g, tree, root: 0, config: (:)) = {
  set page(width: auto, height: auto, margin: 4pt, fill: none)
  set text(size: 10pt)
  context {
    let tree-hedges = tree-hedge-set(tree)
    let depths = subgraph.node-depths(g, tree)
    let styled-nodes = graph.nodes(g)
    let g = graph.map(g, edge: edge => {
      if tree-edge-selected(tree-hedges, edge) {
        (statements: ("source-route-exit": "north", "sink-route-exit": "south"))
      } else {
        non-tree-route-exit-statements(styled-nodes, edge)
      }
    })
    let style = config.at("style", default: none)
    let style = (
      node-label: node-label,
      node-label-style: tree-node-label-style,
      node-style: tree-node-shape-style,
      edge-label: non-tree-edge-label(tree-hedges),
      edge-label-style: tree-edge-label-style,
      unit: diagram-unit,
    ) + if style == none { (:) } else { style }
    let drawing = config.at("draw", default: none)
    let drawing = (
      unit: diagram-unit,
      draw-node: draw-tree-node,
      source-style: source-edge-style(tree-hedges, depths),
      sink-style: sink-edge-style(tree-hedges, depths),
      edge-omega: 0.35,
      padding: 0.4,
    ) + if drawing == none { (:) } else { drawing }
    let layout-defaults = (
      layout-algo: "dot",
      subgraph: tree,
      layout-roots: (root,),
      tree-dx: 0.35,
      tree-dy: 2.2,
      route-label-width-cap: 0.,
      label-steps: 40,
      internal-label-length-scale: 0.35,
      external-label-length-scale: 0.35,
    )
    shared.layout-graph(config + (
      style: style, draw: drawing, layout-defaults: layout-defaults,
    ), g)
  }
}

// Native Python networks arrive as typed values through RenderConfig, exactly
// like the Feynman graph builder. Preserve IDs instead of reconstructing them
// from a DOT parser's insertion order.
#let render-network(network, config: (:)) = {
  let g = graph.build({
    for item in network.nodes {
      graph.node(id: item.id, kind: item.kind, value: item.value,
        label-typst: item.at("label-typst", default: none),
        inspection: (node: item.id) + item.inspection)
    }
    for item in network.edges {
      let endpoints = ()
      if item.source != none {
        endpoints.push(graph.source(item.source.at(0), id: item.source.at(1)))
      }
      if item.sink != none {
        endpoints.push(graph.sink(item.sink.at(0), id: item.sink.at(1)))
      }
      graph.edge(..endpoints, id: item.id, orientation: item.orientation,
        slot: item.slot, label-typst: item.at("label-typst", default: none),
        inspection: (edge: item.id) + item.inspection)
    }
  })
  let tree = subgraph.select(g, edges: network.edges.filter(item => item.tree).map(item => item.id))
  if config.at("title", default: auto) == auto { config += (title: none) }
  render(g, tree, root: network.root, config: config)
}
