// Typst wrapper for the linnest WASM plugin.
//
// Import this file from a directory that also contains `linnest.wasm`:
//
//   #import "linnest.typ": parse-graphs, layout-graph, graph-info
//
// Graph, builder, and subgraph values returned here are archived byte arrays
// intended to be passed back to the functions below.

#let linnest-plugin = plugin("./linnest.wasm")

#let default-layout-config = (
  steps: sys.inputs.at("steps", default: "15"),
  seed: sys.inputs.at("seed", default: "14"),
  step: ".81",
  step_shrink: "0.21",
  temp: ".3",
  beta: "46.1",
  k_spring: "11.",
  g_center: "40.0",
  epochs: "30",
  crossing_penalty: "30",
  gamma_dangling: "40",
  gamma_ee: ".1",
  directional_force: "5",
  label_length_scale: ".6",
  label_spring: "23",
  label_charge: "3",
  label_steps: "20",
  gamma_ev: ".1",
  z_spring: "0.05",
  z_spring_growth: "1.3",
  length_scale: "0.1",
)

// Parse one or more DOT digraphs into archived graph byte arrays.
#let parse-graphs(input) = cbor(linnest-plugin.parse_graph(bytes(input)))

// Build one archived graph byte array from a Typst dictionary:
//
//   (
//     name: "example",
//     nodes: ((name: "a"), (name: "b")),
//     edges: ((source: (node: 0), sink: (node: 1)),),
//   )
#let build-graph(spec) = linnest-plugin.graph_from_spec(cbor.encode(spec))

// Create an archived graph builder. The optional `spec` can set graph-level
// `name`, `statements`, `node_statements`, and `edge_statements`.
#let builder(spec: (:)) = linnest-plugin.graph_builder(cbor.encode(spec))

// Add a node to an archived builder. Returns `(builder: bytes, node: index)`.
#let builder-add-node(builder, name: none, index: none, statements: (:)) = {
  cbor(linnest-plugin.graph_builder_add_node(
    bytes(builder),
    cbor.encode((name: name, index: index, statements: statements)),
  ))
}

// Add a paired or external edge to an archived builder and return the next builder.
// Use `source: none` or `sink: none` for external edges.
#let builder-add-edge(
  builder,
  source: none,
  sink: none,
  orientation: "default",
  flow: none,
  id: none,
  statements: (:),
) = {
  linnest-plugin.graph_builder_add_edge(
    bytes(builder),
    cbor.encode((
      source: source,
      sink: sink,
      orientation: orientation,
      flow: flow,
      id: id,
      statements: statements,
    )),
  )
}

// Finish an archived builder into an archived graph.
#let finish-builder(builder) = linnest-plugin.graph_builder_finish(bytes(builder))

// Layout one parsed graph and return another archived graph byte array.
#let layout-graph(graph, config: default-layout-config) = {
  linnest-plugin.layout_parsed_graph(bytes(graph), cbor.encode(config))
}

// Parse and layout all DOT graphs in `input`.
#let layout-graphs(input, config: default-layout-config) = {
  parse-graphs(input).map(graph => layout-graph(graph, config: config))
}

// Return graph metadata with `name`, `global_statements`,
// `edge_statements`, and `node_statements` fields.
#let graph-info(graph) = cbor(linnest-plugin.graph_info(graph))

// Construct archived subgraph objects.
#let subgraph-from-label(graph, label) = {
  linnest-plugin.graph_archived_subgraph(bytes(graph), cbor.encode(label))
}

#let subgraph-from-bits(graph, bits) = {
  linnest-plugin.graph_archived_subgraph(bytes(graph), cbor.encode(bits))
}

#let subgraph-from-compass(graph, compass) = {
  linnest-plugin.graph_archived_compass_subgraph(bytes(graph), cbor.encode(compass))
}

#let subgraph-label(subgraph) = cbor(linnest-plugin.subgraph_label(bytes(subgraph)))
#let subgraph-hedges(subgraph) = cbor(linnest-plugin.subgraph_hedges(bytes(subgraph)))
#let subgraph-contains-hedge(subgraph, hedge) = {
  cbor(linnest-plugin.subgraph_contains_hedge(bytes(subgraph), cbor.encode(hedge)))
}

// Return node records. With `subgraph`, only nodes incident to that subgraph
// are returned. `subgraph` must be an archived subgraph object.
#let graph-nodes(graph, subgraph: none) = {
  if subgraph == none {
    cbor(linnest-plugin.graph_nodes(graph))
  } else {
    cbor(linnest-plugin.graph_nodes_of_subgraph(bytes(graph), bytes(subgraph)))
  }
}

// Return edge records. With `subgraph`, only edges intersecting that subgraph
// are returned. `subgraph` must be an archived subgraph object.
#let graph-edges(graph, subgraph: none) = {
  if subgraph == none {
    cbor(linnest-plugin.graph_edges(graph))
  } else {
    cbor(linnest-plugin.graph_edges_of_subgraph(bytes(graph), bytes(subgraph)))
  }
}

#let graph-cycle-basis(graph) = cbor(linnest-plugin.graph_cycle_basis(bytes(graph)))
#let graph-spanning-forests(graph) = cbor(linnest-plugin.graph_spanning_forests(bytes(graph)))

#let join-graphs(left, right, key: none) = {
  if key == none {
    panic("join-graphs requires a key")
  }
  linnest-plugin.graph_join_by_hedge_key(bytes(left), bytes(right), cbor.encode((key: key)))
}

// Low-level compatibility helpers that return base62 labels directly.
#let graph-subgraph(graph, subgraph) = cbor(linnest-plugin.graph_subgraph(graph, cbor.encode(subgraph)))
#let graph-compass-subgraph(graph, compass) = {
  cbor(linnest-plugin.graph_compass_subgraph(graph, cbor.encode(compass)))
}
