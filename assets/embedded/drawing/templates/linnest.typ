// Typst wrapper for the linnest WASM plugin.
//
// Import this file from a directory that also contains `linnest.wasm`:
//
//   #import "linnest.typ": parse-graphs, layout-graph, graph-info
//
// Graph values returned by `parse-graphs`, `layout-graph`, and `layout-graphs`
// are archived byte arrays intended to be passed back to the functions below.

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

// Return node records. With `subgraph`, only nodes incident to that subgraph
// are returned. `subgraph` may be a base62 label or a boolean array.
#let graph-nodes(graph, subgraph: none) = if subgraph == none {
  cbor(linnest-plugin.graph_nodes(graph))
} else {
  cbor(linnest-plugin.graph_nodes_of(graph, cbor.encode(subgraph)))
}

// Return edge records. With `subgraph`, only edges intersecting that subgraph
// are returned. `subgraph` may be a base62 label or a boolean array.
#let graph-edges(graph, subgraph: none) = if subgraph == none {
  cbor(linnest-plugin.graph_edges(graph))
} else {
  cbor(linnest-plugin.graph_edges_of(graph, cbor.encode(subgraph)))
}

// Normalize a subgraph spec into the base62 label used by the Rust API.
#let graph-subgraph(graph, subgraph) = {
  cbor(linnest-plugin.graph_subgraph(graph, cbor.encode(subgraph)))
}

// Build a base62 subgraph label from a DOT compass point such as "n", "s",
// "e", "w", "ne", or none.
#let graph-compass-subgraph(graph, compass) = {
  cbor(linnest-plugin.graph_compass_subgraph(graph, cbor.encode(compass)))
}
