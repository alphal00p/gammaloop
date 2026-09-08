// Internal subgraph boundary, deliberately independent of `graph.typ`.
// Only masks and topology metadata cross Wasm; annotations remain Typst values.

#let _plugin = plugin("../../linnest.wasm")

// Validate the object and archive without serializing either annotation field.
#let subgraph-bytes(s) = {
  assert(
    type(s) == dictionary
      and s.at("linnest-kind", default: none) == "linnest-subgraph"
      and type(s.at("bytes", default: none)) == bytes
      and "data" in s
      and type(s.at("hedge-data", default: none)) == dictionary
      and "topology" in s
      and (s.topology == none or type(s.topology) == bytes),
    message: "subgraph: expected a Linnest subgraph object",
  )
  let keys = cbor(_plugin.subgraph_hedges(s.bytes)).map(str)
  assert(
    s.hedge-data.keys().all(key => key in keys),
    message: "subgraph: hedge-data keys must be selected half-edge IDs in str(id) form",
  )
  s.bytes
}

// Native algorithms may return CBOR arrays of bytes. Unbound wrappers support
// graph-free queries, but graph consumers must supply topology: topology(g).
#let wrap(mask, data: none, hedge-data: (:), topology: none) = {
  let mask = if type(mask) == array { bytes(mask) } else { mask }
  assert(
    type(hedge-data) == dictionary,
    message: "subgraph: hedge-data must be a dictionary",
  )
  let annotations = (:)
  for h in cbor(_plugin.subgraph_hedges(mask)) {
    let key = str(h)
    if key in hedge-data {
      annotations.insert(key, hedge-data.at(key))
    }
  }
  let result = (
    linnest-kind: "linnest-subgraph",
    bytes: mask,
    data: data,
    hedge-data: annotations,
    topology: topology,
  )
  let _ = subgraph-bytes(result)
  result
}

// An opaque CBOR signature of IDs, names, and source/sink incidence, not graph
// identity. Drawing labels, orientation, placement, statements, and user data
// are deliberately excluded. Endpoint names are included when exposed natively.
#let topology(g) = {
  assert(
    type(g) == dictionary
      and g.at("linnest-kind", default: none) == "linnest-graph"
      and type(g.at("bytes", default: none)) == bytes,
    message: "subgraph: expected a Linnest graph object",
  )
  let nodes = cbor(_plugin.graph_nodes(g.bytes)).map(node => (
    node: node.node,
    name: node.at("name", default: none),
  ))
  let edges = cbor(_plugin.graph_edges(g.bytes)).map(edge => {
    let name = edge.at("name", default: none)
    let payload = edge.at("data", default: none)
    if payload != none {
      let payload = cbor(if type(payload) == array { bytes(payload) } else {
        payload
      })
      if (
        type(payload) == dictionary
          and payload.at("name", default: none) != none
      ) {
        name = payload.name
      }
    }
    let result = (edge: edge.edge, name: name)
    for flow in ("source", "sink") {
      let half = edge.at(flow, default: none)
      result.insert(flow, if half == none { none } else {
        (
          hedge: half.hedge,
          node: half.node,
          name: half.at("name", default: none),
        )
      })
    }
    result
  })
  cbor.encode((nodes: nodes, edges: edges))
}

// Graph-aware consumers call this before sending subgraph-bytes(s) to Wasm.
// Return s unchanged; an absent signature is not proof of compatibility.
#let validate(g, s) = {
  let mask = subgraph-bytes(s)
  let signature = topology(g)
  assert(
    s.topology == signature,
    message: "subgraph: topology does not match graph",
  )
  let hedges = cbor(signature)
    .edges
    .map(edge => (edge.source, edge.sink))
    .flatten()
    .filter(half => half != none)
    .map(half => half.hedge)
  assert(
    cbor(_plugin.subgraph_size(mask)) == hedges.len(),
    message: "subgraph: mask size does not match graph",
  )
  assert(
    cbor(_plugin.subgraph_hedges(mask)).all(hedge => hedge in hedges),
    message: "subgraph: selected half-edge ID is outside graph",
  )
  s
}
