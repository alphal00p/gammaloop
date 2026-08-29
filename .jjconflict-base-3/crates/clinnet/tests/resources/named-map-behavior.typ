#import "crates/linnest/typst/src/lib.typ": graph, subgraph
#import graph: edge, node, pin, pos, sink, source, start

// Importing this module executes the assertions without rendering any content.
#{
  let content = [native payload]
  let callback = x => x + 7
  let payload = (label: content, callback: callback, nested: (keep: 1, replace: 2))
  let base = graph.build(name: "named-map", data: payload, {
    node(<a>, pos: pos(x: 1, y: 2), statements: (tag: "node statement"), ..payload)
    node(<b>, pos: pos(x: 5, y: 6), label: [unlisted node])
    node(<noop>, label: [no-op node])
    node(id: 3, label: [unnamed node])
    edge(
      <D1>,
      source(<a>, label: [source], callback: callback),
      sink(<b>, label: [sink], callback: callback),
      pos: pos(x: 3, y: 4),
      statements: (tag: "edge statement"),
      ..payload,
    )
    edge(<keep>, source(<b>), label: [unlisted edge])
    edge(<noop>, source(<b>), label: [no-op edge])
    edge(source(3), label: [unnamed edge])
  })
  let original-nodes = graph.nodes(base)
  let original-edges = graph.edges(base)
  let original-node = original-nodes.find(n => n.name == <a>)
  let original-edge = original-edges.find(e => e.name == <D1>)
  assert(original-nodes.last().at("name", default: none) == none)
  assert(original-edges.last().at("name", default: none) == none)

  for mapper in (none, (:), (a: none, b: _ => none, noop: (:))) {
    assert(graph.map(base, node: mapper) == base)
  }
  for mapper in (none, (:), (D1: none, keep: _ => none, noop: (:))) {
    assert(graph.map(base, edge: mapper) == base)
  }
  assert(graph.map(graph.build(), node: (:), edge: (:)) == graph.build())

  let mapped-spring = graph.map(base, edge: (D1: (spring-length: 1.25)))
  let mapped-spring-edge = graph.edges(mapped-spring).find(e => e.name == <D1>)
  assert(mapped-spring-edge.statements.at("spring-length", default: none) == "1.25")
  assert(mapped-spring-edge.data == payload)
  assert((mapped-spring-edge.data.callback)(3) == 10)

  // Constant patches use the same structural/native split and shallow merge.
  let node-patch = (
    pos: pos(z: pin(7)), statements: (mapped: true),
    label: [mapped node], nested: (replace: 9),
  )
  let edge-patch = (
    pos: pos(z: start(-3)), statements: (mapped: true),
    label: [mapped edge], nested: (replace: 10),
  )
  let mapped = graph.map(base, node: (a: node-patch, noop: none), edge: (D1: edge-patch, noop: none))
  let mapped-node = graph.nodes(mapped).find(n => n.name == <a>)
  let mapped-edge = graph.edges(mapped).find(e => e.name == <D1>)
  for (record, original, patch, z, mode) in (
    (mapped-node, original-node, node-patch, "7", "pin"),
    (mapped-edge, original-edge, edge-patch, "-3", "start"),
  ) {
    assert(record.pos == original.pos)
    assert(record.statements == original.statements + (
      mapped: "true", "pos-z": z, "pos-z-mode": mode,
    ))
    assert(record.data == payload + (label: patch.label, nested: patch.nested))
    assert((record.data.callback)(3) == 10)
    assert("keep" not in record.data.nested)
    assert("pos" not in record.data and "statements" not in record.data)
  }
  assert(graph.nodes(mapped).filter(n => n.name != <a>) == original-nodes.filter(n => n.name != <a>))
  assert(graph.edges(mapped).filter(e => e.name != <D1>) == original-edges.filter(e => e.name != <D1>))
  assert(mapped-edge.source == original-edge.source and mapped-edge.sink == original-edge.sink)
  assert(graph.info(mapped).data == payload)

  // A named closure receives the full original record, not just its native data.
  let node-mapper = n => {
    assert(n == original-node + (fields: n.fields))
    assert(n.fields.node == n.node and n.fields.name == <a>)
    assert(n.fields.pos == n.pos and n.fields.data == payload)
    assert(n.fields.tag == "node statement" and n.fields.label == content)
    assert((n.fields.callback)(3) == 10)
    node-patch
  }
  let edge-mapper = e => {
    assert(e == original-edge + (fields: e.fields))
    assert(e.fields.edge == e.edge and e.fields.name == <D1>)
    assert(e.fields.pos == e.pos and e.fields.orientation == e.orientation)
    assert(e.fields.source == e.source and e.fields.sink == e.sink)
    assert(e.fields.tag == "edge statement" and e.fields.data == payload)
    assert(e.fields.label == content and (e.fields.callback)(3) == 10)
    assert(e.source.data.label == [source] and (e.sink.data.callback)(3) == 10)
    edge-patch
  }
  for node-map in ((a: node-mapper), n => if n.name == <a> { node-mapper(n) }) {
    for edge-map in ((D1: edge-mapper), e => if e.name == <D1> { edge-mapper(e) }) {
      assert(graph.map(base, node: node-map, edge: edge-map) == mapped)
    }
  }
  let all = graph.map(base, node: n => (seen: n.node), edge: e => (seen: e.edge))
  assert(graph.nodes(all).all(n => n.data.seen == n.node))
  assert(graph.edges(all).all(e => e.data.seen == e.edge))

  // Explicit data is replacement, including native content and stored functions.
  for data in (content, callback, (label: [replacement], callback: callback)) {
    let replaced = graph.map(base, node: (a: (data: data, ignored: true)), edge: (
      D1: _ => (data: data, ignored: true),
    ))
    assert(graph.node-data(replaced, <a>) == data and graph.edge-data(replaced, <D1>) == data)
    assert(graph.nodes(replaced).first().pos == original-node.pos)
    assert(graph.edges(replaced).first().statements == original-edge.statements)
  }
  assert(graph.map(base, node: (a: (data: none, ignored: true)), edge: (
    D1: _ => (data: none, ignored: true),
  )) == base)
  let deferred = _ => panic("stored functions must not be invoked by graph.map")
  let stored = graph.map(base, node: (a: (deferred: deferred)), edge: (D1: (deferred: deferred)))
  assert(graph.node-data(stored, <a>).deferred == deferred)
  assert(graph.edge-data(stored, <D1>).deferred == deferred)
  assert(graph.node-data(stored, <a>).label == content)
  assert(graph.edge-data(stored, <D1>).label == content)

  // Returning a full record with a z-only change must not reapply its XY fields.
  let depth = graph.map(base,
    node: (a: n => (..n, pos: pos(z: start(-4)))),
    edge: (D1: e => (..e, pos: pos(z: pin(8)))),
  )
  for (record, original, z, mode) in (
    (graph.nodes(depth).first(), original-node, "-4", "start"),
    (graph.edges(depth).first(), original-edge, "8", "pin"),
  ) {
    assert(record.pos == original.pos and record.data == original.data)
    assert(record.statements == original.statements + ("pos-z": z, "pos-z-mode": mode))
  }

  // Named data updates remain replacements even when their data looks structural.
  let replacement = (pos: [native position], statements: (nested: (value: 1)))
  let updated = graph.update-node-data(base, <a>, replacement)
  assert(graph.node-data(updated, <a>) == replacement)
  assert(graph.nodes(updated).first().pos == original-node.pos)
  let updated = graph.update-edge-data(base, <D1>, (data, record) => {
    assert(data == payload and record.name == <D1>)
    replacement
  })
  assert(graph.edge-data(updated, <D1>) == replacement)
  assert(graph.edges(updated).first().statements == original-edge.statements)

  // Cut fragment names are exact selectors; synthetic nodes have no origin and
  // carry boundary metadata.
  let opened = graph.cut(base,
    left: subgraph.with-data(base, subgraph.select(base, source: (<D1>,)), hedge: (winding: 2)),
    right: subgraph.select(base, sink: (<D1>,)),
  )
  let fragments = graph.edges(opened)
    .filter(e => e.origin.edge == original-edge.edge)
    .sorted(key: e => e.origin.segment)
  assert(fragments.map(e => e.name) == (<D1.0>, <D1.1>, <D1.2>))
  let synthetic-nodes = graph.nodes(opened).filter(n => n.origin == none)
  assert(synthetic-nodes.len() == 2 and synthetic-nodes.all(n => n.boundary != none))
  let cut-mapped = graph.map(opened, node: (a: (label: [cut node])), edge: (
    "D1.0": (label: [left fragment]),
    "D1.1": e => {
      assert(e == fragments.at(1) + (fields: e.fields))
      assert(e.origin == (edge: original-edge.edge, name: <D1>, segment: 1, winding: 2))
      assert(e.fields.origin == e.origin and e.fields.name == <D1.1>)
      assert(e.fields.data == payload and e.source != none and e.sink != none)
      (label: [middle fragment])
    },
    "D1.2": none,
  ))
  assert(graph.edge-data(cut-mapped, <D1.0>).label == [left fragment])
  assert(graph.edge-data(cut-mapped, <D1.1>).label == [middle fragment])
  assert((graph.edge-data(cut-mapped, <D1.0>).callback)(3) == 10)
  assert(graph.edge-data(cut-mapped, <D1.2>) == payload)
  assert(
    graph.nodes(cut-mapped).filter(n => n.name != <a>)
      == graph.nodes(opened).filter(n => n.name != <a>),
  )
  assert(
    graph.edges(cut-mapped).filter(e => e.name not in (<D1.0>, <D1.1>))
      == graph.edges(opened).filter(e => e.name not in (<D1.0>, <D1.1>)),
  )
  assert(graph.boundaries(cut-mapped) == graph.boundaries(opened))
  assert(graph.edge-data(opened, <D1.0>) == payload)
  assert(graph.nodes(base) == original-nodes and graph.edges(base) == original-edges)
}
