#import "crates/linnest/typst/src/lib.typ": draw, graph, layout, subgraph
#import graph: edge, node, pos, sink, source
#import "@preview/cetz:0.5.1" as cetz

// Render this value: graph.style measurement and draw callbacks need layout.
// Black-only drawing deliberately avoids the public fixture's SVG color probes.
#let weighted-cut-behavior = context {
  let payload = (content: [native payload], callback: x => x + 7)
  let drawing = none
  for orientation in ("default", "reversed", "undirected") {
    let master = graph.build(name: "weighted", data: payload, {
      node(<a>, pos: pos(x: 0, y: 0), payload: payload)
      node(<b>, pos: pos(x: 8, y: 0), payload: payload)
      node(<isolated>, pos: pos(x: 4, y: 2), payload: payload)
      edge(
        <cut-edge>,
        source(<a>, compass: "e", content: [source], callback: x => x + 10),
        sink(<b>, content: [sink], callback: x => x + 20),
        orientation: orientation,
        payload: payload,
      )
      edge(<keep>, sink(<a>, content: [external], callback: x => x + 30), payload: payload)
    })
    let original-nodes = graph.nodes(master)
    let original-edges = graph.edges(master)
    let original = original-edges.first()
    let empty = subgraph.select(master)
    assert(
      graph.cut(master, left: empty, right: empty, boundary: _ => panic("empty cut callback"))
        == master,
    )
    assert(graph.boundaries(master) == ())
    assert(subgraph.hedges(subgraph.select(master, nodes: (<isolated>,))) == ())

    // Selectors union incidences and deduplicate names, indices, and exact hedges.
    let union = subgraph.select(
      master,
      nodes: (<a>, 0),
      edges: (<keep>, 1),
      source: (<cut-edge>, 0),
      sink: (0,),
      hedges: (original.sink.hedge,),
    )
    assert(subgraph.hedges(union) == (0, 1, 2))
    let selected = subgraph.select(master, source: (<cut-edge>,))
    for s in (
      selected,
      subgraph.bits(master, (true, false, false)),
      subgraph.label(master, subgraph.to-label(selected)),
      subgraph.compass(master, "e"),
    ) {
      assert(type(s) == dictionary and s.linnest-kind == "linnest-subgraph")
      assert(type(s.bytes) == bytes and type(s.topology) == bytes)
      assert(s.data == none and s.hedge-data == (:))
      assert(subgraph.hedges(s) == (original.source.hedge,))
      assert(subgraph.contains(s, original.source.hedge))
      assert(not subgraph.contains(s, original.sink.hedge))
      assert(subgraph.contains-edge(s, original))
      assert(not subgraph.contains-edge(s, original-edges.last()))
      assert(graph.nodes(master, subgraph: s).map(n => n.name) == (<a>,))
      assert(graph.edges(master, subgraph: s) == (original,))
    }

    // Overlapping selections own independent native annotations, not graph data.
    let annotated = subgraph.with-data(master, selected, data: payload, hedge: (tag: [first]))
    let enriched = subgraph.with-data(master, annotated, hedge: h => {
      assert(h.edge == original.edge and h.edge-name == <cut-edge> and h.flow == "source")
      assert(h.node == original.source.node and h.hedge == original.source.hedge)
      assert(h.data == (tag: [first]) and h.graph-data == original.source.data)
      (tag: [second], payload: payload, endpoint: h.graph-data)
    })
    assert(subgraph.with-data(master, enriched) == enriched)
    assert(subgraph.hedge-data(selected, original.source.hedge) == none)
    assert(subgraph.hedge-data(annotated, original.source.hedge) == (tag: [first]))
    assert(enriched.data == payload)
    assert(enriched.hedge-data.keys() == (str(original.source.hedge),))
    assert((subgraph.hedge-data(enriched, original.source.hedge).payload.callback)(3) == 10)
    let cleared = subgraph.with-data(master, enriched, data: none, hedge: _ => none)
    assert(cleared.data == none and subgraph.hedge-data(cleared, original.source.hedge) == none)
    let callable = subgraph.with-data(master, selected, hedge: _ => payload.callback)
    assert(subgraph.hedge-data(callable, original.source.hedge)(3) == 10)
    let inverse = subgraph.complement(master, enriched)
    assert(inverse.data == payload and inverse.hedge-data == (:))
    assert(subgraph.hedges(inverse) == (original.sink.hedge, original-edges.last().sink.hedge))
    let twice = subgraph.complement(master, inverse)
    assert(subgraph.hedges(twice) == subgraph.hedges(selected))
    assert(twice.data == payload and twice.hedge-data == (:))
    let relabeled = subgraph.label(master, subgraph.to-label(enriched))
    assert(relabeled.data == none and relabeled.hedge-data == (:))

    // Structural source/sink order is independent of drawing orientation and L/R.
    // n=1 has no winding field; larger n exercise left-only, right-only, and both.
    for n in (1, 2, 3) {
      for left-flow in ("source", "sink") {
        let right-flow = if left-flow == "source" { "sink" } else { "source" }
        let sides = (:)
        for (side, flow) in (("left", left-flow), ("right", right-flow)) {
          let annotation = (side: side, content: [cut annotation], callback: payload.callback)
          if (
            n > 1
              and (
                orientation == "undirected"
                  or side == if orientation == "default" { "left" } else { "right" }
              )
          ) {
            annotation.winding = n
          }
          sides.insert(side, subgraph.with-data(
            master,
            subgraph.select(master, hedges: (original.at(flow).hedge,)),
            data: (side: side, payload: payload),
            hedge: annotation,
          ))
        }
        let view = graph.cut(master, left: sides.left, right: sides.right)
        let nodes = graph.nodes(view)
        let edges = graph.edges(view)
        let boundaries = graph.boundaries(view)
        assert(nodes.len() == 3 + 2 * (n - 1) and edges.len() == n + 2)
        assert(
          edges.map(e => (e.source, e.sink).filter(h => h != none).len()).sum() == 3 + 2 * (n - 1),
        )
        assert(graph.info(view).name == "weighted" and graph.info(view).data == payload)
        assert((graph.info(view).data.callback)(3) == 10)
        assert(boundaries.len() == 2 * n)
        assert(nodes.filter(v => v.origin != none).map(v => v.origin) == (0, 1, 2))
        for v in nodes {
          if v.origin == none {
            assert(v.boundary != none and v.data == none)
          } else {
            assert(v.data == original-nodes.at(v.origin).data and v.boundary == none)
            assert(
              v.data.payload.content == [native payload] and (v.data.payload.callback)(3) == 10,
            )
          }
        }
        let fragments = edges
          .filter(e => e.origin.edge == original.edge)
          .sorted(key: e => e.origin.segment)
        assert(fragments.map(e => e.name) == range(n + 1).map(i => label("cut-edge." + str(i))))
        for (i, e) in fragments.enumerate() {
          assert(e.origin == (edge: original.edge, name: <cut-edge>, segment: i, winding: n))
          assert(e.orientation == orientation and e.data == original.data)
          assert(e.data.payload.content == [native payload] and (e.data.payload.callback)(3) == 10)
          assert((e.source == none) == (i == n) and (e.sink == none) == (i == 0))
          for flow in ("source", "sink") {
            let h = e.at(flow)
            if h != none {
              // A middle source copies the original sink; its sink copies the source.
              let copy-flow = if i == 0 or i == n { flow } else if flow == "source" {
                "sink"
              } else { "source" }
              let old = original.at(copy-flow)
              assert(h.origin == old.hedge and h.data == old.data)
              assert((h.data.callback)(1) == if copy-flow == "source" { 11 } else { 21 })
              if i == 0 or i == n { assert(nodes.at(h.node).origin == old.node) } else {
                assert(nodes.at(h.node).origin == none)
              }
            }
          }
        }
        let kept = edges.find(e => e.name == <keep>)
        assert(kept.origin == (edge: 1, name: <keep>, segment: none, winding: 0))
        assert(kept.data == original-edges.last().data and kept.boundary == none)
        assert(kept.source == none and kept.sink.data == original-edges.last().sink.data)
        assert(
          kept.sink.origin == original-edges.last().sink.hedge
            and (kept.sink.data.callback)(1) == 31,
        )
        for crossing in range(n) {
          let pair = boundaries.filter(b => b.crossing == crossing)
          assert(pair.len() == 2 and pair.map(b => b.side).sorted() == ("left", "right"))
          assert(pair.map(b => b.hedge).sorted() == (original.source.hedge, original.sink.hedge))
        }
        for b in boundaries {
          let e = edges.at(b.edge)
          let record = if b.node == none { e } else { nodes.at(b.node) }
          let flow = if b.hedge == original.source.hedge { "source" } else { "sink" }
          assert(b.side == if flow == left-flow { "left" } else { "right" })
          assert(b.data == subgraph.hedge-data(sides.at(b.side), b.hedge))
          assert(b.cut-data == sides.at(b.side).data and b.origin == e.origin)
          assert(b.pos == record.pos and record.boundary + (pos: b.pos) == b)
          assert(b.origin.segment == b.crossing + if flow == "source" { 0 } else { 1 })
          assert((b.node == none) == (b.origin.segment == 0 or b.origin.segment == n))
        }
        if n == 3 and orientation == "default" and left-flow == "source" {
          // The callback receives the actual record and returns map-style structural/data patches.
          drawing = graph.cut(master, left: sides.left, right: sides.right, boundary: record => {
            let b = record.boundary
            assert(b != none)
            if b.node == none {
              assert(record.edge == b.edge and record.origin == b.origin)
              assert(record.source == none or record.sink == none)
            } else {
              assert(record.node == b.node and record.origin == none)
            }
            (
              pos: pos(x: 2 * b.crossing + 1, y: if b.side == "left" { -1 } else { 1 }),
              data: (payload: payload, side: b.side),
              statements: (boundary-patched: true),
            )
          })
          for b in graph.boundaries(drawing) {
            let record = if b.node == none { graph.edges(drawing).at(b.edge) } else {
              graph.nodes(drawing).at(b.node)
            }
            assert(record.data == (payload: payload, side: b.side))
            assert(record.statements.at("boundary-patched") == "true")
            assert(b.pos == (x: 2 * b.crossing + 1, y: if b.side == "left" { -1 } else { 1 }))
          }
        }
      }
    }
    assert(graph.nodes(master) == original-nodes and graph.edges(master) == original-edges)
    assert(graph.info(master).data == payload and graph.boundaries(master) == ())
  }

  // Boundary coordinates are live after map/layout, not snapshots taken at cut time.
  let before-map = graph.boundaries(drawing)
  drawing = graph.map(
    drawing,
    node: v => if v.boundary != none { (pos: pos(x: v.pos.x + 1, y: v.pos.y + 2, mode: "start")) },
    edge: e => if e.boundary != none { (pos: pos(x: e.pos.x + 1, y: e.pos.y + 2, mode: "start")) },
  )
  for (old, current) in before-map.zip(graph.boundaries(drawing)) {
    assert(current.pos == (x: old.pos.x + 1, y: old.pos.y + 2))
    assert(current == old + (pos: current.pos))
  }
  let styled = graph.style(
    drawing,
    unit: 8pt,
    node-label: d => {
      assert(d.node.boundary == none)
      [physical]
    },
    node-label-style: d => {
      assert(d.node.boundary == none)
      (padding: 0.1)
    },
    node-style: d => {
      assert(d.node.boundary == none)
      (radius: 0.2, fill: black)
    },
  )
  assert(
    graph
      .nodes(styled)
      .filter(v => v.origin == none)
      .all(v => (
        float(v.statements.at("layout-width")) == 0 and float(v.statements.at("layout-height")) == 0
      )),
  )
  let solved = layout(styled, labels: (steps: 0), solver: (steps: 8, seed: 17))
  assert(graph.boundaries(solved).map(b => b.pos) != graph.boundaries(styled).map(b => b.pos))
  for b in graph.boundaries(solved) {
    let record = if b.node == none { graph.edges(solved).at(b.edge) } else {
      graph.nodes(solved).at(b.node)
    }
    assert(b.pos == record.pos)
    assert(
      b
        == graph.boundaries(styled).find(old => old.edge == b.edge and old.node == b.node)
          + (pos: b.pos),
    )
  }
  draw(solved)
  draw(
    solved,
    node-label: d => {
      assert(d.node.boundary == none)
      [physical]
    },
    node-style: d => {
      assert(d.node.boundary == none)
      (radius: 0.2, fill: black)
    },
    draw-node: (v, box) => {
      assert(v.node.boundary == none)
      cetz.draw.circle(box.center, radius: 0.2, fill: black)
    },
  )

  // K4 momenta use coefficients of (k, p1, p2); superficial arrows do not change flow.
  let momenta = ((1, -1, 0), (1, 0, -1), (0, 1, 0), (0, 0, 1), (-1, 1, 1), (1, 0, 0))
  let master = graph.build(data: payload, {
    for (name, point) in ((<a>, (0, 0)), (<b>, (4, 0)), (<c>, (0, 4)), (<d>, (4, 4))) {
      node(name, pos: pos(x: point.at(0), y: point.at(1)), payload: payload)
    }
    for (i, endpoints) in (
      (<c>, <d>),
      (<a>, <b>),
      (<c>, <a>),
      (<d>, <b>),
      (<a>, <d>),
      (<b>, <c>),
    ).enumerate() {
      edge(
        label("D" + str(i + 1)),
        source(endpoints.at(0)),
        sink(endpoints.at(1)),
        orientation: if i in (1, 2) { "reversed" } else { "default" },
        coefficients: momenta.at(i),
        payload: payload,
      )
    }
  })
  let original-edges = graph.edges(master)
  let original-nodes = graph.nodes(master)
  let views = (
    (source: (), sink: (<D3>, <D4>), counts: (1, 1, 2, 2, 1, 1)),
    (source: (<D1>,), sink: (<D4>, <D6>), counts: (2, 1, 1, 2, 1, 2)),
    (source: (<D2>,), sink: (<D3>, <D6>), counts: (1, 2, 2, 1, 1, 2)),
    (source: (<D1>, <D2>), sink: (<D6>,), counts: (2, 2, 1, 1, 1, 3)),
  )
  for (i, spec) in views.enumerate() {
    let sides = (:)
    for side in ("left", "right") {
      let selected = subgraph.select(
        master,
        source: if side == "left" { spec.source } else { spec.sink },
        sink: if side == "left" { spec.sink } else { spec.source },
      )
      sides.insert(side, subgraph.with-data(
        master,
        selected,
        data: (view: i, side: side, payload: payload),
        hedge: h => (
          winding: if i == 3 and h.edge-name == <D6> { 2 } else { 1 },
          flow: h.flow,
          coefficients: momenta.at(h.edge),
          payload: payload,
        ),
      ))
    }
    let view = graph.cut(master, left: sides.left, right: sides.right)
    let edges = graph.edges(view)
    assert(edges.len() == spec.counts.sum())
    assert(graph.nodes(view).len() == if i == 3 { 6 } else { 4 })
    assert(graph.boundaries(view).len() == 2 * (spec.counts.sum() - 6))
    for (j, count) in spec.counts.enumerate() {
      let fragments = edges.filter(e => e.origin.edge == j).sorted(key: e => e.origin.segment)
      assert(fragments.len() == count)
      assert(fragments.all(e => (
        e.data == original-edges.at(j).data and e.orientation == original-edges.at(j).orientation
      )))
      assert(
        fragments.map(e => e.origin.segment) == if count == 1 { (none,) } else { range(count) },
      )
      assert(fragments.all(e => (
        e.origin.winding == count - 1 and e.origin.name == original-edges.at(j).name
      )))
    }
    // Sum the signed momenta of boundary occurrences; D6 contributes twice in view 4.
    for (side, expected) in (("left", (0, -1, -1)), ("right", (0, 1, 1))) {
      let boundaries = graph.boundaries(view).filter(b => b.side == side)
      let coefficients = range(3).map(axis => boundaries
        .map(b => (
          (if b.data.flow == "source" { 1 } else { -1 }) * b.data.coefficients.at(axis)
        ))
        .sum())
      assert(coefficients == expected)
      assert(boundaries.all(b => b.cut-data == sides.at(side).data))
    }
  }
  assert(graph.edges(master) == original-edges and graph.nodes(master) == original-nodes)
  assert(graph.info(master).data == payload and graph.boundaries(master) == ())

  // Algorithm-produced objects interoperate with public queries, layout, and drawing.
  let cycles = graph.cycles(master)
  let forests = graph.forests(master)
  assert(cycles.len() == 3 and forests.len() == 16)
  let all = subgraph.with-data(
    master,
    subgraph.select(master, edges: range(6)),
    data: payload,
    hedge: payload,
  )
  let solved = layout(master, subgraph: all, labels: (steps: 0), solver: (steps: 2, seed: 17))
  for s in cycles + forests + (all,) {
    assert(
      type(s) == dictionary and s.linnest-kind == "linnest-subgraph" and s.topology == all.topology,
    )
    assert(
      graph.edges(solved, subgraph: s).map(e => e.name)
        == graph.edges(master, subgraph: s).map(e => e.name),
    )
    assert(
      graph.nodes(solved, subgraph: s).map(v => v.name)
        == graph.nodes(master, subgraph: s).map(v => v.name),
    )
    assert(subgraph.with-data(solved, s).bytes == s.bytes)
  }
  assert(forests.all(s => graph.edges(master, subgraph: s).len() == 3))
  let cycle = subgraph.with-data(master, cycles.first(), data: payload, hedge: payload)
  let laid-out = layout(solved, subgraph: cycle, labels: (steps: 0), solver: (steps: 2, seed: 17))
  draw(
    laid-out,
    subgraph: cycle,
    subgraph-edge-style: (stroke: black + 1pt),
    unit: 8pt,
    node-label: none,
  )
  draw(
    laid-out,
    subgraph: (forests.first(), (subgraph: all, edge-style: (stroke: black + 1pt))),
    subgraph-edge-style: (stroke: black + 1pt),
    unit: 8pt,
    node-label: none,
  )

  // Subsequent cuts patch only their new endpoints, not earlier boundary anchors.
  let repeated-master = graph.build({
    node(<a>)
    node(<b>)
    node(<c>)
    edge(<e>, source(<a>, payload: payload), sink(<b>, payload: payload), payload: payload)
    edge(<f>, source(<b>, payload: payload), sink(<c>, payload: payload), payload: payload)
  })
  let first-sides = (:)
  for side in ("left", "right") {
    first-sides.insert(side, subgraph.with-data(
      repeated-master,
      subgraph.select(
        repeated-master,
        source: if side == "left" { (<e>,) } else { () },
        sink: if side == "right" { (<e>,) } else { () },
      ),
      data: (pass: 1, side: side, payload: payload),
      hedge: (winding: 2, side: side, payload: payload),
    ))
  }
  let first = graph.cut(
    repeated-master,
    left: first-sides.left,
    right: first-sides.right,
    boundary: record => {
      let b = record.boundary
      assert(b.origin.name == <e>)
      (
        pos: pos(x: b.crossing + 1, y: if b.side == "left" { -1 } else { 1 }),
        data: (pass: 1, side: b.side, payload: payload),
      )
    },
  )
  let first-nodes = graph.nodes(first)
  let first-edges = graph.edges(first)
  let first-boundaries = graph.boundaries(first)
  assert(first-boundaries.len() == 4 and first-nodes.len() == 5 and first-edges.len() == 4)
  for (target, winding) in ((<f>, 2), (<e.1>, 1)) {
    let sides = (:)
    for side in ("left", "right") {
      sides.insert(side, subgraph.with-data(
        first,
        subgraph.select(
          first,
          source: if side == "left" { (target,) } else { () },
          sink: if side == "right" { (target,) } else { () },
        ),
        data: (pass: 2, side: side, payload: payload),
        hedge: (winding: winding, target: target, side: side, payload: payload),
      ))
    }
    let second = graph.cut(first, left: sides.left, right: sides.right, boundary: record => {
      let b = record.boundary
      assert(b.origin.name == target and b.cut-data.pass == 2)
      if b.node != none { assert(record.origin == none and target == <f>) }
      (
        pos: pos(x: b.crossing + 20, y: if b.side == "left" { -3 } else { 3 }),
        data: (pass: 2, target: target, side: b.side, payload: payload),
      )
    })
    let nodes = graph.nodes(second)
    let edges = graph.edges(second)
    let boundaries = graph.boundaries(second)
    let retained = boundaries.filter(b => b.cut-data.pass == 1)
    let created = boundaries.filter(b => b.cut-data.pass == 2)
    assert(nodes.len() == first-nodes.len() + 2 * (winding - 1))
    assert(edges.len() == first-edges.len() + winding)
    assert(
      boundaries.len() == 4 + 2 * winding and retained.len() == 4 and created.len() == 2 * winding,
    )
    for b in created {
      let record = if b.node == none { edges.at(b.edge) } else { nodes.at(b.node) }
      assert(record.data == (pass: 2, target: target, side: b.side, payload: payload))
      assert(b.data == subgraph.hedge-data(sides.at(b.side), b.hedge))
      assert(b.cut-data == sides.at(b.side).data)
      assert(b.pos == (x: b.crossing + 20, y: if b.side == "left" { -3 } else { 3 }))
      assert(record.boundary + (pos: b.pos) == b)
    }
    // Record origins are immediate-input IDs, including previously generated hedges.
    for e in edges {
      let old = first-edges.at(e.origin.edge)
      assert(e.origin.name == old.name)
      if old.name != target {
        assert(e.origin == (edge: old.edge, name: old.name, segment: none, winding: 0))
      }
      for h in (e.source, e.sink).filter(h => h != none) {
        let old-half = (old.source, old.sink).find(previous => (
          previous != none and previous.hedge == h.origin
        ))
        assert(old-half != none and h.data == old-half.data)
      }
    }
    // Boundary creation provenance stays unchanged; only current anchor IDs may move.
    for old in first-boundaries {
      let matches = retained.filter(b => (
        b.origin == old.origin
          and b.hedge == old.hedge
          and b.side == old.side
          and b.crossing == old.crossing
      ))
      assert(matches.len() == 1)
      let b = matches.first()
      assert(b == old + (edge: b.edge, node: b.node))
      assert(edges.at(b.edge).origin.edge == old.edge)
      let record = if b.node == none { edges.at(b.edge) } else { nodes.at(b.node) }
      let previous = if old.node == none { first-edges.at(old.edge) } else {
        first-nodes.at(old.node)
      }
      assert(record.pos == previous.pos and record.data == previous.data)
      assert(record.boundary + (pos: b.pos) == b)
      if old.node != none { assert(b.node != none and record.origin == old.node) }
    }
    if target == <e.1> {
      // Recutting the middle retains both old boundary nodes on distinct new stubs.
      let anchors = retained.filter(b => b.node != none)
      assert(anchors.len() == 2 and anchors.first().edge != anchors.last().edge)
      assert(created.all(b => b.node == none))
      assert(anchors.all(b => edges.at(b.edge).origin.name == <e.1> and b.origin.name == <e>))
    }
    assert(graph.nodes(first) == first-nodes and graph.edges(first) == first-edges)
    assert(graph.boundaries(first) == first-boundaries and graph.boundaries(repeated-master) == ())
  }
}
