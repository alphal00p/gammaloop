#set page(width: auto, height: auto, margin: 0pt)
#import "crates/linnest/typst/src/lib.typ": graph
#import "crates/linnest/typst/src/draw.typ" as draw
#import "crates/linnest/typst/src/curve.typ" as curve
#import "crates/linnest/typst/src/impl/draw.typ" as drawing
#import graph: node, edge, source, sink, pos

#let xy = p => if type(p) == array { p } else { (p.x, p.y) }
#let knots = path => {
  let segments = curve.segments(path)
  (xy(segments.first().start), ..segments.map(s => xy(s.end)))
}
#let base = graph.build({
  node(<a>, pos: pos(x: 0, y: 0))
  node(<b>, pos: pos(x: 6, y: 0))
  edge(<paired>, source(<a>, id: 3, route-points: ((1, 1), (2, 1))),
    sink(<b>, id: 1, route-points: ((5, -1),)), pos: pos(x: 3, y: 0))
  edge(<out>, source(<a>, route-points: ((-1, 1), (-2, 1))), pos: pos(x: -3, y: 0))
  edge(<in>, sink(<b>, route-points: ((7, 1), (8, 1))), pos: pos(x: 9, y: 0))
  edge(<loop>, source(<b>, route-points: ((5, -2),)),
    sink(<b>, route-points: ((7, -2),)), pos: pos(x: 6, y: -3))
})
#let nodes = graph.nodes(base)
#let edges = graph.edges(base)
#let paired = edges.find(e => e.name == <paired>)
#assert.eq(paired.source.hedge, 3)
#assert.eq(paired.sink.hedge, 1)
#assert.eq(paired.source.route-points.map(xy), ((1, 1), (2, 1)))
#assert.eq(paired.sink.route-points.map(xy), ((5, -1),))
#let points = ((0, 0), (1, 1), (2, 1), (3, 0), (5, -1), (6, 0))
#let halves = draw.edge-halves(paired, nodes)
#assert.eq(curve.segments(halves.curve), curve.segments(curve.hobby-spline(points)))
#assert.eq(knots(halves.source), points.slice(0, 4))
#assert.eq(knots(halves.sink), points.slice(3))
#let loop = draw.edge-halves(edges.find(e => e.name == <loop>), nodes)
#assert.eq(knots(loop.curve), ((6, 0), (5, -2), (6, -3), (7, -2), (6, 0)))

// Endpoint maps change geometry independently of native content and style data.
#let label = [route label]
#let mapped = graph.map(base,
  source: h => if h.hedge == 3 { (route-points: ((x: 1, y: 2),), label: label) },
  sink: h => if h.hedge == 1 { (route-points: ()) },
  edge: (paired: (spring-length: 1.4)),
)
#let changed = graph.edges(mapped).find(e => e.name == <paired>)
#assert.eq(changed.source.route-points.map(xy), ((1, 2),))
#assert.eq(changed.sink.route-points, ())
#assert.eq(changed.source.data.label, label)
#assert(not changed.source.data.keys().contains("route-points"))
#assert.eq(changed.statements.at("spring-length"), "1.4")
#assert.eq(graph.edges(graph.map(mapped, source: h => h)), graph.edges(mapped))
#let cleared = graph.map(base, source: _ => (route-points: ()), sink: _ => (route-points: ()))
#let clear-edge = graph.edges(cleared).find(e => e.name == <paired>)
#let clear-halves = draw.edge-halves(clear-edge, nodes)
#let old = curve.split-through(((0, 0), (3, 0), (6, 0)))
#assert.eq(clear-halves.curve, old.curve)
#assert.eq(clear-halves.source, old.parts.at(0))
#assert.eq(clear-halves.sink, old.parts.at(1))

// Ordinary drawing and its interaction/annotation carrier use the same route.
#let geometry = drawing._edge-geometry-halves(paired, nodes, (), (:), (:), (
  omega: 1.0, source-outset: 0, sink-outset: 0, accuracy: 0.001, label-pos: none,
))
#assert.eq(geometry.source + geometry.sink, curve.segments(halves.curve))
#let ignored = drawing._edge-geometry-halves(paired, nodes, (),
  (route-points: "ignore"), (route-points: "ignore"),
  (omega: 1.0, source-outset: 0, sink-outset: 0, accuracy: 0.001, label-pos: none),
)
#assert.eq(ignored.source + ignored.sink, curve.segments(old.curve))

// Solver-collapsed interior spans do not reach Hobby as repeated knots.
// The geometry remains split at the physical anchor, including a zero half.
#let repeated = paired + (
  source: paired.source + (route-points: ((0, 0), (1, 1), (1, 1), (2, 1), (3, 0))),
  sink: paired.sink + (route-points: ((6, 0), (5, -1), (3, 0))),
)
#let repeated-halves = draw.edge-halves(repeated, nodes)
#assert.eq(curve.segments(repeated-halves.curve), curve.segments(halves.curve))
#assert.eq(curve.segments(repeated-halves.source), curve.segments(halves.source))
#assert.eq(curve.segments(repeated-halves.sink), curve.segments(halves.sink))
#let trimmed-repeated = draw.edge-halves(repeated, nodes, source-outset: 0.2, sink-outset: 0.3)
#assert(calc.abs(curve.length(halves.source) - curve.length(trimmed-repeated.source) - 0.2) < 0.003)
#assert(calc.abs(curve.length(halves.sink) - curve.length(trimmed-repeated.sink) - 0.3) < 0.003)
#let collapsed = paired + (
  pos: (x: 0, y: 0),
  source: paired.source + (route-points: ((0, 0),)),
  sink: paired.sink + (route-points: ()),
)
#assert.eq(curve.segments(draw.edge-halves(collapsed, nodes).source), ())
#let collapsed = collapsed + (sink: collapsed.sink + (node: 0))
#assert.eq(curve.segments(draw.edge-halves(collapsed, nodes).curve), ())
#let empty-geometry = drawing._edge-geometry-halves(collapsed, nodes, (), (:), (:), (
  omega: 1.0, source-outset: 0, sink-outset: 0, accuracy: 0.001, label-pos: none,
))
#assert.eq(empty-geometry.source + empty-geometry.sink, ())

// Both dangling orientations read node-to-anchor route storage. Node outsets
// trim the smooth carrier by arc length; the external endpoint is unchanged.
#for (name, reverse) in ((<out>, false), (<in>, true)) {
  let e = edges.find(e => e.name == name)
  let h = if reverse { e.sink } else { e.source }
  let node = xy(nodes.at(h.node).pos)
  let end = xy(e.pos)
  let points = (node, ..h.route-points.map(xy), end)
  if reverse { points = points.rev() }
  let path = drawing._dangling-path(points.first(), points.last(), none, (:),
    dangling-at-start: reverse, route-points: h.route-points, node-point: node,
  )
  assert.eq(knots(path), points)
  assert.eq(curve.segments(path), curve.segments(curve.hobby-spline(points)))
  let repeated-path = drawing._dangling-path(points.first(), points.last(), none, (:),
    dangling-at-start: reverse, route-points: (node, ..h.route-points, end), node-point: node,
  )
  assert.eq(curve.segments(repeated-path), curve.segments(path))
  let trimmed = drawing._dangling-path(points.first(), points.last(), none, (:),
    dangling-at-start: reverse, route-points: h.route-points, node-point: node, node-outset: 0.2,
  )
  assert(calc.abs(curve.length(path) - curve.length(trimmed) - 0.2) < 0.003)
  let tangent = drawing._dangling-path(points.first(), points.last(), none,
    (dangling-tangent: "horizontal"), dangling-at-start: reverse,
    route-points: h.route-points, node-point: node,
  )
  assert.eq(knots(tangent), points)
  let endpoint = if reverse { curve.segments(tangent).first() } else { curve.segments(tangent).last() }
  assert.eq(xy(if reverse { endpoint.control-start } else { endpoint.control-end }).at(1), end.at(1))
}

#assert.eq(curve.segments(drawing._dangling-path((0, 0), (0, 0), none,
  (dangling-tangent: "horizontal"), route-points: ((0, 0),),
)), ())

// Exercise the public drawing path as well as geometry queries. All logical
// edges retain their source/sink ownership while painted across several cubics.
#draw.draw(base, unit: 12pt, node-radius: 0.08, node-label: none, edge-label: [p])
