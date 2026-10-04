#import "crates/linnest/typst/src/impl/draw.typ" as impl

// Exercise the same native static-cost geometry used by production search.
#let edge-intersection(candidate, lines) = {
  let geometry = (
    placements: ((candidates: ((bounds: candidate.bounds, cost: 0,
      corners: candidate.corners.map(impl._point)),),),),
    obstacles: (), label-padding: 0, obstacle-padding: 0, edge-lines: lines,
  )
  let math = (temperatures: range(84).map(_ => 1.0), exponentials: (), pair-padding: 0)
  cbor(impl._plugin.label_search(cbor.encode(geometry), cbor.encode(math))).costs.first()
}

// Exact text-edge contact has no attraction/repulsion outside the painted
// stroke. Rotation and curve subdivision must not change that rule.
#let candidate(x, cost: 0) = (
  bounds: (left: x, right: x + 1, bottom: 0, top: 1),
  corners: ((x, 1), (x + 1, 1), (x, 0), (x + 1, 0)),
  arrow-bounds: (), cost: cost, position: (x, 0),
)
#let line(a, b, radius: 0) = (start: a, end: b, radius: radius)
#let box = candidate(0)
#assert.eq(edge-intersection(box, (line((-1, 2), (2, 2)),)), 0)
#assert.eq(edge-intersection(box, (line((-1, 0.5), (2, 0.5)),)), 6)
#assert.eq(edge-intersection(box, (
  line((-1, 0.5), (0.5, 0.5)), line((0.5, 0.5), (2, 0.5)),
)), 6)
#assert(edge-intersection(box, (line((-1, 1.02), (2, 1.02), radius: 0.03),)) > 0)
#assert.eq(edge-intersection(box, (line((-1, 1.04), (2, 1.04), radius: 0.03),)), 0)
#let rotated = box + (
  corners: ((0, 1), (1, 2), (1, 0), (2, 1)),
  bounds: (left: 0, right: 2, bottom: 0, top: 2),
)
// This line is inside the AABB but outside the actual rotated text quad.
#assert.eq(edge-intersection(rotated, (line((0.05, 0.05), (0.2, 0.05)),)), 0)
#assert(edge-intersection(rotated, (line((0, 1), (2, 1)),)) > 0)

// Each isolated swap overlaps the other annotation; a coordinated move puts
// both at their preferred positions without an intermediate accepted overlap.
#let labels = (
  (candidates: (candidate(0, cost: 2), candidate(10))),
  (candidates: (candidate(10, cost: 2), candidate(0))),
)
#let relaxed = impl._relax-label-placements(labels, (), label-padding: 0)
#assert.eq(relaxed.map(c => c.position), ((10, 0), (0, 0)))
#assert.eq(impl._relax-label-placements(labels, (), label-padding: 0), relaxed)

// Independent finite segment/polygon distance used by attachment regressions.
// Points follow the perimeter of the measured text quad, in either direction.
#let attachment-distance(lines, points) = {
  let sub = (a, b) => (a.at(0) - b.at(0), a.at(1) - b.at(1))
  let dot = (a, b) => a.at(0) * b.at(0) + a.at(1) * b.at(1)
  let cross = (a, b) => a.at(0) * b.at(1) - a.at(1) * b.at(0)
  let distance = (point, a, b) => {
    let direction = sub(b, a)
    let length = dot(direction, direction)
    let at = if length < 1e-18 { 0 } else { calc.clamp(dot(sub(point, a), direction) / length, 0, 1) }
    let delta = sub(sub(point, a), direction.map(value => value * at))
    dot(delta, delta)
  }
  let edges = range(points.len()).map(i => (points.at(i), points.at(calc.rem(i + 1, points.len()))))
  let result = calc.inf
  for (a, b) in lines {
    let signs = edges.map(((c, d)) => cross(sub(d, c), sub(a, c)))
    if signs.all(sign => sign >= 0) or signs.all(sign => sign <= 0) { return 0 }
    for (c, d) in edges {
      let denominator = cross(sub(b, a), sub(d, c))
      if calc.abs(denominator) > 1e-18 {
        let t = cross(sub(c, a), sub(d, c)) / denominator
        let u = cross(sub(c, a), sub(b, a)) / denominator
        if 0 <= t and t <= 1 and 0 <= u and u <= 1 { return 0 }
      }
      result = calc.min(result, distance(a, c, d), distance(b, c, d), distance(c, a, b), distance(d, a, b))
    }
  }
  calc.sqrt(result)
}

// A wide horizontal label above a short diagonal arrow must clear the finite
// endpoint, not an infinitely extended tangent whose support corner is remote.
#let attachment-corners = ((-4, 0.35), (4, 0.35), (-4, -0.35), (4, -0.35))
#let diagonal = calc.sqrt(0.5)
#let attachment-lines = (((-0.7 * diagonal, -0.7 * diagonal), (0.7 * diagonal, 0.7 * diagonal)),)
#let attachment-normal = (-diagonal, diagonal)
#let attachment-initial = 0.2 + 4.35 * diagonal
#let attachment-offset = impl._label-attachment-offset(
  attachment-lines, attachment-corners, (0, 0), attachment-normal, 0.2, attachment-initial,
)
#assert(calc.abs(attachment-offset - (0.7 + 0.55 / diagonal)) < 1e-10)
#assert(attachment-offset < attachment-initial / 2)
#let translate = (points, offset, normal) => points.map(point => (
  point.at(0) + offset * normal.at(0), point.at(1) + offset * normal.at(1),
))
#let perimeter = points => (points.at(0), points.at(1), points.at(3), points.at(2))
#assert(calc.abs(attachment-distance(attachment-lines,
  perimeter(translate(attachment-corners, attachment-offset, attachment-normal)),
) - 0.2) < 1e-10)

// Rotation and canvas scaling preserve the same actual gap and local offset.
#for angle in (0deg, 30deg, -70deg) {
  for scale in (0.2, 1, 3) {
    let transform = point => (
      scale * (calc.cos(angle) * point.at(0) - calc.sin(angle) * point.at(1)),
      scale * (calc.sin(angle) * point.at(0) + calc.cos(angle) * point.at(1)),
    )
    let lines = attachment-lines.map(line => line.map(transform))
    let corners = attachment-corners.map(transform)
    let normal = transform(attachment-normal)
    let offset = impl._label-attachment-offset(lines, corners, (0, 0), normal, 0.2, attachment-initial)
    assert(calc.abs(offset - attachment-offset) < 1e-10)
    assert(calc.abs(attachment-distance(lines, perimeter(translate(corners, offset, normal))) - 0.2 * scale) < 1e-10)
  }
}

// A remote bend beyond a clear normal interval must not push the label away
// from its local attachment. Zero gap also reaches the boundary, not its center.
#let compact-corners = ((-0.4, 0.2), (0.4, 0.2), (-0.4, -0.2), (0.4, -0.2))
#let bent-lines = (((-1, 0), (1, 0)), ((1, 0), (2, 4)), ((2, 4), (-2, 4)))
#for gap in (0, 0.2) {
  let offset = impl._label-attachment-offset(bent-lines, compact-corners, (0, 0), (0, 1), gap, gap + 0.2)
  assert(calc.abs(offset - (gap + 0.2)) < 1e-10)
}

// The short arrow alone is insufficient when a wide label would cross the
// physical carrier beyond it. Clear the finite union, still at the same frame.
#let physical-line = (-4, 4).map(at => (
  at * diagonal - 0.5 * attachment-normal.at(0),
  at * diagonal - 0.5 * attachment-normal.at(1),
))
#let finite-union = attachment-lines + (physical-line,)
#let union-offset = impl._label-attachment-offset(
  finite-union, attachment-corners, (0, 0), attachment-normal, 0.2, attachment-initial,
)
#assert(calc.abs(union-offset - (4.35 * diagonal - 0.3)) < 1e-10)
#let union-corners = perimeter(translate(attachment-corners, union-offset, attachment-normal))
#assert(calc.abs(attachment-distance(finite-union, union-corners) - 0.2) < 1e-10)
#assert(attachment-distance(attachment-lines, union-corners) > 0.2)

// Empty or collapsed measured content remains finite, including the far cap.
#assert(calc.abs(impl._label-attachment-offset(
  (((-1, -1), (1, 1)),), ((-2, 0), (2, 0), (-2, 0), (2, 0)),
  (0, 0), (0, 1), 0.2, 0.2,
) - 1.2) < 1e-10)
#assert(calc.abs(impl._label-attachment-offset(
  (((-1, 0), (1, 0)),), ((0, 0), (0, 0), (0, 0), (0, 0)),
  (0, 0), (0, 1), 0.2, 0.2,
) - 0.2) < 1e-10)
