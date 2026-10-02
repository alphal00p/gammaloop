#import "crates/linnest/typst/src/curve.typ" as curve
#import "crates/linnest/typst/src/lib.typ": draw, graph
#import "crates/linnest/typst/src/impl/draw.typ" as drawing
#import "@preview/cetz:0.5.1" as cetz
#import "map-style.typ" as feynman
#import "annotation-placement-behavior.typ": attachment-distance, edge-intersection

// Candidate geometry is opaque in production. Decode it only for these
// per-candidate assertions, retaining the original packet for search tests.
#let inspect-layer-label-element(..args) = {
  let placement = drawing._layer-label-element(..args)
  if type(placement) != dictionary or not placement.keys().contains("candidates") {
    return placement
  }
  let packet = placement.candidates
  let footprints = if type(packet.footprints) == bytes { cbor(packet.footprints) }
    else { packet.footprints }
  let groups = ()
  let index = 0
  for batch in packet.batches {
    let candidates = cbor(batch.candidates)
    let records = ()
    for (at, candidate) in batch.positions.zip(candidates) {
      records.push((
        position: candidate.position, at: at, side: batch.side,
        path-index: batch.path-index, path-shift: candidate.path-shift,
        corners: candidate.corners, bounds: candidate.bounds, cost: candidate.cost,
        arrow-bounds: if footprints == none { () } else { footprints.at(index) },
      ))
      index += 1
    }
    groups.push(records)
  }
  let records = if packet.interleave {
    groups.first().zip(groups.last()).flatten()
  } else { groups.flatten() }
  placement + (candidates: records, inspected-candidates: records, packed-candidates: packet)
}

#let relax-inspected-labels(placements, ..args) = drawing._relax-label-placements(
  placements.map(label => if label.keys().contains("packed-candidates")
    and cbor.encode(label.candidates) == cbor.encode(label.inspected-candidates) {
    label + (candidates: label.packed-candidates)
  } else { label }),
  ..args,
)

// Scalar arithmetic retained only as an exact oracle for the native geometry boundary.
#let scalar-label-path-lines(ctx, path, accuracy) = {
  let lines = ()
  if path == none { return lines }
  let pending = curve.segments(path).rev().map(segment => (
    points: (segment.start, segment.control-start, segment.control-end, segment.end).map(
      point => drawing._point(cetz.matrix.mul4x4-vec3(ctx.transform, (..point, 0))),
    ),
    depth: 0,
  ))
  while pending.len() > 0 {
    let current = pending.pop()
    let (a, b, c, d) = current.points
    let chord = drawing._point-sub(d, a)
    let squared = chord.at(0) * chord.at(0) + chord.at(1) * chord.at(1)
    let distance = point => {
      let relative = drawing._point-sub(point, a)
      let at = if squared <= 1e-18 { 0 } else {
        calc.clamp((relative.at(0) * chord.at(0) + relative.at(1) * chord.at(1)) / squared, 0, 1)
      }
      drawing._point-distance(point, drawing._point-lerp(a, d, at))
    }
    let flatness = calc.max(distance(b), distance(c))
    if flatness <= accuracy or current.depth >= 10 {
      lines.push((a, d))
    } else {
      let ab = drawing._point-lerp(a, b, 0.5)
      let bc = drawing._point-lerp(b, c, 0.5)
      let cd = drawing._point-lerp(c, d, 0.5)
      let abc = drawing._point-lerp(ab, bc, 0.5)
      let bcd = drawing._point-lerp(bc, cd, 0.5)
      let mid = drawing._point-lerp(abc, bcd, 0.5)
      pending.push((points: (mid, bcd, cd, d), depth: current.depth + 1))
      pending.push((points: (a, ab, abc, mid), depth: current.depth + 1))
    }
  }
  lines
}

#let scalar-edge-collision-lines(ctx, elements) = {
  let painted = cetz.process.many(ctx, elements.flatten(), compute-bounds: false)
  let lines = ()
  let canvas = ctx + (transform: cetz.matrix.ident(4))
  for drawable in painted.drawables {
    if drawable.type != "path" or cetz.drawable.TAG.hidden in drawable.tags { continue }
    let stroke = cetz.util.resolve-stroke(drawable.at("stroke", default: none))
    if drawable.at("stroke", default: none) == none { continue }
    let radius = cetz.util.resolve-number(ctx, stroke.thickness) / 2
    for (origin, closed, commands) in drawable.segments {
      let start = origin
      for command in commands {
        let pieces = ()
        if command.first() == "l" {
          for end in command.slice(1) {
            pieces.push((drawing._point(start), drawing._point(end)))
            start = end
          }
        } else if command.first() == "c" {
          let segment = (start: drawing._point(start), control-start: drawing._point(command.at(1)),
            control-end: drawing._point(command.at(2)), end: drawing._point(command.at(3)))
          pieces = scalar-label-path-lines(canvas, drawing._segments-path((segment,)), 0.005)
          start = command.last()
        }
        lines += pieces.map(((a, b)) => (start: a, end: b, radius: radius))
      }
      if closed and start != origin {
        lines.push((start: drawing._point(start), end: drawing._point(origin), radius: radius))
      }
    }
  }
  lines
}

#let scalar_label-attachment-offset(lines, corners, frame, outward, gap, initial) = {
  if lines.len() == 0 { return initial }
  let dot = (a, b) => drawing._point-x(a) * drawing._point-x(b) + drawing._point-y(a) * drawing._point-y(b)
  let cross = (a, b) => drawing._point-x(a) * drawing._point-y(b) - drawing._point-y(a) * drawing._point-x(b)
  let direction = drawing._point(outward)
  let scale = drawing._point-length(direction)
  if scale <= 1e-9 { return initial }
  let normal = drawing._point-scale(direction, 1 / scale)
  let tangent = (-normal.at(1), normal.at(0))
  let project = point => (dot(point, tangent), dot(point, normal))
  let quad = (corners.at(0), corners.at(1), corners.at(3), corners.at(2)).map(project)
  let orientation = cross(drawing._point-sub(quad.at(1), quad.at(0)), drawing._point-sub(quad.at(2), quad.at(1)))
  let left = calc.min(..quad.map(point => point.at(0)))
  let right = calc.max(..quad.map(point => point.at(0)))
  let radius = gap * scale
  // Intersect the vertical normal ray with a finite capsule. Its rectangle and
  // endpoint discs give exact clearance to each straight boundary segment.
  let capsule = (a, b) => {
    if calc.min(a.at(0), b.at(0)) > radius or calc.max(a.at(0), b.at(0)) < -radius { return none }
    let low = calc.inf
    let high = -calc.inf
    for point in (a, b) {
      if calc.abs(point.at(0)) <= radius {
        let reach = calc.sqrt(calc.max(0, radius * radius - point.at(0) * point.at(0)))
        low = calc.min(low, point.at(1) - reach)
        high = calc.max(high, point.at(1) + reach)
      }
    }
    let delta = drawing._point-sub(b, a)
    let length = drawing._point-length(delta)
    if length > 1e-12 {
      let along = drawing._point-scale(delta, 1 / length)
      let across = (-along.at(1), along.at(0))
      let bottom = -calc.inf
      let top = calc.inf
      for (axis, minimum, maximum) in ((along, 0, length), (across, -radius, radius)) {
        let intercept = dot(a, axis)
        let slope = axis.at(1)
        if calc.abs(slope) <= 1e-12 {
          if -intercept < minimum or -intercept > maximum {
            bottom = calc.inf
            top = -calc.inf
            break
          }
        } else {
          let first = (minimum + intercept) / slope
          let last = (maximum + intercept) / slope
          bottom = calc.max(bottom, calc.min(first, last))
          top = calc.min(top, calc.max(first, last))
        }
      }
      if bottom <= top {
        low = calc.min(low, bottom)
        high = calc.max(high, top)
      }
    }
    if low <= high { (low, high) } else { none }
  }
  let intervals = ()
  for (start, end) in lines {
    let start = project(drawing._point-sub(start, frame))
    let end = project(drawing._point-sub(end, frame))
    if calc.min(start.at(0), end.at(0)) > right + radius or calc.max(start.at(0), end.at(0)) < left - radius { continue }
    // Sweep the measured quad along this finite segment. Each quad edge is
    // exposed at one endpoint; side transitions add the corner trajectories.
    let first = quad.map(corner => drawing._point-sub(start, corner))
    let last = quad.map(corner => drawing._point-sub(end, corner))
    let delta = drawing._point-sub(end, start)
    let at-end = range(4).map(i => cross(
      drawing._point-sub(first.at(calc.rem(i + 1, 4)), first.at(i)), delta,
    ) * orientation < 0)
    let boundary = ()
    for i in range(4) {
      let next = calc.rem(i + 1, 4)
      let points = if at-end.at(i) { last } else { first }
      boundary.push((points.at(i), points.at(next)))
      if calc.abs(orientation) <= 1e-18 { boundary.push((last.at(i), last.at(next))) }
      if at-end.at(i) != at-end.at(calc.rem(i + 3, 4)) or calc.abs(orientation) <= 1e-18 {
        boundary.push((first.at(i), last.at(i)))
      }
    }
    let low = calc.inf
    let high = -calc.inf
    for (a, b) in boundary {
      let interval = capsule(a, b)
      if interval != none {
        low = calc.min(low, interval.at(0))
        high = calc.max(high, interval.at(1))
      }
    }
    if low <= high { intervals.push((low, high)) }
  }
  // Curves can give disjoint collision intervals. Leave the component touching
  // the attachment frame, rather than jumping over a clear gap to a remote bend.
  let reach = 0
  let attached = false
  for (low, high) in intervals.sorted(key: interval => interval.at(0)) {
    if high < 0 { continue }
    if low > reach + 1e-10 { break }
    attached = true
    reach = calc.max(reach, high)
  }
  // Explicit independent arrow/label shifts may separate the whole finite
  // arrow from this normal ray; retain that deliberate tangential placement.
  if attached { reach / scale } else { initial }
}

#let scalar_label-edge-intersection(candidate, lines) = {
  if lines.len() == 0 { return 0 }
  let corners = candidate.corners
  let origin = drawing._point(corners.at(0))
  let x = drawing._point-sub(corners.at(1), origin)
  let y = drawing._point-sub(corners.at(2), origin)
  let determinant = x.at(0) * y.at(1) - x.at(1) * y.at(0)
  if calc.abs(determinant) < 1e-12 { return 0 }
  let width = drawing._point-length(x)
  let height = drawing._point-length(y)
  let local = point => {
    let delta = drawing._point-sub(point, origin)
    ((delta.at(0) * y.at(1) - delta.at(1) * y.at(0)) / determinant,
      (x.at(0) * delta.at(1) - x.at(1) * delta.at(0)) / determinant)
  }
  let length = 0
  for line in lines {
    let box = candidate.bounds
    if (calc.max(line.start.at(0), line.end.at(0)) + line.radius < box.left
      or calc.min(line.start.at(0), line.end.at(0)) - line.radius > box.right
      or calc.max(line.start.at(1), line.end.at(1)) + line.radius < box.bottom
      or calc.min(line.start.at(1), line.end.at(1)) - line.radius > box.top) { continue }
    let a = local(line.start)
    let b = local(line.end)
    let low = 0
    let high = 1
    for axis in range(2) {
      // The inverse transform's row norm also handles sheared canvases.
      let padding = line.radius * (if axis == 0 { height } else { width }) / calc.abs(determinant)
      let delta = b.at(axis) - a.at(axis)
      if calc.abs(delta) < 1e-12 {
        if a.at(axis) < -padding or a.at(axis) > 1 + padding { high = -1; break }
      } else {
        let from = (-padding - a.at(axis)) / delta
        let to = (1 + padding - a.at(axis)) / delta
        low = calc.max(low, calc.min(from, to))
        high = calc.min(high, calc.max(from, to))
      }
    }
    length += calc.max(0, high - low) * drawing._point-distance(line.start, line.end)
  }
  6 * length / calc.max(1e-9, calc.min(width, height))
}

// Fixed-sample CeTZ measurements retain scalar arithmetic, sample counts and
// integer/float endpoint types; these are distinct from adaptive Kurvst lengths.
#let scalar-cubic-arclen(s, e, c1, c2, samples: 20) = {
  let d = 0
  for i in range(1, samples + 1) {
    let t0 = (i - 1) / samples
    let t1 = i / samples
    d += cetz.vector.dist(
      cetz.path-util.bezier.cubic-point(s, e, c1, c2, t0),
      cetz.path-util.bezier.cubic-point(s, e, c1, c2, t1))
  }
  return d
}

#let scalar-cubic-t-for-distance(s, e, c1, c2, d, samples: 20) = {
  let travel-forwards(s, e, c1, c2, d) = {
    let sum = 0
    for n in range(1, samples + 1) {
      let t0 = (n - 1) / samples
      let t1 = n / samples

      let segment-dist = cetz.vector.dist(cetz.path-util.bezier.cubic-point(s, e, c1, c2, t0),
                                     cetz.path-util.bezier.cubic-point(s, e, c1, c2, t1))
      if sum <= d and d <= sum + segment-dist {
        let lambda = (d - sum) / segment-dist
        return (1 - lambda) * t0 + lambda * t1
      }
      sum += segment-dist
    }
    return 1
  }

  if d == 0 {
    return 0
  }

  if d > 0 {
    return travel-forwards(s, e, c1, c2, d)
  } else {
    return 1 - travel-forwards(e, s, c2, c1, -d)
  }
}

#let curves = (
  ((0, 0), (2, 0), (0, 0), (2, 0)),
  ((-0.0, 0.0), (0.0, -0.0), (-0.0, -0.0), (0.0, 0.0)),
  ((0, 0), (1, 1), (4, -2), (-3, 3)),
  ((0, 0, 0), (1, 1, 3), (4, -2, 7), (-3, 3, 2)),
  ((0, 0), (1, 1, 3), (4, -2), (-3, 3, 2)),
  ((0, 0, 0, 1), (1, 1, 3, 2), (4, -2, 7, 9), (-3, 3, 2, 0)),
  ((1, 2), (1, 2), (1, 2), (1, 2)),
  ((0.000001, -0.000002), (0.000003, 0.000004), (-0.000005, 0.000007), (0.000011, -0.000013)),
)
#for (s, e, c1, c2) in curves {
  for samples in (0, 1, 2, 7, 20, 25, 37) {
    let expected = scalar-cubic-arclen(s, e, c1, c2, samples: samples)
    let actual = cetz.path-util.bezier.cubic-arclen(s, e, c1, c2, samples: samples)
    assert.eq(cbor.encode(actual), cbor.encode(expected),
      message: "sampled length: " + repr((s, e, c1, c2, samples, actual, expected)))
    for d in (0, -0.0, -100, 100, 0.05, 0.1, 0.25, 0.5, 1, 2, 10, -0.05, -0.1, -0.25, -0.5, -1, -2, -10) {
      let expected = scalar-cubic-t-for-distance(s, e, c1, c2, d, samples: samples)
      let actual = cetz.path-util.bezier.cubic-t-for-distance(s, e, c1, c2, d, samples: samples)
      assert.eq(cbor.encode(actual), cbor.encode(expected),
        message: "sampled parameter: " + repr((s, e, c1, c2, samples, d, actual, expected)))
    }
  }
}

// Include this file, or import and render this value: canvas assertions need layout.
#let curved-arrow-behavior = {
  let accuracy = 1e-6
  let epsilon = 1e-8
  let source = curve.cubic((0, 0), (0, 2), (1, 3), (3, 3))
  let sink = curve.cubic((3, 3), (5, 3), (6, 2), (6, 0))
  let path = curve.path(source, sink)
  let geometry-checks = ()

  // Either single move raises the collision cost, but moving both lowers it.
  // The host-certified annealing search must cross this barrier and retain
  // the original candidate records, including their integer cost types.
  let candidate(left, cost) = (
    bounds: (left: left, right: left + 1, bottom: 0, top: 1), cost: cost,
  )
  let labels = (
    (candidates: (candidate(0, 0.001), candidate(11.3495, 0))),
    (candidates: (candidate(10, 0.001), candidate(1.3495, 0))),
  )
  assert.eq(
    cbor.encode(relax-inspected-labels(labels, (), label-padding: 0)),
    cbor.encode(labels.map(label => label.candidates.last())),
    message: "annealing crosses an energy barrier without changing candidate records",
  )

  // The batched sampler must retain the exact arithmetic of scalar prefix
  // trimming, including endpoint clamps, degenerate curves and mixed segments.
  for sampled in (
    path, curve.line((0, 0), (4, 2)), curve.line((0, 0), (0, 0)),
    curve.line((0, 0), (1e-8, 0)),
    curve.path(curve.move-to((1, 2)), curve.quad-to((3, 7), (5, -1)), curve.line-to((7, 3))),
  ) {
    let total = curve.length(sampled, accuracy: accuracy)
    let shifts = (-total, -total / 2, -total / 3, 0, total / 3, total / 2, total)
    let actual = curve.frames(sampled, shifts.map(shift => total / 2 + shift), accuracy: accuracy)
    let expected = shifts.map(shift => drawing._path-mid-frame(sampled, accuracy, shift: shift))
    assert.eq(actual, expected, message: "batched frames preserve scalar sampling")
  }

  // Native obstacle sampling retains the scalar Kurvst/CeTZ arithmetic and
  // ordering, including empty paths, mixed curves and singular transforms.
  for sampled in (
    none, path, curve.line((0, 0), (4, 2)), curve.line((1, -2), (1, -2)),
    curve.path(curve.move-to((1, 2)), curve.quad-to((3, 7), (5, -1)), curve.line-to((7, 3))),
  ) {
    let segments = if sampled == none { () } else { curve.segments(sampled) }
    for transform in (
      ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)),
      ((2, 0.4, 0, 10), (-0.3, 1.5, 0, -4), (0.2, 0.7, 1, 3), (0, 0, 0, 1)),
      ((0, 0, 0, 0), (0, -2, 0, 7), (0, 0, 1, 0), (0, 0, 0, 1)),
    ) {
      for (pad-x, pad-y) in ((0, 0), (0.03, 0.07)) {
        let expected = ()
        for segment in segments {
          let points = range(13).map(step => cetz.matrix.mul4x4-vec3(
            transform, (..curve.cubic-point(segment, step / 12), 0),
          ))
          for (start, end) in points.slice(0, 12).zip(points.slice(1)) {
            expected.push((
              left: calc.min(start.at(0), end.at(0)) - pad-x,
              right: calc.max(start.at(0), end.at(0)) + pad-x,
              bottom: calc.min(start.at(1), end.at(1)) - pad-y,
              top: calc.max(start.at(1), end.at(1)) + pad-y,
            ))
          }
        }
        let actual = cbor(drawing._plugin.label_obstacle_boxes(cbor.encode((
          segments: segments, transform: transform.map(row => row.map(float)),
          pad-x: pad-x, pad-y: pad-y,
        ))))
        assert.eq(cbor.encode(actual), cbor.encode(expected), message: "obstacle boxes preserve scalar sampling exactly")
      }
    }
  }

  // Length, ratio and shift are measured on the offset path, not its centerline.
  for offset in (-0.35, 0.35) {
    let parallel = curve.parallel(path, distance: offset, accuracy: accuracy)
    let total = curve.length(parallel, accuracy: accuracy)
    for ratio in (none, 0.6) {
      for shift in (-0.2, 0.2) {
        let style = (
          offset: offset,
          length: if ratio == none { 2 } else { none },
          ratio: ratio,
          shift: shift,
          accuracy: accuracy,
        )
        let target = if ratio == none { 2 } else { ratio * total }
        let margin = (total - 0.1 - 0.2 - target) / 2
        let expected = curve.trim(
          parallel,
          start-outset: 0.1 + margin + shift,
          end-outset: 0.2 + margin - shift,
          accuracy: accuracy,
        )
        let case = repr((offset: offset, ratio: ratio, shift: shift))
        geometry-checks.push((
          "offset layer " + case,
          drawing._path-layer(path, style, 0.1, 0.2, none, auto),
          expected,
          target,
        ))
        for mode in ("whole", "gap", "opposite") {
          let gap = if mode == "gap" { 0.2 } else { 0 }
          let source-style = style + (split-gap: gap)
          let sink-style = (
            source-style
              + (
                offset: if mode == "opposite" { -offset } else { offset },
              )
          )
          let paired = drawing._split-edge-geometry(
            source,
            sink,
            path,
            source-style,
            sink-style,
            0.1,
            0.2,
            none,
            accuracy,
          )
          assert.eq(paired.split-gap, gap, message: "paired gap " + case)
          if mode == "whole" {
            assert.ne(paired.whole, none, message: "missing whole path " + case)
            geometry-checks.push((
              "paired whole " + case,
              curve.path(..paired.whole.map(curve.from-cubic)),
              expected,
              target,
            ))
          } else {
            assert.eq(
              paired.whole,
              none,
              message: "unexpected whole path " + case,
            )
          }
          for (index, half-style) in (source-style, sink-style).enumerate() {
            let halves = (source, sink).map(part => curve.parallel(
              part,
              distance: half-style.offset,
              accuracy: accuracy,
            ))
            let lengths = halves.map(part => curve.length(
              part,
              accuracy: accuracy,
            ))
            let target = if ratio == none { 2 } else { ratio * lengths.sum() }
            let margin = (lengths.sum() - 0.1 - 0.2 - target) / 2
            let origin = if index == 0 { 0 } else { lengths.first() }
            let start = calc.max(
              if index == 1 { gap / 2 } else { 0 },
              0.1 + margin + shift - origin,
            )
            let end = calc.min(
              lengths.at(index) - if index == 0 { gap / 2 } else { 0 },
              lengths.sum() - 0.2 - margin + shift - origin,
            )
            assert(
              end > start,
              message: "fixture must retain both paired halves",
            )
            geometry-checks.push((
              "paired " + repr((mode, index)) + " " + case,
              curve.path(
                ..(paired.source, paired.sink).at(index).map(curve.from-cubic),
              ),
              curve.trim(
                halves.at(index),
                start-outset: start,
                end-outset: lengths.at(index) - end,
                accuracy: accuracy,
              ),
              end - start,
            ))
          }
        }
      }
    }
  }
  for (case, actual, expected, target) in geometry-checks {
    assert(
      calc.abs(curve.length(actual, accuracy: accuracy) - target)
        < 8 * accuracy,
      message: case + ": visible arc length",
    )
    let actual = curve.segments(actual)
    let expected = curve.segments(expected)
    for (a, b) in (
      (actual.first().start, expected.first().start),
      (actual.last().end, expected.last().end),
    ) {
      assert(
        cetz.vector.dist(a, b) < 8 * accuracy,
        message: case + ": shifted endpoint",
      )
    }
  }

  cetz.canvas({
    cetz.draw.get-ctx(ctx => {
      let segment = (
        start: (0, 0),
        control-start: (0, 2),
        control-end: (2, -2),
        end: (2, 0),
      )
      let curved = curve.from-cubic(segment)
      let straight = curve.from-cubic(curve.line-segment((0, 0), (2, 0)))
      let head = (
        anchor: "center",
        fill: black,
        stroke: black + 0.8pt,
        length: 0.6,
        width: 0.3,
        inset: 0,
        shorten-to: auto,
      )
      let identity = cetz.matrix.ident(4)
      // Native flattening retains transformed coordinates, segment boundaries,
      // zero chords, the depth cap and the old left-to-right subdivision order.
      let flatten-cases = (
        curve.path(),
        curve.path(curve.move-to((1, 2))),
        curve.line((-0.0, 0.0), (2, -1)),
        curve.cubic((0, 0), (0, 2), (3, -1), (4, 0)),
        curve.cubic((1, 2), (1, 2), (1, 2), (1, 2)),
        curve.cubic((0, 0), (2, 3), (-1, 2), (0, 0)),
        curve.path(curve.move-to((0, 0)), curve.quad-to((2, 4), (3, 1)), curve.line-to((5, 2))),
      )
      for transform in (
        identity,
        cetz.matrix.transform-translate(3, -2, 0),
        cetz.matrix.transform-scale((-1.5, 0.6, 1)),
        cetz.matrix.transform-shear-x(0.6),
        cetz.matrix.mul-mat(cetz.matrix.transform-translate(3, 2, 0), cetz.matrix.transform-rotate-z(37deg)),
        cetz.matrix.transform-scale((0, 0, 1)),
      ) {
        let local = ctx + (transform: transform)
        for path in flatten-cases {
          assert.eq(
            cbor.encode(drawing._label-path-lines(local, path, 0.005)),
            cbor.encode(scalar-label-path-lines(local, path, 0.005)),
            message: "native finite-path subdivision retains exact transformed geometry",
          )
        }
      }
      let limit-path = curve.cubic((0, 0), (0, 1), (1, 1), (1, 0))
      let capped = drawing._label-path-lines(ctx, limit-path, 0)
      assert.eq(capped.len(), 1024)
      assert.eq(cbor.encode(capped), cbor.encode(scalar-label-path-lines(ctx, limit-path, 0)))

      // Mixed commands retain straight endpoint number types and close every
      // subpath independently; hidden and unpainted paths stay absent.
      let painted-paths = (
        ((-0.0, 0.0, 0.0), true, (("l", (1, 0), (2.0, 1.0)), ("c", (3, 4), (4, -2), (5, 0)))),
        ((7, 8), false, (("l", (9, 8), (9, 10)),)),
        ((2, 3), true, ()),
      )
      let painted = (local => (ctx: local, drawables: (
        (type: "path", segments: painted-paths, stroke: black + 0.4pt, fill: none, tags: ()),
        (type: "path", segments: painted-paths, stroke: black + 2pt, fill: none, tags: (cetz.drawable.TAG.hidden,)),
        (type: "path", segments: painted-paths, stroke: none, fill: none, tags: ()),
      )),)
      assert.eq(
        cbor.encode(drawing._edge-collision-lines(ctx, painted)),
        cbor.encode(scalar-edge-collision-lines(ctx, painted)),
        message: "batched painted cubics retain line order, closure and radii",
      )

      // Compare numeric results, including signed zeros and scalar early-return
      // types, under reflected, rotated, sheared and degenerate geometry.
      for transform in (
        identity, cetz.matrix.transform-scale((-1.5, 0.6, 1)),
        cetz.matrix.transform-shear-x(0.6), cetz.matrix.transform-rotate-z(37deg),
        cetz.matrix.transform-scale((0, 1, 1)),
      ) {
        let transform-point = point => drawing._point(cetz.matrix.mul4x4-vec3(transform, (..point, 0)))
        let corners = ((-1.0, 0.5), (1.0, 0.5), (-1.0, -0.5), (1.0, -0.5)).map(transform-point)
        let outward = transform-point((0.0, 1.0))
        let box = (
          left: calc.min(..corners.map(p => p.at(0))), right: calc.max(..corners.map(p => p.at(0))),
          bottom: calc.min(..corners.map(p => p.at(1))), top: calc.max(..corners.map(p => p.at(1))),
        )
        let candidate = (bounds: box, corners: corners, cost: 0)
        for points in (
          (), (((-4.0, 0.0), (4.0, 0.0)),),
          (((-4.0, 3.0), (4.0, 3.0)),),
          (((-4.0, 0.0), (4.0, 0.0)), ((-4.0, 3.0), (4.0, 3.0))),
          (((-2.0, -0.4), (3.0, 0.3)), ((3.0, 0.3), (0.0, 2.0))),
          (((0.0, -0.0), (0.0, -0.0)),),
          (((2.0, -4.0), (2.0, 4.0)),),
        ) {
          let lines = points.map(pair => pair.map(transform-point))
          for gap in (0, 0.2, 0.7) {
            for initial in (4, 4.0) {
              assert.eq(
                cbor.encode(drawing._label-attachment-offset(lines, corners, (0, 0), outward, gap, initial)),
                cbor.encode(scalar_label-attachment-offset(lines, corners, (0, 0), outward, gap, initial)),
                message: "native finite attachment retains scalar geometry and number types",
              )
            }
            let strokes = lines.map(pair => (start: pair.first(), end: pair.last(), radius: gap))
            assert.eq(
              cbor.encode(edge-intersection(candidate, strokes)),
              cbor.encode(scalar_label-edge-intersection(candidate, strokes)),
              message: "native finite text-edge intersections retain exact scalar penalties",
            )
          }
        }
      }

      // Final candidate geometry must match the former scalar attachment path,
      // including explicit anchors, fixed labels and dangling endpoint normals.
      let scalar-candidate-final(spec, index) = {
        let frame = spec.frames.at(index)
        let endpoint = spec.endpoint-direction
        let normal = if endpoint != none {
          drawing._point-scale(endpoint, 1 / drawing._point-length(endpoint))
        } else {
          let length = drawing._point-length(frame.tangent)
          if length <= 1e-9 { (0.0, spec.side) } else {
            drawing._point-scale((-frame.tangent.at(1), frame.tangent.at(0)), spec.side / length)
          }
        }
        let outward = cetz.vector.sub(cetz.matrix.mul4x4-vec3(spec.transform, (..normal, 0)), spec.origin)
        let squared = cetz.vector.dot(outward, outward)
        let nearest = if (spec.clear-box or endpoint != none) and squared > 1e-18 {
          calc.min(..spec.corners.map(corner => cetz.vector.dot(corner, outward))) / squared
        } else { 0.0 }
        let offset = spec.gap - nearest
        if spec.attachment != none and spec.clear-box and not spec.fixed {
          let local = ctx + (transform: spec.transform)
          let attachment = spec.attachment
          let lines = (scalar-label-path-lines(local, drawing._segments-path(attachment.paths.at(index)), spec.accuracy)
            + scalar-label-path-lines(local, drawing._segments-path(attachment.carrier), spec.accuracy))
          offset = scalar_label-attachment-offset(lines, spec.corners,
            drawing._point(cetz.matrix.mul4x4-vec3(spec.transform, (..frame.point, 0))),
            outward, spec.gap, offset)
        }
        let position = if spec.fixed and endpoint == none { frame.point } else {
          drawing._point-add(frame.point, drawing._point-scale(normal, offset))
        }
        let center = cetz.matrix.mul4x4-vec3(spec.transform, (..position, 0))
        let corners = spec.corners.map(corner => cetz.vector.add(center, corner))
        (
          normal: normal, outward: outward, nearest: nearest,
          position: position, corners: corners,
          bounds: (
            left: calc.min(..corners.map(p => p.at(0))), right: calc.max(..corners.map(p => p.at(0))),
            bottom: calc.min(..corners.map(p => p.at(1))), top: calc.max(..corners.map(p => p.at(1))),
          ),
        )
      }
      let candidate-frames = (
        (point: (0.0, 0.0), tangent: (1.0, 0.0)),
        (point: (1.2, 0.4), tangent: (0.7, -0.3)),
        (point: (-0.2, -0.1), tangent: (0.0, 0.0)),
      )
      let candidate-arrows = (
        curve.line((-0.4, 0), (0.7, 0)),
        curve.cubic((0.6, 0.2), (0.6, 1.5), (1.8, -1), (1.8, 0.6)),
        curve.line((5, 0), (6, 0)),
      ).map(curve.segments)
      let candidate-carrier = curve.segments(curve.path(
        curve.line((-4, 0), (4, 0)), curve.line((-4, 3), (4, 3)),
      ))
      for (transform-index, transform) in (
        identity,
        cetz.matrix.transform-translate(3, -2, 1),
        cetz.matrix.transform-scale((-1.5, 0.6, 1)),
        cetz.matrix.mul-mat(cetz.matrix.transform-rotate-z(37deg), cetz.matrix.transform-shear-x(0.6)),
        cetz.matrix.transform-scale((0, 1, 1)),
      ).enumerate() {
        let transform = transform.map(row => row.map(float))
        let origin = cetz.matrix.mul4x4-vec3(transform, (0.0, 0.0, 0.0)).map(float)
        let corners = ((-1.0, 0.5), (1.0, 0.5), (-1.0, -0.5), (1.0, -0.5)).map(point =>
          cetz.vector.sub(cetz.matrix.mul4x4-vec3(transform, (..point, 0.0)), origin).map(float))
        for (mode-index, (fixed, clear-box, endpoint)) in (
          (false, true, none), (true, true, none),
          (false, false, none), (true, false, none),
          (false, true, (2.0, -3.0)), (true, true, (2.0, -3.0)),
          (false, false, (-1.0, 0.5)), (true, false, (-1.0, 0.5)),
        ).enumerate() {
          for side in (-1.0, 1.0) {
            for (attachment-index, attachment) in (
              none,
              (carrier: (), paths: ((), (), ())),
              (carrier: candidate-carrier, paths: ((), (), ())),
              (carrier: (), paths: candidate-arrows),
              (carrier: candidate-carrier, paths: candidate-arrows),
            ).enumerate() {
              let spec = (
                frames: candidate-frames, positions: (0.4, 1.2, 2.6),
                transform: transform, origin: origin, corners: corners,
                total: 3.0, preferred: 1.5, accuracy: 0.005,
                side: side, preferred-side: 1.0, path-index: 1, path-shift: -0.25,
                gap: 0.2, clear-box: clear-box, fixed: fixed,
                endpoint-direction: endpoint, attachment: attachment,
              )
              let actual = cbor(drawing._plugin.label_candidates(cbor.encode(spec)))
              assert.eq(actual.len(), candidate-frames.len())
              for (index, candidate) in actual.enumerate() {
                let expected = scalar-candidate-final(spec, index)
                for (key, value) in expected {
                  assert.eq(candidate.at(key), value,
                    message: "batched finite candidate " + repr((transform-index, mode-index, side, attachment-index, index, key)))
                }
              }
            }
          }
        }
      }

      // Batching path drawing preserves CeTZ primitives, including disconnected
      // or empty subpaths and transforms that collapse distinct endpoints.
      for (sampled, unit) in (
        (path, 1), (straight, 1), (curved, 1),
        (curve.path(curve.move-to((0, 0)), curve.line-to((1, 1)), curve.close()), 1),
        (curve.path(curve.move-to((0, 0)), curve.move-to((1, 1)), curve.line-to((2, 3))), 1),
        (curve.path(curve.move-to((0, 0)), curve.line-to((1, 1)), curve.move-to((3, 2)), curve.line-to((4, 1))), 1),
        (curve.path(curve.move-to((-0.0, 0.0)), curve.quad-to((1, 3), (2, 0))), 1.75),
        (curve.path(curve.move-to((0, 0)), curve.line-to((1, 1)), curve.line-to((1, 1)), curve.close()), -0.25),
        (curved, 1pt),
        (curve.path(curve.move-to((0, 0)), curve.cubic-to((1, 2), (2, -1), (3, 0)), curve.close(),
          curve.move-to((3, 0)), curve.cubic-to((4, 2), (5, -2), (6, 0))), 1),
        (curve.path(curve.move-to((0, 0)), curve.line-to((0, 0)),
          curve.move-to((1, 2)), curve.line-to((1, 2))), 1),
      ) {
        let primitives = {
          for (origin, closed, segments) in curve.to-cetz-data(sampled, unit: unit) {
            let current = origin
            for (kind, ..args) in segments {
              if kind == "l" { cetz.draw.line(current, args.last()) }
              else { cetz.draw.bezier(current, args.last(), args.at(0), args.at(1)) }
              current = args.last()
            }
            if closed and current != origin { cetz.draw.line(current, origin) }
          }
        }
        for transform in (
          identity, cetz.matrix.transform-scale((0, 0, 1)),
          cetz.matrix.transform-shear-x(0.6),
          cetz.matrix.transform-rotate-xyz(40deg, 30deg, 0deg),
        ) {
          let local = ctx + (transform: transform)
          let actual = curve.to-cetz(sampled, unit: unit).first()(local)
          let expected = cetz.draw.merge-path(primitives).first()(local)
          assert.eq(actual.drawables, expected.drawables, message: "batched path geometry")
          assert.eq(cbor.encode(actual.drawables.map(d => d.segments)),
            cbor.encode(expected.drawables.map(d => d.segments)), message: "batched path coordinates")
          assert.eq(actual.ctx.prev, expected.ctx.prev, message: "batched path context")
        }
      }
      // CeTZ's close flag affects later primitive joins: once closed,
      // subpath-end is its origin, rather than its final command endpoint.
      let closing-subpaths = curve.path(
        curve.move-to((0, 0)), curve.line-to((1, 0)), curve.line-to((1, 1)),
        curve.move-to((0, 0)), curve.line-to((3, 0)),
      )
      let closing-primitives = {
        cetz.draw.line((0.0, 0.0), (1.0, 0.0))
        cetz.draw.line((1.0, 0.0), (1.0, 1.0))
        cetz.draw.line((0.0, 0.0), (3.0, 0.0))
      }
      for transform in (
        identity, cetz.matrix.transform-scale((-1.5, 0.6, 1)),
        cetz.matrix.transform-scale((0, 0, 1)),
      ) {
        let local = ctx + (transform: transform)
        let actual = curve.to-cetz(closing-subpaths, close: true).first()(local)
        let expected = cetz.draw.merge-path(closing-primitives, close: true).first()(local)
        assert.eq(cbor.encode(actual.drawables.map(d => d.segments)),
          cbor.encode(expected.drawables.map(d => d.segments)),
          message: "closed multi-subpaths retain CeTZ primitive connector semantics")
        assert.eq(actual.ctx.prev, expected.ctx.prev)
      }
      let cases = ()
      let comparisons = ()

      // Native physical sizes, reversal and named anchors survive transforms;
      // boundary heads now move inward instead of extending beyond the carrier.
      for transform in (
        identity,
        cetz.matrix.mul-mat(
          cetz.matrix.transform-translate(3, 2, 0),
          cetz.matrix.transform-rotate-z(37deg),
        ),
        cetz.matrix.transform-scale((-1.5, 0.6, 1)),
        cetz.matrix.transform-shear-x(0.6),
        cetz.matrix.transform-rotate-xyz(40deg, 30deg, 0deg),
      ) {
        for shape in (false, true) {
          for path in (straight, curved) {
            cases.push((path: path, transform: transform, shape: shape))
          }
        }
      }
      for entry in (
        (scale: 1.5),
        (length: 700%, width: 400%, stroke: black + 0.5pt),
      ) {
        cases.push((path: curved, ratio: 0.5, head: entry))
      }
      // Interior heads overlay the entire shaft, even at repeated controls.
      // Exercise the full path, with a head spanning a cubic join, not a terminal half.
      for path in (
        curve.from-cubic(segment + (control-start: segment.start)),
        curve.from-cubic(segment + (control-end: segment.end)),
        curve.path(curved, curve.from-cubic(curve.line-segment(
          (2, 0),
          (2, 0.02),
        ))),
        curve.path(
          curve.cubic((0, 0), (0, 0.5), (0.5, 1), (1, 1)),
          curve.cubic((1, 1), (1.5, 1), (2, 0.5), (2, 0)),
        ),
      ) {
        for ratio in (none, 0, 0.2, 0.5, 0.8, 1) {
          for shift in (-0.1, 0.1) {
            cases.push((path: path, ratio: ratio, shift: shift))
          }
        }
      }
      // A fully consumed shaft still has a compressed mark; only an empty path has none.
      for length in (0, 1e-9, 1e-6, 0.01, 0.1) {
        for ratio in (none, 0.5) {
          cases.push((
            path: curve.from-cubic(curve.line-segment((0, 0), (length, 0))),
            ratio: ratio,
          ))
        }
      }

      // End heads and centered heads have two actual geometric contacts, not tangent alignment.
      for (index, config) in cases.enumerate() {
        let path = config.path
        let head = head + config.at("head", default: (:))
        let transform = config.at("transform", default: identity)
        let shape = config.at("shape", default: false)
        let ratio = config.at("ratio", default: none)
        let shift = config.at("shift", default: 0)
        let local = ctx + (transform: transform)
        let space = ctx + (transform: if shape { identity } else { transform })
        let carrier = curve
          .to-cetz(path, mark: none)
          .first()(space)
          .drawables
          .first()
        if not shape {
          carrier = cetz
            .drawable
            .apply-transform(cetz.matrix.transform-scale((1, 1, 0)), carrier)
            .first()
        }
        let pieces = curve.segments(path)
        let native = if pieces.len() == 1 {
          let piece = pieces.first()
          cetz
            .draw
            .bezier(
              piece.start,
              piece.end,
              piece.control-start,
              piece.control-end,
              mark: none,
            )
            .first()(local)
        } else { curve.to-cetz(path, mark: none).first()(local) }
        let total = cetz.path-util.length(carrier.segments)
        for root in ("start", "end") {
          for symbol in (">", "<", "straight") {
            for reverse in (false, true) {
              for anchor in if ratio == none {
                ("tip", "center", "base")
              } else { ("center",) } {
                let entry = (
                  head + (symbol: symbol, reverse: reverse, anchor: anchor)
                )
                let mark = (transform-shape: shape)
                mark.insert(root, entry)
                let style = (
                  name: "named",
                  stroke: black + 0.8pt,
                  mark: mark,
                  mark-position: if ratio == none { "end" } else { ratio },
                  mark-direction: if root == "start" { "backward" } else {
                    "forward"
                  },
                  mark-shift: shift,
                  accuracy: accuracy,
                )
                let processed = cetz
                  .mark
                  .process-style(
                    ctx,
                    cetz
                      .styles
                      .resolve(
                        ctx.style,
                        merge: drawing._draw-style(style),
                        root: "bezier",
                      )
                      .mark,
                    root,
                    total,
                  )
                  .first()
                let result = cetz.process.many(
                  local,
                  drawing
                    ._derived-path-elements(path, style, auto, true, true)
                    .elements
                    .flatten(),
                  compute-bounds: false,
                )
                let actual = result.drawables

                let marks = actual.filter(item => (
                  cetz.drawable.TAG.mark in item.tags
                ))
                let case = repr((index, root, symbol, reverse, anchor))
                assert.eq(
                  marks.len(),
                  if total == 0 { 0 } else { 1 },
                  message: case + ": head count",
                )
                if total == 0 { continue }
                let (origin, _, commands) = marks.first().segments.first()
                let tip = if symbol == "straight" {
                  commands.first().last()
                } else { origin }
                let wings = if symbol == "straight" {
                  (origin, commands.at(1).last())
                } else { (commands.at(0).last(), commands.at(1).last()) }
                let back = cetz.vector.lerp(..wings, 0.5)
                let reversed = reverse != (symbol == "<")
                let span = calc.min(processed.length, total)
                let tip-at = (tip: 0, center: -span / 2, base: -span).at(anchor)
                let back-at = tip-at + span
                if reversed { (tip-at, back-at) = (back-at, tip-at) }
                let at = if ratio == none { 0 } else {
                  if root == "start" { ratio * total + shift } else {
                    (1 - ratio) * total - shift
                  }
                }
                let low = calc.min(tip-at, back-at) + at
                let inward = calc.clamp(low, 0, total - span) - low
                let stations = (tip-at, back-at).map(s => s + at + inward)
                let contacts = stations.map(s => {
                  cetz
                    .path-util
                    .point-at(
                      carrier.segments,
                      s,
                      reverse: root == "end",
                    )
                    .point
                })
                let axis = cetz.vector.norm(cetz.vector.sub(
                  contacts.last(),
                  contacts.first(),
                ))
                let width = processed.width
                if shape {
                  let normal = (-axis.at(1), axis.at(0), 0)
                  width *= cetz.vector.dist(
                    cetz.matrix.mul4x4-vec3(transform, normal),
                    cetz.matrix.mul4x4-vec3(transform, (0, 0, 0)),
                  )
                  contacts = contacts.map(p => cetz.matrix.mul4x4-vec3(
                    transform,
                    p,
                  ))
                }
                for (actual, expected) in (tip, back).zip(contacts) {
                  assert(
                    cetz.vector.dist(actual, expected) < epsilon,
                    message: case
                      + ": actual tip/back contact "
                      + repr((actual, expected)),
                  )
                }
                assert(
                  calc.abs(cetz.vector.dist(..wings) - width) < epsilon,
                  message: case + ": width must not compress",
                )
                assert.eq(
                  marks.first().stroke.thickness,
                  processed.stroke.thickness,
                  message: case + ": stroke size",
                )
                let shaft = actual.filter(item => (
                  cetz.drawable.TAG.mark not in item.tags
                ))
                if ratio != none {
                  comparisons.push((
                    case + ": continuous interior shaft",
                    shaft.map(d => d.segments),
                    cetz
                      .drawable
                      .apply-transform(
                        if shape { transform } else { identity },
                        carrier,
                      )
                      .map(d => d.segments),
                  ))
                } else if total > processed.length {
                  let contact = if symbol == "straight" or reversed {
                    tip
                  } else { back }
                  let endpoint = if root == "start" {
                    cetz.path-util.first-subpath-start(shaft.first().segments)
                  } else {
                    cetz.path-util.last-subpath-end(shaft.last().segments)
                  }
                  assert(
                    cetz.vector.dist(endpoint, contact) < epsilon,
                    message: case + ": painted shaft contact",
                  )
                }
                // Named anchors continue to describe the unshortened carrier.
                for anchor in (
                  ("start", "end")
                    + if pieces.len() == 1 { ("ctrl-0", "ctrl-1") } else { () }
                ) {
                  let point = (result.ctx.nodes.at("named").anchors)(anchor)
                  let expected = (native.anchors)(anchor)
                  assert(
                    cetz.vector.dist(point, expected) < epsilon,
                    message: case + ": named anchor",
                  )
                }
              }
            }
          }
        }
      }

      // Both ends can shorten the same cubic, or consume the entire short shaft.
      for path in (
        curved,
        curve.from-cubic(curve.line-segment((0, 0), (0.1, 0))),
      ) {
        let style = (mark: (symbol: ">", ..head))
        let actual = cetz
          .process
          .many(
            ctx,
            drawing
              ._derived-path-elements(path, style, auto, true, true)
              .elements
              .flatten(),
            compute-bounds: false,
          )
          .drawables
        assert.eq(actual.len(), 3, message: "two heads retained")
        if actual.first().segments.len() > 0 {
          for (mark, endpoint) in actual
            .slice(1)
            .zip((
              cetz.path-util.first-subpath-start(actual.first().segments),
              cetz.path-util.last-subpath-end(actual.first().segments),
            )) {
            let commands = mark.segments.first().last()
            assert(
              cetz.vector.dist(endpoint, cetz.vector.lerp(
                commands.at(0).last(),
                commands.at(1).last(),
                0.5,
              ))
                < epsilon,
              message: "two-ended shaft contacts",
            )
          }
        } else {
          assert(
            curve.length(path, accuracy: accuracy) < 2 * head.length,
            message: "only overlapping heads consume the shaft",
          )
        }
      }

      // Mnemonics and custom overrides resolve before built-in contact handling.
      let custom = ctx
      custom.marks.marks.insert("straight", cetz.mark-shapes.marks.diamond)
      custom.marks.mnemonics.insert("arrow-alias", "<")
      for (local, symbol, target) in (
        (custom, "straight", "diamond"),
        (custom, "arrow-alias", "<"),
      ) {
        for reverse in (false, true) {
          let variants = (symbol, target).map(symbol => {
            let style = (
              mark: (
                end: head
                  + (symbol: symbol, reverse: reverse, flip: true, slant: 30%),
              ),
              mark-position: "center",
            )
            cetz
              .process
              .many(
                local,
                drawing
                  ._derived-path-elements(curved, style, auto, true, true)
                  .elements
                  .flatten(),
                compute-bounds: false,
              )
              .drawables
          })
          assert.eq(..variants, message: "resolved symbol " + symbol)
        }
      }

      // One prepared carrier must still resolve each caller's marks, inherited
      // stroke, and transform, without retaining a previous caller's context.
      let reusable = drawing._mark-carrier-elements(
        curved,
        (mark: (end: head + (symbol: "straight"))),
        paint: true,
      ).flatten()
      let reference = cetz.process.many(ctx, reusable, compute-bounds: false).drawables
      let styled = cetz.draw.set-style(stroke: red + 1.2pt).first()(ctx).ctx
      let transformed = ctx + (transform: cetz.matrix.transform-scale((2, 0.5, 1)))
      for local in (custom, styled, transformed) {
        let actual = cetz.process.many(local, reusable, compute-bounds: false).drawables
        assert.ne(actual, reference, message: "carrier resolves the current context")
        assert.eq(
          cetz.process.many(ctx, reusable, compute-bounds: false).drawables,
          reference,
          message: "carrier evaluation does not retain another context",
        )
      }

      // Patterns and crossing gaps only change paint, never the single full-path mark.
      for position in ("end", "center") {
        for symbol in (">", "straight") {
          let style = (
            mark: (end: head + (symbol: symbol)),
            mark-position: position,
          )
          let reference = cetz
            .process
            .many(
              ctx,
              drawing
                ._derived-path-elements(curved, style, auto, true, true)
                .elements
                .flatten(),
              compute-bounds: false,
            )
            .drawables
            .filter(d => cetz.drawable.TAG.mark in d.tags)
          for pattern in (none, "wave", "coil") {
            let paint = style + (pattern: pattern, crossing-gap: 0.4)
            let total = curve.length(curved, accuracy: accuracy)
            for elements in (
              drawing
                ._derived-path-elements(curved, paint, auto, true, true)
                .elements,
              drawing._cut-path-elements(
                curved,
                paint,
                (total / 2,),
                mark-style: style,
              ),
            ) {
              let marks = cetz
                .process
                .many(ctx, elements.flatten(), compute-bounds: false)
                .drawables
                .filter(d => cetz.drawable.TAG.mark in d.tags)
              assert.eq(
                marks,
                reference,
                message: repr((position, symbol, pattern))
                  + ": undistorted carrier",
              )
            }
          }
        }
      }

      // Coordinate resolvers apply to the curve, never to the mark's private anchors.
      let local = cetz
        .draw
        .register-coordinate-resolver((ctx, point) => {
          if type(point) == array { point.map(value => value * 2) } else {
            point
          }
        })
        .first()(ctx)
        .ctx
      let style = (mark: (end: head + (symbol: ">")), mark-position: "center")
      let resolved = cetz
        .process
        .many(
          local,
          drawing
            ._derived-path-elements(straight, style, auto, true, true)
            .elements
            .flatten(),
          compute-bounds: false,
        )
        .drawables
      let doubled = curve.from-cubic(curve.line-segment((0, 0), (4, 0)))
      let expected = cetz
        .process
        .many(
          ctx,
          drawing
            ._derived-path-elements(doubled, style, auto, true, true)
            .elements
            .flatten(),
          compute-bounds: false,
        )
        .drawables
      comparisons.push((
        "coordinate resolver",
        resolved.map(d => d.segments),
        expected.map(d => d.segments),
      ))

      // Endpoint clamps remain exact even below the arc-length accuracy.
      let label-accuracy = drawing._style-value((:), "accuracy")
      for length in (0, label-accuracy / 2) {
        let segment = curve.line-segment((0, 0), (length, 0))
        let path = curve.from-cubic(segment)
        for (shift, t) in ((-1, 0), (1, 1)) {
          let frame = drawing._path-mid-frame(path, label-accuracy, shift: shift)
          assert.eq(frame.point, curve.cubic-point(segment, t))
          assert.eq(frame.tangent, curve.cubic-tangent(segment, t))
        }
      }

      // Measure the rendered frame for automatic clearance and the actual CeTZ
      // anchor for explicit placement, never plain Typst text metrics.
      let label-checks = ()
      let halves = (
        source: curve.segments(source),
        sink: curve.segments(sink),
        whole: curve.segments(path),
      )
      for (path-index, (path, paired)) in (
        (straight, none),
        (curved, none),
        (path, (halves: halves, source: true, sink: true)),
        (path, (halves: halves + (whole: none), source: true, sink: true)),
        (source, (halves: halves, source: true, sink: false)),
        (sink, (halves: halves, source: false, sink: true)),
      ).enumerate() {
        for side in ("left", "right") {
          for (label-index, label) in (
            [$k-p_2$],
            [Hgyp],
            [Hg\ yp],
            box(width: 6pt, height: 6pt),
            box(width: 36pt, height: 6pt),
          ).enumerate() {
            for (style-index, label-style) in (
              (:),
              (anchor: auto),
              (anchor: "auto"),
              (anchor: "center"),
              (anchor: "east"),
              (
                wrap: body => box(width: 18pt, text(size: 9pt, body)),
                padding: (left: 0.15, right: 0.25, top: 0.1, bottom: 0.2),
              ),
              (anchor: "east", angle: 13deg, padding: (left: 0.15, top: 0.2)),
              (angle: -23deg, auto-scale: true),
              (angle: 37deg, anchor: "north-east"),
              (angle: -23deg, anchor: "south-west", auto-scale: true),
            ).enumerate() {
              label-checks.push((
                repr((path-index, side, label-index, style-index)),
                path,
                (
                  label: text(size: 6pt, label),
                  label-style: label-style,
                  label-side: side,
                  label-gap: 0.2,
                ),
                paired,
              ))
            }
          }
        }
      }

      // Fixed xbox-opened2 sink geometry: label shifts are absolute on the full
      // offset path. Before annealing, equal and independent shifts must not move
      // the arrow, and its length must not change the label or endpoint clamps.
      let endpoints = (
        (3.7141557137182426, -3.510279312439297),
        (1.5082746017729896, -2.213011602860704),
      )
      let native = ctx + (length: 2.6mm)
      for (side, shift, kind) in (
        ("right", -1, "incoming"), ("right", -0.4, "incoming"),
        ("right", 1, "incoming"), ("left", -0.4, "incoming"),
        ("right", -0.4, "outgoing"), ("right", 0.4, "paired"),
        ("left", -0.4, "paired"),
      ) {
        let path = if kind == "paired" { path } else {
          drawing._dangling-path(
            endpoints.at(if kind == "incoming" { 0 } else { 1 }),
            endpoints.at(if kind == "incoming" { 1 } else { 0 }),
            0,
            (dangling-tangent: "horizontal"),
            dangling-at-start: kind == "incoming",
          )
        }
        let fields = (
          momentum-arrow-side: side,
          momentum-arrow-offset: 0.4,
          momentum-arrow-length: 0.7,
          momentum-arrow-shift: shift,
        )
        let reference = feynman.edge-style((momentum: [], fields: fields)).at(1).label-path
        let reference-path = drawing._path-layer(path, reference, 0, 0, none, auto)
        let reference-paint = cetz
          .process
          .many(
            native,
            drawing
              ._derived-path-elements(
                reference-path,
                reference,
                auto,
                true,
                true,
              )
              .elements
              .flatten(),
            compute-bounds: false,
          )
          .drawables
        let full-style = reference + (
          length: none, ratio: none, resolve-length: "none", shift: 0,
        )
        let full-path = drawing._path-layer(path, full-style, 0, 0, none, auto)
        let full-length = curve.length(
          full-path, accuracy: drawing._style-value(full-style, "accuracy"),
        )
        let paired = if kind == "paired" {
          (
            halves: drawing._split-edge-geometry(
              source, sink, path, full-style, full-style, 0, 0, none,
              drawing._style-value(full-style, "accuracy"),
            ),
            source: true,
            sink: true,
          )
        } else { none }
        for (gap, label, anchor) in (
          (0.2, [$k-p_2$], auto),
          (0.6, [Hg\ yp], "auto"),
          (0, box(width: 6pt, height: 6pt), auto),
          (-0.2, [$k-p_2$], auto),
          (0.2, box(width: 36pt, height: 6pt), "center"),
          (0.2, box(width: 36pt, height: 6pt), "east"),
          (0.2, box(width: 36pt, height: 6pt), "north-east"),
          (0.2, box(width: 36pt, height: 6pt), "south-west"),
          (-0.2, box(width: 36pt, height: 6pt), "center"),
        ) {
          let momentum = text(size: 6pt, label)
          let centers = ()
          for requested in (
            none, shift, shift - 1e-6, shift + 1e-6, -0.1, 0,
            -full-length, full-length,
          ) {
            let label-shift = if requested == none { shift } else { requested }
            let case = "momentum " + repr((side, shift, requested, gap, anchor, kind))
            let fields = fields + (
              momentum-label-gap: gap,
              momentum-label-anchor: anchor,
            ) + if requested == none { (:) } else {
              (momentum-label-shift: requested)
            }
            let layers = feynman.edge-style((momentum: momentum, fields: fields))
            let carrier = layers.last()
            let arrow = carrier.label-path
            assert.eq(layers.len(), 2, message: case + ": one paired momentum annotation")
            assert.eq(
              (carrier.length, carrier.ratio, carrier.resolve-length, carrier.shift),
              (none, none, "none", 0),
              message: case + ": full unshifted label carrier",
            )
            assert.eq(
              carrier.label-shift,
              label-shift,
              message: case + ": requested absolute label shift",
            )
            assert.eq(carrier.stroke, none, message: case + ": invisible carrier")
            assert.eq(carrier.mark, none, message: case + ": no duplicate head")
            assert.eq(
              arrow.at("label", default: none), none,
              message: case + ": no label attached to arrow",
            )
            let arrow-path = drawing._path-layer(path, arrow, 0, 0, none, auto)
            let label-path = drawing._path-layer(path, arrow + (
              length: none, ratio: none, resolve-length: "none", shift: 0,
            ), 0, 0, none, auto)
            assert.eq(
              arrow-path,
              reference-path,
              message: case + ": unchanged arrow path",
            )
            assert.eq(
              cetz
                .process
                .many(
                  native,
                  drawing
                    ._derived-path-elements(arrow-path, arrow, auto, true, true)
                    .elements
                    .flatten(),
                  compute-bounds: false,
                )
                .drawables,
              reference-paint,
              message: case + ": unchanged shaft and chord head",
            )
            assert.eq(label-path, full-path, message: case + ": full offset path")
            let frame = drawing._path-mid-frame(
              label-path, drawing._style-value(carrier, "accuracy"),
              shift: carrier.label-shift,
            )
            let rendered = cetz
              .process
              .many(
                native,
                inspect-layer-label-element(native, path, carrier, none, (:))
                  .flatten(),
                compute-bounds: false,
              )
              .drawables
            let labels = rendered.filter(d => d.type == "content")
            assert.eq(labels.len(), 1, message: case + ": one momentum label")
            centers.push(labels.first().pos)
            label-checks.push((case, label-path, carrier + (label-path: none), paired))
            for length in (0.2, 0.8 * full-length, 2 * full-length) {
              let case = case + " length " + repr(length)
              let layers = feynman.edge-style((
                momentum: momentum,
                fields: fields + (momentum-arrow-length: length),
              ))
              assert.eq(
                layers.at(1).label-path, arrow + (length: length),
                message: case + ": only arrow length changes",
              )
              let carrier = layers.last()
              let label-path = drawing._path-layer(path, carrier.label-path + (
                length: none, ratio: none, resolve-length: "none", shift: 0,
              ), 0, 0, none, auto)
              assert.eq(
                label-path, full-path,
                message: case + ": length-independent path",
              )
              let shifted = drawing._path-mid-frame(
                label-path, drawing._style-value(carrier, "accuracy"),
                shift: carrier.label-shift,
              )
              assert.eq(
                shifted.point, frame.point,
                message: case + ": unchanged label point",
              )
              assert.eq(
                cetz.vector.norm(shifted.tangent), cetz.vector.norm(frame.tangent),
                message: case + ": unchanged label normal",
              )
              let resized = cetz.process.many(
                native,
                inspect-layer-label-element(native, path, carrier, none, (:)).flatten(),
                compute-bounds: false,
              ).drawables
              if anchor in (auto, "auto") {
                let center = resized.filter(d => d.type == "content").first().pos
                let delta = cetz.vector.sub(center, labels.first().pos)
                assert(calc.abs(cetz.vector.dot(delta, (..cetz.vector.norm(frame.tangent), 0))) < epsilon,
                  message: case + ": finite-arrow clearance preserves the label arc coordinate")
              } else {
                assert.eq(resized, rendered, message: case + ": explicit anchor preserves final label placement")
              }
            }
          }
          for center in centers.slice(1, 4) {
            assert(
              cetz.vector.dist(centers.first(), center) < 8 * accuracy,
              message: "equal/adjacent label continuity "
                + repr((side, shift, gap, anchor, kind, centers)),
            )
          }
        }
      }
      for (transform-index, transform) in (
        identity,
        cetz.matrix.mul-mat(
          cetz.matrix.transform-translate(3, 2, 0),
          cetz.matrix.transform-rotate-z(37deg),
        ),
        cetz.matrix.transform-scale((-1.5, 0.6, 1)),
        cetz.matrix.transform-shear-x(0.6),
        cetz.matrix.transform-rotate-xyz(40deg, 30deg, 0deg),
      ).enumerate() {
        let local = native + (transform: transform)
        for (case, path, style, paired) in label-checks {
          let case = "label placement " + repr(transform-index) + " " + case
          let label-style = style.label-style + (name: "checked-label")
          let style = style + (label-style: label-style)
          let anchor = label-style.at("anchor", default: auto)
          let automatic = anchor in (auto, "auto")
          let gap = calc.max(0, style.label-gap)
          let accuracy = drawing._style-value(style, "accuracy")
          let label-shift = style.at("label-shift", default: 0)
          let frame = drawing._path-mid-frame(path, accuracy, shift: label-shift)
          let total = curve.length(path, accuracy: accuracy)
          let at = calc.clamp(total / 2 + label-shift, 0, total)
          let (segment, t) = if at == 0 {
            (curve.segments(path).first(), 0)
          } else if at == total {
            (curve.segments(path).last(), 1)
          } else {
            (curve.segments(curve.trim(
              path, end-outset: total - at, accuracy: accuracy,
            )).last(), 1)
          }
          let tolerance = if at in (0, total) { epsilon } else { 8 * accuracy }
          assert(
            cetz.vector.dist(frame.point, curve.cubic-point(segment, t)) < tolerance,
            message: case + ": shifted arc-length point and full-path endpoint clamp",
          )
          let tangent = cetz.vector.norm(frame.tangent)
          assert(
            cetz.vector.dist(
              tangent, cetz.vector.norm(curve.cubic-tangent(segment, t)),
            ) < tolerance,
            message: case + ": local tangent at shifted label point",
          )
          let normal = cetz.vector.scale(
            (-tangent.at(1), tangent.at(0)),
            if style.label-side == "right" { -1 } else { 1 },
          )
          let origin = cetz.matrix.mul4x4-vec3(transform, (..frame.point, 0))
          let outward = cetz.vector.sub(
            cetz.matrix.mul4x4-vec3(transform, (
              ..cetz.vector.add(frame.point, normal),
              0,
            )),
            origin,
          )
          let element = if paired == none {
            inspect-layer-label-element(local, path, style, none, (:))
          } else {
            drawing._paired-layer-label-element(
              local,
              paired.halves,
              if paired.source { style } else { none },
              if paired.sink { style } else { none },
              none,
              (:),
            )
          }
          let rendered = cetz.process.many(
            local, element.flatten(), compute-bounds: false,
          )
          let point = (rendered.ctx.nodes.at("checked-label").anchors)(
            if automatic { "center" } else { anchor },
          )
          let distance = if automatic {
            (
              cetz.vector.dot(cetz.vector.sub(point, origin), outward)
                / cetz.vector.dot(outward, outward)
            )
          } else { gap }
          assert(
            cetz.vector.dist(
              point,
              cetz.vector.add(origin, cetz.vector.scale(outward, distance)),
            ) < epsilon,
            message: case + if automatic { ": centered legacy position" } else {
              ": explicit anchor at shifted arc-length point plus signed normal gap"
            },
          )
          let frames = rendered.drawables.filter(d => "content-frame" in d.tags)
          assert.eq(frames.len(), 1, message: case + ": one rendered box")
          if automatic {
            let (start, _, commands) = frames.first().segments.first()
            let points = (start, ..commands.map(command => command.last())).dedup()
            let actual = attachment-distance(drawing._label-path-lines(local, path, accuracy), points) / calc.sqrt(outward.at(0) * outward.at(0) + outward.at(1) * outward.at(1))
            assert(
              calc.abs(actual - gap) < 8 * accuracy,
              message: case + ": finite-path rendered clearance " + repr((actual, gap)),
            )
          }
        }
      }

      // Momentum arrows and labels use one collision choice, retaining their
      // independent manual shifts as a fixed relative arc-length displacement.
      for carrier in (curve.line((0, 0), (8, 0)), path) {
        for side in ("left", "right", auto) {
          let annotation = feynman.edge-style((
            momentum: text(size: 6pt, [$k-p_2$]),
            fields: feynman.momentum(
              side: side, offset: 0.4, length: 1.2, shift: 0.35,
              label: (shift: -0.25, gap: 0.2),
            ),
          )).last()
          let placement = inspect-layer-label-element(
            native, carrier, annotation, none, (eid: 7), placement: true,
          )
          assert.eq(placement.edge, 7, message: "momentum pair owns its physical edge")
          assert.eq(placement.paths.len(), if side == auto { 2 } else { 1 }, message: "only automatic momentum sides can flip")
          let chosen = relax-inspected-labels((placement, placement), ())
          assert.eq(chosen, relax-inspected-labels((placement, placement), ()), message: "deterministic momentum pair annealing")
          assert(chosen.any(candidate => (
            candidate.path-index != placement.candidates.first().path-index
              or candidate.path-shift != placement.candidates.first().path-shift
          )), message: "collisions move the momentum arrow together with its label")
          let (a, b) = chosen.map(candidate => candidate.bounds)
          assert(a.right <= b.left or b.right <= a.left or a.top <= b.bottom or b.top <= a.bottom, message: "momentum pairs separate crowded labels")
          for candidate in placement.candidates {
            let full = placement.paths.at(candidate.path-index)
            let accuracy = drawing._style-value(annotation, "accuracy")
            let total = curve.length(full, accuracy: accuracy)
            let arrow-center = total / 2 + candidate.path-shift
            assert(calc.abs(arrow-center - candidate.at - 0.6) < epsilon, message: "arrow and label retain independent relative shifts")
            assert(arrow-center >= 0.6 - 8 * accuracy and arrow-center <= total - 0.6 + 8 * accuracy, message: "annealed arrow remains inside the full offset carrier")
            let arrow-style = placement.path-style + (
              offset: 0, offset-side: none, shift: candidate.path-shift,
            )
            let arrow-path = drawing._path-layer(full, arrow-style, 0, 0, none, auto)
            assert(calc.abs(curve.length(arrow-path, accuracy: accuracy) - 1.2) < 8 * accuracy, message: "sliding does not shorten the momentum arrow")
            let arrow-paint = cetz.process.many(
              native,
              drawing._derived-path-elements(arrow-path, arrow-style, auto, true, true).elements.flatten(),
              compute-bounds: false,
            ).drawables
            assert(arrow-paint.any(item => cetz.drawable.TAG.mark in item.tags), message: "every momentum placement retains its arrowhead")
            let frame = drawing._path-mid-frame(full, accuracy, shift: candidate.at - total / 2)
            let tangent = cetz.vector.norm(frame.tangent)
            let normal = cetz.vector.scale((-tangent.at(1), tangent.at(0)), candidate.side)
            let origin = cetz.matrix.mul4x4-vec3(native.transform, (..frame.point, 0))
            let outward = cetz.vector.sub(
              cetz.matrix.mul4x4-vec3(native.transform, (..cetz.vector.add(frame.point, normal), 0)),
              origin,
            )
            let label = cetz.draw.content(candidate.position, placement.label, padding: 0, ..placement.style).first()(native)
            let corners = ("north-west", "north-east", "south-east", "south-west").map(label.anchors)
            let attachment = drawing._label-path-lines(native, arrow-path, accuracy) + drawing._label-path-lines(native, carrier, accuracy)
            let actual = attachment-distance(attachment, corners) / calc.sqrt(outward.at(0) * outward.at(0) + outward.at(1) * outward.at(1))
            assert(calc.abs(actual - 0.2) < 8 * accuracy, message: "momentum sliding keeps fixed clearance to the finite arrow and physical carrier")
          }
          if side == auto {
            let preferred = placement.candidates.first().side
            let obstacles = placement.candidates.filter(candidate => candidate.side == preferred).map(candidate => candidate.bounds)
            let flipped = relax-inspected-labels((placement,), obstacles).first()
            assert.eq(flipped.side, -preferred, message: "blocked momentum pair flips its arrow and label together")
            assert.eq(flipped, relax-inspected-labels((placement,), obstacles).first(), message: "deterministic momentum side choice")
          }
          for pinned-style in ((label-slide: false), (label-style: (anchor: "east"))) {
            let pinned = inspect-layer-label-element(
              native, carrier, annotation + pinned-style, none, (:), placement: true,
            )
            assert.eq(pinned.candidates.len(), 1, message: "manual momentum placement stays pinned")
            assert(calc.abs(pinned.candidates.first().path-shift - 0.35) < epsilon, message: "pinned momentum preserves its arrow shift")
          }
        }
      }

      // An arrow-only collision must move the pair even when its text is clear.
      // Conversely, empty space between the two pieces is not an obstacle.
      let paired-style = (
        label: box(width: 2pt, height: 2pt, fill: black),
        label-side: "left", label-gap: 1,
        label-path: (
          offset: 0.5, length: 1.4, ratio: 0.5, resolve-length: "min",
          stroke: black + 0.4pt, mark: (end: "straight"), mark-position: 1,
        ),
      )
      let carrier = curve.line((0, 0), (8, 0))
      let pair = inspect-layer-label-element(native, carrier, paired-style, none, (eid: 0), placement: true)
      let first = pair.candidates.first()
      assert(first.arrow-bounds.len() >= 2, message: "shaft and arrowhead have separate collision bounds")
      let blocker = first.arrow-bounds.last()
      assert(first.bounds.bottom > blocker.top + 0.08, message: "arrow-only obstacle is clear of the text")
      let text-only = pair + (candidates: pair.candidates.map(c => c + (arrow-bounds: ())))
      assert.eq(relax-inspected-labels((text-only,), (blocker,), label-padding: 0).first().path-shift, first.path-shift, message: "clear text alone would not move")
      let shifted = relax-inspected-labels((pair,), (blocker,), label-padding: 0).first()
      assert.ne(shifted.path-shift, first.path-shift, message: "arrowhead collision moves the annotation")
      assert(calc.abs(shifted.position.at(0) - first.position.at(0) - shifted.path-shift + first.path-shift) < epsilon, message: "arrow and label have exactly the same tangential displacement")
      assert.eq(shifted.position.at(1), first.position.at(1), message: "shared movement keeps the normal gap")
      let middle = (first.bounds.bottom + calc.max(..first.arrow-bounds.map(b => b.top))) / 2
      let gap-obstacle = (left: 3.95, right: 4.05, bottom: middle - 0.05, top: middle + 0.05)
      assert.eq(relax-inspected-labels((pair,), (gap-obstacle,), label-padding: 0).first().path-shift, first.path-shift, message: "unoccupied arrow/text gap is not one enclosing collision box")
      // A fixed standalone arrow has the same clearance as a pinned paired
      // arrow, and its actual marker remains part of the obstacle footprint.
      let fixed-paint = drawing._derived-path-elements(
        curve.line((4, -0.2), (4, 1.2)),
        (stroke: black + 0.4pt, mark: (end: "straight"), mark-position: 1),
        auto, true, true,
      ).elements
      let fixed-bounds = drawing._annotation-bounds(native, (fixed-paint,)).first()
      assert(fixed-bounds.len() >= 2, message: "fixed arrows retain shaft and marker footprints")
      let shifted-paint = (cetz.draw.translate((10, 0)), fixed-paint)
      let grouped-bounds = drawing._annotation-bounds(native, (shifted-paint, fixed-paint, cetz.draw.hide(fixed-paint.flatten())))
      assert.eq(grouped-bounds.at(1), fixed-bounds, message: "candidate groups do not inherit another arrow's transform")
      assert.eq(grouped-bounds.at(2), (), message: "hidden candidate groups remain empty")
      assert.eq(grouped-bounds.first(), drawing._annotation-bounds(native, (shifted-paint,)).first(),
        message: "batched and individual candidate footprints agree")

      // Prepared candidates share native placement while retaining the same
      // independent contexts and actual custom marker geometry as rendering.
      let prepared-fixed = drawing._derived-path(
        curve.line((4, -0.2), (4, 1.2)),
        (stroke: black + 0.4pt, mark: (end: "straight"), mark-position: 1),
        auto, true, true,
      )
      let prepared-shifted = prepared-fixed + (elements: (cetz.draw.translate((10, 0)), prepared-fixed.elements))
      let prepared-empty = drawing._derived-path(curve.path(), (:), auto, true, true)
      let custom-native = native
      custom-native.marks.marks.insert("straight", cetz.mark-shapes.marks.diamond)
      for local in (native, custom-native, native + (transform: cetz.matrix.transform-scale((-2, 0.5, 1)))) {
        let actual = drawing._annotation-path-bounds(local, (prepared-shifted, prepared-fixed, prepared-empty))
        let expected = drawing._annotation-bounds(local, (shifted-paint, fixed-paint, ()))
        assert.eq(cbor.encode(actual), cbor.encode(expected),
          message: "prepared candidate batches preserve custom marks, transforms and empty carriers")
        assert.eq(actual.at(1), drawing._annotation-path-bounds(local, (prepared-fixed,)).first(),
          message: "an earlier prepared candidate cannot alter a later candidate's geometry")
      }

      let avoiding-fixed = relax-inspected-labels((pair,), (), label-padding: 0, fixed-arrows: fixed-bounds).first()
      assert.ne(avoiding-fixed.path-shift, first.path-shift, message: "a fixed standalone arrow moves the paired annotation")
      let fixed-pair = (candidates: ((bounds: (left: 100, right: 101, bottom: 100, top: 101), arrow-bounds: fixed-bounds, cost: 0),))
      assert.eq(
        avoiding-fixed.position,
        relax-inspected-labels((pair, fixed-pair), (), label-padding: 0).first().position,
        message: "fixed and pinned annotation arrows have identical collision clearance",
      )
      assert.eq(drawing._annotation-bounds(native, (cetz.draw.hide(fixed-paint.flatten()),)).first(), (), message: "hidden arrow geometry is not an obstacle")

      let blocker-label = inspect-layer-label-element(
        native, carrier, (label: box(width: 2pt, height: 2pt, fill: black)), none, (eid: 1),
        fixed-position: (4, 0.5), placement: true,
      )
      assert.ne(relax-inspected-labels((pair, blocker-label), (), label-padding: 0).first().path-shift, first.path-shift, message: "other labels repel the arrow even when both text boxes are clear")
      let opposing = (-2, 2).map(shift => inspect-layer-label-element(
        native, carrier, paired-style + (label-shift: shift), none, (:), placement: true,
      ))
      let separated = relax-inspected-labels(opposing, (), label-padding: 0)
      assert(separated.zip(opposing).any(((choice, original)) => choice.path-shift != original.candidates.first().path-shift), message: "overlapping arrows separate even when their text is apart")
      for (length, expected) in ((8, 1.4), (2, 1), (0.8, 0.4)) {
        let capped = drawing._path-layer(curve.line((0, 0), (length, 0)), paired-style.label-path, 0, 0, none, auto)
        assert(calc.abs(curve.length(capped) - expected) < 1e-6, message: "arrow length is capped by its carrier percentage")
      }

      // Crowded labels may slide or switch sides on a straight or curved carrier
      // at fixed clearance. Check the rendered text bounds after either move.
      for carrier in (curve.line((0, 0), (8, 0)), path) {
        for side in ("left", "right") {
          let style = (label: [AB], label-gap: 0.3, label-side: side)
          let label = inspect-layer-label-element(native, carrier, style, none, (eid: 0), placement: true)
          assert.eq(label.edge, 0, message: "path labels retain carrier ownership")
          let total = curve.length(carrier, accuracy: 0.001)
          let chosen = relax-inspected-labels((label, label), ())
          assert.eq(chosen, relax-inspected-labels((label, label), ()), message: "deterministic label annealing")
          assert(chosen.any(candidate => candidate.position != label.candidates.first().position), message: "crowded labels slide")
          let (a, b) = chosen.map(candidate => candidate.bounds)
          assert(a.right <= b.left or b.right <= a.left or a.top <= b.bottom or b.top <= a.bottom, message: "sliding separates overlapping labels")
          assert.eq(relax-inspected-labels((label,), ()).first().position, label.candidates.first().position, message: "uncrowded labels keep their preferred position")
          let sign = if side == "left" { 1 } else { -1 }
          assert(label.candidates.all(candidate => candidate.side == sign), message: "explicit label sides stay fixed")
          let frame = drawing._path-mid-frame(carrier, 0.001)
          let tangent = cetz.vector.norm(frame.tangent)
          let preferred = cetz.vector.add(frame.point, cetz.vector.scale((-tangent.at(1), tangent.at(0)), sign))
          let automatic = inspect-layer-label-element(
            native, carrier, style + (label-side: auto), preferred,
            (edge: (pos: frame.point)), placement: true,
          )
          assert.eq(relax-inspected-labels((automatic,), ()).first().side, sign, message: "uncrowded labels keep their preferred side")
          let obstacles = automatic.candidates.filter(candidate => candidate.side == sign).map(candidate => candidate.bounds)
          let flipped = relax-inspected-labels((automatic,), obstacles).first()
          assert.eq(flipped.side, -sign, message: "labels flip when the preferred side is blocked")
          assert.eq(flipped, relax-inspected-labels((automatic,), obstacles).first(), message: "deterministic side flipping")
          for candidate in (..chosen, flipped) {
            assert(candidate.at >= 0 and candidate.at <= total, message: "labels stay on the carrier")
            let frame = drawing._path-mid-frame(carrier, 0.001, shift: candidate.at - total / 2)
            let tangent = cetz.vector.norm(frame.tangent)
            let normal = cetz.vector.scale((-tangent.at(1), tangent.at(0)), candidate.side)
            let origin = cetz.matrix.mul4x4-vec3(native.transform, (..frame.point, 0))
            let outward = cetz.vector.sub(cetz.matrix.mul4x4-vec3(native.transform, (..cetz.vector.add(frame.point, normal), 0)), origin)
            let element = cetz.draw.content(candidate.position, label.label, padding: 0, ..label.style)
            let measured = element.first()(native)
            let corners = ("north-west", "north-east", "south-east", "south-west").map(measured.anchors)
            let actual = attachment-distance(drawing._label-path-lines(native, carrier, 0.001), corners)
            assert(calc.abs(actual - 0.3) < 0.008, message: "sliding preserves actual clearance to the finite carrier")
          }
          let blocked = relax-inspected-labels((label,), (label.candidates.first().bounds,))
          assert.ne(blocked.first().position, label.candidates.first().position, message: "labels slide away from node obstacles")
          for pinned-style in ((label-slide: false, label-side: auto), (label-style: (anchor: "center"), label-side: auto)) {
            let pinned = inspect-layer-label-element(native, carrier, style + pinned-style + (label-shift: 1000), none, (:), placement: true)
            assert.eq(pinned.candidates.len(), 1, message: "explicit placements stay fixed")
            assert.eq(pinned.candidates.first().at, total, message: "explicit shifts clamp to the carrier endpoint")
          }
        }
      }

      // A near miss at base clearance becomes a collision when label boxes grow.
      let small = (bounds: (left: 0, right: 1, bottom: 0, top: 1), arrow-bounds: (), cost: 0)
      let moved = (bounds: (left: -2, right: -1, bottom: 0, top: 1), arrow-bounds: (), cost: 0.001)
      let obstacle = (left: 1.2, right: 2.2, bottom: 0, top: 1)
      let movable = (candidates: (small, moved))
      assert.eq(relax-inspected-labels((movable,), (obstacle,), label-padding: 0).first(), small, message: "zero padding keeps a clear label in place")
      assert.eq(relax-inspected-labels((movable,), (obstacle,), label-padding: 0.3).first().bounds, moved.bounds, message: "padding repels labels from nearby obstacles")
      // Inspect a fixed candidate too: the enclosing label must not enlarge the obstacle.
      let fixed = (candidates: (small,))
      let scored = relax-inspected-labels((fixed, movable), ((left: 0.3, right: 0.7, bottom: 0.3, top: 0.7),), label-padding: 0.3)
      assert(calc.abs(scored.first().cost - calc.pow(0.48, 2)) < epsilon, message: "collision cost matches the displayed contained obstacle box")
      // Padding must not move a label away from ordinary edge carriers.
      let owned = movable + (edge: 0)
      for padding in (0, 0.25, 0.6) {
        let own-edge = small.bounds + (edge: 0, self-loop: false)
        for carrier in (own-edge, own-edge + (edge: 1)) {
          assert.eq(relax-inspected-labels((owned,), (carrier,), label-padding: padding).first(), small, message: "ordinary edges do not cause label drift")
        }
        for blocking in (own-edge + (self-loop: true), small.bounds) {
          assert.eq(relax-inspected-labels((owned,), (blocking,), label-padding: padding).first().bounds, moved.bounds, message: "self-loops and nodes still repel labels")
        }
      }
      let same-edge = relax-inspected-labels((owned, owned), (), label-padding: 0.25)
      assert.ne(same-edge.first().bounds, same-edge.last().bounds, message: "labels on the same edge still repel each other")
      let neighbour = (candidates: ((bounds: obstacle, arrow-bounds: (), cost: 0),))
      // Separate the labels beyond the base pair padding, but inside added padding.
      neighbour.candidates.at(0).bounds.left = 1.5
      neighbour.candidates.at(0).bounds.right = 2.5
      assert.eq(relax-inspected-labels((movable, neighbour), (), label-padding: 0).first(), small, message: "base pair clearance is preserved")
      assert.eq(relax-inspected-labels((movable, neighbour), (), label-padding: 0.3).first().bounds, moved.bounds, message: "padding also repels labels from each other")

      let offset-label = inspect-layer-label-element(
        native, curve.line((0, 0), (8, 0)),
        (label: [AB], offset: -0.2), (4, 0), (edge: (pos: (4, 0))), placement: true,
      )
      assert.eq(offset-label.candidates.first().side, -1, message: "signed offsets choose the preferred automatic side")
      assert(offset-label.candidates.any(candidate => candidate.side == 1), message: "inferred sides remain free to flip")

      // Compare drawable structure and points in CeTZ's resolved 3D coordinates.
      for (case, actual, expected) in comparisons {
        assert.eq(
          actual,
          expected,
          message: case + ": unchanged shaft geometry",
        )
      }
      ()
    })
  })
  // Exercise fixed-arrow and actual text-edge contact through the complete
  // drawing path. A clear ordinary carrier and invisible decoration stay inert.
  for mode in ("arrow", "shaft", "ordinary", "clear-ordinary", "hidden") {
    let incoming-x = if mode == "clear-ordinary" { 4.4 } else { 4 }
    let scene = graph.build({
      graph.node(<a>, pos: graph.pos(x: 0, y: 0))
      graph.node(<b>, pos: graph.pos(x: 8, y: 0))
      graph.node(<c>, pos: graph.pos(x: incoming-x, y: 4))
      graph.edge(<pair>, graph.source(<a>), graph.sink(<b>), pos: graph.pos(x: 4, y: 0))
      graph.edge(<incoming>, graph.sink(<c>), pos: graph.pos(x: incoming-x, y: -4))
    })
    draw(scene, unit: 10pt, node-label: none, edge-label: none,
      node-style: (radius: 0, stroke: none, fill: none), node-outset: 0,
      label-collision-padding: 0,
      edge-style: edge => if edge.eid == 0 { (
        (stroke: gray + 0.3pt),
        (label-only: true, label: box(width: 2pt, height: 2pt, fill: red),
          label-side: "left", label-gap: 1, label-style: (name: "moving-pair"),
          label-path: (offset: 0.5, length: 1.4, ratio: 0.5,
            stroke: blue + 0.4pt, mark: (end: "straight"), mark-position: 1)),
      ) } else if mode in ("ordinary", "clear-ordinary") { (stroke: gray + 0.3pt) } else { (
        offset: 0.5, length: 1.4,
        stroke: if mode == "hidden" { none } else { black + 0.4pt },
        mark: if mode == "arrow" { (end: "straight") } else { none },
      ) },
      draw-after: (g, bounds) => cetz.draw.get-ctx(ctx => {
        let label = (ctx.nodes.at("moving-pair").anchors)("center")
        if mode in ("arrow", "shaft", "ordinary") {
          assert(calc.abs(label.at(0) - 4) > 0.1, message: "standalone arrows and actual text-edge contact move the shared annotation")
        } else {
          assert(calc.abs(label.at(0) - 4) < 1e-6, message: "clear ordinary and invisible carriers do not move the annotation")
        }
        ()
      }),
    )
  }
  // Automatic external labels clear their measured text boxes along the final
  // tangent. Check both flows, diagonal boxes and a curved dangling carrier.
  for incoming in (false, true) {
    for (endpoint, tangent, outward) in (
      ((4, 0), auto, (1, 0)),
      ((-4, 0), auto, (-1, 0)),
      ((3, 3), auto, (1 / calc.sqrt(2), 1 / calc.sqrt(2))),
      ((3, -3), "horizontal", (1, 0)),
    ) {
      for gap in (0, 0.45, 0.8) {
        let scene = graph.build({
          graph.node(<blob>, pos: graph.pos(x: 0, y: 0))
          graph.edge(<external>, if incoming { graph.sink(<blob>) } else { graph.source(<blob>) },
            pos: graph.pos(x: endpoint.at(0), y: endpoint.at(1)))
        })
        draw(scene, unit: 10pt, node-label: none,
          edge-label: box(width: 15pt, height: 7pt, fill: red),
          edge-label-style: (name: "external-label", rotate: 15deg),
          edge-dangling-tangent: tangent, external-label-gap: gap,
          draw-after: (g, bounds) => cetz.draw.get-ctx(ctx => {
            let anchors = ctx.nodes.at("external-label").anchors
            let nearest = calc.min(..("north-west", "north-east", "south-west", "south-east").map(name => {
              let corner = anchors(name)
              ((corner.at(0) - endpoint.at(0)) * outward.at(0)
                + (corner.at(1) - endpoint.at(1)) * outward.at(1))
            }))
            assert(calc.abs(nearest - gap) < 1e-6, message: "external labels keep the requested text-box clearance")
            ()
          }),
        )
      }
    }
  }
  // Authored positions remain exact, even with a larger automatic gap.
  draw(graph.build({
    graph.node(<manual>, pos: graph.pos(x: 0, y: 0))
    graph.edge(<label>, graph.sink(<manual>), pos: graph.pos(x: 3, y: 0),
      label-pos: (3, 2))
  }), node-label: none, edge-label: [manual], external-label-gap: 2,
    edge-label-style: (anchor: "center", name: "manual-label"),
    draw-after: (g, bounds) => cetz.draw.get-ctx(ctx => {
      let center = (ctx.nodes.at("manual-label").anchors)("center")
      assert(calc.abs(center.at(0) - 3) < 1e-6 and calc.abs(center.at(1) - 2) < 1e-6)
      ()
    }),
  )
  // Hidden edges provide neither label carriers nor stroke obstacles.
  draw(graph.build({
    graph.node(<hidden-a>, pos: graph.pos(x: 0, y: 0))
    graph.node(<hidden-b>, pos: graph.pos(x: 2, y: 0))
    graph.edge(graph.source(<hidden-a>), graph.sink(<hidden-b>))
  }), source-style: none, sink-style: none, node-label: none)
}

#curved-arrow-behavior

#cetz.canvas({
  cetz.draw.get-ctx(ctx => {
    // Preserve exact coordinates and primitive callback types when the native
    // pattern operation emits CeTZ tuples directly, including implicit origins.
    let typed-coordinate(ctx, point) = if type(point) == array and type(point.at(0)) in (int, float) {
      (point.at(0) + if type(point.at(0)) == int { 0.125 } else { 0.0 }, point.at(1))
    } else { point }
    let zero = curve.from-elements((curve.move-to((0, 0)), curve.close(), curve.line-to((0, 0))))
    let cases = (
      curve.path(), zero,
      curve.quad((0, 0), (1e-18, 2e-18), (3e-18, 0)),
      curve.line((-0.0, 0.0), (1.5, -0.5)),
      curve.cubic((0, 0), (0.5, 1), (1, -0.5), (2, 0)),
    )
    for path in cases {
      for unit in (1, -1, -0.0, 1.75, 1pt) {
        for options in (
          (pattern: "wave", samples-per-period: 4),
          (pattern: "zigzag", samples-per-period: 4),
          (pattern: "coil", samples-per-period: 4, anchor-start: false, anchor-end: false),
        ) {
          let generated = curve.pattern(path, wavelength: 2.5, amplitude: 0.04, ..options)
          for variant in range(4) {
            let local = ctx
            let style = (stroke: black + 0.3pt)
            if variant == 1 {
              local.transform = cetz.matrix.transform-shear-x(0.3)
              local.resolve-coordinate = (typed-coordinate,)
            } else if variant == 2 { local.debug = true }
            else if variant == 3 {
              local.style.line.mark = (end: ">",)
              local.style.bezier.mark = (end: ">",)
            }
            let actual = curve.pattern-to-cetz(path, unit: unit, style: style,
              wavelength: 2.5, amplitude: 0.04, ..options)
            let expected = curve.to-cetz(generated, unit: unit, ..style)
            if expected.len() == 0 { assert.eq(actual, ()) }
            else {
              let actual = actual.first()(local)
              let expected = expected.first()(local)
              assert.eq(actual.drawables, expected.drawables, message: "fused pattern drawables")
              assert.eq(cbor.encode(actual.drawables.map(d => d.segments)),
                cbor.encode(expected.drawables.map(d => d.segments)), message: "fused pattern coordinates")
              assert.eq(cbor.encode(actual.ctx.prev), cbor.encode(expected.ctx.prev), message: "fused pattern context")
            }
          }
        }
      }
    }
    // A `unit` in an edge's drawing style remains a coordinate scale, rather
    // than being accidentally forwarded as a CeTZ stroke-style property.
    let path = curve.line((0, 0), (1.5, -0.5))
    for unit in (-1, 1.75, 1pt, 9007199254740993) {
      let style = (pattern: "wave", pattern-wavelength: 2.5, pattern-amplitude: 0.04,
        pattern-samples-per-period: 4, unit: unit)
      let actual = drawing._segments-elements(curve.segments(path), style, auto, true, true).elements
      let generated = curve.pattern(drawing._segments-path(curve.segments(path)), pattern: "wave", wavelength: 2.5, amplitude: 0.04, samples-per-period: 4)
      let expected = curve.to-cetz(generated, ..drawing._draw-style(style))
      let actual = cetz.process.many(ctx, actual.flatten())
      let expected = cetz.process.many(ctx, expected)
      assert.eq(actual.drawables, expected.drawables, message: "patterned edge unit override")
      assert.eq(cbor.encode(actual.drawables.map(d => d.segments)),
        cbor.encode(expected.drawables.map(d => d.segments)), message: "patterned edge scaled coordinates")
    }
    ()
  })
})

// Native candidate footprints retain exact painted geometry.
#cetz.canvas(length: 1pt, {
  (ctx => {
    let paths = (
      curve.line((0, 0), (4, 0)),
      curve.from-cubic((start: (0, 0), control-start: (1, 2), control-end: (3, -1), end: (4, 0))),
      curve.path(curve.line((0, 0), (1, 1)), curve.line((3, 1), (5, 0))),
      curve.path(curve.move-to((0, 0)), curve.line-to((1, 1)), curve.close(),
        curve.move-to((3, 0)), curve.cubic-to((4, 1), (5, -2), (6, 0))),
      curve.path(), curve.line((2, 2), (2, 2)),
    )
    let styles = (
      (mark: (end: (symbol: "straight", length: 0.8, width: 0.4)), mark-position: 1),
      (mark: (end: (symbol: ">", length: 0.8, width: 0.4)), mark-position: "center"),
      (mark: (start: (symbol: "<", pos: 20%, offset: -5%), end: (symbol: "straight", pos: 35%, offset: 10%))),
      (mark: (end: ((symbol: ">", sep: 0.15), (symbol: "straight", length: 75%, width: 50%)))),
      (mark: (end: (symbol: "straight", reverse: true, flip: true, harpoon: true, slant: 30%)), mark-position: 0.6, mark-shift: -0.1),
      (mark: (end: (symbol: ">", fill: red, stroke: none, length: 0.5, width: 0.2))),
      (mark: (end: (symbol: ">", stroke: blue + 0.7pt), transform-shape: true)),
      (mark: (end: (symbol: "straight", length: 0.8, width: 0.4)), unit: 2),
    )
    let transforms = (
      cetz.matrix.ident(4),
      cetz.matrix.transform-scale((-1.7, 0.6, 1)),
      cetz.matrix.mul-mat(cetz.matrix.transform-translate(2, -3, 1), cetz.matrix.transform-shear-x(0.4)),
      cetz.matrix.transform-rotate-xyz(35deg, 20deg, 0deg),
      cetz.matrix.transform-scale((0, 0, 1)),
    )
    for (transform-index, transform) in transforms.enumerate() {
      for (style-index, options) in styles.enumerate() {
        let local = ctx + (transform: transform)
        local.style.bezier.stroke = red + 0.5pt
        let style = (stroke: black + 0.4pt) + options
        if transform-index != 4 {
          assert.ne(drawing._mark-carrier-footprint-spec(local, paths, style), none,
            message: "ordinary matrix case must exercise the native footprint kernel")
        }
        let actual = drawing._annotation-candidate-bounds(local, paths, style)
        let expected = drawing._annotation-path-bounds(local, paths.map(path => drawing._derived-path(path, style, auto, true, true)))
        assert.eq(cbor.encode(actual), cbor.encode(expected), message: "native candidate footprints " + repr((transform-index, style-index)))
      }
    }
    let empty-prefix = curve.path(curve.move-to((9, 9)), curve.move-to((0, 0)),
      curve.line-to((1, 0)), curve.line-to((2, 1)))
    let empty-style = (stroke: black + 0.4pt) + styles.first()
    assert.eq(drawing._mark-carrier-footprint-spec(ctx, (empty-prefix,), empty-style), none)
    assert.eq(cbor.encode(drawing._annotation-candidate-bounds(ctx, (empty-prefix,), empty-style)),
      cbor.encode(drawing._annotation-path-bounds(ctx, (drawing._derived-path(empty-prefix, empty-style, auto, true, true),))))
    (ctx: ctx)
  },)
})

#let footprint-custom-mark(entry) = {
  cetz.draw.anchor("tip", (0, 0))
  cetz.draw.anchor("base", (entry.length, 0))
  cetz.draw.line((0, 0), (entry.length, 0.3), stroke: entry.stroke)
  cetz.draw.content((entry.length / 2, 0), [X], padding: 0)
  cetz.draw.hide({ cetz.draw.line((-20, -20), (20, 20), stroke: red + 2pt) })
}
#cetz.canvas(length: 1pt, {
  (ctx => {
    let paths = (
      curve.line((0, 0), (4, 0)),
      curve.from-cubic((start: (0, 0), control-start: (1, 2), control-end: (3, -1), end: (4, 0))),
      curve.path(curve.line((0, 0), (1, 1)), curve.line((3, 1), (5, 0))),
    )
    let custom = ctx
    custom.marks.marks.insert("custom", footprint-custom-mark)
    custom.marks.mnemonics.insert("my-arrow", "custom")
    for transform in (
      cetz.matrix.ident(4),
      cetz.matrix.transform-scale((-1.7, 0.6, 1)),
      cetz.matrix.mul-mat(cetz.matrix.transform-translate(2, -3, 1), cetz.matrix.transform-shear-x(0.4)),
    ) {
      for shape in (false, true) {
        let local = custom + (transform: transform)
        let style = (stroke: black + 0.4pt,
          mark: (end: (symbol: "my-arrow", length: 0.8, width: 0.4, pos: 0.7), transform-shape: shape),
          mark-position: 1)
        assert.ne(drawing._mark-carrier-footprint-spec(local, paths, style), none)
        let actual = drawing._annotation-candidate-bounds(local, paths, style)
        let expected = drawing._annotation-path-bounds(local, paths.map(path => drawing._derived-path(path, style, auto, true, true)))
        assert.eq(cbor.encode(actual), cbor.encode(expected), message: "custom content, hidden pieces and transforms")
      }
    }
    let relative = (stroke: black + 0.4pt,
      mark: (end: (symbol: "my-arrow", length: 0.8, width: 0.4, pos: 25%)))
    assert.eq(drawing._mark-carrier-footprint-spec(custom, paths, relative), none)
    assert.eq(cbor.encode(drawing._annotation-candidate-bounds(custom, paths, relative)),
      cbor.encode(drawing._annotation-path-bounds(custom, paths.map(path => drawing._derived-path(path, relative, auto, true, true)))))
    let unevaluated = ctx
    unevaluated.marks.marks.insert("unevaluated", _ => panic("zero carriers must not evaluate custom marks"))
    let zero-style = (stroke: black + 0.4pt, mark: (end: "unevaluated"))
    let collapsed = unevaluated + (transform: cetz.matrix.transform-scale((0, 0, 1)))
    assert.eq(drawing._mark-carrier-footprint-spec(collapsed, paths, zero-style), none)
    assert.eq(cbor.encode(drawing._annotation-candidate-bounds(collapsed, paths, zero-style)),
      cbor.encode(drawing._annotation-path-bounds(collapsed, paths.map(path => drawing._derived-path(path, zero-style, auto, true, true)))))
    let zero-unit = zero-style + (unit: 0)
    let multi = (paths.last(),)
    assert.eq(drawing._mark-carrier-footprint-spec(unevaluated, multi, zero-unit), none)
    assert.eq(cbor.encode(drawing._annotation-candidate-bounds(unevaluated, multi, zero-unit)),
      cbor.encode(drawing._annotation-path-bounds(unevaluated, multi.map(path => drawing._derived-path(path, zero-unit, auto, true, true)))))
    (ctx: ctx)
  },)
})
