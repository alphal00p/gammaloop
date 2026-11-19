#import "@preview/cetz:0.5.1" as cetz

// Exact geometry oracle retained from the former Typst contact/trim owner.
// Native results use f64 coordinates; integer/float CBOR representation is not
// part of this internal API. Numerical equality below has no tolerance.
#let reference(carrier) = {
  let total = carrier.total
  let cuts = (none, none)
  let marks = ()
  for entry in carrier.marks {
    let (tip, back, length, at, side) = (entry.tip, entry.back, entry.length, entry.at, entry.side)
    let axis = cetz.vector.sub(back, tip)
    let compression = if length == 0 { 1 } else { calc.min(1, total / length) }
    let direction = if axis.first() < 0 { -1 } else { 1 }
    let stations = ()
    for point in (tip, back) {
      stations.push(at + if length == 0 { 0 } else {
        cetz.vector.dot(point, axis) / length * direction * compression
      })
    }
    let low = calc.min(..stations)
    let span = calc.abs(stations.last() - stations.first())
    let inward = calc.clamp(low, 0, calc.max(0, total - span)) - low
    stations = (
      calc.clamp(stations.first() + inward, 0, total),
      calc.clamp(stations.last() + inward, 0, total),
    )
    let contacts = ()
    for station in stations {
      contacts.push(cetz.path-util.point-at(carrier.segments, if side == 1 { total - station } else { station }))
    }
    let chord = cetz.vector.sub(contacts.last().point, contacts.first().point)
    let angle = if cetz.vector.len(chord) == 0 { 0deg } else { calc.atan2(..chord.slice(0, 2)) }
    let reference-angle = if length == 0 { 0deg } else { calc.atan2(..axis.slice(0, 2)) }
    let alignment = cetz.matrix.mul-mat(
      cetz.matrix.transform-translate(..contacts.first().point),
      cetz.matrix.transform-rotate-z(angle),
      cetz.matrix.transform-scale((if length == 0 { 1 } else { cetz.vector.len(chord) / length }, 1, 1)),
      cetz.matrix.transform-rotate-z(-reference-angle),
      cetz.matrix.transform-translate(..cetz.vector.scale(tip, -1)),
    )
    marks.push((alignment: alignment))
    if entry.shorten {
      let contact = if entry.straight { 0 } else { if stations.first() > stations.last() { 0 } else { 1 } }
      cuts.at(side) = contacts.at(contact)
    }
  }
  let painted = carrier.segments
  if carrier.paint and cuts != (none, none) {
    let (origin, closed, commands) = carrier.segments.first()
    let bezier = cetz.path-util.bezier
    let bounds = cuts.enumerate().map(((side, cut)) => {
      if cut == none { return side * commands.len() }
      let (kind, ..args) = commands.at(cut.segment-index)
      let t = if kind == "l" {
        cut.distance / cetz.vector.dist(cut.previous-point, args.last())
      } else {
        bezier.cubic-t-for-distance(cut.previous-point, args.last(), ..args.slice(0, 2), cut.distance, samples: cetz.path-util.number-of-samples(auto))
      }
      cut.segment-index + calc.clamp(t, 0, 1)
    })
    painted = ()
    if bounds.first() < bounds.last() {
      let first = calc.floor(bounds.first())
      let last = calc.ceil(bounds.last()) - 1
      if first > 0 { origin = commands.at(first - 1).last() }
      commands = commands.slice(first, last + 1)
      let end = bounds.last() - last
      let start = (bounds.first() - first) / if first == last { end } else { 1 }
      for (side, t) in (end, start).enumerate() {
        let index = if side == 0 { commands.len() - 1 } else { 0 }
        let previous = if index == 0 { origin } else { commands.at(index - 1).last() }
        let (kind, ..args) = commands.at(index)
        if kind == "c" {
          let (s, e, c1, c2) = bezier.split(previous, args.last(), ..args.slice(0, 2), t).at(side)
          if side == 1 { origin = s }
          args = (c1, c2, e)
        } else if side == 0 {
          args = (cetz.vector.lerp(previous, args.last(), t),)
        } else { origin = cetz.vector.lerp(previous, args.last(), t) }
        commands.at(index) = (kind, ..args)
      }
      painted = ((origin, closed, commands),)
    }
  }
  (segments: painted, marks: marks)
}

#let mark(tip: (0, 0, 0), back: (0.4, 0, 0), at: 0, side: 0, shorten: true, straight: false) = (
  tip: tip, back: back, length: cetz.vector.len(cetz.vector.sub(back, tip)), at: at,
  side: side, shorten: shorten, straight: straight,
)
#let paths = (
  ((0,0,0), false, (("l", (3,0,0)),)),
  ((-0.0,0.0,0.0), false, (("l", (-3.0,1.7,0.0)),)),
  ((1,2,0), false, (("l", (1,-2,0)),)),
  ((0,0,0), false, (("c", (0,2,0), (3,-1,0), (3,0,0)),)),
  ((0,0,0), false, (("c", (0,0,0), (2,0,0), (2,0,0)),)),
  ((0,0,1), false, (("c", (1,2,3), (2,-1,-1), (3,0,2)),)),
  ((3,0,0), false, (("c", (3,-1,0), (0,2,0), (0,0,0)),)),
  ((0,0,0), false, (("l", (1,0,0)), ("c", (1,2,0), (2,-1,0), (3,0,0)), ("l", (4,1,0)))),
  ((0,0,0), false, (("c", (0,1,0), (1,1,0), (1,0,0)), ("c", (2,-1,0), (3,-1,0), (3,0,0)))),
  ((0,0,0), false, (("l", (0.01,0.02,0.03)),)),
  ((0,0,0), true, (("l", (2,0,0)), ("l", (0,2,0)), ("l", (0,0,0)))),
  ((0,0), false, (("l", (3,1,0)),)),
  ((0,0,0), false, (("c", (0,2), (3,-1,0), (3,0)),)),
)
#let carriers = ()
#let names = ()
#for (path-index, path) in paths.enumerate() {
  let segments = (path,)
  let total = cetz.path-util.length(segments)
  let cases = (
    ("empty", ()),
    ("start", (mark(),)),
    ("end", (mark(side: 1),)),
    ("start-straight", (mark(straight: true),)),
    ("end-straight", (mark(side: 1, straight: true),)),
    ("reversed", (mark(tip: (0.4,0,0), back: (0,0,0)), mark(tip: (0,0,0), back: (-0.4,0,0), side: 1))),
    ("offset-axis", (mark(tip: (-0.3,0.1,0), back: (0.4,-0.1,0), at: 0.2),)),
    ("negative-station", (mark(at: -0.7),)),
    ("beyond-end", (mark(at: total + 0.2),)),
    ("both-cuts", (mark(), mark(side: 1))),
    ("multiple-last-wins", (mark(), mark(at: 0.6), mark(side: 1), mark(side: 1, at: 0.8))),
    ("overlapping-cuts", (mark(back: (total,0,0)), mark(back: (total,0,0), side: 1))),
    ("zero-mark", (mark(back: (0,0,0), at: 0.2), mark(tip: (0.3,0.2,0.1), back: (0.3,0.2,0.1), side: 1))),
    ("three-dimensional-mark", (mark(tip: (0.1,-0.2,0.3), back: (0.6,0.2,-0.1)), mark(tip: (0,0,0), back: (0,0,0.4), side: 1))),
    ("interior-overlay", (mark(at: total * 0.4, shorten: false), mark(side: 1, at: total * 0.3, shorten: false))),
    ("compressed", (mark(back: (total * 3,0,0)),)),
    ("exact-segment-boundary", (mark(back: (0,0,0), at: 1),)),
  )
  for (name, marks) in cases {
    for paint in (false, true) {
      names.push(str(path-index) + "/" + name + "/paint=" + repr(paint))
      carriers.push((segments: segments, total: total, paint: paint, marks: marks))
    }
  }
}
// Empty/no-mark inputs must not need contact sampling or nonzero path lengths.
#for segments in ((), (((0, 0, 0), false, ()),), (((0, 0, 0), false, (("l", (0, 0, 0)),)),)) {
  names.push("empty/no-mark")
  carriers.push((segments: segments, total: 0, paint: true, marks: ()))
}
#let kernel = cetz.path-util.bezier.cetz-core
#assert.eq(cbor(kernel.marker_geometry_func(cbor.encode((carriers: ())))), (carriers: ()))
#let actual = cbor(kernel.marker_geometry_func(cbor.encode((carriers: carriers))))
#assert.eq(actual.carriers.len(), carriers.len())
#for (name, carrier, result) in names.zip(carriers, actual.carriers) {
  assert.eq(result, reference(carrier), message: "marker geometry: " + name)
}
