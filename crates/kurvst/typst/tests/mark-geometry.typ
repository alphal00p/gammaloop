#import "../src/lib.typ" as kurvst
#import kurvst: mark

#let ctx = (units-per-pt: 1, line-thickness: 1)
#let constructors = (
  mark.triangle,
  mark.straight,
  mark.stealth,
  mark.round,
  mark.tikz,
  mark.barb,
  mark.hooks,
  mark.bar,
  mark.bracket,
  mark.circle,
  mark.square,
  mark.diamond,
  mark.rays,
)
#let templates = constructors.map(constructor => (
  mark: constructor(),
  "context": ctx,
))
#let line = kurvst.line((0, 0), (100, 0))
#let placements = range(constructors.len()).map(index => (
  template: index,
  carrier: 0,
  station: (kind: "end"),
))
#let candidates = mark.geometry(templates, (line,), placements)
#assert.eq(candidates.marks.len(), 13)
#assert.eq(candidates.shafts, ())
#assert.eq(candidates.marks.at(0).end, 7.5)
#assert.eq(candidates.marks.at(2).end, 4.5)
#let packets = cbor(mark.geometry(templates, (line,), placements, format: "cbor"))
#let same-path(native, packet) = {
  assert.eq(type(packet), bytes)
  assert.eq(cbor.encode(cbor(packet)), cbor.encode(native))
  assert.eq(kurvst.elements(packet), kurvst.elements(native))
  assert.eq(kurvst.points(packet), kurvst.points(native))
  assert.eq(kurvst.segments(packet), kurvst.segments(native))
  assert.eq(kurvst.to-cetz-data(packet), kurvst.to-cetz-data(native))
  assert.eq(kurvst.to-native(packet, unit: 1pt), kurvst.to-native(native, unit: 1pt))
}
#for (native, packet) in candidates.marks.zip(packets.marks) {
  same-path(native.shaft, packet.shaft)
  same-path(native.shaft-outline, packet.shaft-outline)
  same-path(native.footprint, packet.footprint)
  assert.eq(cbor.encode(native.footprint-bounds), cbor.encode(packet.footprint-bounds))
  for (a, b) in native.paths.zip(packet.paths) {
    same-path(a.path, b.path)
    same-path(a.outline, b.outline)
    assert.eq(cbor.encode(a.outline-bounds), cbor.encode(b.outline-bounds))
  }
}
#for result in candidates.marks {
  assert.eq(result.tip, (100.0, 0.0))
  for point in kurvst.points(result.footprint) {
    assert(point.all(value => value == value and calc.abs(value) < calc.inf))
  }
}

// Selected endpoints share one shortened shaft, rather than repainting
// an independent candidate shaft under each head.
#let selected = mark.geometry(
  templates,
  (line,),
  (
    (template: 0, carrier: 0, station: (kind: "start"), direction: "backward"),
    (template: 0, carrier: 0, station: (kind: "end")),
    (template: 3, carrier: 0, station: (kind: "ratio", value: 0.5)),
  ),
  mode: "selected",
)
#assert.eq(selected.marks.len(), 3)
#assert.eq(selected.shafts.len(), 1)
#assert.eq(kurvst.points(selected.shafts.first().shaft).first(), (7.5, 0.0))
#assert.eq(kurvst.points(selected.shafts.first().shaft).last(), (92.5, 0.0))
#for head in selected.marks {
  assert.eq(kurvst.elements(head.shaft), ())
  assert.eq(kurvst.elements(head.shaft-outline), ())
}
#let selected-packets = cbor(mark.geometry(
  templates, (line,),
  (
    (template: 0, carrier: 0, station: (kind: "start"), direction: "backward"),
    (template: 0, carrier: 0, station: (kind: "end")),
    (template: 3, carrier: 0, station: (kind: "ratio", value: 0.5)),
  ),
  mode: "selected", format: "cbor",
))
#for (a, b) in selected.shafts.zip(selected-packets.shafts) {
  same-path(a.shaft, b.shaft)
  same-path(a.shaft-outline, b.shaft-outline)
  same-path(a.footprint, b.footprint)
}
#for (a, b) in selected.marks.zip(selected-packets.marks) {
  same-path(a.shaft, b.shaft)
  same-path(a.shaft-outline, b.shaft-outline)
  same-path(a.footprint, b.footprint)
}

// Context belongs to each template; mixed line widths still use one batch.
#let mixed = mark.geometry(
  (
    (mark: mark.triangle(), "context": ctx),
    (mark: mark.triangle(), "context": (units-per-pt: 1, line-thickness: 2)),
    (
      mark: mark.combine(mark.bar(), 2pt, mark.bar(), fit: "bend"),
      "context": ctx,
    ),
  ),
  (line,),
  (
    (template: 0, carrier: 0, station: (kind: "end")),
    (template: 1, carrier: 0, station: (kind: "end")),
    (template: 2, carrier: 0, station: (kind: "end")),
  ),
)
#assert.eq(mixed.marks.at(0).end, 7.5)
#assert.eq(mixed.marks.at(1).end, 12.0)

// Independent pre-P4 traversal oracle. Compare encoded coordinates as well as
// values so equal signed zeros cannot conceal an ordering regression.
#let old-outline-bounds(path) = {
  let boxes = ()
  let points = ()
  let flush(points) = {
    if points.len() == 0 { return () }
    let xs = points.map(point => point.at(0))
    let ys = points.map(point => point.at(1))
    ((left: calc.min(..xs), right: calc.max(..xs),
      bottom: calc.min(..ys), top: calc.max(..ys)),)
  }
  for element in kurvst.elements(path) {
    if element.kind == "move" {
      boxes += flush(points)
      points = ()
    }
    for key in ("start", "end", "control", "control-start", "control-end") {
      if key in element { points.push(element.at(key)) }
    }
  }
  boxes + flush(points)
}
#let check-bounds(result) = {
  let compare(actual, expected) = {
    assert.eq(actual, expected)
    for (box, old) in actual.zip(expected) {
      for key in ("left", "right", "bottom", "top") {
        assert.eq(cbor.encode(box.at(key)), cbor.encode(old.at(key)))
      }
    }
  }
  compare(result.footprint-bounds, old-outline-bounds(result.footprint))
  for path in result.paths {
    compare(path.outline-bounds, old-outline-bounds(path.outline))
  }
}
#for result in candidates.marks + selected.marks + mixed.marks {
  check-bounds(result)
}
#let bounds-specs = constructors.map(constructor => constructor()) + (
  mark.triangle(rev: true, stroke: black),
  mark.hooks(rev: true),
  mark.combine(mark.circle(), -3pt,
    mark.combine(mark.bar(), 25%, mark.triangle(stroke: black))),
)
#let bounds-carriers = (
  kurvst.cubic((12, -17), (12, 23), (72, 23), (72, -17)),
  kurvst.line((100, -20), (-30, 40)),
  kurvst.line((-0.0, 0.0), (0.0, -0.0)),
  kurvst.path(),
)
#for thickness in (0, 1) {
  let templates = bounds-specs.map(spec => (
    mark: spec,
    "context": (units-per-pt: 1, line-thickness: thickness),
  ))
  for mode in ("candidates", "selected") {
    for direction in ("forward", "backward") {
      let placements = ()
      for template in range(templates.len()) {
        for carrier in range(bounds-carriers.len()) {
          placements.push((
            template: template, carrier: carrier,
            station: (kind: "ratio", value: 0.5),
            direction: direction, shift: -3,
          ))
        }
      }
      let result = mark.geometry(templates, bounds-carriers, placements, mode: mode)
      for head in result.marks { check-bounds(head) }
    }
  }
}

Mark geometry passed.
