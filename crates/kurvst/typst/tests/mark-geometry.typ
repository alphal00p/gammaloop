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

Mark geometry passed.
