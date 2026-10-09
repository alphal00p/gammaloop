#import "crates/kurvst/typst/src/lib.typ" as curve

// Public Wasm integration replaces the obsolete patched-CeTZ core oracle.
// Native catalogue/contact coverage lives with the shared Kurvst engine.
#let paths = (
  curve.line((0, 0), (100, 0)),
  curve.cubic((0, 0), (0, 40), (100, -20), (100, 0)),
  curve.path(curve.line((0, 0), (50, 0)), curve.line((50, 0), (100, 20))),
  curve.line((0, 0), (0.1, 0)),
)
#let templates = (
  curve.mark.triangle(length: 9pt, width: 6pt),
  curve.mark.straight(length: 9pt, width: 6pt),
  curve.mark.stealth(length: 9pt, width: 6pt, inset: 40%),
  curve.mark.combine(curve.mark.bar(width: 6pt), 2pt, curve.mark.bar(width: 6pt)),
).map(mark => (mark: mark, "context": (units-per-pt: 1, line-thickness: 0.5)))
#for path in paths {
  for template in range(templates.len()) {
    for station in ((kind: "start"), (kind: "end"), (kind: "ratio", value: 0.5)) {
      let placements = ((template: template, carrier: 0, station: station),)
      let candidate = curve.mark.geometry(templates, (path,), placements)
      let selected = curve.mark.geometry(templates, (path,), placements, mode: "selected")
      assert.eq(candidate.marks.len(), 1)
      assert.eq(candidate.shafts, ())
      assert.eq(selected.marks.len(), 1)
      assert.eq(selected.shafts.len(), 1)
      assert.eq(candidate.marks.first().tip, selected.marks.first().tip)
      for result in (candidate, selected) {
        for point in curve.points(result.marks.first().footprint) {
          assert(point.all(value => value == value and calc.abs(value) < calc.inf),
            message: "public mark geometry remains finite on short and joined carriers")
        }
      }
    }
  }
}
#assert.eq(curve.mark.geometry((), (), ()).marks, ())
