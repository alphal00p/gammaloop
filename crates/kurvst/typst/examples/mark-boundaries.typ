#import "../src/lib.typ" as kurvst
#import kurvst: mark
#import "@preview/cetz:0.5.2" as cetz
#let dark = sys.inputs.at("theme", default: "light") == "dark"
#let ink = if dark { rgb("#edf0ff") } else { rgb("#24314c") }
#set page(width: auto, height: auto, margin: 12pt, fill: if dark { rgb("#101827") } else { white })
#set text(size: 8pt, fill: ink)
#let samples = (
  ("tiny positive triangle (1e-20 pt)", mark.triangle(length: float("1e-20") * 1pt, width: float("1e-20") * 1pt, fill: false, stroke: true)),
  ("tiny positive triangle (1e-200 pt)", mark.triangle(length: float("1e-200") * 1pt, width: float("1e-200") * 1pt, fill: false, stroke: true)),
  ("near-flat stealth · grown full extent", mark.stealth(length: 10pt, width: 8pt, inset: 99.9999%, fill: false, stroke: true)),
  ("thin triangle · grown full extent", mark.triangle(length: 10pt, width: 0.00001pt, fill: false, stroke: true)),
  ("barb · full positive turn", mark.barb(width: 10pt, arc: 360deg)),
  ("barb · open neighboring turn", mark.barb(width: 10pt, arc: 359.999999deg)),
  ("barb · full negative turn", mark.barb(width: 10pt, arc: -360deg)),
  ("barb · reversed full turn", mark.barb(width: 10pt, arc: 360deg, rev: true)),
  ("full-turn composite", mark.combine(mark.barb(width: 10pt, arc: 360deg), 2pt, mark.bar(width: 6pt))),
  ("full turn on a curve", mark.barb(width: 10pt, arc: 360deg)),
  ("near-flat stealth · clipped tip view", mark.stealth(length: 10pt, width: 8pt, inset: 99.9999%, fill: false, stroke: true), true),
  ("thin triangle · clipped tip view", mark.triangle(length: 10pt, width: 0.00001pt, fill: false, stroke: true), true),
  ("reviewed thinner triangle · full extent", mark.triangle(length: 10pt, width: 0.0000001pt, fill: false, stroke: true)),
  ("reviewed deeper stealth · full extent", mark.stealth(length: 10pt, width: 8pt, inset: float("0.99999999") * 100%, fill: false, stroke: true)),
  ("reviewed thinner triangle · clipped tip", mark.triangle(length: 10pt, width: 0.0000001pt, fill: false, stroke: true), true),
  ("reviewed deeper stealth · clipped tip", mark.stealth(length: 10pt, width: 8pt, inset: float("0.99999999") * 100%, fill: false, stroke: true), true),
)
#let picture(spec, curved: false, detail: false) = {
  let carrier = if curved { kurvst.cubic((0, 0), (0, 20), (20, 20), (20, 0)) } else { kurvst.line((0, 0), (20, 0)) }
  let head = mark.geometry(
    ((mark: spec, "context": (units-per-pt: 1, line-thickness: 1)),), (carrier,),
    ((template: 0, carrier: 0, station: (kind: "end")),),
  ).marks.first()
  let boxes = head.footprint-bounds
  let left = calc.min(..boxes.map(box => box.left))
  let right-bound = calc.max(..boxes.map(box => box.right))
  let top = calc.max(..boxes.map(box => box.top))
  let bottom = calc.min(..boxes.map(box => box.bottom))
  // Use the same per-path cap/join/miter paint contract as the ordinary
  // gallery. The footprint is a coverage region, not one filled CeTZ contour.
  let zoom = if detail { 10 } else { 110 / calc.max(right-bound - left, top - bottom, 20) }
  let commands = {
    kurvst.to-cetz(head.shaft, unit: zoom, stroke: ink + zoom * 1pt)
    for part in head.paths {
      kurvst.to-cetz(part.path, unit: zoom,
        fill: if part.fill { ink } else { none },
        stroke: if part.stroke { (paint: ink, thickness: zoom * 1pt, join: part.join, cap: part.cap, miter-limit: part.miter-limit) } else { none },
      )
    }
  }
  let painting = cetz.canvas(length: 1pt, {
    if detail {
      // A hidden viewport fixes display bounds in the engine's coordinates;
      // floating changes only layout bounds, never the source geometry/paint.
      cetz.draw.hide(bounds: true, {
        cetz.draw.rect((12 * zoom, -3.5 * zoom), (23 * zoom, 3.5 * zoom), fill: ink, stroke: none)
      })
      cetz.draw.floating(commands)
    } else { commands }
  })
  let image = if detail { box(clip: true, painting) } else { painting }
  stack(spacing: 4pt, image,
    text(size: 6pt)[paint x: (#left, #right-bound); y: (#bottom, #top)],
    text(size: 6pt)[shaft end: #calc.round(head.end, digits: 6); display zoom: #zoom],
  )
}
#grid(columns: (auto, auto), gutter: 14pt, ..samples.enumerate().map(((index, sample)) => {
  stack(spacing: 6pt, [#sample.first()], picture(sample.at(1), curved: index == 9, detail: sample.at(2, default: false)))
}))
