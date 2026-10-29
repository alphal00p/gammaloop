#import "../src/lib.typ" as kurvst
#import kurvst: mark
#import "@preview/cetz:0.5.2" as cetz

#let dark = sys.inputs.at("theme", default: "light") == "dark"
#let ink = if dark { rgb("#edf0ff") } else { rgb("#24314c") }
#set page(width: auto, height: auto, margin: 12pt, fill: if dark { rgb("#101827") } else { white })
#set text(size: 8pt, fill: ink)

#let samples = (
  ("triangle · sharp miter", mark.triangle(length: 18pt, width: 4pt, stroke: true)),
  ("straight · cap tangents", mark.straight(length: 18pt, width: 12pt)),
  ("stealth · drawn rear", mark.stealth(length: 18pt, width: 12pt, stroke: true)),
  ("stealth · deep notch growth", mark.stealth(length: 8pt, width: 6pt, inset: 90%, stroke: true)),
  ("round · inner box", mark.round(length: 18pt, width: 12pt)),
  ("tikz · quarter circles", mark.tikz(width: 18pt)),
  ("barb", mark.barb(width: 18pt)),
  ("hooks", mark.hooks(width: 18pt)),
  ("square", mark.square(length: 10pt)),
  ("diamond · rotated square", mark.diamond(length: 14.142135623730951pt)),
  ("rays · 6", mark.rays(n: 6, length: 10pt)),
  ("rays · 8", mark.rays(n: 8, length: 10pt)),
  ("bracket", mark.bracket(width: 16pt)),
  ("reversed stealth", mark.stealth(length: 18pt, width: 12pt, rev: true, stroke: true)),
  ("stroked composite", mark.combine(mark.stealth(length: 14pt, width: 10pt, stroke: true), 2pt, mark.bar(width: 16pt))),
  ("bend composite", mark.combine(mark.straight(length: 14pt, width: 10pt), 2pt, mark.circle(length: 6pt), fit: "bend")),
)
#let draw-sample(spec, curved: false) = {
  let carrier = if curved { kurvst.cubic((0, 0), (30, 28), (70, -22), (100, 0)) } else { kurvst.line((0, 0), (100, 0)) }
  let width = 1.2
  let geometry = mark.geometry(
    ((mark: spec, "context": (units-per-pt: 1, line-thickness: width)),),
    (carrier,), ((template: 0, carrier: 0, station: (kind: "end")),),
    mode: "selected",
  )
  cetz.canvas(length: 1pt, {
    kurvst.to-cetz(geometry.shafts.first().shaft, stroke: ink + width * 1pt)
    for part in geometry.marks.first().paths {
      kurvst.to-cetz(part.path,
        fill: if part.fill { ink } else { none },
        stroke: if part.stroke { (paint: ink, thickness: width * 1pt, join: part.join, cap: part.cap, miter-limit: part.miter-limit) } else { none },
      )
    }
  })
}
#grid(columns: (auto, auto), gutter: 12pt, ..samples.enumerate().map(((index, sample)) => {
  stack(spacing: 4pt, [#sample.first()], draw-sample(sample.last(), curved: calc.odd(index)))
}))
