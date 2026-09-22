#import "../notation.typ" as notation

#set page(width: auto, height: auto, margin: 8pt)

#let one = (kind: "number", source: "1")
#let two = (kind: "number", source: "2")
#let x = (kind: "variable", source: "x", symbol: (name: "x"))
#let y = (kind: "variable", source: "y", symbol: (name: "y"))
#let z = (kind: "variable", source: "z", symbol: (name: "z"))
#let xy = (kind: "sum", terms: (x, y))
#let yz = (kind: "sum", terms: (y, z))
#let payload = bytes("original ordered product")
#let ordered = (
  kind: "function", symbol: (name: "spenso::bracket"),
  arguments: (xy, yz), atom: payload,
)

// Each additive factor is grouped; no function head or argument comma remains.
#let grouped = notation.render(ordered)
#let expected = math.lr($ (#notation.render(xy)) $) + h(0.12em) + math.lr($ (#notation.render(yz)) $)
#assert.eq(grouped, expected)

// Nested products preserve their supplied order without redundant parentheses.
#let zx = (kind: "function", symbol: (name: "spenso::bracket"), arguments: (z, x))
#let nested = (kind: "function", symbol: (name: "spenso::bracket"), arguments: (zx, y))
#let nested-visual = notation.render(nested)
#assert.eq(nested-visual, notation.render(z) + h(0.12em) + notation.render(x) + h(0.12em) + notation.render(y))

// A power groups the entire product; exponent one and reciprocals need no extra group.
#let squared = notation.render((kind: "power", base: ordered, exponent: two))
#assert.eq(squared, math.attach(math.lr($ (#grouped) $), t: notation.render(two)))
#assert.eq(notation.render((kind: "power", base: ordered, exponent: one)), grouped)
#let inverse = notation.render((kind: "power", base: ordered, exponent: (kind: "number", source: "-1")))
#assert.eq(inverse, math.frac(notation.render(one), grouped))

// The custom visual retains the exact whole-call payload and document overrides.
#let annotated = notation.render(ordered, annotate: (atom, visual, node) => {
  assert.eq(atom, payload)
  assert.eq(node, ordered)
  visual + metadata(atom)
})
#assert.eq(annotated, grouped + metadata(payload))
#let custom = notation.notation(calls: ("spenso::bracket": ctx => $cal(P)$))
#assert.eq(notation.render(ordered, notation: custom), $cal(P)$.body)

$ #grouped quad #nested-visual quad #squared quad #inverse $
