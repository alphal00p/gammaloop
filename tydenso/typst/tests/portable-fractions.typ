#import "../notation.typ" as notation

#set page(width: auto, height: auto, margin: 8pt)

#let one = (kind: "number", source: "1")
#let two = (kind: "number", source: "2")
#let x = (kind: "variable", source: "x", symbol: (name: "x"))
#let y = (kind: "variable", source: "y", symbol: (name: "y"))
#let inverse-x = (kind: "power", base: x, exponent: (kind: "number", source: "-1"))
#let inverse-y = (kind: "power", base: y, exponent: (kind: "number", source: "-1"))

// Independent reciprocals share one fraction bar, including rational coefficients.
#let collected = notation.render((kind: "product", factors: (inverse-x, inverse-y)))
#assert.eq(collected.func(), math.frac)
#assert.eq(collected.num, notation.render(one))
#assert.eq(collected.denom, notation.render((kind: "product", factors: (x, y))))
#let scaled = notation.render((kind: "product", factors: ((kind: "number", source: "1/8"), inverse-x, inverse-y)))
#assert.eq(scaled.func(), math.frac)
#assert.eq(scaled.num, notation.render(one))
#assert.eq(scaled.denom, notation.render((kind: "product", factors: ((kind: "number", source: "8"), x, y))))
#let mixed = notation.render((kind: "product", factors: (x, inverse-y)))
#assert.eq(mixed.num, notation.render(x))
#assert.eq(mixed.denom, notation.render(y))

// Standalone inverses, inverse powers, and sums retain their mathematical grouping.
#let squared = notation.render((kind: "power", base: x, exponent: (kind: "number", source: "-2")))
#assert.eq(squared.denom, notation.render((kind: "power", base: x, exponent: two)))
#let sum = (kind: "sum", terms: (x, y))
#let inverse-sum = notation.render((kind: "power", base: sum, exponent: (kind: "number", source: "-1")))
#assert.eq(inverse-sum.denom, notation.render(sum))
#let scaled-sum = notation.render((kind: "product", factors: ((kind: "number", source: "1/8"), (kind: "power", base: sum, exponent: (kind: "number", source: "-1")))))
#assert.eq(scaled-sum.denom, notation.render((kind: "product", factors: ((kind: "number", source: "8"), sum))))
#assert(notation.render((kind: "product", factors: (x, y))).func() != math.frac)

// A leading minus in a complex exponent is not a reciprocal of its unsigned text.
#let complex-power = (kind: "power", base: x, exponent: (kind: "number", source: "-1+2𝑖"))
#assert(notation.render(complex-power).func() != math.frac)
#assert.eq(notation.render((kind: "product", factors: (complex-power, y))).func(), notation.render((kind: "product", factors: (x, y))).func())

// Only the ordinary endpoint slot is implicit. Higher-spin and dummy slots
// stay distinct, and every label uses subscripts.
#for label in ("s", "t", "e", "v") {
  for slot in ("0", "1", "2") {
    let node = (
      kind: "function",
      symbol: (name: "indices::" + label, tags: ("spenso::index", "spenso::index-label:" + label)),
      arguments: ((kind: "number", source: "7"), (kind: "number", source: slot)),
    )
    let rendered = notation.render(node)
    let expected = if label in ("s", "t") and slot == "1" { [7] } else { [7] + [.] + text(slot) }
    assert.eq(rendered, math.attach(eval(label, mode: "math").body, b: text(size: 0.75em, expected)))
  }
}

$ #collected quad #scaled quad #mixed quad #squared quad #inverse-sum $
