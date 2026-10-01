#import "../notation.typ" as notation

#let number(source) = (kind: "number", source: source)
#let x = (kind: "variable", source: "x", symbol: (name: "x"))
#let visual-x = notation.render(x)
#let flatten-sequences(visual) = {
  if repr(visual.func()) == "sequence" {
    visual.children.map(flatten-sequences).join()
  } else { visual }
}

// Complex coefficients are single numeric nodes even when their visible
// notation is a sum. Preserve that sum as one factor.
#for source in ("1+𝑖", "-1+2𝑖", "1-𝑖", "1/2+𝑖/3", "1+1𝑖", "1.2 dot 10^(-3)+2.4 dot 10^(-5)𝑖") {
  let coefficient = number(source)
  let visual = eval(source, mode: "math").body
  let product = (kind: "product", factors: (coefficient, x))
  assert.eq(notation.render(product), $(#visual) #visual-x$.body)

  // A lone numerator or denominator already has the fraction bar as grouping.
  let inverse = (kind: "power", base: x, exponent: number("-1"))
  assert.eq(notation.render((kind: "product", factors: (coefficient, inverse))), math.frac(visual, visual-x))
}

// Pure imaginary factors need no product parentheses, including signed values
// and scientific notation whose exponent contains a minus sign.
#for source in ("𝑖", "-𝑖", "2𝑖", "-2𝑖", "𝑖/2", "1.2 dot 10^(-3)𝑖") {
  let visual = eval(source, mode: "math").body
  assert.eq(
    flatten-sequences(notation.render((kind: "product", factors: (number(source), x)))),
    flatten-sequences($#visual #visual-x$.body),
  )
}

// Both mixed and pure imaginary bases must be grouped before exponentiation.
#for source in ("1+𝑖", "-1+2𝑖", "𝑖", "-2𝑖", "1/2+𝑖/3") {
  let visual = eval(source, mode: "math").body
  assert.eq(notation.render((kind: "power", base: number(source), exponent: x)), math.attach($(#visual)$.body, t: visual-x))
}
