#import "../notation.typ" as notation

#let number(source) = (kind: "number", source: source)
#let x = (kind: "variable", source: "x", symbol: (name: "x"))
#let power(base, exponent) = (kind: "power", base: base, exponent: number(exponent))
#let visual = notation.render(x)

// Standalone rational powers use roots, including negative powers and large
// exact exponents that must not pass through Typst's bounded integers.
#assert.eq(notation.render(power(x, "1/2")), math.sqrt(visual))
#assert.eq(notation.render(power(x, "2/3")), math.attach(math.root($3$.body, visual), t: $2$.body))
#assert.eq(notation.render(power(x, "-1/2")), math.frac($1$.body, math.sqrt(visual)))
#assert.eq(notation.render(power(x, "1/18446744073709551617")), math.root($18446744073709551617$.body, visual))

// In a product's denominator Symbolica retains a positive fractional exponent.
#let fraction = notation.render((kind: "product", factors: (x, power(x, "-1/2"))))
#assert.eq(fraction.num, visual)
#assert.eq(fraction.denom, math.attach(visual, t: $(1)/(2)$.body))

// Negative and already-powered bases need parentheses before another power.
#let negative = notation.render(number("-2"))
#assert.eq(notation.render(power(number("-2"), "2")), math.attach($(#negative)$.body, t: $2$.body))
#let inner = notation.render(power(x, "2"))
#assert.eq(notation.render(power(power(x, "2"), "3")), math.attach($(#inner)$.body, t: $3$.body))
