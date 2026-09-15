#import "../lib.typ": *

#set page(width: auto, height: auto, margin: 12pt)

#let V = mink(4)
#let mu-slot = slot(V, 1)
#let nu-slot = slot(V, 2)
#let p = vector("p")

#let expression = mul(metric(V, mu-slot, nu-slot), p(nu-slot))
#let contracted = simplify-metrics(expression)
#let W = lor(4)
#let T = tensor("T")
#let mixed = T(slot(W, 1), slot(W, 2, dual: true))

$
  #to-typst(expression)
  quad arrow.r.long
  #to-typst(contracted)
$

#let M = representation(
  "M",
  4,
  namespace: "palette_example",
  self-dual: true,
  indices: ($mu$, $nu$),
)
#let X = tensor("X", namespace: "palette_example")
#let named-indices = X(slot(M, 1), slot(M, 2), slot(M, 3), slot(M, $rho_2$))
#let shown-named-indices = to-typst(named-indices)

#let R = representation("R", 4, namespace: "dual_example")
#let R-dual = dual-representation(R)
#let i-name = $i$
#let j-name = $j$
#let i = slot(R, i-name)
#let j = slot(R, j-name)
#let j-dual = slot(R-dual, j-name)
#let q = vector("q", namespace: "dual_example")
#let custom-dual-contraction = simplify-metrics(mul(metric(R, i, j-dual), q(j)))

// Raw Typst math is one portable display label in tensor and vector calls.
// Wrap the same input in `math` when it should be Symbolica algebra instead.
#let display-label = p($arrow(x + y)$, V)
#let algebraic-argument = p(math($x + y$), V)

// A string is a semantic identifier; math content owns its written notation.
#let identifier-mu = p("mu", V)
#let notation-mu = p($mu$, V)

#grid(
  columns: 2,
  gutter: 0.8em,
  [mixed variance], $ #to-typst(mixed) $,
  [custom metric contraction], $ #to-typst(custom-dual-contraction) $,
  [automatic and manual indices], $ #shown-named-indices $,
  [indices with representation data], $ #to-typst(
    named-indices,
    notation: notation(with-dim: true),
  ) $,
  [arbitrary math as a label], $ #to-typst(display-label) $,
  [the same input as algebra], $ #to-typst(algebraic-argument) $,
  [the identifier `"mu"`], $ #to-typst(identifier-mu) $,
  [the written symbol $mu$], $ #to-typst(notation-mu) $,
)
