#import "../lib.typ": *

#set page(width: auto, height: auto, margin: 8pt)

#let M = mink(4)
#let B = bis(4)
#let mu = slot(M, $u$)
#let nu = slot(M, 2)
#let a = slot(B, 1)
#let b = slot(B, 2)

#let p = vector("p")
#let q = vector("q")
#let u = vector("u")
#let v = vector("v")
#let chi = vector("chi")
#let r = vector("r")
#let s = vector("s")
#let T = tensor("T")
#let A = tensor("A")
#let H = tensor("H")
#let C = tensor("C")
#let D = tensor("D")
#let K = tensor("K")
#let cin = symbol("in")
#let cout = symbol("out")
#let color-t = function("t")
#let ordinary-f = function("ordinary", namespace: "example")
#let ordinary-x = symbol("x", namespace: "example")

#let product = dot(p(1, M), q(2, M))
#let gamma-tensor = gamma(mu, a, b)
#let open-chain = chain(
  v(p(2)),
  u(p(1)),
  gamma(mu),
  gamma(p(1, M)),
  gamma(nu),
)
#let explicit-slot-chain = chain(
  slot(B, $std.sym.mu$),
  slot(B, "blou"),
  gamma(mu),
  gamma(p(1, M)),
  gamma(nu),
)
#let closed-chain = trace(
  B,
  cyclic(gamma(mu), gamma(p(1, M)), gamma(nu)),
)
#let interleaved = T(mu, a, nu, b)
#let two-mink-ports = A(mu, p(1, M), nu, q(2, M), slot(M, 3))
#let heterogeneous-bra = H(mu, p(1, M), a, chi(2, B))
#let L = lor(4)
#let compact-layout = A(
  slot(L, 1),
  p(1, L),
  slot(L, 2),
  slot(L, 3, dual: true),
)
#let heterogeneous-ket = C(
  slot(M, "i"), r(1, L), slot(M, "j"), s(2, L),
)
#let mixed-polarity = D(
  slot(M, "i"), p(1, L), slot(M, "j"), q(2, dual-representation(L)),
)
#let nested-chain = chain(
  u(1, B),
  v(2, B),
  gamma(slot(M, $alpha$)),
  A("in", "out", mu, p(3, M), nu),
  gamma(slot(M, $beta$)),
)
#let middle-representation-product = p(1, q(2, M))
#let reversed-marker-factor = K(mu, cout, a, cin, nu)
#let gamma0-factor = gamma0()
#let gamma5-factor = gamma5()
#let projector = projp(a, b)
#let Adj = coad("Na")
#let Fund = cof("Nc")
#let color-generator = color-t(
  slot(Adj, 3),
  slot(Fund, 1),
  slot(dual-representation(Fund), 2),
)

#let ordinary-call = ordinary-f(add(ordinary-x, 1))

#let labelled-momentum = p($arrow(x + y)$, M)

// Exact heads, complete calls, tags, and tensor classes can all be styled in Typst.
#let special = function(
  "weighted",
  namespace: "example",
  tags: ("display::operator",),
)
#let F = tensor("F", namespace: "example", antisymmetric: true)
#let custom-display = notation(
  heads: ("example::x": $xi$),
  calls: (
    "example::ordinary": ctx => {
      let (argument,) = ctx.visual-arguments
      $cal(F)[#argument]$
    },
  ),
  tags: (
    "display::operator": ctx => {
      let (argument,) = ctx.visual-arguments
      $bold(W)(#argument)$
    },
  ),
  classes: (
    "antisymmetric": ctx => {
      let arguments = ctx.visual-arguments.join($,$)
      $cal(A)[#arguments]$
    },
  ),
)
#let custom-call = to-typst(ordinary-call, notation: custom-display)
#let custom-tag = to-typst(special(ordinary-x), notation: custom-display)
#let custom-class = to-typst(F(mu, nu), notation: custom-display)

#grid(
  columns: 2,
  gutter: 1em,
  [explicit rows], $ #to-typst(interleaved) $,
  [gamma tensor], $ #to-typst(gamma-tensor) $,
  [dot product], $ #to-typst(product) $,
  [open chain], $ #to-typst(open-chain) $,
  [explicit-slot chain], $ #to-typst(explicit-slot-chain) $,
  [closed trace], $ #to-typst(closed-chain) $,
  [two Mink ports], $ #to-typst(two-mink-ports) $,
  [Schoonschip scripts], $ #to-typst(
    compact-layout,
    notation: notation(tensor-layout: "schoonschip"),
  ) $,
  [row-grouped call], $ #to-typst(
    compact-layout,
    notation: notation(tensor-layout: "call"),
  ) $,
  [nested factor], $ #to-typst(nested-chain) $,
  [heterogeneous bra], $ #to-typst(heterogeneous-bra) $,
  [heterogeneous ket], $ #to-typst(heterogeneous-ket) $,
  [mixed polarity], $ #to-typst(mixed-polarity) $,
  [typed ports], $ #to-typst(
    heterogeneous-bra,
    notation: notation(with-dim: true),
  ) $,
  [arbitrary momentum label], $ #to-typst(labelled-momentum) $,
  [middle rep], $ #to-typst(middle-representation-product) $,
  [reversed factor], $ #to-typst(reversed-marker-factor) $,
  [Idenso heads], $ #to-typst(projector) #to-typst(gamma0-factor) #to-typst(gamma5-factor) $,
  [color rows], $ #to-typst(color-generator) $,
  [custom exact call], $ #custom-call $,
  [custom tag], $ #custom-tag $,
  [custom class], $ #custom-class $,
)
