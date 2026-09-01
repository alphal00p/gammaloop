#import "../lib.typ": *
#import "@local/tymbolica:0.1.0" as algebra

#set page(width: auto, height: auto, margin: 0pt)

#let M = mink(4)
#let B = bis(4)
#let L = lor(4)
#let mu = slot(M, 1)
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

#let atom-shape(node) = {
  if node.kind == "function" {
    node.short-name + "(" + node.arguments.map(atom-shape).join(",") + ")"
  } else if node.kind == "symbol" {
    node.short-name
  } else if node.kind == "number" {
    node.value
  } else {
    panic("unexpected Atom node in rendering test: " + repr(node))
  }
}

#let product = dot(p(1, M), q(2, M))
#let gamma-factor = gamma(mu)
#let gamma-tensor = gamma(mu, a, b)
#let open-chain = chain(
  u(1, B),
  v(2, B),
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
#let middle-representation-product = dot(p(1, M, 2), q(3, M, 4))
#let marker-factor = K(mu, cin, a, cout, nu)
#let reversed-marker-factor = K(mu, cout, a, cin, nu)
#let gamma0-factor = gamma0()
#let reversed-gamma0-factor = gamma0(cout, cin)
#let gamma5-factor = gamma5()
#let reversed-gamma5-factor = gamma5(cout, cin)
#let projector = projp(a, b)
#let Adj = coad("Na")
#let Fund = cof("Nc")
#let color-generator = color-t(
  slot(Adj, 3),
  slot(Fund, 1),
  slot(dual-representation(Fund), 2),
)

// A dualizable representation supplies a genuine lower row for the two
// alternate compact-vector layouts.
#let layout-L = lor(4)
#let layout-L-dual = dual-representation(layout-L)
#let layout-x = symbol("x", namespace: "layout_test")
#let layout-expression = A(
  slot(layout-L, 1),
  p(1, layout-L),
  slot(layout-L, 2),
  slot(layout-L, 3, dual: true),
)
#let lower-compact-expression = A(
  slot(layout-L, 1),
  q(2, layout-L-dual),
  slot(layout-L-dual, 3),
)
#let lower-only-expression = A(slot(layout-L-dual, 3))

#assert(atom-shape(inspect(product)) == "dot(p(1,mink(4)),q(2,mink(4)))")
#assert(atom-shape(inspect(gamma-factor)) == "gamma(in,out,mink(4,1))")
#assert(
  atom-shape(inspect(gamma-tensor))
    == "gamma(bis(4,1),bis(4,2),mink(4,1))",
)
#assert(
  atom-shape(inspect(open-chain))
    == "chain(u(1,bis(4)),v(2,bis(4)),gamma(in,out,mink(4,1)),gamma(in,out,p(1,mink(4))),gamma(in,out,mink(4,2)))",
)
#assert(
  atom-shape(inspect(mixed-polarity))
    == "D(mink(4,i),p(1,lor(4)),mink(4,j),q(2,dind(lor(4))))",
)

#let shown-open-chain = to-typst(open-chain)
#assert.eq(inspect(math($#shown-open-chain$)), inspect(open-chain))

#let shown-product = to-typst(product)
#let composed = math($#shown-open-chain + #shown-product$)
#assert.eq(inspect(composed), inspect(add(open-chain, product)))

#let shown-block = to-typst(open-chain, block: true)
#assert(shown-block.fields().block)
#assert.eq(inspect(math(shown-block)), inspect(open-chain))

#let labelled-momentum = p($arrow(x + y)$, M)
#assert.eq(
  inspect(math($#to-typst(labelled-momentum)$)),
  inspect(labelled-momentum),
)

#let closed-tree = inspect(closed-chain)
#assert(closed-tree.kind == "function")
#assert(closed-tree.short-name == "trace")
#assert(closed-tree.arguments.len() == 2)
#assert(closed-tree.arguments.at(0).short-name == "bis")
#assert(closed-tree.arguments.at(1).short-name == "cyclic")
#assert(closed-tree.arguments.at(1).cycle-symmetric)
#assert(closed-tree.arguments.at(1).arguments.len() == 3)

#let rendered-cases = (
  product,
  middle-representation-product,
  interleaved,
  gamma(p(1, M)),
  open-chain,
  nested-chain,
  closed-chain,
  explicit-slot-chain,
  two-mink-ports,
  mixed-polarity,
  heterogeneous-bra,
  heterogeneous-ket,
  dot(p(1, M), q(2, mink(5))),
  dot(p(1, M), q(2, B)),
  dot(p(1, L), q(2, L)),
  marker-factor,
  reversed-marker-factor,
  gamma0-factor,
  reversed-gamma0-factor,
  gamma5-factor,
  reversed-gamma5-factor,
  projector,
  color-generator,
)
#for expression in rendered-cases {
  let shown = to-typst(expression)
  assert.eq(inspect(math($#shown$)), inspect(expression))
}

#for layout in ("schoonschip", "call") {
  for expression in (
    layout-expression,
    lower-compact-expression,
    lower-only-expression,
  ) {
    let shown = algebra.to-typst(
      expression,
      notation: notation(tensor-layout: layout),
    )
    assert.eq(inspect(math($#shown$)), inspect(expression))
  }
}

#let call-without-commas = algebra.to-typst(
  layout-expression,
  notation: notation(tensor-layout: "call", commas: false),
)
#let qualified-schoonschip = algebra.to-typst(
  layout-expression,
  notation: notation(tensor-layout: "schoonschip", with-dim: true),
)
#let call-from-settings = algebra.to-typst(
  layout-expression,
  notation: notation(settings: (tensor-layout: "call")),
)
#assert.eq(inspect(math($#call-without-commas$)), inspect(layout-expression))
#assert.eq(inspect(math($#qualified-schoonschip$)), inspect(layout-expression))
#assert.eq(inspect(math($#call-from-settings$)), inspect(layout-expression))

// Turning symbol scripts off keeps ordinary tensor arguments in their own
// parentheses instead of assigning them an index row in call layout.
#let layout-ordinary = A(layout-x, slot(layout-L, 1), slot(layout-L-dual, 3))
#let call-with-ordinary = algebra.to-typst(
  layout-ordinary,
  notation: notation(tensor-layout: "call", symbol-scripts: false),
)
#assert.eq(inspect(math($#call-with-ordinary$)), inspect(layout-ordinary))

// Specialized compact-vector uses remain dots, slashes, and chain endpoints
// under both alternate generic tensor layouts.
#for expression in (product, gamma(p(1, M)), open-chain) {
  for layout in ("schoonschip", "call") {
    let shown = algebra.to-typst(expression, notation: notation(tensor-layout: layout))
    assert.eq(inspect(math($#shown$)), inspect(expression))
  }
}

#let ordinary = function("ordinary", namespace: "render_test")
#let x = symbol("x", namespace: "render_test")
#let ordinary-call = ordinary(add(x, 1))
#assert.eq(inspect(math($#to-typst(ordinary-call)$)), inspect(ordinary-call))

#let special = function(
  "weighted",
  namespace: "render_test",
  tags: ("display::operator",),
)
#let F = tensor("F", namespace: "render_test", antisymmetric: true)
#let custom-display = notation(
  heads: ("render_test::x": $xi$),
  calls: (
    "render_test::ordinary": ctx => {
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
      let visual = ctx.visual-arguments
      assert.eq(visual.len(), 2)
      assert.eq(repr(visual.at(0)), repr($std.sym.mu$))
      assert.eq(repr(visual.at(1)), repr($std.sym.nu$))
      let arguments = visual.join($,$)
      $cal(A)[#arguments]$
    },
  ),
)
#let custom-call = to-typst(ordinary-call, notation: custom-display)
#let custom-tag = to-typst(special(x), notation: custom-display)
#let custom-class = to-typst(F(mu, nu), notation: custom-display)
#assert.eq(inspect(math($#custom-call$)), inspect(ordinary-call))
#assert.eq(inspect(math($#custom-tag$)), inspect(special(x)))
#assert.eq(inspect(math($#custom-class$)), inspect(F(mu, nu)))
