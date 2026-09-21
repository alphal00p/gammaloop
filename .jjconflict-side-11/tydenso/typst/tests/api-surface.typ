#import "../lib.typ": *

#set page(width: auto, height: auto, margin: 0pt)

#let V = mink(4)
#let mu = slot(V, 1)
#let nu = slot(V, 2)
#let p = vector("p")
#let expression = mul(metric(V, mu, nu), p(nu))
#let contracted = simplify-metrics(expression)
#let parsed = math($#metric(V, mu, nu) #p(nu)$)

#assert.eq(mu.kind, "slot")
#assert.eq(mu.representation.name, "mink")
#assert.eq(mu.representation.dimension, 4)
#assert.eq(inspect(simplify-metrics(parsed)), inspect(p(mu)))
#assert(type(to-typst(contracted)) == content)
#assert(type(to-string(contracted)) == str)
#assert.eq(inspect(math($#to-typst(contracted)$)), inspect(contracted))

#let namespaced = init(namespace: "representation_test")
#let W = (namespaced.representation)("W", 2, self-dual: true)
#assert.eq(W.namespace, "representation_test")
#assert.eq(
  ((namespaced.inspect)((namespaced.construct)(W))).name,
  "representation_test::W",
)

#let M = representation(
  "M",
  4,
  namespace: "palette_test",
  self-dual: true,
  indices: ($std.sym.mu$, $std.sym.nu$),
)
#let M-dual = dual-representation(M)
#assert.eq(M-dual.indices, M.indices)
#assert.eq(M-dual.index-start, 1)
#assert.eq(M-dual.name, M.name)

#let X = tensor("X", namespace: "palette_test")
#let named = X(slot(M, 1), slot(M, 2), slot(M, 3), slot(M, $rho_2$))
#assert.eq(inspect(math($#to-typst(named)$)), inspect(named))

#let R = representation("R", 4, namespace: "dual_test")
#let i = slot(R, $i$)
#let j = slot(R, $j$)
#let j-dual = slot(dual-representation(R), $j$)
#let q = vector("q", namespace: "dual_test")
#let dual-contraction = simplify-metrics(mul(metric(R, i, j-dual), q(j)))
#assert.eq(inspect(dual-contraction), inspect(q(i)))

#let exact-index = symbol("i", namespace: "exact_index_test")
#assert.eq(
  inspect(p(slot(V, exact-index))),
  inspect(p(slot(V, atom(exact-index)))),
)

#let display-label = p($arrow(x + y)$, V)
#let algebraic-argument = p(math($x + y$), V)
#assert(inspect(display-label) != inspect(algebraic-argument))
#assert.eq(inspect(math($#to-typst(display-label)$)), inspect(display-label))

#let identifier-mu = p("mu", V)
#let notation-mu = p($std.sym.mu$, V)
#assert(inspect(identifier-mu) != inspect(notation-mu))

#let color-dual = dual-representation(cof(3))
#assert(color-dual.is-dual)
#assert(slot(color-dual, 1).dual)

#let cases = (
  (rep: mink(4), labels: ("mu", "nu", "rho", "sigma"), self-dual: true, row: "top"),
  (rep: euc(4), labels: ("i", "j", "k", "l"), self-dual: true, row: "top"),
  (rep: lor(4), labels: ("mu", "nu", "rho", "sigma"), self-dual: false, row: "top"),
  (rep: bis(4), labels: ("a", "b", "c", "d"), self-dual: true, row: "bottom"),
  (rep: spf(4), labels: ("alpha", "beta", "gamma", "delta"), self-dual: false, row: "top"),
  (rep: cof(3), labels: ("i", "j", "k", "l"), self-dual: false, row: "top"),
  (rep: coad(8), labels: ("a", "b", "c", "d"), self-dual: true, row: "top"),
  (rep: cos(6), labels: ("I", "J", "K", "L"), self-dual: false, row: "top"),
)

#for case in cases {
  assert.eq(case.rep.indices, case.labels)
  assert.eq(case.rep.index-start, 1)
  assert.eq(case.rep.self-dual, case.self-dual)
  assert.eq(case.rep.index-row, case.row)

  let dual = dual-representation(case.rep)
  assert.eq(dual.indices, case.labels)
  assert.eq(dual.index-start, 1)
  assert.eq(dual.is-dual, if case.self-dual { false } else { true })
}

#let generic-mink = representation("mink", 4, self-dual: true)
#assert.eq(generic-mink.indices, ("mu", "nu", "rho", "sigma"))

#let F = tensor("F", antisymmetric: true)
#assert.eq(to-string(add(F(mu, nu), F(nu, mu))), "0")
#assert.eq(to-string(math($#F(mu, nu) + #F(nu, mu)$)), "0")
