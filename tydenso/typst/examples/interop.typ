#import "../lib.typ" as tensors
#import "@local/tymbolica:0.1.0" as algebra
#import "@local/tymbolica-rubi:0.1.0" as calculus

#set page(width: auto, height: auto, margin: 12pt)

#let V = tensors.mink(4)
#let p = tensors.vector("p")
#let mu = tensors.slot(V, 1)
#let nu = tensors.slot(V, 2)
#let expression = tensors.mul(tensors.metric(V, mu, nu), p(nu))

Tydenso prints its own tensor notation:

$ #tensors.to-typst(expression) $

The same native Atom payload can be inspected or transformed by Tymbolica:

#let expanded = algebra.expand(expression)
#raw(algebra.canonical(expanded))

// Antisymmetry survives a trip through Tymbolica.
#let W = tensors.euc(3)
#let a = tensors.slot(W, 1)
#let b = tensors.slot(W, 2)
#let F = tensors.tensor("Finterop", antisymmetric: true)
#let carried = algebra.expand(F(a, b))
#let cancelled = tensors.add(carried, F(b, a))

// The other direction retains a namespace that is not visible in the glyph.
#let mass = algebra.symbol("m", namespace: "model")
#let roundtrip = tensors.add(mass, 0)

// Nested function heads exercise Symbolica's partial-import remapping in both
// directions, not merely leaf symbols.
#let nested = algebra.math($f(2*g(r-x))$)
#let nested-roundtrip = algebra.expand(tensors.add(nested, 0))

// Spinor representations make the same round trip.
#let B = tensors.bis(4)
#let u = tensors.vector("u")
#let spinor = u(tensors.slot(B, 1))
#let spinor-roundtrip = algebra.expand(spinor)

// Structured Typst notation is preserved too.
#let routed = p($arrow(x + y)$, V)
#let routed-roundtrip = algebra.expand(routed)

// A custom representation keeps its index palette.
#let M = tensors.representation(
  "M",
  3,
  namespace: "interop_representation",
  self-dual: true,
  indices: ($std.sym.mu$, $std.sym.nu$),
)
#let q = tensors.vector("q", namespace: "interop_representation")
#let custom = q(tensors.slot(M, 1))
#let custom-roundtrip = algebra.expand(custom)

// Rubi accepts the same Atom payload and returns something Tydenso can print.
#let x = algebra.math($x$)
#let custom-primitive = calculus.integrate(tensors.atom(custom), x)

#grid(
  columns: 2,
  gutter: 0.8em,
  [antisymmetry after Tymbolica], $ #tensors.to-typst(cancelled) $,
  [namespace after Tymbolica], raw(algebra.canonical(roundtrip, namespaces: true)),
  [nested calls after both plugins], raw(algebra.canonical(nested-roundtrip)),
  [spinor after Tymbolica], $ #tensors.to-typst(spinor-roundtrip) $,
  [display metadata after Tymbolica], $ #tensors.to-typst(routed-roundtrip) $,
  [representation metadata after Tymbolica], $ #tensors.to-typst(custom-roundtrip) $,
  [representation metadata after Rubi], $ #tensors.to-typst(custom-primitive) $,
)
