#import "../notation.typ" as notation

#set page(width: auto, height: auto, margin: 8pt)

#let number(value) = (kind: "number", text: str(value), source: str(value))
#let call(name, arguments, tags: (), atom: none) = (
  kind: "function", source: repr(name),
  symbol: (name: name, tags: tags), arguments: arguments, atom: atom,
)
#let coordinates = call("spenso::cind", (number(0), number(12)))
#let component(atom: none) = call("A", (coordinates,), tags: ("spenso::tensor",), atom: atom)
#let tree(root) = (kind: "atom-render-tree", root: root, attachments: ())
#let actual = notation.render(tree(component()))
#let expected = notation.render(tree(component()), notation: notation.notation(
  calls: ("A": ctx => math.attach([A], t: ([0], [12]).join($,$))),
))
#assert.eq(actual, expected)

// Components are coordinates, independently of abstract-index display layouts.
#for layout in ("ports", "schoonschip", "call") {
  assert.eq(actual, notation.render(tree(component()), notation: notation.notation(tensor-layout: layout)))
}

// Ordinary parameters stay in parentheses before the coordinate superscript.
#let x = (kind: "variable", symbol: (name: "x",), source: "x")
#let parameterized = call("A", (x, number(7), coordinates), tags: ("spenso::tensor",))
#let parameters = (notation.render(x), [7]).join($,$)
#let parameterized-expected = math.attach(
  ([A], math.lr($ (#parameters) $)).join(),
  t: ([0], [12]).join($,$),
)
#for layout in ("ports", "schoonschip", "call") {
  for scripts in (false, true) {
    assert.eq(parameterized-expected, notation.render(tree(parameterized), notation: notation.notation(
      tensor-layout: layout, symbol-scripts: scripts,
    )))
  }
}

// Array coordinates share parameter placement and preserve powers without extra grouping.
#let array-notation = notation.notation(component-style: "array")
#let array-bare = notation.render(tree(component()), notation: array-notation)
#let array-expected = ([A], math.lr($ [#([0], [12]).join($,$)] $)).join()
#assert.eq(array-bare, array-expected)
#let array-parameterized = notation.render(tree(parameterized), notation: array-notation)
#assert.eq(array-parameterized,
  ([A], math.lr($ (#parameters) $), math.lr($ [#([0], [12]).join($,$)] $)).join())
#assert.eq(
  notation.render(tree((kind: "power", base: parameterized, exponent: number(2))), notation: array-notation),
  math.attach(array-parameterized, t: notation.render(number(2))),
)
#let scalar = call("A", (call("spenso::cind", ()),), tags: ("spenso::tensor",))
#assert.eq(notation.render(tree(scalar), notation: array-notation), [A])

// A power must not merge with the component superscript.
#let squared = notation.render(tree((kind: "power", base: component(), exponent: number(2))))
#assert.eq(squared, math.attach(math.lr($ (#actual) $), t: notation.render(number(2))))

// Printing preserves the original component Atom and its coordinate payload.
#let payload = bytes("A(spenso::cind(0,12))")
#let annotated = notation.render(tree(component(atom: payload)), annotate: (atom, visual, node) => {
  assert.eq(atom, payload)
  assert.eq(node.arguments.first(), coordinates)
  visual
})
#assert.eq(actual, annotated)

// An unrelated function named cind must not acquire coordinate semantics.
#let unrelated = call("A", (call("other::cind", (number(0), number(12))),), tags: ("spenso::tensor",))
#assert.ne(actual, notation.render(tree(unrelated)))

$ #actual quad #parameterized-expected quad #array-parameterized $
