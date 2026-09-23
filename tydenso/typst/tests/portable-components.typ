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

$ #actual $
