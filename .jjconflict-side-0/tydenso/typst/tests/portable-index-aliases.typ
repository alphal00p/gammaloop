#import "../notation.typ" as notation

#set page(width: auto, height: auto, margin: 8pt)

#let number(value) = (kind: "number", text: str(value), source: str(value))
#let call(name, arguments, tags: (), atom: none) = (
  kind: "function", source: repr(name),
  symbol: (name: name, tags: tags), arguments: arguments, atom: atom,
)
#let index = call("gammalooprs::hedge", (number(4), number(1)), tags: ("spenso::index",))
#let slot(index) = call("spenso::mink", (number(4), index), tags: ("spenso::representation",))
#let tensor(index, atom: none) = call("T", (slot(index),), tags: ("spenso::tensor",), atom: atom)
#let tree(root) = (
  kind: "atom-render-tree",
  root: root,
  attachments: ((
    schema: "spenso.representation", version: 1, identity: bytes("spenso::mink"),
    data: cbor.encode(("inline-metric", ("cyclic", 1, (("symbol", "mu"), ("symbol", "nu"))), "top")),
  ),),
)
#let alphabet = notation.notation(settings: (
  index-aliases: (("spenso::mink", "gammalooprs::hedge", (4, 1), 1),),
))
#let actual = notation.render(tree(tensor(index)), notation: alphabet)
#assert.eq(actual, notation.render(tree(tensor(number(1)))))

// A same-numbered edge has a distinct identity and must not borrow the alias.
#let edge = call("gammalooprs::edge", (number(4), number(1)), tags: ("spenso::index",))
#assert.ne(notation.render(tree(tensor(edge)), notation: alphabet), actual)

// Graph labels can provide safe structured math directly.
#let linked = notation.notation(settings: (
  index-aliases: (("spenso::mink", "gammalooprs::hedge", (4, 1), $mu_(upright("h4"))$),),
))
#let linked-visual = notation.render(tree(tensor(index)), notation: linked)
#assert.ne(linked-visual, actual)

// Visual aliases must leave the complete tensor's exact Atom payload intact.
#let payload = bytes("original tensor with hedge(4,1)")
#let annotated = notation.render(
  tree(tensor(index, atom: payload)), notation: alphabet,
  annotate: (atom, visual, node) => {
    assert.eq(atom, payload)
    assert.eq(node.arguments.first().arguments.at(1), index)
    visual
  },
)
#assert.eq(annotated, actual)

$ #actual quad #linked-visual $
