#import "@preview/parsely:0.1.0"
#import "../lib.typ": *

#set page(width: auto, height: auto, margin: 0pt)

#let V = mink(4)
#let mu = slot(V, 1)
#let nu = slot(V, 2)
#let T = tensor("T", tags: ("physics::field-strength",))
#let annotation = parsely.parse(
  $#T(mu, nu)$,
  (
    semantic-metadata: (postfix: metadata, prec: 5),
    mul: (infix: $$, prec: 2.5, assoc: true),
  ),
).tree.slots.value

#assert.eq(annotation.protocol, "tymbolica")
#assert.eq(annotation.kind, "atom")
#assert(type(annotation.atom) == bytes)
#assert.eq(annotation.semantic.kind, "tensor")
#assert.eq(annotation.semantic.name, "T")
#assert.eq(annotation.semantic.tags, ("physics::field-strength",))
#assert.eq(annotation.semantic.arguments.len(), 2)
#assert.eq(annotation.semantic.arguments.at(0).index, 1)
#assert.eq(annotation.semantic.arguments.at(1).index, 2)
#assert.eq(annotation.semantic.arguments.at(0).representation.name, "mink")
#assert.eq(annotation.semantic.arguments.at(0).representation.dimension, 4)
