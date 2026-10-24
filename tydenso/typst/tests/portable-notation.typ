#import "../notation.typ" as notation
#set page(width: auto, height: auto, margin: 8pt)
#let name = "portable::tensor"
#let tree = (
  kind: "atom-render-tree",
  root: (kind: "variable", source: "tensor", symbol: (name: name, tags: ())),
  attachments: ((schema: "spenso.notation", version: 1, identity: bytes(name),
    data: cbor.encode((head: ("symbol", "eta"), "index-label": ("symbol", "mu"))),
  ),),
)
#let actual = notation.render(tree)
#assert.eq(actual, notation.render(tree, notation: notation.notation(heads: ((name): $eta$))))
#let index-tree = (
  kind: "atom-render-tree",
  root: (kind: "function", source: "tensor(2)",
    symbol: (name: name, tags: ("spenso::index",)),
    arguments: ((kind: "number", source: "2", text: "2"),)),
  attachments: tree.attachments,
)
#assert.eq(notation.render(index-tree), notation.render(index-tree, notation: notation.notation(settings: (
  index-labels: ((name): (label: "tensor", visual: $mu$)),
))))
$ #actual quad #notation.render(index-tree) $

// Exercise the real plugin import boundary as well as the pure render tree.
#import "../lib.typ" as tensors
#let engine = tensors.init()
#let imported-name = "portable_import::Head"
#let base = (engine.atom)((engine.symbol)(imported-name))
// Append one attachment to the version-2 envelope (big-endian lengths).
#assert.eq(base.slice(0, 8), bytes("SYMATOM\u{0}"))
#assert.eq(base.slice(8, 14), bytes((0, 2, 0, 0, 0, 0)))
#let be16(value) = bytes((calc.rem(calc.quo(value, 256), 256), calc.rem(value, 256)))
#let be32(value) = bytes((0, 0)) + be16(value)
#let schema = bytes("spenso.notation")
#let identity = bytes(imported-name)
#let data = cbor.encode((head: ("symbol", "eta"),))
#let payload = base.slice(0, 12) + be16(1) + base.slice(14) + (
  be16(schema.len()) + be32(1) + be32(identity.len()) + be32(data.len())
  + schema + identity + data
)
#let imported-tree = cbor(engine.plugin.render_tree(payload))
#assert.eq(notation.render(imported-tree), notation.render(imported-tree, notation: notation.notation(
  heads: ((imported-name): $eta$),
)))
$ #notation.render(imported-tree) $

// A sequence containing only s still has the native endpoint label s.
#let endpoint-tree = (
  ..index-tree,
  root: (..index-tree.root, arguments: (
    (kind: "number", source: "4", text: "4"),
    (kind: "number", source: "1", text: "1"),
  )),
  attachments: ((schema: "spenso.notation", version: 1, identity: bytes(name),
    data: cbor.encode((head: ("symbol", "eta"), "index-label": ("sequence", (("symbol", "s"),)))),
  ),),
)
#assert.eq(notation.render(endpoint-tree), notation.render(endpoint-tree, notation: notation.notation(settings: (
  index-labels: ((name): (label: "s", visual: math.italic("s"))),
))))
