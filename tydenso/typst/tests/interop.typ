// Cross-plugin payload-contract checks for the cleaned interoperability example.
#import "../examples/interop.typ" as example
#import "../lib.typ" as tensors
#import "@local/tymbolica:0.1.0" as algebra

#assert(type(example.expanded) == bytes)
#assert(type(tensors.inspect(example.expanded)) == dictionary)
#assert.eq(tensors.to-string(example.cancelled), "0")
#assert.eq(
  algebra.canonical(example.roundtrip, namespaces: true),
  algebra.canonical(example.mass, namespaces: true),
)
#assert.eq(
  algebra.canonical(example.nested-roundtrip, namespaces: true),
  algebra.canonical(example.nested, namespaces: true),
)
#assert.eq(
  tensors.to-string(example.spinor-roundtrip),
  tensors.to-string(example.spinor),
)
#assert.eq(tensors.inspect(example.routed-roundtrip), tensors.inspect(example.routed))
#assert.eq(
  tensors.inspect(tensors.math($#tensors.to-typst(example.routed-roundtrip)$)),
  tensors.inspect(example.routed),
)
#assert.eq(tensors.inspect(example.custom-roundtrip), tensors.inspect(example.custom))
#assert(type(tensors.to-typst(example.custom-primitive)) == content)
