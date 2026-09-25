#import "../src/lib.typ" as kurvst

#let base = kurvst.line((0, 0), (4, 0))
#for (side, expected) in (((2, 1), 0.3), ((2, -1), -0.3)) {
  let layer = kurvst.layer(base, offset: 0.3, side-point: side, length: 2)
  assert.eq(layer.offset, expected)
  assert.eq(layer.path, kurvst.layer(base, offset: expected, length: 2).path)
}

// This curved carrier and its first quarter choose different local signs.
#let base = kurvst.cubic((0, 0), (1, 3), (4, -2), (5, 1))
#let L = kurvst.length(base)
#let first = kurvst.trim(base, end-outset: L * 3 / 4)
#let whole = kurvst.layer(base, offset: 0.3, side-point: (0, 1))
#let local = kurvst.layer(first, offset: 0.3, side-point: (0, 1))
#assert.eq(whole.offset, -0.3)
#assert.eq(local.offset, 0.3)
#let shared = kurvst.layer(first, offset: whole.offset)
#assert.eq(shared.offset, whole.offset)
#assert.eq(shared.path, kurvst.parallel(first, distance: -0.3).path)

#for offset in (0, 0.3, -0.3) {
  let empty = kurvst.layer(kurvst.path(), offset: offset, side-point: (0, 1))
  assert.eq(empty.offset, offset)
  assert.eq(kurvst.elements(empty), ())
}
