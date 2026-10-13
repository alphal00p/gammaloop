#import "crates/kurvst/typst/src/lib.typ" as kurvst

#set page(width: auto, height: auto, margin: 5mm)

// Closure has the same straight geometry in Rust, Typst, and CeTZ.
#let closed = kurvst.path(
  kurvst.line((0, 0), (3, 0)),
  kurvst.line((3, 0), (3, 4)),
  kurvst.close(),
)
#assert.eq(kurvst.elements(closed).last(), (kind: "close"))
#assert.eq(kurvst.segments(closed).len(), 3)
#assert.eq(kurvst.segments(closed).last().end, (0, 0))
#assert(calc.abs(kurvst.length(closed) - 12) < 1e-6)
#assert.eq(kurvst.to-cetz-data(closed).first().at(1), true)
#let outline = kurvst.outline(closed, width: 0.2)
#assert.eq(kurvst.elements(outline).filter(e => e.kind == "close").len(), 2)

// Positive point outsets are capped; negative outsets remain signed.
#assert.eq(kurvst.outset-point((0, 0), (10, 0), distance: 2), (2, 0))
#assert.eq(kurvst.outset-point((0, 0), (10, 0), distance: 8), (8, 0))
#assert.eq(kurvst.outset-point((0, 0), (10, 0), distance: 12), (10, 0))
#assert.eq(kurvst.outset-point((0, 0), (10, 0), distance: -2), (-2, 0))
#assert.eq(kurvst.outset-point((1, 2), (1, 2), distance: 8), (1, 2))

// Every policy and alias has explicit behavior, including with a sole limit.
#for (method, expected) in (
  ("min", 2),
  ("shorter", 2),
  ("max", 3),
  ("longer", 3),
  ("length", 3),
  ("fixed", 3),
  ("ratio", 2),
  ("relative", 2),
  ("none", none),
  ("full", none),
) {
  assert.eq(
    kurvst.resolve-length(10, length: 3, ratio: 0.2, method: method),
    expected,
  )
  if method in ("none", "full") {
    assert.eq(kurvst.resolve-length(10, length: 3, method: method), none)
    assert.eq(kurvst.resolve-length(10, ratio: 0.2, method: method), none)
  } else {
    assert.eq(kurvst.resolve-length(10, length: 3, method: method), 3)
    assert.eq(kurvst.resolve-length(10, ratio: 0.2, method: method), 2)
  }
  assert.eq(kurvst.resolve-length(10, method: method), none)
}
#assert.eq(kurvst.resolve-length(10, length: 0, ratio: -1), none)
#assert.eq(
  kurvst.resolve-length(10, length: 3, ratio: 0.2, method: limits => {
    assert.eq(limits, (base-length: 10, length: 3, ratio: 2))
    limits.length + limits.ratio
  }),
  5,
)
#assert.eq(
  kurvst.resolve-length(10, length: 0, ratio: -1, method: limits => {
    assert.eq(limits, (base-length: 10, length: none, ratio: none))
    none
  }),
  none,
)
#assert.eq(
  kurvst.resolve-length(0, ratio: 0.5, method: limits => {
    assert.eq(limits, (base-length: 0, length: none, ratio: none))
    none
  }),
  none,
)

// Layer defaults can be unpacked directly, without unsupported options.
#let base = kurvst.line((0, 0), (10, 0))
#assert(not ("optimize" in kurvst.layer-defaults))
#assert(
  calc.abs(kurvst.length(kurvst.layer(base, ..kurvst.layer-defaults)) - 10)
    < 1e-6,
)
#assert(
  calc.abs(kurvst.length(kurvst.parallel(base, distance: 0.2)) - 10) < 1e-6,
)
#for method in ("none", "full") {
  assert.eq(kurvst.center-outset(10, length: 3, resolve-length: method), 0)
  assert(
    calc.abs(
      kurvst.length(kurvst.layer(base, ratio: 0.2, resolve-length: method))
        - 10,
    )
      < 1e-6,
  )
}
#let empty = kurvst.from-elements(())
#for offset in (0, 1) {
  assert.eq(kurvst.elements(kurvst.layer(empty, offset: offset)), ())
}

#kurvst.to-native(closed, unit: 10pt, fill: none, stroke: black + 1pt)
