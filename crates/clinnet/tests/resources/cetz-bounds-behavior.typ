#import "@preview/cetz:0.5.2" as cetz

// Stock CeTZ owns public path bounds, including zero-length segments,
// repeated roots and varying coordinate sizes.
#let paths = (
  (((-0.0, 0.0, 0.0), false, ()),),
  (((0, 0), true, (("l", (0, 0)), ("l", (2, -3)), ("l", (0, 0)))),),
  (((1, 2, 3), false, (("c", (1, 2, 3), (1, 2, 3), (1, 2, 3)),)),),
  (((0, 0, 0), false, (("c", (1, 3, -2), (2, -4, 5), (3, 0, 1)),)),),
  (((0, 0), false, (
    ("c", (1, 1e-10), (2, -1e-10), (3, 0)),
    ("l", (4, 1, -2)),
    ("c", (3, 1, -3), (2, -1, 4), (0, 0, 0)),
  )),),
  (
    ((-2, 3, 0, 7), false, (("l", (1, -4, 2, 9)),)),
    ((5, 6, -1), true, (("c", (2, 9, 3), (-1, 0, -5), (5, 6, -1)),)),
  ),
)

#let aabb = cetz.path-util.bezier.aabb.aabb
#for path in paths {
  let expected = aabb(cetz.path-util.bounds(path))
  cetz.canvas({
    cetz.draw.get-ctx(ctx => {
      let actual = cetz.process.element(ctx, ctx => (
        ctx: ctx, drawables: (cetz.drawable.path(path),),
      )).bounds
      assert.eq(cbor.encode(actual), cbor.encode(expected))
      ()
    })
  })
}

// Aggregating per-path corners must retain the first occurrence of equal
// extrema, including signed zeros, as well as the order of mixed path bounds.
#for subset in (paths, paths.rev()) {
  let points = subset.map(cetz.path-util.bounds).join()
  let corners = subset.map(path => {
    let bounds = aabb(cetz.path-util.bounds(path))
    (bounds.low, bounds.high)
  }).join()
  assert.eq(cbor.encode(aabb(corners)), cbor.encode(aabb(points)))
}

// Exercise the canvas and debug bounds consumers with gradients and content.
#for debug in (false, true) {
  cetz.canvas(debug: debug, {
    cetz.draw.get-ctx(ctx => {
      let nonfinite = (((calc.inf - calc.inf, 0, 0), false,
        (("l", (-3, -7, 2)),)),)
      let combined = paths.first() + nonfinite
      let actual = cetz.process.element(ctx, ctx => (
        ctx: ctx,
        drawables: (cetz.drawable.path(paths.first()), cetz.drawable.path(nonfinite)),
      )).bounds
      assert.eq(cbor.encode(actual), cbor.encode(aabb(cetz.path-util.bounds(combined))))
      ()
    })
    cetz.draw.scale(x: 0.7, y: 1.2)
    cetz.draw.bezier((0, 0), (3, 0), (1, 3), (2, -4),
      fill: gradient.linear(red, blue), stroke: black)
    cetz.draw.content((1, 1), [bounds], frame: "rect", fill: white)
  })
}

// These no-ink paths still own bounds and named anchors. Only their final
// painting can be omitted; content links and visible custom paint remain.
#let no-ink-scene = context {
  for debug in (false, true) {
    let scene = cetz.canvas(length: 1cm, padding: 0, debug: debug, {
      cetz.draw.line((-3, -2), (5, 4), name: "ghost", stroke: none)
      cetz.draw.get-ctx(ctx => {
        let (_, point) = cetz.coordinate.resolve(ctx, "ghost.end")
        assert.eq(point, (5, 4, 0), message: "unpainted named path retains its endpoint anchor")
        ()
      })
      cetz.draw.line("ghost.end", (6, 4), stroke: rgb("#8732a8") + 0.5pt)
      cetz.draw.rect((0, 0), (2, 1), stroke: none,
        fill: gradient.linear(rgb("#08a17b"), rgb("#e98224")))
      cetz.draw.line((0, 2), (2, 2), stroke: rgb("#137c45") + 0.7pt)
      cetz.draw.content((1, 3), link("https://example.com/cetz-no-ink",
        box(width: 8pt, height: 8pt)), padding: 0)
    })
    let size = measure(scene)
    assert(calc.abs(size.width - 9cm) < 1e-8pt,
      message: "unpainted path retains the canvas width")
    assert(calc.abs(size.height - 6cm) < 1e-8pt,
      message: "unpainted path retains the canvas height")
    scene
  }
}
#no-ink-scene
