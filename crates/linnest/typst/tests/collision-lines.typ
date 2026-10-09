#import "@preview/cetz:0.5.2" as cetz
#import "../src/impl/draw.typ" as draw

#cetz.canvas({
  cetz.draw.line((0, 0), (1, 1))
  cetz.draw.get-ctx(ctx => {
    let path(points, thickness: 2pt, tags: ()) = cetz.drawable.path(
      ((points.first(), false, (("l", ..points.slice(1)),)),),
      stroke: thickness, tags: tags,
    )
    let drawables = (
      path(((0, 0, 0), (1, 0, 0))),
      path(((2, 0, 0), (3, 0, 0)), thickness: 4pt),
      path(((0, 1, 0), (1, 1, 0)), tags: (cetz.drawable.TAG.hidden,)),
      cetz.drawable.content((0, 0, 0), 1, 1, (), []),
      cetz.drawable.path((((0, 0, 0), false, (("l", (9, 9, 0)),)),), stroke: none),
    )
    let actual = draw._edge-collision-lines(ctx, drawables, painted: true)
    assert.eq(actual, (
      (start: (0.0, 0.0), end: (1.0, 0.0), radius: 1pt / ctx.length),
      (start: (2.0, 0.0), end: (3.0, 0.0), radius: 2pt / ctx.length),
    ))
    let transformed = cetz.drawable.apply-transform(
      cetz.matrix.transform-scale((-2, 3, 1)), drawables,
    )
    let shifted = draw._edge-collision-lines(ctx, transformed, painted: true)
    assert.eq(shifted.first().end, (-2.0, 0.0))
    assert.eq(shifted.last().start, (-4.0, 0.0))
    assert.eq(shifted.last().end, (-6.0, 0.0))
    ()
  })
})
