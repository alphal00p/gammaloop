#import "@preview/cetz:0.5.2" as cetz
#import "../src/impl/identity-targets.typ" as targets
#set page(width: auto, height: auto, margin: 0pt)

#cetz.canvas({
  cetz.draw.line((-1, -1), (2, 1), stroke: 0pt)
  cetz.draw.get-ctx(ctx => {
    let body = link("#linnet-edge-7?source", box(width: 8pt, height: 6pt))
    for transform in (
      cetz.matrix.ident(4),
      cetz.matrix.transform-scale((-1.5, 0.6, 1)),
      cetz.matrix.transform-rotate-xyz(40deg, 30deg, 12deg),
    ) {
      let local = ctx
      local.transform = transform
      let input = (
        (position: (1, 2), body: body),
        (position: (rel: (0.2, 0.4)), body: body),
      )
      let actual = targets.draw(local, input)
      let expected = local
      for target in input {
        let (next, point) = cetz.coordinate.resolve(expected, target.position)
        expected = next
      }
      assert.eq(actual.ctx.prev, expected.prev)
      assert.eq(actual.drawables.len(), input.len())
      for drawable in actual.drawables {
        assert.eq(drawable.body, body)
        assert.eq(drawable.width, 8pt / ctx.length)
        assert.eq(drawable.height, 6pt / ctx.length)
        assert.eq(drawable.tags, (cetz.drawable.TAG.no-bounds,))
      }
      assert.eq(cetz.process.element(local, ctx => targets.draw(ctx, input)).bounds, none)
    }
    // Emit one prepared element, retaining both independent source identities.
    (ctx => targets.draw(ctx, (
      (position: (0, 0), body: body),
      (position: (1, 0), body: link("#linnet-node-3", box(width: 8pt, height: 8pt))),
    )),)
  })
})
