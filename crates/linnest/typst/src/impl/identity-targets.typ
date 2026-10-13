// Fixed-size identity links are not measured content or styled annotations.
#import "@preview/cetz:0.5.2" as cetz

#let draw(ctx, targets) = {
  let drawables = ()
  for target in targets {
    let body = target.body
    assert(body.func() == link and body.body.func() == box,
      message: "identity targets require a link containing a fixed-size box")
    let width = body.body.width
    let height = body.body.height
    assert(width != auto and height != auto and width.ratio == 0% and height.ratio == 0%,
      message: "identity target dimensions must be absolute lengths")
    width = width.length
    height = height.length
    let (next, point) = cetz.coordinate.resolve(ctx, target.position)
    ctx = next
    let center = cetz.matrix.mul4x4-vec3(ctx.transform, point)
    // Like ordinary CeTZ content, only the center follows the canvas
    // transform. The transparent link rectangle retains its physical size.
    let drawable = cetz.drawable.content(
      center, width / ctx.length, height / ctx.length, (), body,
    )
    drawable.tags = (cetz.drawable.TAG.no-bounds,)
    drawables.push(drawable)
  }
  (ctx: ctx, drawables: drawables)
}
