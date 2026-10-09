#import "@preview/cetz:0.5.2" as cetz
#import "crates/linnest/typst/src/impl/identity-targets.typ" as identity-targets

// Identity targets retain fixed-size link rectangles and coordinate resolver
// context; styled annotation content remains the ordinary CeTZ owner's job.
#cetz.canvas({
  cetz.draw.get-ctx(ctx => {
    // A line strip normalizes exactly like a general path, including repeated
    // points, explicit closing points and degenerate one-point paths.
    for points in (
      ((0, 0, 0),),
      ((0, 0, 0), (0, 0, 0)),
      ((0, 0, 0), (1, 2, 3), (1, 2, 3), (4, 2, 1)),
      ((0, 0, 0), (1, 2, 3), (0, 0, 0), (0, 0, 0)),
    ) {
      for close in (false, true) {
        let style = (stroke: red + 0.4pt, fill: blue, tags: ("test",))
        assert.eq(
          cetz.drawable.line-strip(points, close: close, ..style),
          cetz.drawable.path(((points.first(), close, points.slice(1).map(p => ("l", p))),), ..style),
        )
      }
    }

    ()
  })
})

#{
  cetz.canvas({
    cetz.draw.get-ctx(ctx => {
      let check(local, targets, ..style) = {
        let expected = cetz.process.many(local, targets.map(target => {
          cetz.draw.floating(cetz.draw.content(target.position, target.body,
            tags: (cetz.drawable.TAG.no-bounds,), ..style))
        }).flatten(), compute-bounds: false)
        let actual = identity-targets.draw(local, targets)
        // Geometry can contain NaNs, so compare its encoded representation.
        assert.eq(cbor.encode(actual.ctx.prev), cbor.encode(expected.ctx.prev))
        // Ordinary content also emits an empty frame; identity targets paint
        // only the fixed-size linked content, not that annotation decoration.
        let content = expected.drawables.filter(drawable => drawable.type == "content")
        assert.eq(actual.drawables.len(), content.len())
        for (a, b, target) in actual.drawables.zip(content, targets) {
          assert.eq(cbor.encode(a.pos), cbor.encode(b.pos), message: "identity content center")
          assert.eq(a.width, target.body.body.width.length / local.length,
            message: "identity width is the authored physical width, not measured text")
          assert.eq(a.height, target.body.body.height.length / local.length,
            message: "identity height is the authored physical height, not measured text")
          assert.eq(a.body, target.body, message: "identity preserves each independent link")
          assert.eq(a.segments, (), message: "identity targets have no annotation frame geometry")
          assert.eq(a.tags, (cetz.drawable.TAG.no-bounds,))
        }
      }
      let body = link("https://example.com/identity", box(width: 8pt, height: 8pt))
      let coords = ((0, 0), (-0.0, 0.0), (1.25, -4e-8), (9007199254740995, -9007199254740995))
      for transform in (
        cetz.matrix.ident(4),
        cetz.matrix.transform-scale((0, 0, 1)),
        cetz.matrix.transform-scale((-1.5, 0.6, 1)),
        cetz.matrix.transform-rotate-xyz(40deg, 30deg, 12deg),
        ((1.2, 0.75, 0, 0.1), (-0.1, -1.2, 0, -0.75), (0, 0, 1, 0), (0, 0, 0, 1)),
      ) {
        let local = ctx
        local.transform = transform
        for target-body in (body, link("https://example.com/empty", box(width: 0pt, height: 0pt))) {
          let targets = coords.map(position => (position: position, body: target-body))
          check(local, targets, padding: 0)
        }
      }
      for position in ((1cm, -2cm), (1, 2, 3), (x: 1.25, y: -2.5), (30deg, 2), (rel: (0.2, 0.4)), ()) {
        check(ctx, ((position: position, body: body), (position: (1, -1), body: body)), padding: 0)
      }
      // Overflowing geometry uses scalar comparison semantics.
      for position in ((1e308, 1e308),) {
        let local = ctx
        local.transform = cetz.matrix.transform-scale((2.0, -2.0, 1.0))
        check(local, ((position: position, body: body),), padding: 0)
      }
      ()
    })
  })
}

// Use the public sampling owner directly, without the unrelated Linnest wrapper.
#import "crates/kurvst/typst/src/lib.typ" as region-curve
#let region-reference(parts, regions: 4, unit: 1pt, step: 4pt) = {
  let distance(a, b) = {
    let dx = a.at(0) - b.at(0)
    let dy = a.at(1) - b.at(1)
    calc.sqrt(dx * dx + dy * dy)
  }
  let paths = parts.map(part => region-curve.path(..part.segments.map(region-curve.from-cubic)))
  let lengths = paths.map(region-curve.length)
  let total = lengths.sum()
  let offset = 0
  let result = ()
  for (part, path, length) in parts.zip(paths, lengths) {
    if part.visible and length > 0 {
      for region in range(regions) {
        let start = calc.max(0, total * region / regions - offset)
        let end = calc.min(length, total * (region + 1) / regions - offset)
        if start >= end { continue }
        let trimmed = region-curve.trim(path, start-outset: start, end-outset: length - end)
        let points = ()
        for segment in region-curve.segments(trimmed) {
          let speed = 3 * calc.max(
            distance(segment.start, segment.control-start),
            distance(segment.control-start, segment.control-end),
            distance(segment.control-end, segment.end),
          )
          let steps = calc.max(1, int(calc.ceil(speed * calc.abs(unit / step))))
          for index in range(steps + 1) {
            points.push(region-curve.cubic-point(segment, index / steps))
          }
        }
        result.push((region, points))
      }
    }
    offset += length
  }
  result
}
#let check-regions(parts, regions: 4, unit: 1pt, step: 4pt) = {
  let expected = region-reference(parts, regions: regions, unit: unit, step: step)
  let actual = region-curve.region-samples(parts, regions: regions, unit: unit / 1pt, step: step / 1pt)
  assert.eq(cbor.encode(actual), cbor.encode(expected))
}

// Native region sampling must preserve the scalar arc-length partitions,
// both endpoint samples and de Casteljau coordinates bit for bit.
#let region-line(a, b) = region-curve.line-segment(a, b)
#let region-parts = (
  (segments: (region-line((-0.0, 0.0), (0, 0)),), visible: true),
  (segments: (region-line((0, 0), (1, 0)), region-line((1, 0), (2, 0))), visible: true),
  (segments: (region-line((2, 0), (3, 0)),), visible: false),
  (segments: (
    (start: (3, 0), control-start: (3, 4), control-end: (-1, 2), end: (4, 0)),
    (start: (7, 2), control-start: (7, 2), control-end: (9, -1), end: (9, -1)),
  ), visible: true),
)
#for regions in (1, 2, 4, 7) {
  for unit in (0pt, 1pt, -1pt, 1mm, -1mm) {
    check-regions(region-parts, regions: regions, unit: unit)
  }
}
#for segments in (
  (region-line((1e-20, -1e-20), (2e-20, -2e-20)),),
  (region-line((9007199254740995, 9007199254740993), (9007199255740995, 9007199255740993)),),
) {
  check-regions(((segments: segments, visible: true),), unit: 1e-6pt)
}
#assert.eq(region-curve.region-samples(()), ())
