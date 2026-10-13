#import "crates/linnest/typst/src/impl/draw.typ" as drawing
#import "crates/linnest/typst/src/lib.typ": graph, draw
#import "crates/kurvst/typst/src/lib.typ" as curve
#import "@preview/cetz:0.5.2" as cetz

#set page(width: auto, height: auto, margin: 0pt)
#let length = 0.45 * 3.688720
#let accuracy = 0.001
#let style = (
  stroke: black + 0.5pt,
  pattern: "coil",
  pattern-amplitude: 0.15,
  pattern-wavelength: 0.45,
  pattern-fit: true,
  pattern-accuracy: accuracy,
  pattern-phase: calc.pi / 2,
  pattern-coil-longitudinal-scale: 1.4,
)

// Resolve label-directed offsets on the full carrier. Choosing a side on each
// half independently can put its marks on the opposite side from the coil.
#{
  let path = curve.cubic((0, 0), (1, 3), (4, -2), (5, 1))
  let accuracy = 1e-8
  let total = curve.length(path, accuracy: accuracy)
  let source = curve.trim(path, end-outset: total * 3 / 4, accuracy: accuracy)
  let sink = curve.trim(path, start-outset: total / 4, accuracy: accuracy)
  for label-pos in ((0, 1), (2, 3), (4, -1)) {
    let offset = curve.layer(path, offset: 0.3, side-point: label-pos, accuracy: accuracy).offset
    for gap in (0, 0.4) {
      let automatic = style + (
        offset: 0.3, offset-side: "label", split-gap: gap,
        accuracy: accuracy, pattern-accuracy: accuracy,
      )
      let fixed = automatic + (offset: offset, offset-side: none)
      let halves(style) = drawing._split-edge-geometry(
        source, sink, path, style, style, 0, 0, label-pos, accuracy,
      )
      let automatic = halves(automatic)
      let fixed = halves(fixed)
      assert.eq(automatic.source, fixed.source)
      assert.eq(automatic.sink, fixed.sink)
      assert.eq(automatic.pattern-whole, fixed.pattern-whole)
    }
  }
}

#for (length, wavelength, periods) in ((2.0, 0.4, 4.5), (2.0, 0.43, 4.5), (0.1, 1.0, 0.5)) {
  for amplitude in (-0.15, 0, 0.15) {
    let fitted = curve.coil(fit-length: length, amplitude: amplitude, wavelength: wavelength, longitudinal-scale: 1.4)
    assert(not fitted.endpoint-ramp)
    assert.eq(fitted.points.len(), int(periods * 16) + 1)
    assert.eq(fitted.points.first(), (at: 0, x: 0, y: 0))
    assert.eq(fitted.points.last(), (at: 1, x: 0, y: 0))
    let rendered = curve.points(curve.pattern(
      curve.line((0, 0), (length, 0)),
      pattern: fitted, amplitude: amplitude, wavelength: length,
      samples-per-period: fitted.points.len() - 1,
      anchor-start: false, anchor-end: false,
    ))
    assert.eq(rendered.len(), fitted.points.len())
    for (p, actual) in fitted.points.zip(rendered) {
      let theta = calc.pi + 2 * calc.pi * periods * p.at
      let longitudinal = calc.abs(amplitude) * 1.4
      let x = (p.at * length + longitudinal * (1 + calc.cos(theta))) * length / (length + 2 * longitudinal)
      assert(calc.abs(actual.at(0) - x) < 1e-8)
      assert(calc.abs(actual.at(1) - amplitude * calc.sin(theta)) < 1e-8)
    }
  }
}

// A short first segment must not shorten the taper or independently fit a turn.
#let path = curve.path(
  curve.line((0, 0), (0.1, 0)),
  curve.line((0.1, 0), (length, 0)),
)
#let fitted = curve.coil(fit-length: length, amplitude: 0.15, wavelength: 0.45, longitudinal-scale: 1.4)
#cetz.canvas({
  (ctx => {
    assert(not drawing._same-pattern-geometry(style, style + (pattern-endpoint-slope: 1)))
    assert(not drawing._same-pattern-geometry(style, style + (pattern-natural-endpoints: true)))
    for (fit, start, end, wavelength, slope, natural) in (
      (true, true, true, length / 4, 0, false),
      (false, true, true, 0.45, 0, false),
      (true, true, false, 0.45, 0, false),
      (true, false, true, 0.45, 0, false),
      (true, true, true, length / 4, 1, false),
      (true, true, true, length / 4, 3, false),
      (true, true, false, 0.45, 1, false),
      (true, false, true, 0.45, 1, false),
      (true, true, true, 0.45, 0, true),
      (false, true, true, 0.45, 0, true),
      (true, true, true, 0.45, 3, true),
      (true, true, false, 0.45, 1, true),
      (true, false, true, 0.45, 1, true),
      (true, false, false, 0.45, 1, true),
    ) {
      let actual = drawing._segments-elements(
        curve.segments(path),
        style + (pattern-fit: fit, pattern-endpoint-slope: slope, pattern-natural-endpoints: natural),
        auto, start, end,
      )
      let natural = natural and start and end
      let expected = curve.pattern(
        curve.line((0, 0), (length, 0)),
        pattern: if natural { fitted } else { "coil" },
        amplitude: 0.15,
        wavelength: if natural { length } else { wavelength },
        phase: if natural { 0 } else { calc.pi / 2 },
        samples-per-period: if natural { fitted.points.len() - 1 } else { 16 },
        coil-longitudinal-scale: 1.4,
        anchor-start: start,
        anchor-end: end,
        endpoint-slope: slope,
        accuracy: accuracy,
      )
      assert.eq(type(actual), array)
      let rendered = cetz.process.many(ctx, actual.flatten(), compute-bounds: true)
      let reference = cetz.process.many(ctx, curve.to-cetz(expected, stroke: style.stroke).flatten(), compute-bounds: true)
      // Arc-length inversion may differ within its accuracy on split baselines;
      // drawing structure and styles must still match exactly.
      assert.eq(rendered.drawables.len(), reference.drawables.len())
      for (actual, expected) in rendered.drawables.zip(reference.drawables) {
        let actual-paths = actual.remove("segments")
        let expected-paths = expected.remove("segments")
        assert.eq(actual, expected)
        assert.eq(actual-paths.len(), expected-paths.len())
        for ((a-origin, a-closed, a-curves), (b-origin, b-closed, b-curves)) in actual-paths.zip(expected-paths) {
          assert(cetz.vector.dist(a-origin, b-origin) <= accuracy)
          assert.eq(a-closed, b-closed)
          assert.eq(a-curves.len(), b-curves.len())
          for (a, b) in a-curves.zip(b-curves) {
            assert.eq(a.first(), b.first())
            assert.eq(a.len(), b.len())
            for (a-point, b-point) in a.slice(1).zip(b.slice(1)) {
              assert(cetz.vector.dist(a-point, b-point) <= accuracy)
            }
          }
        }
      }
    }
    (ctx: ctx)
  },)
})

// The public crossing path must also accept two paints on one continuous coil,
// retaining the carrier arrow while the crossing removes part of its stroke.
#for natural in (false, true) {
  let g = graph.build({
    graph.node(<a>, pos: graph.pos(x: 0, y: 0))
    graph.node(<b>, pos: graph.pos(x: 4, y: 0))
    graph.node(<c>, pos: graph.pos(x: 2, y: -1))
    graph.node(<d>, pos: graph.pos(x: 2, y: 1))
    graph.edge(graph.source(<a>), <coil>, graph.sink(<b>), pos: graph.pos(x: 1, y: 0))
    graph.edge(graph.source(<c>), <over>, graph.sink(<d>), pos: graph.pos(x: 2, y: 0))
  })
  let coil = style + (
    pattern-natural-endpoints: natural,
    crossing-under: <over>,
    crossing-gap: 0.4,
    mark: curve.mark.straight(length: 4.5pt, width: 3.375pt),
    mark-position: 0.8,
  )
  draw(graph.style(g, node-label: none, edge-label: none,
      node-style: (radius: 0, stroke: none, fill: none)),
    title: none, node-outset: 0,
    source-style: edge => if edge.eid == 0 { coil + (stroke: red + 0.5pt) } else { (stroke: black + 0.5pt) },
    sink-style: edge => if edge.eid == 0 { coil + (stroke: blue + 0.5pt) } else { (stroke: black + 0.5pt) },
  )
}

// Colour and gap boundaries must cut the already fitted coil. The reference
// slices known sample knots, independently of the carrier-window implementation.
#cetz.canvas({
  (ctx => {
    let same-geometry(actual, expected) = {
      // Arc-length inversion may put a boundary just beyond a sample knot.
      // Ignore only curve fragments whose entire control polygon is below the
      // same coordinate tolerance used for the geometry comparison.
      let nondegenerate(origin, curves) = {
        let kept = ()
        for segment in curves {
          if segment.slice(1).any(p => cetz.vector.dist(origin, p) >= 1e-6) {
            kept.push(segment)
          }
          origin = segment.last()
        }
        kept
      }
      let paths(elements) = cetz.process.many(ctx, elements.flatten(), compute-bounds: true)
        .drawables.fold((), (paths, drawable) => paths + drawable.segments)
      let actual = paths(actual)
      let expected = paths(expected)
      assert.eq(actual.len(), expected.len())
      for ((a-origin, a-closed, a-curves), (b-origin, b-closed, b-curves)) in actual.zip(expected) {
        assert(cetz.vector.dist(a-origin, b-origin) < 1e-6)
        assert.eq(a-closed, b-closed)
        let a-curves = nondegenerate(a-origin, a-curves)
        let b-curves = nondegenerate(b-origin, b-curves)
        assert.eq(a-curves.len(), b-curves.len())
        for (a, b) in a-curves.zip(b-curves) {
          assert.eq(a.first(), b.first())
          assert.eq(a.len(), b.len())
          for (p, q) in a.slice(1).zip(b.slice(1)) {
            assert(cetz.vector.dist(p, q) < 1e-6)
          }
        }
      }
    }
    for path in (
      curve.line((0, 0), (4, 0)),
      curve.cubic((0, 0), (1, 3), (4, -2), (5, 1)),
      curve.cubic((0, 0), (3, 3), (-3, 3), (0, 0)),
    ) {
      let accuracy = 1e-8
      let total = curve.length(path, accuracy: accuracy)
      let source = curve.trim(path, end-outset: total * 3 / 4, accuracy: accuracy)
      let sink = curve.trim(path, start-outset: total / 4, accuracy: accuracy)
      for natural in (false, true) {
        let style = style + (
          pattern-wavelength: total / 4,
          pattern-natural-endpoints: natural,
          accuracy: accuracy,
          pattern-accuracy: accuracy,
        )
        let fitted = curve.coil(fit-length: total, amplitude: 0.15,
          wavelength: total / 4, longitudinal-scale: 1.4)
        let reference = curve.pattern(path,
          pattern: if natural { fitted } else { "coil" },
          amplitude: 0.15,
          wavelength: if natural { total } else { total / 4 },
          phase: if natural { 0 } else { calc.pi / 2 },
          samples-per-period: if natural { fitted.points.len() - 1 } else { 16 },
          coil-longitudinal-scale: 1.4,
          accuracy: accuracy,
        )
        let segments = curve.segments(reference)
        assert.eq(calc.rem(segments.len(), 8), 0)
        let eighth = int(segments.len() / 8)
        let painted(start, end, paint) = curve.to-cetz(
          curve.path(..segments.slice(start, end).map(curve.from-cubic)),
          stroke: paint + 0.5pt,
        )
        for gap in (0, total / 4) {
          let source-style = style + (stroke: red + 0.5pt, split-gap: gap)
          let sink-style = style + (stroke: blue + 0.5pt, split-gap: gap)
          let halves = drawing._split-edge-geometry(
            source, sink, path, source-style, sink-style, 0, 0, none, accuracy,
          )
          let inset = if gap == 0 { 0 } else { eighth }
          same-geometry(
            drawing._pattern-edge-halves(halves, source-style, sink-style),
            painted(0, 2 * eighth - inset, red)
              + painted(2 * eighth + inset, segments.len(), blue),
          )
        }
        same-geometry(
          drawing._cut-path-elements(path, style + (crossing-gap: total / 4), (total / 2,)),
          painted(0, 3 * eighth, black) + painted(5 * eighth, segments.len(), black),
        )
        same-geometry(
          drawing._cut-path-elements(path, style + (crossing-gap: total / 4), (total / 2,),
            paint-windows: (
              (start: 0, end: total / 4, style: style + (stroke: red + 0.5pt)),
              (start: total / 4, end: total, style: style + (stroke: blue + 0.5pt)),
            ),
          ),
          painted(0, 2 * eighth, red) + painted(2 * eighth, 3 * eighth, blue)
            + painted(5 * eighth, segments.len(), blue),
        )
      }
    }
    (ctx: ctx)
  },)
})
