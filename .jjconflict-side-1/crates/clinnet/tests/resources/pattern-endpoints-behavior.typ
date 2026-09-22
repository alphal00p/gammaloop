#import "crates/linnest/typst/src/impl/draw.typ" as drawing
#import "crates/kurvst/typst/src/lib.typ" as curve
#import "@preview/cetz:0.5.1" as cetz

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
      assert(calc.abs(actual.length - length) < 1e-8)
      let rendered = cetz.process.many(ctx, actual.elements.flatten(), compute-bounds: true)
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
