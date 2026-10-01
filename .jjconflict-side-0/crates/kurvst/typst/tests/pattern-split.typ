#import "../src/lib.typ" as kurvst

// Exercise the public Typst -> CBOR -> WASM boundary with both coil modes.
#for base in (
  kurvst.path(kurvst.line((0, 0), (3, 0))),
  kurvst.hobby-spline(((0, 0), (0.8, 1), (2, -0.4), (3, 0))),
) {
  let L = kurvst.length(base)
  let natural = kurvst.coil(fit-length: L, amplitude: 0.12, wavelength: 0.55)
  for options in (
    (pattern: "coil", amplitude: 0.12, wavelength: 0.55),
    (
      pattern: natural,
      amplitude: 0.12,
      wavelength: L,
      samples-per-period: natural.points.len() - 1,
    ),
    (pattern: "wave", amplitude: 0.12, wavelength: 0.55),
    (pattern: "zigzag", amplitude: 0.12, wavelength: 0.55),
  ) {
    let whole = kurvst.pattern(base, ..options)
    let split = kurvst.pattern(
      base,
      ..options,
      split-at: (L * 0.413, L * 0.627),
    )
    assert.eq(split.path, whole.path)
    assert.eq(whole.parts.at(0).path, whole.path)
    assert.eq(split.parts.len(), 3)
    for index in range(2) {
      assert.eq(
        kurvst.points(split.parts.at(index)).last(),
        kurvst.points(split.parts.at(index + 1)).first(),
      )
    }
    // Integer cuts and repeated/end cuts must retain stable part indices.
    let boundary = kurvst.pattern(base, ..options, split-at: (0, L, 0))
    assert.eq(boundary.parts.len(), 4)
    for index in (0, 1, 3) {
      assert.eq(kurvst.elements(boundary.parts.at(index)), ())
    }
    assert.eq(boundary.parts.at(2).path, whole.path)
  }
}

// Compare the numeric owner with the original Typst callback at its first boundary.
// CBOR equality also checks float/signed-zero values and integer endpoint coordinates.
#for length in (0.01, 0.1, 1.375, 3.0, 3.141592653589793) {
  for amplitude in (-0.12, -0.0, 0.0, 0.12) {
    for samples-per-period in (-2.5, -2, 0, 0.5, 1.0, 1, 3, 16, 31) {
      let wavelength = 0.55
      let longitudinal-scale = 1.25
      let span = length + calc.abs(amplitude) * longitudinal-scale * 2
      let periods = calc.max(1, calc.round(length / wavelength)) - 0.5
      let samples = calc.max(2, int(calc.ceil(periods * calc.max(1, samples-per-period))))
      let scale = longitudinal-scale * (length / span) * if amplitude < 0 { -1 } else { 1 }
      let expected = range(0, samples + 1).map(index => {
        let at = index / samples
        let theta = 2 * calc.pi * at
        let position = theta / (2 * calc.pi)
        let phase = calc.pi + periods * theta
        let offset = if position == 0 or position == 1 { (0, 0) } else {
          (scale * (1 + calc.cos(phase) - 2 * position), calc.sin(phase))
        }
        (at: at, x: offset.at(0), y: offset.at(1))
      })
      let actual = kurvst.coil(fit-length: length, amplitude: amplitude,
        wavelength: wavelength, samples-per-period: samples-per-period)
      assert.eq(cbor.encode(actual.points), cbor.encode(expected))
    }
  }
}

#for base in (
  kurvst.path(kurvst.line((0, 0), (3, 0))),
  kurvst.hobby-spline(((0, 0), (0.8, 1), (2, -0.4), (3, 0))),
) {
  let L = kurvst.length(base)
  let natural = kurvst.coil(fit-length: L, amplitude: 0.12, wavelength: 0.55)
  let public = kurvst.pattern(base, pattern: natural, amplitude: 0.12, wavelength: L,
    samples-per-period: natural.points.len() - 1, split-at: (L * 0.413, L * 0.627))
  let direct = kurvst.pattern(base, pattern: (
    kind: "fitted-coil", fit-length: L, amplitude: 0.12, wavelength: 0.55,
    samples-per-period: 16, longitudinal-scale: 1.25,
  ), amplitude: 0.12, wavelength: L, split-at: (L * 0.413, L * 0.627))
  assert.eq(public, direct)
}
