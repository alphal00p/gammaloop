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
