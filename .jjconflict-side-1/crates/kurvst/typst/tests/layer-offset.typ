#import "../src/lib.typ" as kurvst

#let base = kurvst.line((0, 0), (4, 0))
#for (side, expected) in (((2, 1), 0.3), ((2, -1), -0.3)) {
  let layer = kurvst.layer(base, offset: 0.3, side-point: side, length: 2)
  assert.eq(layer.offset, expected)
  assert.eq(layer.path, kurvst.layer(base, offset: expected, length: 2).path)
}

// This curved carrier and its first quarter choose different local signs.
#let base = kurvst.cubic((0, 0), (1, 3), (4, -2), (5, 1))
#let L = kurvst.length(base)
#let first = kurvst.trim(base, end-outset: L * 3 / 4)
#let whole = kurvst.layer(base, offset: 0.3, side-point: (0, 1))
#let local = kurvst.layer(first, offset: 0.3, side-point: (0, 1))
#assert.eq(whole.offset, -0.3)
#assert.eq(local.offset, 0.3)
#let shared = kurvst.layer(first, offset: whole.offset)
#assert.eq(shared.offset, whole.offset)
#assert.eq(shared.path, kurvst.parallel(first, distance: -0.3).path)

#for offset in (0, 0.3, -0.3) {
  let empty = kurvst.layer(kurvst.path(), offset: offset, side-point: (0, 1))
  assert.eq(empty.offset, offset)
  assert.eq(kurvst.elements(empty), ())
}

// Retain the scalar orchestration as an independent oracle: delegating layer()
// back to layers() would not check the batch's clamps, outsets, or cubics.
#let shifts = (-100, -0.3, -0.0, 0, 0.2, 100)
#for carrier in (
  kurvst.path(), kurvst.path(kurvst.move-to((2, 3))),
  kurvst.line((0, 0), (0, 0)), kurvst.line((0, 0), (5, 1)),
  kurvst.quad((0, 0), (1, 3), (4, 0)), base,
  kurvst.path(kurvst.move-to((0, 0)), kurvst.line-to((2, 1)), kurvst.close()),
  kurvst.path(kurvst.move-to((0, 0)), kurvst.line-to((2, 1)),
    kurvst.move-to((4, 1)), kurvst.quad-to((5, 3), (7, 0))),
) {
  for options in (
    (:), (length: 1.2), (ratio: 0.35),
    (length: 1.2, ratio: 0.5, resolve-length: "max"),
    (length: 1.2, ratio: 0.5, resolve-length: "none"),
    (length: 1.2, ratio: 0.5, resolve-length: args => (args.length + args.ratio) / 3),
    (start-outset: 0.2, end-outset: 0.3, length: 1.0),
    (start-outset: 20, end-outset: 30, ratio: 0.2),
  ) {
    let slices = kurvst.layers(carrier, shifts, ..options)
    assert.eq(slices.len(), shifts.len())
    assert.eq(kurvst.layers(carrier, (), ..options), ())
    let packet = kurvst.layers(carrier, shifts, format: "cbor", unit: -2, ..options)
    assert.eq(packet.count, shifts.len())
    let restored = cbor(packet.layers).map(slice => (
      path: slice.path + (offset: packet.offset), segments: slice.segments,
    ))
    assert.eq(cbor.encode(restored), cbor.encode(slices))
    let footprints = slices.map(slice => {
      let single = slice.segments.len() == 1
      let path = if single { kurvst.from-cubic(slice.segments.first()) } else { slice.path }
      (segments: kurvst.to-cetz-data(path, unit: if single { 1 } else { -2 }), single: single)
    })
    let nonzero = slices.any(slice => slice.segments.any(segment => (
      segment.start != segment.control-start or segment.start != segment.control-end
        or segment.start != segment.end
    )))
    assert.eq(packet.nonzero, nonzero)
    assert.eq(packet.all-single, footprints.all(path => path.single))
    if packet.supported {
      assert.eq(packet.footprints, cbor.encode(footprints))
    } else {
      assert(footprints.all(path => path.segments.len() == 0
        or path.segments.any(subpath => subpath.last().len() == 0)))
    }
    assert.eq(kurvst.frames(carrier, shifts, format: "cbor"),
      cbor.encode(kurvst.frames(carrier, shifts)))
    for (shift, slice) in shifts.zip(slices) {
      let expected = (path: carrier.at("path", default: carrier), offset: 0)
      if kurvst.segments(carrier).len() > 0 {
        let center = kurvst.center-outset(kurvst.length(carrier), ..options)
        let at = calc.max(-center, calc.min(center, shift))
        expected = kurvst.trim(carrier,
          start-outset: options.at("start-outset", default: 0) + center + at,
          end-outset: options.at("end-outset", default: 0) + center - at,
        ) + (offset: 0)
      }
      assert.eq(cbor.encode(slice.path), cbor.encode(expected))
      assert.eq(cbor.encode(slice.segments), cbor.encode(kurvst.segments(expected)))
    }
  }
}
#for side in ((0, 1), (0, -1)) {
  let offset = kurvst.layer(base, offset: 0.3, side-point: side).offset
  let carrier = kurvst.parallel(base, distance: offset)
  let expected = kurvst.layers(carrier, shifts, length: 1.2, ratio: 0.5)
  let actual = kurvst.layers(base, shifts, offset: 0.3, side-point: side,
    length: 1.2, ratio: 0.5)
  for (slice, reference) in actual.zip(expected) {
    assert.eq(slice.path.offset, offset)
    assert.eq(cbor.encode(slice.path.path), cbor.encode(reference.path.path))
    assert.eq(cbor.encode(slice.segments), cbor.encode(reference.segments))
  }
}

// The single Bezier drawing path ignores an attached multi-path unit.
#let single = kurvst.layers(base, (0,), format: "cbor", unit: 0)
#assert(single.supported and single.all-single)
#let multi = kurvst.path(kurvst.move-to((0, 0)), kurvst.line-to((1, 2)), kurvst.line-to((3, 0)))
#let collapsed = kurvst.layers(multi, (0,), format: "cbor", unit: 0)
#assert(not collapsed.supported and not collapsed.all-single)
#let empty-batch = kurvst.layers(base, (), format: "cbor")
#assert.eq(empty-batch.count, 0)
#assert.eq(cbor(empty-batch.layers), ())
#assert.eq(cbor(empty-batch.footprints), ())
#assert(not empty-batch.all-single and not empty-batch.nonzero)
