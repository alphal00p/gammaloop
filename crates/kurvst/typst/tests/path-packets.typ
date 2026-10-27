#import "../src/lib.typ" as kurvst

// Run with --input contract=close-mode, unknown-field, or geometry-format
// to verify the corresponding strict rejection through the packet reader.
#let contract = sys.inputs.at("contract", default: "valid")
#if contract == "geometry-format" {
  kurvst.mark.geometry((), (), (), format: "unknown")
} else if contract == "close-mode" {
  kurvst.elements(cbor.encode((path: (elements: (
    (kind: "close", mode: "straight"),
  ),))))
} else if contract == "unknown-field" {
  kurvst.elements(cbor.encode((path: (elements: (
    (kind: "move", start: (0, 0), extra: true),
  ),))))
} else {
  assert.eq(contract, "valid")
  for native in (
    kurvst.path(),
    kurvst.path(kurvst.move-to((-0.0, 0.0)), kurvst.line-to((3, 2)), kurvst.close()),
    kurvst.quad((0, 0), (2, 4), (5, 1)),
    kurvst.cubic((0, 0), (1, 3), (4, 2), (6, 0)),
  ) {
    let packet = cbor.encode(native)
    assert.eq(kurvst.elements(packet), kurvst.elements(native))
    assert.eq(kurvst.points(packet), kurvst.points(native))
    assert.eq(kurvst.segments(packet), kurvst.segments(native))
    assert.eq(kurvst.to-cetz-data(packet), kurvst.to-cetz-data(native))
    assert.eq(kurvst.to-native(packet, unit: 1pt), kurvst.to-native(native, unit: 1pt))
  }
}
