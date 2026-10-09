#import "../src/mark.typ" as mark

#let mode = sys.inputs.at("contract", default: "valid")
#if mode == "valid" {
  let constructors = (
    mark.triangle,
    mark.straight,
    mark.stealth,
    mark.round,
    mark.tikz,
    mark.barb,
    mark.hooks,
    mark.bar,
    mark.bracket,
    mark.circle,
    mark.square,
    mark.diamond,
    mark.rays,
  )
  for constructor in constructors {
    let spec = constructor()
    assert(type(spec) == dictionary and spec.kind == "kurvst-mark")
    assert(spec.len() == 2)
    assert(mark.prepare(spec) == (shape: spec.shape))
  }
  for constructor in constructors.filter(constructor => (
    constructor != mark.rays
  )) {
    assert(not ("width" in constructor(width: auto)))
  }
  assert(not ("phase" in mark.rays(phase: auto)))
  assert(not ("length" in mark.bracket(length: auto)))
  assert(mark.prepare(mark.bar(stroke: stroke(paint: red))).stroke)
  assert(mark.triangle(fill: false).fill == false)
  let spec = mark.triangle(
    length: 3pt + 450%,
    width: 2pt,
    inset: 25%,
    shorten: 100%,
    fill: red,
    stroke: none,
    rev: true,
  )
  let prepared = mark.prepare(spec)
  assert(spec.fill == red and spec.stroke == none)
  assert(prepared.length == (points: 3.0, ratio: 4.5))
  assert(prepared.width == (points: 2.0, ratio: 0.0))
  assert(prepared.inset == 0.25 and prepared.shorten == 1.0)
  assert(prepared.fill and not prepared.stroke and prepared.rev)
  assert(not ("kind" in prepared))
  assert(
    mark.prepare(mark.barb(width: 200%, arc: 180deg)).width
      == (points: 0.0, ratio: 2.0),
  )
  assert(mark.prepare(mark.barb(arc: 180deg)).arc == calc.pi)
  assert(
    mark.prepare(mark.rays(n: 3, phase: 90deg, align: "end")).phase
      == calc.pi / 2,
  )
  assert(
    mark.prepare(mark.circle(fill: auto, stroke: auto)) == (shape: "circle"),
  )
  let composite = mark.combine(
    mark.triangle(),
    1pt,
    25%,
    mark.combine(mark.bar(), 2pt + 50%, mark.circle()),
    fit: "bend",
  )
  let numeric = mark.prepare(composite)
  assert(numeric.parts.at(1) == (gap: (points: 1.0, ratio: 0.0)))
  assert(numeric.parts.at(2) == (gap: (points: 0.0, ratio: 0.25)))
  assert(numeric.parts.at(3).parts.at(1) == (gap: (points: 2.0, ratio: 0.5)))
  assert(type(cbor.encode(numeric)) == bytes)
  assert(type(cbor.encode(prepared)) == bytes)
} else {
  let failures = (
    positional: () => mark.triangle(3pt),
    unsupported: () => mark.straight(fill: red),
    tikz-arc: () => mark.tikz(arc: 30deg),
    symbol: () => mark.triangle(symbol: "triangle"),
    scale: () => mark.triangle(scale: 2),
    anchor: () => mark.triangle(anchor: "tip"),
    shorten-to: () => mark.triangle(shorten-to: 100%),
    cetz: () => mark.prepare((symbol: "triangle")),
    n: () => mark.rays(n: 0),
    n-fraction: () => mark.rays(n: 1.5),
    align: () => mark.bar(align: "start"),
    fit: () => mark.triangle(fit: "curve"),
    nonfinite: () => mark.triangle(width: float("inf") * 1pt),
    bad-context: () => mark.geometry(
      ((mark: mark.bar(), "context": (units-per-pt: 0, line-thickness: 1)),),
      (),
      (),
    ),
    thickness: () => mark.geometry(
      (
        (
          mark: mark.bar(),
          "context": (units-per-pt: 1, line-thickness: float("inf")),
        ),
      ),
      (),
      (),
    ),
    child-fit: () => mark.combine(mark.bar(fit: "bend")),
    missing-parts: () => mark.combine(),
    named-parts: () => mark.combine(parts: (mark.bar(),)),
    fill-content: () => mark.triangle(fill: [paint]),
    fill-function: () => mark.triangle(fill: () => red),
    stroke-dictionary: () => mark.bar(stroke: (meaningless: true)),
    stroke-thickness: () => mark.bar(stroke: red + 2pt),
    stroke-cap: () => mark.bar(stroke: stroke(paint: red, cap: "round")),
    stroke-join: () => mark.bar(stroke: stroke(paint: red, join: "bevel")),
    stroke-dash: () => mark.bar(stroke: stroke(paint: red, dash: "dashed")),
    stroke-miter: () => mark.bar(stroke: stroke(paint: red, miter-limit: 3)),
  )
  assert(mode in failures, message: "unknown mark contract mode")
  failures.at(mode)()
}
Mark contract passed.
